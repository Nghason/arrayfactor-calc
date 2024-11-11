% clear;
% load"data_m4_c6";
%%
n_links    = fabric.numFabricPorts;     % nb of fabric ports (ingress+egress)    
n_machines = n_links/2;                 % nb of machines on the fabric
n_coflows  = length(coflows);           % nb of coflows
n_flows    = [coflows.numFlows];        % nb of flows of each coflow
n_flows_all = sum(n_flows);             % total nb of flows

portCapacity = [[fabric.machinesPorts.ingress] [fabric.machinesPorts.egress]];
portCapacity = [portCapacity.linkCapacity];
port_cumsum = [];
for i = 1:n_links 
    if i == 1
        port_cumsum(i) = 1;
    elseif i == 2
        port_cumsum(i) = port_cumsum(i-1) -1 + n_flows_all;
    else 
        port_cumsum(i) = port_cumsum(i-1) + n_flows_all;
    end
end
flows_cumsum = cumsum(n_flows);
time_instants = [0];
for c = coflows
   time_instants = [time_instants, c.deadline];
end
time_instants = unique(time_instants);
time_instants = sort(time_instants);
% SMOOTH
num_subsegments = 2;
time_instants_subsegment = [];
for i = 1:length(time_instants) - 1;
    subsegment = linspace(time_instants(i), time_instants(i + 1), num_subsegments +1);
    if i ~=length(time_instants) - 1
        subsegment = subsegment(1:end - 1);
    end
    time_instants_subsegment = [time_instants_subsegment, subsegment];
end
n_slots =  size(time_instants_subsegment, 2) - 1;
time_slots = cell(1, n_slots); 
slot_duration = zeros(1, n_slots);


for i = 1:n_slots
   time_slots{i} = [time_instants_subsegment(i), time_instants_subsegment(i+1)];
   slot_duration(i) = time_instants_subsegment(i + 1) - time_instants_subsegment(i);
end

n_slots_of_coflow = zeros(1, n_coflows);

for c = coflows
	c.addParam.k_arrival = 1;    % all coflows arrive at t=0
    for i = 1:n_slots
        if c.deadline == time_slots{i}(2)
            c.addParam.k_deadline = i;
        end
    end
    n_slots_of_coflow(c.id) = c.addParam.k_deadline; 
    c.addParam.slot_id = (c.addParam.k_arrival:c.addParam.k_deadline); % slot indexes of coflow c
end
slot_duration_per_coflows = zeros(n_links * n_flows_all,n_slots);
zn = binvar(n_coflows,1);
send_rate = sdpvar(n_links * n_flows_all,n_slots);
cct_var = binvar(n_coflows,size(time_slots,2));

constraints = [ ];
for c = coflows
    cid = c.id;
    for f = c.flows
        fid = f.id;
        if cid ~= 1
            fid = fid + flows_cumsum(cid - 1);
        end
        for i = 1:n_slots_of_coflow(cid)% Chạy từng timeslot
%             for i = 1:n_slots
            for j = 1:n_links % Chạy từng port
                if ismember(j,f.links) % Kiểm tra flow có thuộc port ko
                    if j == 1
                        slot_duration_per_coflows(fid,i) = slot_duration(i);
                    else
                        slot_duration_per_coflows(fid + port_cumsum(j),i) = slot_duration(i);
                    end
                end
            end
        end
    end
end
%ma trận C1 sẽ cho ra dữ liệu ở đúng vị trí của các flow ứng với port và
%time slot của nó còn những ô khác thì sẽ có giá trị 0
C1 = send_rate .* slot_duration_per_coflows;
%% Constraint 1
count = 0;
for c = coflows
    cid = c.id;
    for f = c.flows
        fid = f.id;
        if cid ~= 1
            fid = fid + flows_cumsum(cid - 1);
        end
        f_vol = f.volume;
        for i = 1:n_links
            if ismember(i,f.links)
                if count == 1
                    constraints = [constraints,zn(cid)*f_vol - sum(C1(fid + port_cumsum(i), :)) == 0];
                end
                if count == 0
                if i == 1
                    constraints = [constraints,zn(cid)*f_vol - sum(C1(fid, :)) == 0];
                else
                    constraints = [constraints,zn(cid)*f_vol - sum(C1(fid + port_cumsum(i), :)) == 0];
                end
                count = count + 1;
                end
            end
        end
        count = 0;
    end
end
%% C2
send_data_location = zeros(n_links * n_flows_all,n_slots);
for c = coflows
    cid = c.id;
    for f = c.flows
        fid = f.id;
        if cid ~= 1
            fid = fid + flows_cumsum(cid - 1);
        end
        for i = 1:n_slots_of_coflow(cid)
            for j = 1:n_links
                if ismember(j,f.links)
                    if j == 1
                        send_data_location(fid,i) = 1;
                    else
                        send_data_location(fid + port_cumsum(j),i) = 1;
                    end
                end
            end
        end
    end
end

C2 = send_rate .* send_data_location;
for p = 1:n_links
    for slot_id = 1:n_slots
        if p == 1 
            start_idx = port_cumsum(p);
            end_idx = port_cumsum(p + 1);
        else if p ==n_links
            start_idx = port_cumsum(p) + 1;
            end_idx = port_cumsum(p) + n_flows_all;
        else 
            start_idx = port_cumsum(p) + 1;
            end_idx = port_cumsum(p + 1);
        end
        end
        constraints = [constraints, sum(C2(start_idx:end_idx, slot_id)) <= portCapacity(p)];
    end
end
%%
%% C3
%ma trận C3 sẽ chứa các giá trị send_rate theo thứ tự 
C3 = sdpvar;
count = 1;
for c = coflows
    cid = c.id;
    for f =c.flows
        fid = f.id;
        if cid ~= 1
            fid = fid + flows_cumsum(cid - 1);
        end
        for i = 1:n_links
            for j = 1:n_slots_of_coflow(1,cid)
            if ismember(i,f.links)
                if i == 1
                    C3(count,j) = send_rate(fid,j);
                else
                    C3(count,j) = send_rate(fid + port_cumsum(i),j);
                end
            end
            end
            if ismember(i,f.links)
                count = count + 1;
            end
        end
    end
end
%biến nhớ tạm cho phần xử lý
cache_mem_for_C3 = 0;
%biến Xj theo định dạng yalmip
C3_1 = sdpvar;
count = 1;
for c = coflows
    cid = c.id;
    c_vol = c.addParam.volumes(1,:);
    c_vol_value = sum(c_vol(:,:));
    c_vol_value = c_vol_value * 2;
    if cid == 1
        start_idx_C3 = 1;
        end_idx_C3 = n_flows(cid)*2 + start_idx_C3 - 1;
        cache_mem_for_C3 = end_idx_C3; 
    else
        start_idx_C3 =  cache_mem_for_C3 + 1;
        end_idx_C3 = n_flows(cid)*2 + start_idx_C3 - 1;
        cache_mem_for_C3 = end_idx_C3; 
    end
    for i = 1:1:n_slots
        if i <= n_slots_of_coflow(cid)
            if i ==1
                data_cumsum_var = sum(C3(start_idx_C3:end_idx_C3, i)) * slot_duration(1,i);
            else
                data_cumsum_var = data_cumsum_var + sum(C3(start_idx_C3:end_idx_C3, i)) * slot_duration(1,i);
            end
            C3_1(cid, i) = data_cumsum_var/c_vol_value;
        else
            C3_1(cid, i) = 1 * zn(cid);
        end
         constraints = [constraints, C3_1(cid, i) - cct_var(cid, i) >= 0];
%          constraints = [constraints, C3_1(cid, i) <= 1 - 0.01 + cct_var(cid, i) * 0.01];
    end
end


%%
constraints = [constraints, send_rate >= 0];
constraints = [constraints, C2 >= 0];
constraints = [constraints, C1 >= 0];
ops = sdpsettings();
cct_var_for_obj = sdpvar;
for i = 1:n_coflows
    cct_var_for_obj(1,i) = sum(cct_var(i,:));
end
 obj = -sum(zn(:,1))-(1/(n_links))*sum(cct_var_for_obj(1,:));
result = optimize(constraints,obj,ops);  
x = value(cct_var);
Xj = value(C3_1);
Zn = value(zn);
c2 = value(C2);
c1 = value(C1);

%% CCT_average
zXj = value(C3_1);
epsilon = 1e-6;
cct_avgsm = zeros(2,n_coflows+1);
for n = 1:n_coflows
    one = find(abs(zXj(n,:) - 1) < epsilon);
    if isempty(one)
    continue;
    else
        cct_avgsm(1,n) = one(1);
        cct_avgsm(2,n) = time_instants_subsegment(one(1)+1);
    end
end
% calc cct at the end
% for i = 1:n_coflows
%    if cct_avg3(1,i) ~= 0
%        indx = cct_avg3(1,i);
%        cct_avg3(2,i) = time_instants(indx+1);
%    end
% end
cnt = 0;
for i = 1:size(cct_avgsm,2)
    if cct_avgsm(1,i) ~= 0
        cnt = cnt + 1;
    end
end
cct_avgsm(2,n_coflows+1) = sum(cct_avgsm(2,1:n_coflows))/cnt;
cct_avgsm(2,n_coflows+1);


