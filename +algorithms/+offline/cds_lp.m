% clear;
% load"data_m10_c20";
%
function outputs = cds_lp(fabric,coflows)
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
n_slots = length(time_instants) - 1; % n_coflows; % nb of time slots
time_slots = cell(1, n_slots); 
slot_duration = zeros(1, n_slots); 
for i = 1:n_slots
   time_slots{i} = [time_instants(i), time_instants(i+1)];
   slot_duration(i) = time_instants(i+1) - time_instants(i);
end

n_slots_of_coflow = zeros(1, n_coflows);

for c = coflows
	c.addParam.k_arrival = 1;    % all coflows arrive at t=0
    for i = 1:n_slots
        A = time_slots{i}(2);
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
cct_var = intvar(n_coflows,size(time_slots,2));

constraints = [ ];
for c = coflows
    cid = c.id;
    for f = c.flows
        fid = f.id;
        if cid ~= 1
            fid = fid + flows_cumsum(cid - 1);
        end
        for i = 1:n_slots_of_coflow(cid)% Chạy từng timeslot
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
C1 = send_rate .* slot_duration_per_coflows;
%% Constraint 1
count = 0;
for c = coflows % coflow
    cid = c.id;
    for f = c.flows % coflow/flow
        fid = f.id;
        if cid ~= 1
            fid = fid + flows_cumsum(cid - 1);
        end
        f_vol = f.volume;
        for i = 1:n_links % coflow/flow/port
            if ismember(i,f.links)
                if count == 1
                    constraints = [constraints,zn(cid)*f_vol - sum(C1(fid + port_cumsum(i), :)) == 0];
                end
                if count == 0
                if i == 1 % coflow/flow/port/data sum
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
% C2
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

%
constraints = [constraints, send_rate >= 0];
constraints = [constraints, C2 >= 0];
constraints = [constraints, C1 >= 0];
addpath('E:\Study_corner\Do_an_Tot_Nghiep\Coflow\Solver\CPLEX_Studio\cplex\matlab\x64_win64');
ops = sdpsettings('solver','cplex','verbose',0);
% cct_var_for_obj = sdpvar;
% for i = 1:n_coflows
%     cct_var_for_obj(1,i) = sum(cct_var(i,:));
% end
% obj = -sum(zn(:,1))-(1/(n_links))*sum(cct_var_for_obj(:,1));
obj = -sum(zn(:,1));
result = optimize(constraints,obj,ops);  
z_x = value(cct_var);
% x_1 = value(cct_var_1);
% x_2 = value(cct_var_2);
% z_Xj = value(C3_1);
Zn = value(zn);
z_c2 = value(C2);
z_c1 = value(C1);
Z_sr = value(send_rate);
% successfully solved
% if result.problem == 0 
%     value(send_rate)
%     % value(C1)
%     %Value(C2)
%     value(zn)
%     -value(obj)   
% else
%     disp('Wrong!');
% end

%% CCT average:
cct_avg= zeros(2,n_coflows+1);
Zn = value(zn);
Tslot = 0;
for c = coflows 
    cid = c.id;
   % Kiểm tra xem coflow có được accept trong Zn hay không?
   if Zn(cid) ~= 0
   for f = c.flows
       fid = f.id;
       % Kiểm tra flow đó nằm tại port nào, xác định vị trí trong Z_c2
       % port (in and out):
       for i = 1:2
       port_fl = f.links(i);
       if port_fl == 1 % Lưu ý khi xử lý port và coflow 1 
          if cid == 1
              row = fid;
          else
              row = flows_cumsum(cid -1) + fid;
          end
       else
          if cid == 1
          row = port_cumsum(port_fl) + fid;
          else
          row = port_cumsum(port_fl)+flows_cumsum(cid-1) + fid;
          end
       end     
       Tslot = find(z_c2(row, :) ~= 0, 1, 'last' );
       if Tslot > cct_avg(1,cid)
          cct_avg(1,cid) = Tslot; 
       end
       end           
   end 
   end
end
% Calc avg_cct at the end
for i = 1:n_coflows
   if cct_avg(1,i) ~= 0
       indx = cct_avg(1,i);
       cct_avg(2,i) = time_instants(indx+1);
   end
end
cnt = 0;
for i = 1:size(cct_avg,2)
    if cct_avg(1,i) ~= 0
        cnt = cnt + 1;
    end
end
cct_avg(2,n_coflows+1) = sum(cct_avg(2,1:n_coflows))/cnt;
cct_avg(2,n_coflows+1);
Zn = Zn';
% Output
    outputs.zn = Zn;
    outputs.accepts = find(Zn >= 0.99); % id of accepted coflows
    outputs.rejects = find(Zn < 0.1);   % id of rejected coflows   
    outputs.cct = cct_avg(:, 1:n_coflows);
    outputs.ccts = cct_avg(2,find(cct_avg(2, 1:n_coflows) > 0));
    outputs.avg_cct = cct_avg(2, n_coflows+1);
end