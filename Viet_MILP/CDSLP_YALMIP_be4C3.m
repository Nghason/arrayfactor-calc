
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
    else if i == 2
        port_cumsum(i) = port_cumsum(i-1) -1 + n_flows_all;
    else 
        port_cumsum(i) = port_cumsum(i-1) + n_flows_all;
    end
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
            % fprintf("coflow %d: slot %d\n", c.id, i);
        end
    end
    n_slots_of_coflow(c.id) = c.addParam.k_deadline; 
    c.addParam.slot_id = (c.addParam.k_arrival:c.addParam.k_deadline); % slot indexes of coflow c
    % fprintf("coflow %d:\n", c.id); disp(c.addParam.slot_id);
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
                    %C1_2 = sum(C1(fid + port_cumsum(i), :));
                end
                if count == 0
                if i == 1
                    constraints = [constraints,zn(cid)*f_vol - sum(C1(fid, :)) == 0];
                    % C1_1 = sum(C1(fid, :));
                else
                    constraints = [constraints,zn(cid)*f_vol - sum(C1(fid + port_cumsum(i), :)) == 0];
                    % C1_1 = sum(C1(fid + port_cumsum(i), :));
                end
                count = count + 1;
                end
            end
        end
        count = 0;
%         constraints = [constraints,zn(cid)*f_vol - C1_1 == 0];
%         constraints = [constraints,zn(cid)*f_vol - C1_2 == 0];
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
C3 = sdpvar;
count = 0;
for c = coflows
    cid = c.id;
    for i = 1:n_slots
        for j = 1:n_links
        for f = c.flows
        fid = f.id;
        if cid ~= 1
            fid = fid + flows_cumsum(cid - 1);
        end
        if ismember(j,f.links)
                if count == 1
                    C3(2,fid) = send_rate(fid + port_cumsum(j),i);
                end
                if count == 0
                if j == 1
                    C3(1,fid) = send_rate(fid, i);
                else
                    C3(1,fid) = send_rate(fid + port_cumsum(j), i);
                end
                if fid + 1 > c.numFlows
                count = count + 1;
                end
                end
        end
        end
        end
        constraints = [constraints, sum(C3(1,:))/sum(c.addParam.volumes(1,:)) - cct_var(cid,i) >= 0];
        constraints = [constraints, sum(C3(2,:))/sum(c.addParam.volumes(1,:)) - cct_var(cid,i) >= 0];
        C3 = sdpvar;
        count = 0;
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
obj = -sum(zn(:,1))
result = optimize(constraints,obj,ops);  
% successfully solved
if result.problem == 0 
    sr = value(send_rate)
    % value(C1)
    %Value(C2)
    z = value(zn)
    -value(obj)   
else
    disp('Wrong!');
end
