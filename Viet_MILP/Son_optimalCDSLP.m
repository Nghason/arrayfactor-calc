n_links    = fabric.numFabricPorts;     % nb of fabric ports (ingress+egress)    
n_machines = n_links/2;                 % nb of machines on the fabric
n_coflows  = length(coflows);           % nb of coflows
n_flows    = [coflows.numFlows];        % nb of flows of each coflow
n_flows_all = sum(n_flows);             % total nb of flows

portCapacity = [[fabric.machinesPorts.ingress] [fabric.machinesPorts.egress]];
portCapacity = [portCapacity.linkCapacity];



time_instants = [0];
%   Create time instant array [0, T1, T2, ... Tn]
for c = coflows
   time_instants = [time_instants, c.deadline];
   c.addParam.cstatus = "alive";
end
%   Clear duplication and sort the order (min-max)
time_instants = unique(time_instants);
time_instants = sort(time_instants);
n_slots = length(time_instants) - 1; % n_coflows; % nb of time slots
time_slots = cell(1, n_slots); 
slot_duration = zeros(1, n_slots); 

%   Tạo 1 cell chứa thông tin điểm đầu, cuối của timeslot.
%   Tạo 1 mảng chứa khoảng thời gian delta của các timeslot
for i = 1:n_slots
   time_slots{i} = [time_instants(i), time_instants(i+1)];
   slot_duration(i) = time_instants(i+1) - time_instants(i);
end

%   Tạo 1 struct chứa 2 trường time instance và time slot
LP.time_instants = time_instants;
LP.time_slots = time_slots;
n_slots_of_coflow = zeros(1, n_coflows);
%   Tạo mảng chứa timeslot gặp deadline của mỗi coflow, tạo thêm trường
%   timeslot gặp deadline (k_deadline), tạo thêm trường với kiểu dữ liệu là
%   1 mảng (slot_id) gồm các timeslot mà coflow sẽ đi qua để tới deadline.
for c = coflows
 c.addParam.k_arrival = 1;    % all coflows arrive at t=0
    for i = 1:n_slots
        A = time_slots{i}(2);
        if c.deadline == time_slots{i}(2)
            c.addParam.k_deadline = i;  % Thêm trường k_deadline: coflow gặp deadline tại timeslot thứ mấy
            % fprintf("coflow %d: slot %d\n", c.id, i);
        end
    end
    n_slots_of_coflow(c.id) = c.addParam.k_deadline; 
    c.addParam.slot_id = (c.addParam.k_arrival:c.addParam.k_deadline); % slot indexes of coflow c
    % fprintf("coflow %d:\n", c.id); disp(c.addParam.slot_id);
end

%   Tạo id đánh dấu index min và max của các slot (chưa hiểu)
id_slot_mn = 1;
id_slot_mx = 1;
for c = coflows
    if id_slot_mn > min(c.addParam.slot_id)
        id_slot_mn = min(c.addParam.slot_id);
    end
    if id_slot_mx < max(c.addParam.slot_id)
        id_slot_mx = max(c.addParam.slot_id);
    end 
end

% Nb of variables tổng các biến
nvar_bijm = sum(n_flows .* n_slots_of_coflow);   % total nb of variables bijm
nvar_all = n_coflows + nvar_bijm;           % total nb of variables

% Tạo 1 cell gồm 20 coflows, mỗi coflow chứa n cell flow, mỗi flow chứa
% mảng zero có kích cỡ (1, n_slot_of_coflows) n_slot... là index mà coflow
% đó gặp deadline.
bijm_cell = cell(1, n_coflows);
for ci = 1:n_coflows
    bijm_cell{ci} = cell(1, n_flows(ci));
    for fi = 1:n_flows(ci)
        bijm_cell{ci}{fi} = zeros(1, n_slots_of_coflow(ci));
    end
end
%   Option: Weight coflow đều = 1 nên obj là như nhau, khác nhau zn_type
mod_name = 'w_CdsRelax';
switch mod_name
    case {'w_CdsOptim','w_CdsRelax'}    %  weighted objective: sum_{k in C} w_k*z_k
        obj = [-([coflows.weight]'); zeros(nvar_bijm, 1)];
        
    case {'CdsOptim','CdsRelax'}        % non-weighted objective: sum_{k in C} z_k
        obj = [-ones(n_coflows, 1); zeros(nvar_bijm, 1)];
end
switch mod_name
    case {'w_CdsOptim','CdsOptim'}    
        zn_type = 'B';
    case {'w_CdsRelax','CdsRelax'}    
        zn_type = 'C';
end
%%   TẠO MẢNG C1 = [ ma trận đường chéo volume| ma trận "đường chéo" slot
%   duration của các flow thuộc coflow 
% any_coflow_dead = 1;
%while any_coflow_dead 
C1 = [];
for c = coflows 
    cid = c.id;
    for f = c.flows
        fid = f.id;
        % zn: volume of flow f_j; %   Mỗi vòng lặp tạo 1 mảng zn mới
        zn = zeros(1, n_coflows);
        zn(cid) = f.volume;
        if c.addParam.cstatus == "dead" % Note: flow bị bỏ
            zn(cid) = 0;
        end
        % bijm: -slot_duration
        bijm = bijm_cell;   
        for mi = 1:n_slots_of_coflow(cid)
            slot_id = c.addParam.slot_id(mi);
            % Gán slot duration ứng với các timeslot của các flow
            bijm{cid}{fid}(mi) = -slot_duration(slot_id);
        end
        
        % Transform bijm to a row vector, Tạo mới bijm mỗi loop
        bijm_transform = [];
        for id = 1:n_coflows
            bijm_transform = [bijm_transform, [bijm{id}{:}]];
        end

        % Append zn and bijm to C1
        C1 = [C1; zn, bijm_transform];
        
    end % end loop of fi
end % end loop of ci
%%
C1_R = zeros(size(C1,1), 1);

C2 = [];
C2_R = [];

for p = 1:n_links   % foreach fabric port (ingress or egress link)
    for mi = 1:n_slots % foreach slot in the time horizon [0,T]
         
        bijm = bijm_cell;
        
        for c = coflows
            ci = c.id;
            for f = c.flows
                fi = f.id;
                % Kiểm tra coflow c có đi qua timeslot mi và kiểm tra flow
                % có đi qua port p:
                if ismember(mi, c.addParam.slot_id) && ismember(p, f.links)
                    bijm{ci}{fi}(mi) = 1;   % Nếu có gán giá trị bijm{c}{f} tại timeslot đó = 1 
                end
            end
        end
        
        % Transform bijm to a row vector
        bijm_transform = [];
        for id = 1:n_coflows
            bijm_transform = [bijm_transform, [bijm{id}{:}]];
        end
        
        % Append zn and bijm to C2
        if any(bijm_transform) % only append if bijm_transform is not null
            C2 = [C2; zeros(1, n_coflows), bijm_transform];
            
            % Right-hand of C2
            C2_R = [C2_R; portCapacity(p)];
        end
    end % end loop of mi
end % end loop of e
% Vai trò Ctype?
ctype = strcat( char(ones([1, n_coflows])*(zn_type)), ... % type of zn
                char(ones([1, nvar_bijm])*('S')));
switch mod_name
    case {'w_CdsOptim','CdsOptim'}   
        lb = [-inf*ones(1, n_coflows), zeros(1, nvar_bijm)];
        ub = [inf*ones(1, n_coflows), inf*ones(1, nvar_bijm)];
        
    case {'w_CdsRelax','CdsRelax'}       
        lb = zeros(1, nvar_all);
        ub = [ones(1, n_coflows), inf*ones(1, nvar_bijm)];
end
%%
[Aineq, bineq, Aeq, beq] = deal(C2, C2_R, C1, C1_R);     
%[x,fval,exitflag,output] = cplexmilp(obj,Aineq,bineq,Aeq,beq,[],[],[],lb,ub,ctype,x0,options);
[x, fval, exitflag, output] = linprog(obj,Aineq,bineq,Aeq,beq,lb,ub);
any_coflow_dead = 0;
for c_reject = coflows  % Check Zn, nhỏ hơn 0.99 thì loại bỏ
    c_reject_id = c_reject.id;
    if x(c_reject_id,1) < 0.99
        c_reject.addParam.cstatus = "dead"; % status coflow loại bỏ
        any_coflow_dead = 1;    % Chưa hiểu tại sao lại set về 1, tại sao lại phải để x = 1
        break;
    end 
end
%%
%end
%%
if exitflag > 0 %isempty(x) == 0
    % Rearrange result
    zn = x(1:n_coflows)';
    bijm_vect = x(n_coflows+1 : end)';
    bijm = bijm_cell;

    k = 1;
    for ci = 1:n_coflows
        for fi = 1:n_flows(ci)
            for mi = 1:n_slots_of_coflow(ci)
                bijm{ci}{fi}(mi) = bijm_vect(k); % Gán kết quả xj(delta_m) vào bijm 
                k = k + 1;
            end
        end
    end
end
cct_vect = cell(1, n_coflows);  % slots that coflow traverses
    cct_slot = zeros(1, n_coflows); % slot at which coflow finishes
    ccts = zeros(1, n_coflows);     % cct of coflow
    
    for ci = 1:n_coflows
        cct_vect{ci} = [];
        if zn(ci) >= 0.99 % for safety
            for fi = 1:n_flows(ci)
                fi_ct = find(bijm{ci}{fi} > 0, 1, 'last'); % completion time of flow fi
                cct_vect{ci} = [cct_vect{ci}, fi_ct];
            end
            if isempty(fi_ct) 
                continue;
            end
            cct_slot(ci) = max(cct_vect{ci});
            ccts(ci) = time_slots{cct_slot(ci)}(2);
        end
    end
        