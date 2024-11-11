% (c) 2022 Quang-Trung Luu (LAAS-CNRS)
% Implementation of linear programs CDS-LP and CDS-LPA
% Reference:    S.-H. Tseng and A. Tang, “Coflow deadline scheduling via network-aware optimization,” 
%               in Proc. Annu. Allert. Conf. Commun. Control Comput., 2018, pp. 829–833.
% Remarks:
% - [0,T] is divided into M disjoint time slots such that their endpoints are either a_j or d_j if not 0 or T
% Inputs:
% - fabric: the Fabric ojbect representing the network
% - coflows: an array of Coflow objects within the fabric
% - zn_type: the type of variable z_{n}: integer for CDS-LP, and continuous for CDS-LPA
% Outputs:
% - resource allocation scheme for the given coflows


function outputs = optim_cds(inputs, fabric, coflows, mod_name)


global params;

% Type of zn ('B' for CDS-LP and 'C' for CDS-LPA)
switch mod_name
    case {'w_CdsOptim','CdsOptim'}    
        zn_type = 'B';
    case {'w_CdsRelax','CdsRelax'}    
        zn_type = 'C';
end


%% Initialization
n_links    = fabric.numFabricPorts;     % nb of fabric ports (ingress+egress)    
n_machines = n_links/2;                 % nb of machines on the fabric
n_coflows  = length(coflows);           % nb of coflows
n_flows    = [coflows.numFlows];        % nb of flows of each coflow
n_flows_all = sum(n_flows);             % total nb of flows

portCapacity = [[fabric.machinesPorts.ingress] [fabric.machinesPorts.egress]];
portCapacity = [portCapacity.linkCapacity];

% Find bottlenecks:
% port_load = zeros(n_links, n_coflows); % full rate processing time of each coflow on each link
% for c = coflows
%     port_load(:,c.id) = (c.indicator*c.getFlowsVolume')./portCapacity'; % load on each port of each coflow 
% end
% cum_port_load = sum(port_load,2);
% bottleneck = max(cum_port_load);
% bottle_links = find(cum_port_load == bottleneck);
% epsilon = 0;
% T_bottleneck = (1 + epsilon)*max(bottleneck./portCapacity(bottle_links));
% T_bottleneck = round(T_bottleneck, 3);


%% Temporal parameters of coflows
% Time instants [0, T_1,..., T_N] (T = T_N for simplicity)
time_instants = [0];
for c = coflows
   time_instants = [time_instants, c.deadline];
end
time_instants = unique(time_instants);
time_instants = sort(time_instants); % ascending sort

% Time slots
n_slots = length(time_instants) - 1; % n_coflows; % nb of time slots
time_slots = cell(1, n_slots); 
slot_duration = zeros(1, n_slots);   % duration of each time slot (not fixed!)

for i = 1:n_slots
   time_slots{i} = [time_instants(i), time_instants(i+1)];
   slot_duration(i) = time_instants(i+1) - time_instants(i);
end

LP.time_instants = time_instants;
LP.time_slots = time_slots;

% Slot assignment for coflows
n_slots_of_coflow = zeros(1, n_coflows); % nb of slots of each coflow
for c = coflows
	c.addParam.k_arrival = 1;    % all coflows arrive at t=0
    for i = 1:n_slots
        if c.deadline == time_slots{i}(2)
            c.addParam.k_deadline = i;
            % fprintf("coflow %d: slot %d\n", c.id, i);
        end
    end
    n_slots_of_coflow(c.id) = c.addParam.k_deadline; 
    c.addParam.slot_id = (c.addParam.k_arrival:c.addParam.k_deadline); % slot indexes of coflow c
    % fprintf("coflow %d:\n", c.id); disp(c.addParam.slot_id);
end


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

% Nb of variables
nvar_bijm = sum(n_flows .* n_slots_of_coflow);   % total nb of variables bijm
nvar_all = n_coflows + nvar_bijm;           % total nb of variables

%% CONSTRUCTION OF THE LP   

% An empty cell for all ci-fi-mi to be used widely
bijm_cell = cell(1, n_coflows);
for ci = 1:n_coflows
    bijm_cell{ci} = cell(1, n_flows(ci));
    for fi = 1:n_flows(ci)
        bijm_cell{ci}{fi} = zeros(1, n_slots_of_coflow(ci));
    end
end


%% Objective function (column vector)

% Append zn and bijm to f: max(sum of zn) -> min(-sum of zn)
switch mod_name
    case {'w_CdsOptim','w_CdsRelax'}    %  weighted objective: sum_{k in C} w_k*z_k
        obj = [-([coflows.weight]'); zeros(nvar_bijm, 1)];
        
    case {'CdsOptim','CdsRelax'}        % non-weighted objective: sum_{k in C} z_k
        obj = [-ones(n_coflows, 1); zeros(nvar_bijm, 1)];
end

%% C1: Demand constraint for each flow
% sum_{all slots m of flows fj in coflow i} (x_j(Delta_m)*|Delta_m|) = s_j*z^n
% volume(fj)*zn - sum(bijm(i,j,m) * slot_duration) = 0, forall i, j
% expected nb of rows of C1: n_flows_all
C1 = [];
for c = coflows 
    cid = c.id;
    
    for f = c.flows
        fid = f.id;
        
        % zn: volume of flow f_j
        zn = zeros(1, n_coflows);
        zn(cid) = f.volume;
    
        % bijm: -slot_duration
        bijm = bijm_cell;
        for mi = 1:n_slots_of_coflow(cid)
            slot_id = c.addParam.slot_id(mi);
            bijm{cid}{fid}(mi) = -slot_duration(slot_id);
        end
        
        % Transform bijm to a row vector
        bijm_transform = [];
        for id = 1:n_coflows
            bijm_transform = [bijm_transform, [bijm{id}{:}]];
        end

        % Append zn and bijm to C1
        C1 = [C1; zn, bijm_transform];
        
    end % end loop of fi
end % end loop of ci

% Right-hand of C1
C1_R = zeros(size(C1,1), 1);


%% C2: Link capacity limit constraint
% for each link, all flows traversing that link should not exceed the link capacity
% sum_(all fi of all ci that use link e) bijm <= bandwidth(e), forall link e, slot m

C2 = [];
C2_R = [];

for p = 1:n_links   % foreach fabric port (ingress or egress link)
    for mi = 1:n_slots % foreach slot in the time horizon [0,T]
         
        bijm = bijm_cell;
        
        for c = coflows
            ci = c.id;
            for f = c.flows
                fi = f.id;
                % Check if slot_id of flow fi contains mi and flow fi uses port e
                if ismember(mi, c.addParam.slot_id) && ismember(p, f.links)
                    bijm{ci}{fi}(mi) = 1;
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


%% Optimizer setup

% Types of variables
ctype = strcat( char(ones([1, n_coflows])*(zn_type)), ... % type of zn
                char(ones([1, nvar_bijm])*('S')));        % type of bijm
            
% Bound: zn in {0,1}, bijm in [0, inf]
switch mod_name
    case {'w_CdsOptim','CdsOptim'}   
        lb = [-inf*ones(1, n_coflows), zeros(1, nvar_bijm)];
        ub = [inf*ones(1, n_coflows), inf*ones(1, nvar_bijm)];
        
    case {'w_CdsRelax','CdsRelax'}       
        lb = zeros(1, nvar_all);
        ub = [ones(1, n_coflows), inf*ones(1, nvar_bijm)];
end


[Aineq, bineq, Aeq, beq] = deal(C2, C2_R, C1, C1_R);       

LP.Aineq = Aineq;         
LP.bineq = bineq; 
LP.Aeq = Aeq;
LP.beq = beq;
LP.x0 = [];

% Save parameters of optimizer
LP.solver   = params.solver;
LP.f        = obj;
LP.C1       = C1;
LP.C2       = C2;
LP.C1_R     = C1_R;
LP.C2_R     = C2_R;
LP.lb       = lb;
LP.ub       = ub;
LP.ctype    = ctype;

%======================================
% Solving
%======================================

t1 = cputime;

if strcmp(LP.solver, 'cplexmilp')  
    %[x,fval,exitflag,output] = cplexmilp(f,Aineq,bineq,Aeq,beq,[],[],[],lb,ub,ctype,x0,options);
    
    LP.options = []; %opt_options.options; 
    
    % options = cplexoptimset('cplex');
    % LP.timelimit = opt_options.MaxTime;  
    % LP.mip.tolerances.mipgap = 5e-3; % 0-1 default 1e-4
    % LP.mip.limits.nodes = opt_options.MaxNodes;   
    % options.exportmodel = '+ResourceAllocationAlgos/+Trung/ifip/cplex_model.lp';
    % LP.options = options;

    if strcmp(zn_type, 'B') || strcmp(zn_type, 'I') || strcmp(zn_type, 'N') % CDS-LP
        [x, fval, exitflag, output] = cplexmilp(LP);
    elseif strcmp(zn_type, 'C') || strcmp(zn_type, 'S')% CDS-LPA
        [x, fval, exitflag, output] = cplexmilp(LP); %cplexlp(LP);
    end
    outputs.output = output;
    
elseif strcmp(LP.solver, 'gurobi')  
    model.obj = obj;
    model.A   = sparse([Aineq; Aeq]);
    model.rhs = [bineq; beq];
    model.lb  = lb;
    model.ub  = ub;
    model.sense = [repmat('<', size(Aineq,1), 1); repmat('=', size(Aeq,1), 1)];
    model.vtype = ctype;
    model.modelsense = 'min';
    
    params.OutputFlag = 0; % disable solver output log
    %params.LogToConsole = 0; % disable solver output log
    
      
    names = {};
    for ci = 1:n_coflows
        names{end+1} = {strcat('z', num2str(ci))}; 
    end

    for ci = 1:n_coflows
        for fi = 1:n_flows(ci)
            bijm_cell{ci}{fi} = zeros(1, n_slots_of_coflow(ci));
            for si = 1: n_slots_of_coflow(ci)
                names{end+1} = {strcat('x(', num2str(ci), ',', num2str(fi), ',', num2str(si), ')')}; 
            end
        end
    end
    for i=1:length(names); names{i} = names{i}{1}; end
    
    model.varnames = names';
    % gurobi_write(model, '+ResourceAllocationAlgos/+Trung/ifip/gurobi_model.lp'); % save model
    
    % Solve
    results = gurobi(model, params);
    [x, fval] = deal(results.x, results.objval);
    if strcmp(results.status, 'OPTIMAL')
        exitflag = 1; 
    else
        exitflag = 0; 
    end
    outputs.output = results;

elseif strcmp(LP.solver, 'intlinprog')   
	
	% LP.options = [];
    % LP.options.MaxTime = 20;      % max compute time   
    LP.options = inputs.options;
    %[x,fval,exitflag,output] = intlinprog(obj,intcon,Aineq,bineq,Aeq,beq,lb,ub,options);
    
    if strcmp(zn_type, 'B') || strcmp(zn_type, 'I') || strcmp(zn_type, 'N') % CDS-LP
        LP.intcon = find(ctype == 'B' | ctype == 'I' | ctype == 'N'); % indexes of integer numbers
        [x, fval, exitflag, output] = intlinprog(LP);
        outputs.output = output;

    elseif strcmp(zn_type, 'C') || strcmp(zn_type, 'S')% CDS-LPA
        %[x, fval, exitflag, output] = linprog(LP);
        options = optimoptions('linprog');
        options.Display = 'off';
        [x, fval, exitflag, output] = linprog(obj,Aineq,bineq,Aeq,beq,lb,ub,options);
        outputs.output = output;
    end
end

t2 = cputime;
runtime = t2 - t1;
outputs.runtime = runtime;

clearvars C1 C2 C1_R C2_R lb ub ctype x0;

outputs.LP = LP;
outputs.exitflag = exitflag;
outputs.x = x;
outputs.fval = fval;

% flag indicating if at least one coflow is accepted
outputs.isAllocSuccess = false; 
    
if exitflag > 0 %isempty(x) == 0

    % Rearrange result
    zn = x(1:n_coflows)';
    bijm_vect = x(n_coflows+1 : end)';
    bijm = bijm_cell;

    k = 1;
    for ci = 1:n_coflows
        for fi = 1:n_flows(ci)
            for mi = 1:n_slots_of_coflow(ci)
                bijm{ci}{fi}(mi) = bijm_vect(k);
                k = k + 1;
            end
        end
    end
    
    %%
    bijm_synthesis = cell(1, n_links);
    
    for p = 1:n_links    
        bijm_synthesis{p} = zeros(n_flows_all, id_slot_mx);
        k = 1;
        for c = coflows
            ci = c.id;
            mn = c.addParam.slot_id(1);
            mx = c.addParam.slot_id(end);
            for f = c.flows
                fi = f.id;
                if ismember(p, f.links)
                    bijm_synthesis{p}(k, mn:mx) = bijm{ci}{fi}; 
                end
                k = k + 1;
            end
        end
    end
    
    
    %%
    % Calculate CCT
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
            cct_slot(ci) = max(cct_vect{ci});
            ccts(ci) = time_slots{cct_slot(ci)}(2);
        end
    end
    
    % cct_avg = mean(cct);
    % cct_min = min(cct);
    % cct_max = max(cct);

    % Save results to OptimOut
    if sum(zn) > 0 
        outputs.isAllocSuccess = true;
    end
    outputs.zn = (zn >= 0.99);
    outputs.accepts = find(zn >= 0.99); % id of accepted coflows
    outputs.rejects = find(zn < 0.1);   % id of rejected coflows
    outputs.bijm = bijm;
    outputs.bijm_synthesis = bijm_synthesis;
    outputs.nac = length(outputs.accepts); % for safety, instead of sum(zn)
    
    outputs.cct_vect = cct_vect;
    outputs.cct_slot = cct_slot;
    outputs.ccts = ccts(outputs.accepts); % only keep ccts of accepted coflows
    
    
    % OptimOut.cct_avg = cct_avg;
    % OptimOut.cct_min = cct_min;
    % OptimOut.cct_max = cct_max;
    
else
    disp('feasible solution not found');
    % OptimOut = {};
end


end