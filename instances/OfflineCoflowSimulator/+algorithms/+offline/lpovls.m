% (c)   2023 Quang-Trung Luu
%       School of Electrical and Electronic Engineering
% Implementation of LP-OV-LS algorithm
% Reference:    M. Shafiee and J. Ghaderi, "An Improved Bound for Minimizing
%               the Total Weighted Completion Time of Coflows in Datacenters," 
%               in IEEE/ACM Transactions on Networking, vol. 26, no. 4, 
%               pp. 1674-1687, Aug. 2018
% Inputs:
% - fabric: the Fabric ojbect representing the network
% - coflows: an array of Coflow objects within the fabric
% Outputs:
% - resource allocation scheme for the given coflows


function order = lpovls(fabric, coflows)

global params;

%% Initialization
n_links    = fabric.numFabricPorts;     % nb of fabric ports (ingress+egress)    
n_machines = n_links/2;                 % nb of machines on the fabric
n_coflows  = length(coflows);           % nb of coflows
n_flows    = [coflows.numFlows];        % nb of flows of each coflow
n_flows_all = sum(n_flows);             % total nb of flows

portCapacity = [[fabric.machinesPorts.ingress] [fabric.machinesPorts.egress]];
portCapacity = [portCapacity.linkCapacity];

coflows_ids = [coflows.id];

% Aggregate volume of coflows at source i and destination j: d^{k}_{i} and d^{k}_{j}
for c = coflows
    cid = c.id;
    agg_src = zeros(1, n_machines);
    agg_dst = zeros(1, n_machines);
    
    tmp1 = [c.flows.source];
    tmp2 = [c.flows.destination];
    
    f_src = [tmp1.id];
    f_dst = [tmp2.id];
    
    for i = 1:n_machines
        j = i + n_machines;         % destination port number
        for f = c.flows
            if f.links(1) == i      % flow f uses source i
                agg_src(i) = agg_src(i) + f.volume;
            end
            if f.links(2) == j      % f uses destination j
                agg_dst(i) = agg_dst(i) + f.volume;
            end
        end
    end
    c.addParam.f_src = f_src;
    c.addParam.f_dst = f_dst;
    c.addParam.agg_src = agg_src;
    c.addParam.agg_dst = agg_dst;
    
    % Effective volume of coflows: W(k)= max{max d^k_i , max d^k_j}
    % c.addParam.eff_vol = max(max(agg_src), max(agg_dst)); % exactly the CCT0!
    
end
clearvars tmp1 tmp2 agg_src agg_dst;

% Nb of variables
nvar_fk = n_coflows;
% nvar_delta = nchoosek(n_coflows, 2); % nb of combinations (C^K_2)
nvar_delta = n_coflows^2;



%% 1. Solve the linear program
%{
min     sum_{k=1}^{K} w_k * f_k
s.t.    (c1) f_k >= d^{k}_{i} + sum_{k' in K} d^{k'}_{i} * delta_{k'k},
                    forall i,k
        (c2) f_k >= d^{k}_{j} + sum_{k' in K} d^{k'}_{j} * delta_{k'k},
                    forall j,k
        (c3) f_k >= W(k) + r_k,                 forall k
        (c4) delta_{k'k} + delta_{kk'} = 1,     forall k,k'
        (c5) f_k in R, delta_{kk'} in {0,1},              forall k,k'
%}


%% Objective function (column vector)
% minimize sum_{k=1}^{K} w_k * f_k
obj = [[coflows.weight]'; zeros(nvar_delta, 1)];


%% C1: CCT bound at ingress port i
% -f_k + sum_{k' in K} d^{k'}_{i} * delta_{k'k} <= -d^{k}_{i}, forall i,k

C1 = [];
C1_R = [];
for ci = 1:n_coflows
    not_ci = setdiff(coflows_ids, ci); % id of other coflows
    c = coflows(ci);
    
    % Fill fk
    fk = zeros(1, n_coflows);
    fk(ci) = -1;
    
    % Fill delta
    for i = 1:n_machines
        delta = zeros(n_coflows, n_coflows);
        for nci = not_ci
            delta(nci, ci) = coflows(nci).addParam.agg_src(i); % d^{k'}_{i}
        end
        rhs = -c.addParam.agg_src(i); % -d^{k}_{i}
        
        % Append fk and delta to C1
        C1 = [C1; fk, reshape(delta, [1, nvar_delta])];
        C1_R = [C1_R; rhs];
    end
end


%% C2: CCT bound at egress port j
% -f_k + sum_{k' in K} d^{k'}_{j} * delta_{k'k} <= -d^{k}_{j}, forall j,k

C2 = [];
C2_R = [];
for ci = 1:n_coflows
    not_ci = setdiff(coflows_ids, ci); % id of other coflows
    c = coflows(ci);
    
    % Fill fk
    fk = zeros(1, n_coflows);
    fk(ci) = -1;
    
    % Fill delta
    for j = 1:n_machines
        delta = zeros(n_coflows, n_coflows);
        for nci = not_ci
            delta(nci, ci) = coflows(nci).addParam.agg_dst(j); % d^{k'}_{j}
        end
        rhs = -c.addParam.agg_dst(j); % -d^{k}_{j}
        
        % Append fk and delta to C2
        C2 = [C2; fk, reshape(delta, [1, nvar_delta])];
        C2_R = [C2_R; rhs];
    end
end


%% C3: CCT bound wrt effective volume
% -f_k <= W(k) + r_k = CCT0(k) + k.arrival, forall k

C3 = [];
C3_R = [];
for c = coflows
    cid = c.id;
    
    % Fill fk
    fk = zeros(1, n_coflows);
    fk(cid) = -1;
    
    rhs = c.addParam.CCT0 + c.arrival;
    
    % Append fk and delta to C3
    C3 = [C3; fk, zeros(1, nvar_delta)];
    C3_R = [C3_R; rhs];
end


%% C4: For each coflow pari, one coflow precedes the other
% delta_{k'k} + delta_{kk'} = 1, forall k,k'

C4 = [];
C4_R = [];
for ci = 1:n_coflows
    not_ci = setdiff(coflows_ids, ci); % id of other coflows
    
    % Fill delta
    delta = zeros(n_coflows, n_coflows);
    for nci = not_ci
        delta(ci, nci) = 1;
        delta(nci, ci) = 1; 
    end
    C4 = [C4; zeros(1, n_coflows), reshape(delta, [1, nvar_delta])];
    C4_R = [C4_R; 1];
end


%% Optimizer setup

% Types of variables
ctype = strcat( char(ones([1, nvar_fk])*('S')), ...  % type of f_k
                char(ones([1, nvar_delta])*('B')));  % type of delta
            
% Bound: f_k in [0, inf], delta in {0,1} 
lb = [zeros(1, nvar_fk), zeros(1, nvar_delta)];
ub = [inf*ones(1, nvar_fk), ones(1, nvar_delta)];
     
Aineq = [C1; C2; C3];         
bineq = [C1_R; C2_R; C3_R]; 
Aeq = C4;
beq = C4_R;

% Save parameters of optimizer
LP.solver   = 'intlinprog';
LP.f        = obj;
LP.Aineq    = Aineq;         
LP.bineq    = bineq; 
LP.Aeq      = Aeq;
LP.beq      = beq;
LP.C1       = C1;
LP.C2       = C2;
LP.C1_R     = C1_R;
LP.C2_R     = C2_R;
LP.lb       = lb;
LP.ub       = ub;
LP.ctype    = ctype;
LP.x0       = [];

%======================================
% Solving
%======================================

t1 = cputime;

if strcmp(LP.solver, 'cplexmilp')  
    %[x,fval,exitflag,output] = cplexmilp(f,Aineq,bineq,Aeq,beq,[],[],[],lb,ub,ctype,x0,options);
    % options = cplexoptimset('cplex');
    % LP.timelimit = opt_options.MaxTime;  
    % LP.mip.tolerances.mipgap = 5e-3; % 0-1 default 1e-4
    % LP.mip.limits.nodes = opt_options.MaxNodes;   
    % LP.options = options;
    if params.export_lp
        options.exportmodel = 'logs/cplex_model.lp';
    end
    
    LP.options = []; %opt_options.options; 
    [x, fval, exitflag, output] = cplexmilp(LP);
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
    
    % Set variable names
    names = {};
    for ci = 1:n_coflows
        names{end+1} = {strcat('f', num2str(ci))}; 
    end
    
    for i = 1:n_coflows
       for j = 1:n_coflows
           names{end+1} = {strcat('delta(', num2str(i), ',', num2str(j), ')')}; 
       end
    end
    
    for i=1:length(names)
        names{i} = names{i}{1}; 
    end
    
    model.varnames = names';
    
    if params.export_lp
        gurobi_write(model, 'logs/gurobi_model.lp'); % save model
    end
    
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
    LP.options = optimoptions('intlinprog');
    LP.options.Display = 'off';
    % LP.options.MaxTime = 10;      % max compute time         
    % LP.options.MaxNodes = 1e6;  % max compute nodes
    
    %[x,fval,exitflag,output] = intlinprog(obj,intcon,Aineq,bineq,Aeq,beq,lb,ub,options);
    LP.intcon = find(ctype == 'B' | ctype == 'I' | ctype == 'N'); % indexes of integer numbers
        [x, fval, exitflag, output] = intlinprog(LP);
        outputs.output = output;
        
end

t2 = cputime;
runtime = t2 - t1;
outputs.runtime = runtime;

% clearvars C1 C2 C1_R C2_R lb ub ctype x0;

outputs.LP = LP;
outputs.exitflag = exitflag;
outputs.x = x;
outputs.fval = fval;

% flag indicating if at least one coflow is accepted
outputs.isAllocSuccess = false; 

if exitflag > 0
    % Rearrange result
    fk = x(1:nvar_fk)';
    delta = x(nvar_fk+1 : end)';
    delta = reshape(delta, [n_coflows, n_coflows]);
    outputs.fk = fk;
    outputs.delta = delta;
    
    if params.debug_mode
        fprintf('Results of solving the MILP:\n');
        fprintf(strcat('\t fk = [',  strjoin(string(fk), ', '), ']\n'));
        fprintf('\t delta = \n');
        disp(delta);
        
    end
    
    % Get coflow scheduling order
    [~, order] = sort(fk);
    
    
    
    
%     Acceptance rate synthesis
%     outputs.zn = (zn >= 0.99);
%     outputs.accepts = find(zn >= 0.99); % id of accepted coflows
%     outputs.rejects = find(zn < 0.1);   % id of rejected coflows
%     outputs.bijm = bijm;
%     outputs.bijm_synthesis = bijm_synthesis;
%     outputs.nac = length(outputs.accepts); % for safety, instead of sum(zn)
    
%     CCT synthesis
%     outputs.cct_vect = cct_vect;
%     outputs.cct_slot = cct_slot;
%     outputs.ccts = ccts(outputs.accepts);  
    
    
else
    disp('feasible solution not found');
end










end % END FUNCTION