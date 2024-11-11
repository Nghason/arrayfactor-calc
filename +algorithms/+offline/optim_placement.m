% (c) 2022 Quang-Trung Luu (LAAS-CNRS)
% Linear program to solve the optimal reducer placement problem
% Inputs:
% - fabric: the Fabric ojbect representing the network
% - coflows: an array of Coflow objects within the fabric
% Outputs:
% - order: scheduling order of coflows




function outputs = placement_optim(inputs, fabric, coflows)


%% Initialization
n_links    = fabric.numFabricPorts;     % nb of fabric ports (ingress+egress)    
n_machines = n_links/2;                 % nb of machines on the fabric
n_coflows  = length(coflows);           % nb of coflows
n_flows    = [coflows.numFlows];        % nb of flows of each coflow
n_flows_all = sum(n_flows);             % total nb of flows

portCapacity = [[fabric.machinesPorts.ingress] [fabric.machinesPorts.egress]];
portCapacity = [portCapacity.linkCapacity];

ingresses = (1:n_machines);
egresses  = (n_machines+1 : n_links);

% Nb of variables [zk-delta-phi-y]
nvar_zk    = n_coflows;                              % z_{k}
nvar_delta = n_coflows^2;                            % delta_{k,k'}: redundant elements when k = k' (diagonal values)
nvar_phi   = n_flows_all * n_links;                  % phi^{k,j}_{ell}
nvar_y     = n_coflows^2 * n_flows_all * n_links;    % y_{k,k',j,ell}: redundant elements when k = k'
nvar_y_ci  = n_coflows * n_flows_all * n_links;      % nb of variables in each y_{k}{k'}
nvar_all   = nvar_zk + nvar_delta + nvar_phi + nvar_y;   % total nb of variables

% An empty cell for all ci-fi-ell to be used widely
phi_cell = cell(1, n_coflows);
for ci = 1:n_coflows
    phi_cell{ci} = zeros(n_flows(ci), n_links);
end

% Empty cell for y_{k,k',j,ell}
y_cell = cell(1, n_coflows);
for c1 = 1:n_coflows
    y_cell{c1} = cell(1, n_coflows-1); % all k' != k in C
    for c2 = 1:n_coflows
        y_cell{c1}{c2} = zeros(n_flows(c2), n_links); % each k' has matrix of F_k' * L
    end
end

% Trick to transform phi and y to vector
% phi_vec = [];   for ii = 1:n_coflows;   phi_vec = [phi_vec, [phi{ii}]];         end
% y_c1_vec = [];  for ii = 1:n_coflows;   y_c1_vec = [y_c1_vec, [y_c1{ii}{:}]];   end


%====================================
%                                   |
% CONSTRUCTION OF THE LP            |
%                                   |
%==================================== 


%% Objective function (column vector)

LP = struct; % we gonna save all setup into LP struct

% Append zk to f: max(sum of zk) -> min(-sum of zk)
obj = [-ones(nvar_zk, 1); zeros(nvar_all - nvar_zk, 1)];


%% C1: Ordering of two coflows
% delta_{k,k'} + delta_{k',k} = 1 for all k, k'
% nb of rows: combination C_2^{n_coflows}

C1 = [];
combs = combnk(1:n_coflows, 2); % combination of 2 in n_coflows
for i = 1:size(combs,1)
    cmb = combs(i,:);
    delta = zeros(n_coflows, n_coflows);    
    delta(cmb(1), cmb(2)) = 1;
    delta(cmb(2), cmb(1)) = 1;
    delta = reshape(delta, [1, nvar_delta]);
    
    % Append to C1 [zk-delta-phi-y]
    C1 = [C1; zeros(1, nvar_zk), delta, zeros(1, nvar_phi + nvar_y)];
end

% Right-hand of C1
C1_R = ones(size(C1,1), 1);


%% C2: Ordering of three coflows
% delta_{k,k'} + delta_{k',k"} + delta_{k",k} <= 2 for all k, k', k"

C2 = [];
combs = combnk(1:n_coflows, 3); % combination of 3 in n_coflows
for i = 1:size(combs,1)
    cmb = combs(i,:);
    permus = perms(cmb); % permutation of each triplet (k,k',k")
    
    for j = 1:size(permus,1)
        pm = permus(j,:);
        delta = zeros(n_coflows, n_coflows);    
        delta(pm(1), pm(2)) = 1;
        delta(pm(2), pm(3)) = 1;
        delta(pm(3), pm(1)) = 1;
        delta = reshape(delta, [1, nvar_delta]);

        % Append to C2 [zk-delta-phi-y]
        C2 = [C2; zeros(1, nvar_zk), delta, zeros(1, nvar_phi + nvar_y)];
    end
end
% Remove duplicated rows in C2
C2 = unique(C2, 'row', 'stable');

% Right-hand of C2
C2_R = ones(size(C2,1), 1);


%% C3: Completion time c_{ell, k} 
% 1/B_ell * sum_{j in F_k} [v_{k,j}*delta_{k,j,ell}] + 
% 1/B_ell * sum_{k'! = k} sum_{j in F_k'} y_{k,k',j,ell} <= T_k forall ell, k
% nb of rows: n_coflows * n_flows_all

C3 = [];
C3_R = [];

for c1 = 1:n_coflows
    for ell = 1:n_links
        
        % Each {c1, ell} has a line of constraint
        phi = phi_cell;
        y_c1 = y_cell; % y{c1}
        
        % Proc time of k: p_{ell,k} = 1/B_ell * sum_{j in F_k} [v_{k,j}*delta_{k,j,ell}]
        for f = coflows(c1).flows
            f1 = f.id;
            phi{c1}(f1, ell) = coflows(c1).addParam.volumes(f1) / portCapacity(ell); 
        end
        
        % Proc time of all k' != k: 1/B_ell * sum_{k' != k} sum_{j in F_k'} y_{k,k',j,ell}
        for c2 = 1:n_coflows
            if c2 ~= c1
                for f = coflows(c2).flows
                    f2 = f.id;
                    y_c1{c1}{c2}(f2, ell) = coflows(c2).addParam.volumes(f2) / portCapacity(ell); 
                end
            end
        end
        
        % Transform phi and y to vector
        phi_vec = [];   for ii = 1:n_coflows;   phi_vec = [phi_vec, [phi{ii}]];         end
        y_c1_vec = [];  for ii = 1:n_coflows;   y_c1_vec = [y_c1_vec, [y_c1{ii}{:}]];   end
        
        % Append to C3 [zk-delta-phi-y]
        before = (c1-1)*nvar_y_ci;              
        after  = nvar_y - before - nvar_y_ci;   
        C3 = [  C3; ...
                zeros(1, nvar_zk), ...
                zeros(1, nvar_delta), ...
                phi_vec, ...        
                zeros(1, before), ...   % all y{ci} before y{c1}
                y_c1_vec, ...       
                zeros(1, after)...      % all y{ci} after y{c1}
        ];
    
        % Right-hand of C3
        C3_R = [C3_R; coflows(c1).deadline];
    
    end % end loop of ell
end % end loop of c1


%% C4: Definition of phi{k}(j, ell}
% phi{k}(j, ell) = z_{k} if ell = Input(j) and 0 if ell != Input(j)

C4 = [];
for ci = 1:n_coflows
    for ell = 1:n_links
        for f = coflows(ci).flows
            fi = f.id;
            
            % Each {ci, fi, ell} has a line of constraint
            phi = phi_cell;
            zk = zeros(1, n_coflows);
            
            if ell == f.links(1) 
                % if ell is ingress of f then z_{k} - phi{k}(j, ell) = 0 
                zk(ci) = 1;
                phi{ci}(fi, ell) = -1; 
            else
                % if ell is not ingress of f then phi{k}(j, ell) = 0
                phi{ci}(fi, ell) = 1; 
            end
            % Transform phi to vector
            phi_vec = []; for ii = 1:n_coflows; phi_vec = [phi_vec, [phi{ii}]]; end

            % Append to C4 [zk-delta-phi-y]
            C4 = [C4; zk, zeros(1, nvar_delta), phi_vec, zeros(1, nvar_y)];   
            
        end % end loop of f
    end % end loop of ell
end % end loop of c1

% Right-hand of C4
C4_R = zeros(size(C4,1), 1);

%% C5: Mapped only once of phi{k}(j, ell} with ell is egress
% sum_{ell in L_out} phi{k}(j, ell) = 1 forall j in F_k and k in C

C5 = [];
for ci = 1:n_coflows
    for ell = n_machines+1 : n_links % ell in L_out (egress)
        for f = coflows(ci).flows
            fi = f.id;
            
            % Each {ci, fi, ell} has a line of constraint
            phi = phi_cell;
            phi{ci}(fi, ell) = 1; 
            
            % Transform phi to vector
            phi_vec = []; for ii = 1:n_coflows; phi_vec = [phi_vec, [phi{ii}]]; end
            
            % Append to C5 [zk-delta-phi-y]
            C5 = [C5; zeros(1, nvar_zk), zeros(1, nvar_delta), phi_vec, zeros(1, nvar_y)]; 

        end % end loop of f
    end % end loop of ell
end % end loop of c1

% Right-hand of C5
C5_R = ones(size(C5,1), 1);


%% C6: Supplemental constraints for y_{k,k',j,ell}
% y_{k,k',j,ell} <= phi^{k',j}_{ell}                     => - phi^{k',j}_{ell} + y_{k,k',j,ell}  <= 0  
% y_{k,k',j,ell} <= delta_{k,k'}                         => - delta_{k,k'} + y_{k,k',j,ell} <= 0         
% y_{k,k',j,ell} >= phi^{k',j}_{ell} + delta_{k,k'} - 1  => delta_{k,k'} + phi^{k',j}_{ell} - y_{k,k',j,ell} <= 1   

% HERE !!!!
C6 = [];
C6_R = [];

for c1 = 1:n_coflows
    for c2 = 1:n_coflows
        if c2 ~= c1
            for ell = 1:n_links
        
                % Each {c1, c2, f2, ell} has a line of constraint
                delta = zeros(n_coflows, n_coflows);
                phi = phi_cell;
                y_c1 = y_cell; % y{c1}

                for f = coflows(c2).flows
                    f2 = f.id;
                    y_c1{c1}{c2}(f2, ell) = coflows(c2).addParam.volumes(f2) / portCapacity(ell); 
                    
                    % Transform phi and y to vector
                    phi_vec = [];   for ii = 1:n_coflows;   phi_vec = [phi_vec, [phi{ii}]];         end
                    y_c1_vec = [];  for ii = 1:n_coflows;   y_c1_vec = [y_c1_vec, [y_c1{ii}{:}]];   end

                    % Append to C6 [zk-delta-phi-y]
                    before = (c1-1)*nvar_y_ci;              
                    after  = nvar_y - before - nvar_y_ci;   
                    C6 = [  C6; ...
                            zeros(1, nvar_zk), ...
                            zeros(1, nvar_delta), ...
                            phi_vec, ...        
                            zeros(1, before), ...   % all y{ci} before y{c1}
                            y_c1_vec, ...       
                            zeros(1, after)...      % all y{ci} after y{c1}
                    ];

                    % Right-hand of C3
                    C6_R = [C6_R; coflows(c1).deadline];
        
                end % end loop of f2
            end % end loop of ell 
        end % end if c2 != c1
    end % end loop of c2
end % end loop of c1





%% Optimizer setup

% Bound: zn in {0,1}, bijm in [0, inf]
lb = zeros(1, nvar_all);
ub = [ones(1, n_coflows), inf*ones(1, nvar_bijm)];
ctype = strcat( char(ones([1, n_coflows])*(zn_type)), ... % type of zn
                char(ones([1, nvar_bijm])*('S')));        % type of bijm

[Aineq, bineq, Aeq, beq] = deal(C2, C2_R, C1, C1_R);       

LP.Aineq = Aineq;         
LP.bineq = bineq; 
LP.Aeq = Aeq;
LP.beq = beq;

% LP.Aineq = [C1; C2];         
% LP.bineq = [C1_R; C2_R];     
% LP.Aeq = [];
% LP.beq = [];

LP.x0 = [];

% Save parameters of optimizer
LP.solver   = inputs.solver;
LP.f        = obj;
LP.C1       = C1;
LP.C2       = C2;
LP.C1_R     = C1_R;
LP.C2_R     = C2_R;
LP.lb       = lb;
LP.ub       = ub;
LP.ctype    = ctype;


t1 = cputime;

if strcmp(LP.solver, 'cplexmilp')  
    LP.options = []; %opt_options.options; 
    %options = cplexoptimset('cplex');
    %options.timelimit = opt_options.MaxTime;  
    %options.mip.tolerances.mipgap = 5e-3; % 0-1 default 1e-4
    %options.mip.limits.nodes = opt_options.MaxNodes;
    %options.exportmodel = 'dyn.lp';
	%[x,fval,exitflag,output] = cplexmilp(f,Aineq,bineq,Aeq,beq,[],[],[],lb,ub,ctype,x0,options);
    if strcmp(zn_type, 'N') % CDS-LP
        [x, fval, exitflag, output] = cplexmilp(LP);
    elseif strcmp(zn_type, 'S') % CDS-LPA
        [x, fval, exitflag, output] = cplexlp(LP);
    end
    
elseif strcmp(LP.solver, 'intlinprog')   
	
	% LP.options = [];
    % LP.options.MaxTime = 20;      % max compute time   
    LP.options = inputs.options;
    %[x,fval,exitflag,output] = intlinprog(obj,intcon,Aineq,bineq,Aeq,beq,lb,ub,options);
    
    if strcmp(zn_type, 'N') % CDS-LP
        LP.intcon = find(ctype == 'B' | ctype == 'I' | ctype == 'N'); % indexes of integer numbers
        [x, fval, exitflag, output] = intlinprog(LP);
    elseif strcmp(zn_type, 'S') % CDS-LPA
        %[x, fval, exitflag, output] = linprog(LP);
        options = optimoptions('linprog');
        options.Display = 'off';
        [x, fval, exitflag, output] = linprog(obj,Aineq,bineq,Aeq,beq,lb,ub,options);
    end
   
end

t2 = cputime;
runtime = t2 - t1;
outputs.runtime = runtime;

clearvars C1 C2 C1_R C2_R lb ub ctype x0;

outputs.LP = LP;
outputs.exitflag = exitflag;
outputs.output = output;
outputs.x = x;
outputs.fval = fval;

% flag indicating if at least one coflow is accepted
outputs.isAllocSuccess = false; 
    
if exitflag > 0 %isempty(x) == 0

    % Rearrange result
    zn = x(1:n_coflows)';
    bijm_vect = x(n_coflows+1 : end)';
    bijm = phi_cell;

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
end


end