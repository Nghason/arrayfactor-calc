% (c) 2022 Quang-Trung Luu (LAAS-CNRS)
% Determine the scheduling of coflows with DCoflow algorithm
% Reference: Q-T. Luu et al., "DCoflow: Deadline-Aware Scheduling Algorithm for Coflows in Datacenter Networks," 2022.
% Inputs:
% - method: name of DCoflow variant: 'heu_v2_min_link', 'heu_v2_min_sum_negative', 'heu_v2_min_sum_congested'
% - fabric: the Fabric ojbect representing the network
% - coflows: an array of Coflow objects within the fabric
% Outputs:
% - order: scheduling order of coflows



function outputs = dcoflow(mod_name, fabric, coflows)


outputs.mod_name = mod_name;
global params;


%% Initializations:

n_links = fabric.numFabricPorts;    % nb of links
n_coflows = length(coflows);        % nb of coflows
deadlines = [coflows.deadline];     % deadline of coflows
order = zeros(1,n_coflows);         % scheduling order

% Capacity of each link
portCapacity = [[fabric.machinesPorts.ingress] [fabric.machinesPorts.egress]];
portCapacity = [portCapacity.linkCapacity];

% Unscheduled coflows IDs:
S = [coflows.id];
unsch_coflow_ids = S; % S = [n]

% Matrix of processing time of each coflow on each link:
D = zeros(n_links, n_coflows);
for c = coflows
    D(:, c.id) = c.indicator*[c.flows.volume]' ./ portCapacity';
end

% Debug
if params.debug_mode
    fid = fopen(params.debug_path, 'a+');
    fprintf(fid, "\nEXECUTION OF DCOFLOW ALGORITHM (variant: '%s')\n", mod_name);
    fprintf(fid, 'Processing time matrix D: [n_links * n_coflow]\n');
    %writematrix(D, fname, 'Delimiter', 'tab') % only available from R2019a
    dlmwrite(params.debug_path, D, 'delimiter', '\t', '-append', 'precision', '%.5f')
    % fclose(fid);
    fclose('all');
end


k = n_coflows;
order_star = []; % keep track of rejected coflows


%=================================================================================
% PART 1: COFLOW ORDERING 
%=================================================================================

while ~isempty(unsch_coflow_ids)
    
    % Debug
    if params.debug_mode
        fid = fopen(params.debug_path, 'a+'); 
        fprintf(fid, strcat('\n-----------------Unscheduled coflows = [',strjoin(string(unsch_coflow_ids),','),']-----------------\n'));
        fclose('all');
    end
    
    % Find the most bottlenecked links:
    cumulD = sum(D,2);
    b_candidates = find(cumulD == max(cumulD));
    
    % Randomly pick one such link:
    r_ind = randi(length(b_candidates));
    b = b_candidates(r_ind);   

    % Coflows that use bottleneck
    Sb = D(b,:) > 0; % coflows using b (S_b)
    Sb_id = S(Sb);
    
    if params.debug_mode
        fid = fopen(params.debug_path, 'a+'); 
        fprintf(fid, 'Bottleneck is %d with completion time %.4f\n', b, max(cumulD));
        fprintf(fid, strcat('Coflows using b: Sb = [',strjoin(string(Sb_id),','),']\n'));
        fclose('all');
    end
    
    
    % Coflows that have sum_{j in S_b} p_{jb} Vb/Bb <= DL_j
    % (any coflow j that can finish before deadline when scheduled last)
    C_d = Sb_id(deadlines(Sb) >= cumulD(b)/portCapacity(b));
    
    if ~isempty(C_d)
        
        % Take only the coflow with largest deadline in C_d
        j_star = C_d(deadlines(C_d) == max(deadlines(C_d)));    
        
        % Schedule j_star last
        order(k) = j_star; 

        % Remove j_star from S: S = S\{j_star}
        unsch_coflow_ids = setdiff(unsch_coflow_ids, order(k)); 
        
        % Set demand of last scheduled coflow to zero
        D(:, order(k)) = 0;     
        k = k - 1; % update k
        
        % Debug
        if params.debug_mode
            fid = fopen(params.debug_path, 'a+'); 
            fprintf(fid, 'Coflow %d can be scheduled as the last one on link %d\n', j_star, b);
            % fclose(fid);
            fclose('all');
        end
        
    else
        % x_star = argmin psi_{ell,x}
        psi = []; 
        
        %-----------------------------------------------------------------
        % SELECT DCFOLOW VERSION
        %-----------------------------------------------------------------
        % Considered coflows on the bottleneck
        switch mod_name
            
            %-----------------------------------------------------------------
            % STANDARD DCOFLOW & WEIGHTED DCOFLOW
            %-----------------------------------------------------------------
            case {'dcoflow_min_all','dcoflow_min_negative','dcoflow_min_congested',...
                  'w_dcoflow_min_all','w_dcoflow_min_negative','w_dcoflow_min_congested'} 
                
                considered_coflows = Sb_id;
                for cid = considered_coflows
                    c = coflows(cid);
                    D_ci = D(:, cid); % column of coflow cid in matrix D
                    psi_c_all = D_ci(c.addParam.indicator, :) ./ portCapacity(c.addParam.used_links)' .* ...
                            (deadlines(cid) * ones(length(c.addParam.used_links), 1) - ...
                            cumulD(c.addParam.used_links) ./ portCapacity(c.addParam.used_links)');  
                    
                    if strcmp(mod_name, 'dcoflow_min_all') || strcmp(mod_name, 'w_dcoflow_min_all')
                        psi_c_row = min(psi_c_all);     % take directly the minimum element of each row
                        
                    elseif strcmp(mod_name, 'dcoflow_min_negative') || strcmp(mod_name, 'w_dcoflow_min_negative')
                        psi_c_all(psi_c_all > 0) = 0;   % only summing negative elements
                        psi_c_row = sum(psi_c_all);
                        
                    elseif strcmp(mod_name, 'dcoflow_min_congested') || strcmp(mod_name, 'w_dcoflow_min_congested')
                        % Check if links used by coflow c are congested ~bottleneck
                        is_congested = cumulD(c.addParam.indicator, :) >= 0.7*cumulD(b);
                        psi_c_row = psi_c_all(is_congested, :);
                        psi_c_row = sum(psi_c_row);
                    end
                    
                    psi = [psi; psi_c_row];
                    
                end % end for cid
            
            %-----------------------------------------------------------------
            % DCOFLOW WITH DYNAMIC PROGRAMMING (min_negative only)
            %-----------------------------------------------------------------
            case {'dcoflow_dp'}
                % Run DP with the coflow set S_b on bottleneck link b
                [accepts, rejects] = algorithms.offline.DP(fabric, coflows(Sb), b);
                
                considered_coflows = rejects;
                for cid = considered_coflows
                    c = coflows(cid);
                    D_ci = D(:, cid); % column of coflow cid in matrix D
                    
                    % Different to DCoflow standard: Psi is divided by coflow weights
                    psi_c_all = D_ci(c.addParam.indicator, :) ./ portCapacity(c.addParam.used_links)' .* ...
                                (deadlines(cid) * ones(length(c.addParam.used_links), 1) - ...
                                cumulD(c.addParam.used_links) ./ portCapacity(c.addParam.used_links)')/c.weight;
                            
                    psi_c_all(psi_c_all > 0) = 0;   % only summing negative elements
                    psi_c_row = sum(psi_c_all);
                    psi = [psi; psi_c_row];
                   
                end % end for cid
                
            %-----------------------------------------------------------------
            % DCOFLOW WITH MOORE-HODGSON
            %-----------------------------------------------------------------
            case {'dcoflow_mh'}
                
                % Run Moore-Hodgson algorithm on bottleneck link b
                [S_ell, E_ell] = algorithms.offline.moore_hodgson(fabric, coflows(Sb), b);
                
                considered_coflows = E_ell;
                for cid = considered_coflows
                    c = coflows(cid);
                    D_ci = D(:, cid); % column of coflow cid in matrix D
                    
                    % psi similar to standard Dcoflow Different to DCoflow standard: Psi is divided by coflow weights
                    psi_c_all = D_ci(c.addParam.indicator, :) ./ portCapacity(c.addParam.used_links)' .* ...
                                (deadlines(cid) * ones(length(c.addParam.used_links), 1) - ...
                                cumulD(c.addParam.used_links) ./ portCapacity(c.addParam.used_links)'); %/c.weight;
                            
                    psi_c_all(psi_c_all > 0) = 0;   % only summing negative elements
                    psi_c_row = sum(psi_c_all);
                    psi = [psi; psi_c_row];
                   
                end % end for cid
              
                
            otherwise
                warning('No dcoflow option specified!');   
                
        end % end switch methods            
        
        % Select coflow x_star to reject
        x_star = considered_coflows((psi == min(psi)));
        x_star = x_star(1); % take the first in x_star (update 07/10/2021)
        
        %-----------------------------------------------------------------
        % Debug
        if params.debug_mode
            fid = fopen(params.debug_path, 'a+');
            fprintf(fid, '\nRejectCoflow: considering coflow %d with deadline %.4f\n', c.id, c.deadline);
            fprintf(fid, "Calculate psi_(l,k') = p_{l,k'}*(T_k'-sum_{k in S_b} p_{l,k})\n");
            for ii = 1:length(c.addParam.used_links)
                lnk = c.addParam.used_links(ii);
                fprintf(fid, "\ton link %d: psi = %.4f*(%.4f-%.4f) = %.4f\n", ...
                    lnk, D_ci(lnk), c.deadline, cumulD(lnk)/portCapacity(lnk), psi_c_all(ii));
            end
            fprintf(fid, '\nDecision: reject coflow %d among those using bottleneck %d\n', x_star, b);
            % fclose(fid);
            fclose('all');
        end
        %-----------------------------------------------------------------
        order_star = [x_star(1), order_star]; 
        order(k) = x_star;

        % Remove x_star from S: S = S\{x_star}
        unsch_coflow_ids = setdiff(unsch_coflow_ids, x_star);
        
        % Set demand of last scheduled (removed) coflow to zero
        D(:, x_star) = 0; 
        k = k - 1; % update k   
        
    end % end if ~isempty(C_d)
    
end % end while 

outputs.ini_order = order;

% Debug
if params.debug_mode
    fid = fopen(params.debug_path, 'a+'); 
    fprintf(fid, strcat('sigma = [',strjoin(string(order),','),']\n'));
    fprintf(fid, strcat('sigma* = [',strjoin(string(order_star),','),']\n'));
    % fclose(fid);
    fclose('all');
end

%=================================================================================
% PART 2: REMOVE LATE COFLOWS FROM THE ORDER 
%=================================================================================

% Initial guess, with only CCTs, without backtracking
ini_ccts = utils.evalCCTfromOrder(fabric, coflows, order); 
ini_zn = (ini_ccts <= deadlines(order));
ini_nac = sum(ini_zn);           
ini_order_star = order_star;
pred_rejects = [];  % final sigma^star



% Backtrack to remove rejected coflows from the order
while ~isempty(order_star)
    kstar = order_star(1);              % first coflow in sigma^star
    kstar_id = find(order == kstar);    % id of k^star in sigma

    % Calculate the CCT of sigma(k^star)
    temp = utils.evalCCTfromOrder(fabric, coflows, order(1:kstar_id));
    kstar_cct = temp(kstar_id); % note: the order in temp is 1:n_coflows
    kstar_deadline = coflows(kstar).deadline;

    %fprintf("c%d: cct = %.4f\n", kstar, kstar_cct);
    
    % Debug
    if params.debug_mode
        fid = fopen(params.debug_path, 'a+');
        fprintf(fid, '\nRemoveLateCoflow: consider coflow %d: cct = %.4f, deadline = %.4f\n', kstar, kstar_cct, kstar_deadline);
        if kstar_cct > kstar_deadline
            fprintf(fid, 'Coflow %d is definitely rejected\n', kstar);
        end
        % fclose(fid);
        fclose('all');
    end

    if kstar_cct > kstar_deadline
        % Remove k^star from sigma
        order(kstar_id) = []; 
        pred_rejects = [pred_rejects, kstar];
    end
    % Remove k^star from sigma^star
    order_star(1) = []; 
end   


% Save
outputs.ini_ccts = ini_ccts; 
outputs.deadlines = deadlines(order);
outputs.ini_zn = ini_zn;
outputs.ini_nac = ini_nac;           
outputs.ini_order_star = ini_order_star; % sigma^star: coflows that do not meet their deadlines
outputs.order_star = order_star;
outputs.pred_rejects = pred_rejects;

% Final predicted results
outputs.pred_zn      = ~ismember(outputs.ini_order, outputs.pred_rejects);
outputs.pred_accepts = outputs.ini_order(outputs.pred_zn);
outputs.pred_ccts    = outputs.ini_ccts(outputs.pred_zn); % predicted ccts (of only accepted coflows)
outputs.pred_nac     = length(order);
outputs.order 		 = order; % final scheduling order


end