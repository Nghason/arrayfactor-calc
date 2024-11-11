% (c) 2022 Quang-Trung Luu (LAAS-CNRS)
% Modification of DCoflow algorithm by adding coflow weights (WDCoflow)
% Inputs:
% - method: name of WDCoflow variant: 'wDcoflow_min_link', 'wDcoflow_min_sum_negative', 'wDcoflow_min_sum_congested'
% - fabric: the Fabric ojbect representing the network
% - coflows: an array of Coflow objects within the fabric
% Outputs:
% - order: scheduling order of coflows



function outputs = weighted_dcoflow(method, fabric, coflows)


outputs.method = method;


%% Initializations:
% Number of links
n_links = fabric.numFabricPorts;

% Capacity of each link
portCapacity = [[fabric.machinesPorts.ingress] [fabric.machinesPorts.egress]];
portCapacity = [portCapacity.linkCapacity];

% Number of coflows
n_coflows = length(coflows);

% Deadline of coflows
deadlines = [coflows.deadline];

% Unscheduled coflows IDs:
S = [coflows.id];
unsch_coflow_ids = S; % S = [n]

% Demand of each coflow on each link:
D_ini = zeros(n_links, n_coflows);

% Scheduling order
order = zeros(1,n_coflows);

% Compute D:
for c = coflows
    D_ini(:, c.id) = c.indicator*[c.flows.volume]';
end
D = D_ini; % keep track of D_ini

outputs.D_ini = D_ini;
outputs.bottleneck = zeros(1, n_coflows);

k = n_coflows;
order_star = []; % keep track of rejected coflows
last_ok = []; % coflows in the IF of the algorithm


%=================================================================================
% PART 1: COFLOW ORDERING 
%=================================================================================

while ~isempty(unsch_coflow_ids)
    
    % Find the most bottlenecked links:
    cumulD = sum(D,2);
    b_candidates = find(cumulD == max(cumulD));
    
    % Randomly pick one such link:
    r_ind = randi(length(b_candidates));
    b = b_candidates(r_ind);   
 
    
    %% TEST: to match implementation of example in sincronia paper:
    % b = max(b_canditates(b_canditates<=4));
    % end TEST
    %%
    
    % Coflows that use bottleneck
    Sb = D(b,:) > 0; % coflows using b (S_b)
    Sb_id = S(Sb);
    
    % Coflows that have sum_{j in S_b} p_{jb} Vb/Bb <= DL_j
    % (any coflow j that can finish before deadline when scheduled last)
    C_d = Sb_id(deadlines(Sb) >= cumulD(b)/portCapacity(b));
    
    
    % fprintf("k = %d, b = %d with vol = %.4f\n", k, b, max(cumulD));
    % fprintf(['\t C_d = [' repmat('%d,', 1, numel(C_d)) ']\n'], C_d);
    
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
        
    else
        % x_star = argmin psi_{ell,x}
        psi = []; 
        
        for cid = Sb_id
            c = coflows(cid);
            D_ci = D(:, cid); % column of coflow cid in matrix D
            % psi = p_{l,k'}*(T_{k'} - sum_{k in S} p_{l,k})
            psi_c_all = D_ci(c.addParam.indicator, :) ./ portCapacity(c.addParam.used_links)' .* ...
                        (deadlines(cid) * ones(length(c.addParam.used_links), 1) - ...
                        cumulD(c.addParam.used_links) ./ portCapacity(c.addParam.used_links)')/c.weight;
            
            switch method
                case 'wDcoflow_min_link'
                    psi_c_row = min(psi_c_all);     % take directly the minimum element of each row
                case 'wDcoflow_min_sum_negative'
                    psi_c_all(psi_c_all > 0) = 0;   % only summing negative elements
                    psi_c_row = sum(psi_c_all);
                case 'wDcoflow_min_sum_congested'
                    % Check if links used by coflow c are congested ~bottleneck
                    is_congested = cumulD(c.addParam.indicator, :) >= 0.7*cumulD(b);
                    psi_c_row = psi_c_all(is_congested, :);
                    psi_c_row = sum(psi_c_row);
                otherwise
                    warning('No WDcoflow option specified!');   
            end
            
      
            psi = [psi; psi_c_row];
        end
        
        % Select coflow x_star to reject
        x_star = Sb_id((psi == min(psi)));
        x_star = x_star(1); % take the first in x_star (update 07/10/2021)
        order_star = [order_star, x_star(1)]; 
        order(k) = x_star;

        % Remove x_star from S: S = S\{x_star}
        unsch_coflow_ids = setdiff(unsch_coflow_ids, x_star);
        
        % Set demand of last scheduled (removed) coflow to zero
        D(:, x_star) = 0; 
        k = k - 1; % update k   
    end
    
end % end while 

outputs.ini_order = order;

%=================================================================================
% PART 2: REMOVE LATE COFLOWS FROM THE ORDER 
%=================================================================================

% Initial guess, with only CCTs, without backtracking
ini_ccts = utils.evalCCTfromOrder(fabric, coflows, order); 
ini_zn = (ini_ccts <= deadlines(order));
ini_nac = sum(ini_zn);           
ini_order_star = order_star; % sigma^star: coflows that do not meet their deadlines
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