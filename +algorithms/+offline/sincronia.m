% (c) 2022 Cedric Richier (LIA-Avignon)
% Implementation of Sincronia ordering: Bottleneck-Select-Scale-Iterate Algorithm
% Reference: S. Agarwal et al., “Sincronia: Near-optimal network design for coflows,” in Proc. ACM SIGCOMM, 2018, pp. 16–29.
% Inputs:
% - fabric: the Fabric ojbect representing the network
% - coflows: an array of Coflow objects within the fabric
% Outputs:
% - order: scheduling order of coflows


function sincronia_order = sincronia(fabric,coflows)



%% Initializations:
% Number of links
n_links = fabric.numFabricPorts;

% Number of coflows
n_coflows = length(coflows);

% Unscheduled coflows IDs:
unsch_coflow_ids = [coflows.id];

% Demand of each coflow on each link:
D = zeros(n_links, n_coflows);

% Permutation
sincronia_order = zeros(1,n_coflows);

% Compute D:
for c = coflows
    D(:,c.id) = c.indicator*[c.flows.volume]';
end

% weights (initialized to one)
W = ones(1,n_coflows);

% Last index in permutation
K = n_coflows;

% index in permutation:
k = K;

%% Main loop: starting by finding the last coflow to schedule
while k > 0
    
    % Find the most bottlenecked links:
    cumulD = sum(D,2);
    b_canditates = find(cumulD == max(cumulD));
    % Randomly pick one such link:
    r_ind = randi(length(b_canditates));
    b = b_canditates(r_ind);   
    %% TEST: to match implementation of example in sincronia paper:
%    b = max(b_canditates(b_canditates<=4));
    % end TEST
    %%
    
    % Select weighted largest coflow to schedule last
    % argmin w_c/D_c_b: the min is computed among coflows that use link b
    set_idx = find(D(b,:)>0);    
    c_candidates = find( W./D(b,:) == min(W(set_idx)./D(b,set_idx)));
    % Randomly pick one such coflow:
    r_ind = randi(length(c_candidates));
    sincronia_order(k) = c_candidates(r_ind);
    %perm(k) = c_id;
    
    % Scale the weights:
    unsch_coflow_ids = setdiff(unsch_coflow_ids,sincronia_order(k));
    for c = coflows(unsch_coflow_ids)
        W(c.id) = W(c.id) - W(sincronia_order(k))*D(b,c.id)/D(b,sincronia_order(k));
    end
    
    % Set demand of last scheduled coflow to zero:
    D(:,sincronia_order(k)) = 0;
    W(sincronia_order(k)) = 1;
    
    % Set k:
    k = k-1;
end


end




