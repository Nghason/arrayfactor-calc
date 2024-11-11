function outputs = varys(fabric,coflows)
% Simulates a workload of coflows where the resource management is done
% according to the method used in Varys system in an offline scenario
% The function returns an outputs variable that stores the average 
% completion time of coflows, a list of coflow completion times and the 
% rates allocated to all flows at each iteration of the process
% Parameters:
% - fabric: the Fabric ojbect representing the network
% - coflows: an array of Coflow objects within the fabric
% Cedric Richier, LIA
% (c) LIA, 2021

%% Formating inputs

% A parameter to manage roundings
precision = 5;

% Capacities of links
tmp = [[fabric.machinesPorts.ingress] [fabric.machinesPorts.egress]];
Cap = [tmp.linkCapacity];

% Remaining volumes of flows
tmp = [coflows.flows];
f_rvol = [tmp.volume];

% List of active coflows to schedule
coflows_to_schedule = coflows;

% Number of flows by coflows
n_flows_by_coflows = [coflows.numFlows];

% Number of coflows
n_coflows = length(coflows);

% Array of index bases used to retrieve flow indexes according to their
% coflow
index_bases = [0 cumsum(n_flows_by_coflows(1:n_coflows-1))];

% Array for coflow completion times
cct = zeros(1,n_coflows);

% Structure used to store the rate of each flow at each step of the main
% loop
rates = cell(n_coflows,n_coflows);

prio_order = cell(n_coflows,1);

% Iteration index
s=0;


%% Main Loop
% Main loop in Varys' method: at each iteration, one coflow is prioritized
% and is allocated some resources in order to finish during this iteration
while(~isempty(coflows_to_schedule))
    % Incrementing the number of iterations
    s = s+1;
    
    % Number of coflows to schedule during this iteration
    n_coflows_to_schedule = length(coflows_to_schedule);
    
    % Initializing the values of gamma to zero for each coflow in the list
    gammas = zeros(1,n_coflows_to_schedule);
    
    % Setting the remaining bandwidth to full capacity for each link in the
    % fabric
    Rem = Cap;
    
    % array used to store ids of coflows that will end during this
    % iteration 
    finished_coflows_ids = zeros(1,n_coflows);
    
    % Iterator over the list of coflows to schedule
    k = 0;
    
    % Computing the values of gamma for each remaining coflow to schedule
    for c = coflows_to_schedule        
        k = k+1;
        c_flows = c.flows;
        c_rvol = f_rvol(index_bases(c.id)+[c_flows.id]);
        gammas(k) = utils.computeCCT0(c,c_rvol,Rem);
    end
    
    % Ordering the coflows according to ascendent gamma order
    [Ordered_gammas, c_Index] = sort(gammas);
    coflows_to_schedule = coflows_to_schedule(c_Index);
    
    % Storing ids of active coflows
    active_coflows_ids = [coflows_to_schedule.id];
    prio_order{s} = active_coflows_ids;
    
    % Storing the minimum value of all gammas (representing the completion
    % time of all flows that will finish during this iteration)
    min_fct = Ordered_gammas(1);
    
    % Applying MADD algorithm to the ordered list of coflows to schedule
    % MADD = Minimum Allocation for Desired Duration
    % It aims at allocating the minimum bandwidth to flows of a coflow such
    % that all the flows finish in the same time as the slowest one
    for co = coflows_to_schedule
        % Storing the id of the current coflow
        co_id = co.id;
        
        % Getting the index base used to retrieve flows
        index_base = index_bases(co_id);
        
        % Initializing rates to zero for all the flows of this coflow
        rates{s,co_id} = zeros(1,co.numFlows);
        
        % Computing the current value of gamma of this coflow according to the remaining bandwidth
        co_flows = co.flows;
        co_rvol = f_rvol(index_bases(co.id)+[co_flows.id]);
        gamma_co = utils.computeCCT0(co,co_rvol,Rem);
        
        % Assuming the coflow will finish during this iteration
        finished = 1;
        
        % Allocating bandwidth to each flow of the coflow
        for fl = co.flows
            % Storing the flow id
            fl_id = fl.id;
            
            % Computing the rate according to the remaining volume of the
            % flow and the value of gamma of the corresponding coflow
            rates{s,co_id}(fl_id) = round(f_rvol(index_base+fl_id)/gamma_co,precision);
            
            % Getting the input and output links of this flow
            l_ind = find(co.indicator(:,fl_id));
            
            % Updating remaining bandwidth according to the allocated
            % resource
            Rem(l_ind) = max(round(Rem(l_ind)-rates{s,co_id}(fl_id),precision),0);
            
            % Updating the volume of the current flow
            f_rvol(index_base+fl_id) = round(f_rvol(index_base+fl_id) - min_fct*rates{s,co_id}(fl_id),precision-2);

            if f_rvol(index_base+fl_id) > 0
                % If there's at least one flow with non zero remaining volume, then the
                % coflow will not finish during this iteration
                finished = 0;
            end
        end
        
        if finished
            % Updating the array that stores ids of finished coflows
            finished_coflows_ids(co_id) = 1;
        end
    end
    
    % Updating the cct of active coflows
    cct(active_coflows_ids) = cct(active_coflows_ids) + round(min_fct,precision);
    
    % Removing finished coflows from the list of coflows to schedule
    coflows_to_schedule(ismember(active_coflows_ids,find(finished_coflows_ids))) = [];
end

%% Formating outputs

% Computing the average coflow completion time
avg_cct = mean(cct);

% Array of completion time of each coflow
outputs.cct = cct;

% Average coflow completion time
outputs.avg_cct = avg_cct;

% Rates allocated to flows at each iteration
outputs.rates = rates;

% Priority order of coflows in each step
%outputs.prio_order = prio_order;

end
