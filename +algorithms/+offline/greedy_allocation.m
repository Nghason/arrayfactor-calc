function outputs = greedy_allocation(fabric, coflows, prio_order)
% Process the Greedy Rate Allocation Algorithm presented in sincronia paper
%
% Cedric Richier, LIA
% (c) LIA, 2020

%% Initialisations:

% Store the priorities order in the output
%outputs.prio_order = prio_order;

% Store indexes of unfinished coflows sorted by ascendant priority
unfinished = prio_order;

% Number of coflows
n_coflows = length(coflows);

% Set of indexes of finished coflows
finished=[];

% Link capacities
tmp = [[fabric.machinesPorts.ingress] [fabric.machinesPorts.egress]];
B = [tmp.linkCapacity];

% Initialising the final cct of each coflow to zero
final_ccts = zeros(1,n_coflows);

% total number of flows
n_flows = sum([coflows.numFlows]);

% Initializing a (n_flows x n_flows) matrix to store rates of all flows at each step:
All_rates_by_step = zeros(n_flows,n_flows);

% Store flow completion time of all flows:
all_fcts = zeros(n_flows,1);

% Set d_rates to zero for all flows:
for c = coflows
    for f = c.flows
        f.d_rate = 0;
    end
end

% Store the number of flows by coflows:
n_flows_by_coflows = [coflows.numFlows];

% Array that stores the offset of flows according to their coflows
index_bases = [0 cumsum(n_flows_by_coflows(1:n_coflows-1))];

% Counter for number of steps:
n_steps = 0;

instants = [];

%% Main loop:

n_links = fabric.numFabricPorts;

while ~isempty(unfinished)
    
    % Incrementing number of steps:
    n_steps = n_steps+1;
    
    % Greedy flow allocation
    for k = unfinished
        c = coflows(k);
        for f = c.flows
            % Allocation to flows that are not finished
            if f.state_f % 0: finished, 1: unfinished
                % Check bandwidth avalability on flow's path
                if sum(B(f.links)) == 2 % NOTE: capacity of ALL links is 1
                    % Allocate all bandwidth to flow
                    f.d_rate = 1;
                    % Update bandwidth on path
                    B(f.links) = 0;
                end
            end
        end
    end
 
    % Get flows' remaining volume and rate in flow objects
    tmp = [coflows.flows];
    f_r_vols = [tmp.remainingVolume];
    f_rates  = [tmp.d_rate];
    
    % Store rates of flows for current step:
    All_rates_by_step(:,n_steps) = f_rates';     
    
    % Store only remaining volume of flows with non-zero rate (active volumes):
    active_f_vols = f_r_vols.*(f_rates>0);
    
    % Compute the minimum duration to end at least one of all active flows:
    %  - Divide remaining volumes of active flows by their positive rate
    durations = active_f_vols(f_rates>0)./f_rates(f_rates>0);
    %  - Get the minimum duration:
    min_duration = min(durations);
    
    %fprintf('step %d: \n', n_steps);
    instants = [instants, min_duration];
    
    % Update remaining volumes of active flows:
    f_r_vols(f_rates>0) = f_r_vols(f_rates>0)-min_duration*f_rates(f_rates>0);
    
    % Store the number finished coflows:
    len_finished = length(finished);
    
    % Update fcts, remaining volume, d_rate and state_f of all flow objects
    % that are not finished
    for k = unfinished        
        c = coflows(k);
        if c.state_c
            for f = c.flows
                if f.state_f
                    % Reset rate of flow to zero for next step
                    f.d_rate = 0;
                    
                    % Update flow completion time
                    f.fct = f.fct + min_duration;
                    
                    %f.fct c.deadline
                    
                    % Update remaining volume of flow
                    f.remainingVolume = f_r_vols(index_bases(c.id)+f.id);
                    if ~f.remainingVolume % flow ends now
                        % Update the unfinished flow status
                        f.state_f = 0;
                        % Update the array of fcts:
                        all_fcts(index_bases(c.id)+f.id) = f.fct;
                    end
                end
            end            
            % Check if coflow is finished
            tmp = [c.flows];
            c_f_states = [tmp.state_f];
            if ~sum(c_f_states) % all flows of coflow are done
                % Update the active status of coflow
                c.state_c = 0;
                % Store the final CCT of coflow
                final_ccts(c.id) = max([tmp.fct]);
                % Update the set of indexes for finished coflows
                finished = union(finished,c.id);
            end            
        end        
    end
    
    % Update the set of unfinished coflows if necessary
    if len_finished < length(finished)
        unfinished = setdiff(unfinished,finished,'stable');
    end
    
    % Reset available bandwidth on all links for next step (i.e., allow preemption <18/11/2021>)
    B(:) = 1; % NOTE: capacity of ALL links is 1
    
end % end while ~isempty(unfinished)

% Reset properties of coflow objects and their flow objects:
% utils.resetCoflows(coflows);
for c = coflows
    c.state_c = 1;
    %c.priority = 0;
    c.prices = zeros(1,length(c.prices));
    c.maxPrices = c.prices;
    c.prices_prev = c.prices;
    c.stability = 0;
    c.current_w = 0;
    c.max_diff_prices = 0;
    c.ts_counter = 0;
    c.remaning_vol = 0;
    c.departure = -1;
    for f = c.flows
        f.remainingVolume = f.volume;
        f.d_rate = 1;
        f.fct = 0;
        f.state_f = 1;
    end
end



%% Formating outputs:
instants = [0, instants];
outputs.instants = instants;
outputs.avg_cct = mean(final_ccts);
outputs.ccts = final_ccts;
outputs.fcts_all = all_fcts;
if n_steps < n_flows
    All_rates_by_step(:,n_steps+1:n_flows) = [];
end
outputs.rates = All_rates_by_step;
outputs.total_steps = n_steps;

end