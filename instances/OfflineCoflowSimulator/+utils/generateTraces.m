% Random generator of fabric and coflow objects

function outputs = generateTraces(trace_type, inputs, n_machines, n_arrivals)

global params;

scenario_onoff = inputs.scenario_onoff;

ARRIVAL_RATE = inputs.ARRIVAL_RATE;         % poisson arrival rate (for online simulation only)
FLOW_RANGE   = inputs.FLOW_RANGE;           % range of nb of flows (for Facebook traces only)
CLASS_ONE_PROBA = inputs.CLASS_ONE_PROBA;   % probability to decide coflow type, default: [0.6, 0.4]
PRIO_CLASS_NUM = inputs.PRIO_CLASS_NUM;     % nb of priority classes
PRIO_PROBAS = inputs.PRIO_PROBAS;           % probabilities of choosing priority class for coflows
PRIO_WEIGHTS = inputs.PRIO_WEIGHTS;         % weight associated to each priority class     




%% Nb of coflows is derived from n_arrivals and n_coflows_per_arrival
switch scenario_onoff 
    case 'offline'
        n_coflows_per_arrival = ones(1, n_arrivals); 
        % use no batch distribution in offline simulation
        inputs.batch_dist = '';
    case 'online' 
        n_coflows_per_arrival = zeros(1, n_arrivals); 
        for i = 1:n_arrivals
            % if no specific batch distribution is declared -> one coflow per arrival
            if strcmp(inputs.batch_dist, '') 
                n_coflows_per_arrival(i) = 1;
            elseif strcmp(inputs.batch_dist, 'Uniform') 
                n_coflows_per_arrival(i) = inputs.n_coflows_per_arrival(i);
            else
                n_coflows_per_arrival(i) = inputs.pd.random;
            end
        end
end

% total nb of coflows
n_coflows = sum(n_coflows_per_arrival);


%% Run generator: {'random', 'facebook'}

% Fabric setup
link_caps = [1, 1, 1];
config.NumMachines = n_machines;
config.minNumMachines = n_machines; 
config.maxNumMachines = n_machines;
config.linkCapacitiesAvailable = link_caps; % possible link capacities (should not exceed 1x3 vector)

% Coflow setup
switch trace_type
    case 'random'   % Synthetic traces
        
        config.minNumCoflows = n_coflows;
        config.maxNumCoflows = n_coflows;
        
        % Setting parameters:
        [vol_avg, vol_std] = deal(1, 0.2);  % default: [1, 0.2]
        
        % Parameters
        config.architecture_type = 'twoTypeCoflows_architectureVol';
        config.typeCoflowProba = [CLASS_ONE_PROBA, 1-CLASS_ONE_PROBA];   % probability to decide of coflows type, default: [0.6, 0.4]
        % type 1: coflows with one flow
        % type 2: coflows with multiple flows
        config.avgFlowVolume = vol_avg;       % average size of flow*
        config.standardDivVolume = vol_std;   % standard deviation (since flows' volumes will be generated according to normal dist.       
        config.ratioClass = 0.8;              % volume ratio of class 2/class 1, default: 0.8
    case 'facebook' % Facebook traces
        
        % File path for exported CSV
        csv_out = params.csv_out;                   
        
        % Load Facebook raw coflow data
        load(inputs.fb_raw); % contains struct 'coflows_raw'
        
        %---------------------------------------------------------------------------------------
        % Construct CSV file from coflow raw data of Facebook
        % Content of function buildFBTracesCSVFile_maxVol(raw_struct, max_n_flows, csv_out)
        % usage example: buildFBTracesCSVFile_maxVol(coflows, 200, 'ouput.csv');

        % Filter: min and max nb of flows per coflow
        filter = coflows_raw.n_map(:).*coflows_raw.n_red(:) >= FLOW_RANGE(1) & ...
                 coflows_raw.n_map(:).*coflows_raw.n_red(:) <= FLOW_RANGE(2);

        filter_struct.n_machines    = coflows_raw.n_machines;
        filter_struct.n_coflows     = sum(filter);
        filter_struct.id            = (1:filter_struct.n_coflows)';
        filter_struct.arrival_time  = coflows_raw.arrival_time(filter);
        filter_struct.n_map         = coflows_raw.n_map(filter);
        filter_struct.map_list      = coflows_raw.map_list(filter);
        filter_struct.n_red         = coflows_raw.n_red(filter);
        filter_struct.red_list      = coflows_raw.red_list(filter);
        filter_struct.red_shuffle   = coflows_raw.red_shuffle(filter);

        %---------------------------------------------------------------------------------------
        % Content of function buildFBTracesCSVFile(filter_struct, csv_out)

        % nb of coflows
        n_coflows_original = filter_struct.n_coflows;

        % total nb of flows
        tot_n_flows = sum(filter_struct.n_map.*filter_struct.n_red);

        % initialize table
        flowID    = zeros(tot_n_flows,1);
        coflowID  = zeros(tot_n_flows,1);
        src       = zeros(tot_n_flows,1);
        dst       = zeros(tot_n_flows,1);
        bandwidth = zeros(tot_n_flows,1);
        weight    = zeros(tot_n_flows,1);

        T = table(flowID,coflowID,src,dst,bandwidth,weight);
        n_line = 0;
        for k = 1:n_coflows_original
            c_ID = filter_struct.id(k)-1;
            map_list = filter_struct.map_list{k};
            red_list = filter_struct.red_list{k};
            red_shuffle = filter_struct.red_shuffle{k};
            for r = 1:filter_struct.n_red(k)
                for m = 1:filter_struct.n_map(k)
                    n_line = n_line+1;
                    T.flowID(n_line)    = n_line-1;
                    T.coflowID(n_line)  = c_ID;
                    T.src(n_line)       = map_list(m);
                    T.dst(n_line)       = red_list(r);
                    T.bandwidth(n_line) = red_shuffle(r)/filter_struct.n_map(k);
                    T.weight(n_line)    = 1;
                end
            end
        end

        % modify table T to adapt with the desired n_machines (using modulo). Ref: function utils.reduceFB
        for i = 1:size(T, 1)
            T.src(i) = mod(T.src(i) - 1, n_machines) + 1; % normalized ingress port
            T.dst(i) = mod(T.dst(i) - n_machines - 1, n_machines) + n_machines + 1; % normalized egress port
        end

        % table of each coflow
        for i = 1:n_coflows_original
            row_id = T.coflowID == i-1;
            tmp.flowID = T.flowID(row_id); % in T, id of coflow and flow starts from 0
            tmp.coflowID = T.coflowID(row_id);
            tmp.src = T.src(row_id);
            tmp.dst = T.dst(row_id);
            tmp.bandwidth = T.bandwidth(row_id);
            tmp.weight = T.weight(row_id);
            c_tab = table(tmp.flowID, tmp.coflowID, tmp.src, tmp.dst, tmp.bandwidth, tmp.weight);
            c_tab.Properties.VariableNames = {'flowID', 'coflowID', 'src', 'dst', 'bandwidth', 'weight'};
            CoflowTable{i} = c_tab;
        end
        
        % each coflow take a random id from coflows_original of facebook
        gen_ids = randi(n_coflows_original, [1, n_coflows]);

        T_with_n_coflows = table;
        for i = 1:n_coflows
            % coflows(i) = coflows_original(tmp(i));
            tmp = CoflowTable{gen_ids(i)};
            % reset id of coflow (starts from 0)
            tmp.coflowID = (i)*ones(size(tmp,1), 1);

            T_with_n_coflows = vertcat(T_with_n_coflows, tmp);
        end

        % reset if of flows (starts from 0)
        T_with_n_coflows.flowID = (1 : size(T_with_n_coflows,1))';

        % export CSV file
        writetable(T_with_n_coflows, csv_out);

        % Input parameters to build coflows
        config.architecture_type = 'csv_architecture';
        config.filename = csv_out;
        config.format = 0; % 0 by default
end

% Construct fabric and coflows objects
config = utils.fabric_generation(config);
config = utils.coflows_generation(config);



% Generating coflows
% CONFIG.coflows = [];
% CONFIG.NumCoflows = n_coflows; 
% for ii = 1:CONFIG.NumCoflows
%     if (ii == 1) % insure that there is at least 1 coflow of class 1
%         flag_type = 0;
%         CONFIG.coflows = [CONFIG.coflows, network_elements.Coflow2classesVol(ii, CONFIG, flag_type)];
%     elseif (ii == 2) % insure that there is at least 1 coflow of class 2
%         flag_type = 1;
%         CONFIG.coflows = [CONFIG.coflows, network_elements.Coflow2classesVol(ii, CONFIG, flag_type)];
%     else % decide of the coflow class randomly
%         flag_type = 2;
%         CONFIG.coflows = [CONFIG.coflows, network_elements.Coflow2classesVol(ii, CONFIG, flag_type)];
%     end
% end       



% Get main variables
fabric      = config.fabric;
coflows     = config.coflows;
n_links     = fabric.numFabricPorts;     % nb of fabric ports (ingress+egress)    
n_coflows   = length(coflows);           % nb of coflows
n_flows     = [coflows.numFlows];        % nb of flows of each coflow
n_flows_all = sum(n_flows);             % total nb of flows


%% Find bottleneck

portCapacity = [[fabric.machinesPorts.ingress] [fabric.machinesPorts.egress]];
portCapacity = [portCapacity.linkCapacity];

% Full rate processing time of each coflow on each link:
port_load = zeros(n_links, n_coflows);

% Compute P (load on each port of each coflow)
for c = coflows
    port_load(:,c.id) = (c.indicator*c.getFlowsVolume')./portCapacity';
    
    % save volumes of each flows in coflow c
    c.addParam.volumes = [c.flows.volume_initial];
     
    % CCT isolated of coflow c    
    c.addParam.CCT0 = max(port_load(:, c.id));
    
end



%% Arrival times and deadlines

% Arrival time of each coflow
arrival_times = zeros(1, n_coflows);
switch scenario_onoff 
    case 'offline'
        arrival_times = zeros(1, n_coflows);
    case 'online'  
        % Arrival times are exponentially distributed with rate Lambda
        % Lambda = 1.5;          % lambda   (events per unit time)
        NumEvents = n_arrivals;  % number of events to generate
        InterEventTimes = exprnd(1/ARRIVAL_RATE, NumEvents, 1);
        EventTimes = cumsum(InterEventTimes);
        k = 1;
        for i = 1:n_arrivals
            for j = 1:n_coflows_per_arrival(i)
                arrival_times(k) = EventTimes(i);
                k = k + 1;
            end
        end
end

% Priority class and corresponding weights
class_draw = randsample(PRIO_CLASS_NUM, n_coflows, true, PRIO_PROBAS);
    
% Assign coflow parameters
for c = coflows      
    % Indicator of used links of the whole coflow
    c.addParam.indicator = (sum(c.indicator, 2) > 0);
    c.addParam.used_links = find(c.addParam.indicator == 1);
    
    % Arrival time
    c.arrival = arrival_times(c.id); % 0; arrival time
    c.arrival = round(c.arrival, 5);
    
    % Deadline: a + (b-a).*rand(N,1)
	a = inputs.DEADLINE_RANGE(1) * c.addParam.CCT0;
	b = inputs.DEADLINE_RANGE(2) * c.addParam.CCT0;
    temp = a + rand(1)*(b-a);    
    c.deadline = round(c.arrival + temp, 5);
    
    % Assign priority class
    c.priority = class_draw(c.id);
    
    % Assign weight according to priority class
    c.weight = PRIO_WEIGHTS(c.priority);
    
end


outputs.coflows = coflows;
outputs.fabric = fabric;
outputs.n_coflows_per_arrival = n_coflows_per_arrival;

end

