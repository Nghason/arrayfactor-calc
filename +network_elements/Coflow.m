classdef Coflow < handle
    
    % AUTHOR: Afaf & Cedric
    % LAST MODIFIED: 18/05/2021 by Trung: added 'addParam' parameter
    
    properties
        id;                 % ID of the coflow
        numFlows;           % total # of flows forming the coflow
        flows;              % the flows forming the coflow (array of objects)
        indicator;          % indicates whether flow j of coflow k uses link i
        weight = 1;         % weight of the coflows. Default value is one
        priority = 0;       % will change according to the order that maybe given to coflows
        state_c = 1;        % (initially all flows are active) if there is at least one flow active state_c = 1 if no flow is active state_c = 0
        prices;             % vector of prices on each link
        maxPrices;          % maximal value reached
        arrival = 0;        % slot of arrival of coflow (dynamic case) (offline case arrival = 0)
        departure = -1;     % slot of departure of the LAST FLOW
        %% price convergence
        stability = 0;      % vector of prices converges stability = 1 otherwise stability = 0
        current_w = 0;      % counter of the number of slots where the vector of prices converges (set to 0 whenever the change is bigger than epsi)
        prices_prev;        % vector of prices of the previous slot
        max_diff_prices = 0;    % Max_diff_prices = max(|prices - prices_prev|)
        ts_counter = 0;         % total number of slots (during all simulation) where the vector of prices is stable
        remaning_vol;   % remaining volume of coflow (sum of remaining volumes of its flows)
        deadline = -1; % deadline for the coflow
        
        %% additional parameters, gathered into the struct addParam
        % Examples:
        % c.addParam.CCT0: isolated CCT of coflow c
        % c.addParam.t_arrival: arrival time
        % c.addParam.t_deadline: deadline
        % c.addParam.k_arrival: arrival time slot
        % c.addParam.k_deadline: deadline time slot
        % c.addParam.slot_id: array of time slots between [k_arrival, k_deadline]
        addParam = {};

    end
    
    methods
        function obj = Coflow(coflowID) % Coflow Construct an instance of this class
            obj.id = coflowID;
            obj.remaning_vol = 0;
        end
        
        function flowsVolume = getFlowsVolume(obj) % !!!! Cédric !!!
            tmp = [obj.flows];
            flowsVolume = [tmp.volume];
        end
    end
end

