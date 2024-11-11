classdef Flow < handle
    
    % AUTHOR: Afaf
    % LAST MODIFIED: 18/11/2020
    
    properties
        id;                      % id of flow
        idCoflow;                % the id of the coflow to which the flow belongs
        volume;                  % flow size in Mbit
        volume_initial;          % saving initial volume value
        remainingVolume;         % remaining size in Mbit
        source;                  % the source of the flow (one of the ingress ports of the fabric)
        destination;             % the destination of the flow (one of the egress ports of the fabric)
        d_rate = 1;              % dynamic rate
        d_rate_old = 1,          % for updating prices !!!
        ad_rate = 1;             % adjusted rate
        fct = 0;                 % flow completion time (initializing to zero (unit = second)) --> 1Paris FCTs
        state_f = 1;             % 1 if remainingVolume>0 0 if remainingVolume = 0
        price;                   % overall price on the path used by the flow
        arrival = 0;             % slot of arrival of COFLOW (dynamic case) (offline case arrival = 0)
        departure = -1;          % slot of departure of the FLOW
        links;                   % links used by the flow
        fct_pricing = 0;         % flow completion time (initializing to zero (unit = second)) --> pricing FCTs (with volume update)
        fct_pricing_2 = 0;       % flow completion time (initializing to zero (unit = second)) --> pricing FCTs (without volume update)
        
        % source and destination are chosen such that they do not correspond to the same machine
    end
    
    methods
        function obj = Flow(id, coflowID)
            % FLOW Construct an instance of this class
            obj.id = id;
            obj.idCoflow = coflowID;
        end
    end
end

