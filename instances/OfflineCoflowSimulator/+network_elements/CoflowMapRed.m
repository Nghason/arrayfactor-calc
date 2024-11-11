classdef CoflowMapRed < network_elements.Coflow
    %COFLOWMAPRED Summary of this class goes here
    %   Detailed explanation goes here
    
    % LAST MODIFIED: 26/01/2021
    
    methods
        
        function obj = CoflowMapRed(coflowID, architecture_config)
            
            obj = obj@network_elements.Coflow(coflowID);
            
            minMaxMap              = architecture_config.minMaxMap;
            minMaxRed              = architecture_config.minMaxRed;
            avgFlowVolume          = architecture_config.avgFlowVolume;
            fabricMachinesPorts    = architecture_config.fabric.machinesPorts;
            numMachines = architecture_config.NumMachines;
            
            numMap = randi(minMaxMap,1,1); % # of mappers assigned to the current coflow
            numRed = randi(minMaxRed,1,1); % # of reducers assigned to the current coflow
            
            fprintf('Number of mappers -> %d, ', numMap);
            fprintf('Number of reducers -> %d, ', numRed);
            
            sources = [fabricMachinesPorts(randperm(numMachines, numMap)).ingress]; % random vector of different ingresses
            destinations = [fabricMachinesPorts(randperm(numMachines, numRed)).egress]; % random vector of different egresses
            
            ii = 0;
            
            % creating flows :
            %   flows source and destination should not be on the same machine
            %   2 flows of the same coflow should not share the exact same path
            for src = sources
                for dst = destinations
                    if(src.id + numMachines ~= dst.id)
                        ii = ii+1;
                        obj.flows = [obj.flows network_elements.FlowMapRed(ii,avgFlowVolume,coflowID,src,dst)];
                    else % if there is one mapper and one reducer and they are on the same machine
                        if (numMap == 1 && numRed == 1)
                            ii = ii+1;
                            a = fabricMachinesPorts(randperm(numMachines, numRed)).egress;
                            while a.id == dst.id
                                a = fabricMachinesPorts(randperm(numMachines, numRed)).egress;
                            end
                            obj.flows = [obj.flows network_elements.FlowMapRed(ii,avgFlowVolume,coflowID,src,a)];
                        end
                    end
                end
            end
            obj.numFlows = length(obj.flows);
            
            fprintf('Number of flows -> %d \n', obj.numFlows)
            
            
            numPorts = 2*numMachines;
            
            obj.indicator = zeros(numPorts,obj.numFlows);
            
            %                     for jj = 1:numPorts
            %                         for ii = 1:obj.numFlows
            %                             if (obj.flows(ii).source.id == jj)
            %                                 obj.indicator(jj,ii) = 1;
            %                             elseif (obj.flows(ii).destination.id == jj + numPorts/2)
            %                                 obj.indicator(jj + numPorts/2, ii) = 1;
            %                             end
            %                         end
            %                     end
            for flow = obj.flows
                obj.indicator(flow.source.id,flow.id) = 1;
                obj.indicator(flow.destination.id,flow.id) = 1;
            end
            
            obj.prices = zeros(1,2*numMachines);
            %obj.prices(find(sum(obj.indicator,2)>=1)) = 1;
            obj.maxPrices = obj.prices;
            obj.prices_prev = obj.prices;
            
        end
        
    end
    
end

