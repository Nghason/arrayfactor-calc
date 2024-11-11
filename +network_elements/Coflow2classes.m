classdef Coflow2classes < network_elements.Coflow
    % DEPRECATED
    % AUTHOR: Afaf
    % LAST MODIFIED: 14/12/2020
    % DEPRECATED: 26/01/2021
    
    methods
        function obj = Coflow2classes(coflowID, architecture_config)
            
            obj = obj@network_elements.Coflow(coflowID);

            set = [];
            for i = 1:length(architecture_config.typeCoflowProba)
                set = [set,i*ones(1,architecture_config.typeCoflowProba(i)*10)];
            end
            index = randi([1,length(set)]);
            
            switch set(index)   % randomly pack one of the traffic models according to the aPrioriPdf
                case 1
                    obj.numFlows = 1;
                case 2
                    obj.numFlows = randi([ceil(architecture_config.NumMachines/2) architecture_config.NumMachines]);
            end
            
%             if rand <= architecture_config.typeCoflowProba % coflow with 1 flow
%                 obj.numFlows = 1;
%             else % coflow with multiple flows
%                 obj.numFlows = randi([ceil(architecture_config.NumMachines/2) architecture_config.NumMachines]);
%             end
            
            fprintf('Number of flows %d \n', obj.numFlows)
            
            %avgFlowVolume = architecture_config.avgFlowVolume / obj.numFlows;
            avgFlowVolume = architecture_config.avgFlowVolume ;
            
            fabricMachinesPorts = architecture_config.fabric.machinesPorts;
            
            
            for ii = 1:obj.numFlows
                obj.flows = [obj.flows network_elements.Flow2classes(ii,avgFlowVolume, fabricMachinesPorts, coflowID, architecture_config.standardDivVolume)]; % creating the flows and assigning a size and source/destination to each one
            end
            
            numPorts = 2*length(fabricMachinesPorts);
            
            obj.indicator = zeros(numPorts,obj.numFlows);
            
%             for jj = 1:numPorts
%                 for ii = 1:obj.numFlows
%                     if (obj.flows(ii).source.id == jj)
%                         obj.indicator(jj,ii) = 1;
%                     elseif (obj.flows(ii).destination.id == jj + numPorts/2)
%                         obj.indicator(jj + numPorts/2, ii) = 1;
%                     end
%                 end
%             end
            
            for f = obj.flows
                ii = f.source.id;
                jj = f.destination.id;
                obj.indicator(ii,f.id) = 1;
                obj.indicator(jj,f.id) = 1;
            end
            
            obj.prices = zeros(1,numPorts);
            obj.maxPrices = obj.prices;
            obj.prices_prev = obj.prices;
        end
    end
end
