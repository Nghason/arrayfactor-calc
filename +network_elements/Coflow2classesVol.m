classdef Coflow2classesVol < network_elements.Coflow
    
    % AUTHOR: Afaf
    % LAST MODIFIED: 07/01/2021
    
    properties
        coflowClass;
    end
    
    methods
        function obj = Coflow2classesVol(coflowID, architecture_config,flag_type)
            
            obj = obj@network_elements.Coflow(coflowID);

            switch flag_type
                case 0 % the coflow belongs to class 1 (insure that there is at least 1 coflow of class 1)
                    obj.numFlows = 1;
                    obj.coflowClass = 1;
                case 1 % the coflow belongs to class 2 (insure that there is at least 1 coflow of class 2)
                    obj.numFlows = randi([ceil(2*architecture_config.NumMachines/3) architecture_config.NumMachines]);
                    obj.coflowClass = 2;
                case 2 % decide of the coflow class randomly
                    set = [];
                    for i = 1:length(architecture_config.typeCoflowProba)
                        set = [set,i*ones(1, uint8(architecture_config.typeCoflowProba(i)*10))];
                    end
                    index = randi([1,length(set)]);
                    
                    switch set(index)   % randomly pack one of the traffic models according to proba dist.
                        case 1
                            obj.numFlows = 1;
                            obj.coflowClass = 1;
                        case 2
                            obj.numFlows = randi([ceil(2*architecture_config.NumMachines/3) architecture_config.NumMachines]);
                            obj.coflowClass = 2;
                    end
            end
            
            %fprintf('Number of flows %d \n', obj.numFlows)
            
            %avgFlowVolume = architecture_config.avgFlowVolume / obj.numFlows;
             if (obj.coflowClass == 1)  
                avgFlowVolume1 = architecture_config.avgFlowVolume;
             else 
                 avgFlowVolume2 = architecture_config.avgFlowVolume*architecture_config.ratioClass;
             end
            
            fabricMachinesPorts = architecture_config.fabric.machinesPorts;
            
            %% find all possible permutations (flows of the same coflow will not share the same path)
            ingress_id = [];
            egress_id = [];
            for ii = 1:length(fabricMachinesPorts)
                ingress_id = [ingress_id fabricMachinesPorts(ii).ingress.id];
                egress_id = [egress_id fabricMachinesPorts(ii).egress.id];
            end
            
            [A,B] = meshgrid(ingress_id,egress_id);
            c=cat(2,A',B');
            permutation_mat_ports = reshape(c,[],2); 
            %%
            
            
            for ii = 1:obj.numFlows
                % source and destination of each coflow
                numFabricMachines = length(fabricMachinesPorts); % # of machines
                p = randi([1 size(permutation_mat_ports,1)],1,1);
                source = fabricMachinesPorts(permutation_mat_ports(p,1)).ingress;
                if (mod(permutation_mat_ports(p,2),numFabricMachines) == 0)
                    dest = fabricMachinesPorts(numFabricMachines).egress;
                else
                    dest = fabricMachinesPorts(mod(permutation_mat_ports(p,2),numFabricMachines)).egress;
                end
                if (obj.coflowClass == 1)
                    obj.flows = [obj.flows network_elements.Flow2classes(ii,avgFlowVolume1, coflowID, architecture_config.standardDivVolume, source, dest)]; % creating the flows and assigning a size and source/destination to each one
                else
                    obj.flows = [obj.flows network_elements.Flow2classes(ii,avgFlowVolume2, coflowID, architecture_config.standardDivVolume, source, dest)]; % creating the flows and assigning a size and source/destination to each one 
                end
                permutation_mat_ports(p,:) = [];
            end
            
            numPorts = 2*length(fabricMachinesPorts);
            
            obj.indicator = zeros(numPorts,obj.numFlows);
            
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

