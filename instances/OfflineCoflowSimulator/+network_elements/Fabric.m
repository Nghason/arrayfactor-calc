classdef Fabric < handle
    
    % AUTHOR: Afaf
    % LAST MODIFIED: 17/12/2020
    
    properties
        numFabricPorts;          % the total # of DC fabric's ports
        machinesPorts;           % machines of the DC fabric  
                                 % we assume that a machine has only 2 ports (ingress / egress)
        sumRatesOnLink;          % vector of sum of rates of flows using a link 
        sumAdjustedRatesOnLink;  % vector of sum of adjusted rates of flows using a link 
    end
    
    methods
        
        
        function obj = Fabric(NumMachines)
            
            % INPUTS:
            % NumMachines: number of machines in the fabric
            
            obj.numFabricPorts = 2*NumMachines;
           % fprintf('* Number of Fabric ports %d \n', obj.numFabricPorts)
            
            
            %% Used in pricing
            obj.sumRatesOnLink = zeros(obj.numFabricPorts,1);
            obj.sumAdjustedRatesOnLink = zeros(obj.numFabricPorts,1);
        end
        
%         function setLinkCapacity(obj,port_id,cap)
%             obj.machinesPorts(port_id).ingress.linkCapacity = cap;
%             obj.machinesPorts(port_id).egress.linkCapacity = cap;
%         end
%         
%         function setLinkCapacities(obj,cap)
%             numMachines = obj.numFabricPorts/2;          
%             for i = 1:numMachines
%                 obj.setLinkCapacity(i,cap);
%             end
%         end
        
    end
end

