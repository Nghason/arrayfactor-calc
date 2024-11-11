classdef Fabric2Classes < network_elements.Fabric
    
    % LAST MODIFIED: 28/12/2020
    
    methods
        function obj = Fabric2Classes(NumMachines,linkCapacitiesAvailable)
            
            obj = obj@network_elements.Fabric(NumMachines);
            
            for ii = 1:NumMachines
                obj.machinesPorts(end+1).ingress.id = ii; % id of ingress port
                obj.machinesPorts(end).ingress.linkCapacity = utils.getLinkCapacity(linkCapacitiesAvailable); % capacity of ingress port link in bps
                obj.machinesPorts(end).egress.id = ii + obj.numFabricPorts/2; % id of egress port
                obj.machinesPorts(end).egress.linkCapacity = utils.getLinkCapacity(linkCapacitiesAvailable); % capacity of egress port link in bps
            end
            
            
        end

    end
end

