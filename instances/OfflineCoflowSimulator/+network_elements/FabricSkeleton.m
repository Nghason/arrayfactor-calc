classdef FabricSkeleton < network_elements.Fabric

    % LAST MODIFIED: 17/12/2020    
    
    methods
        function obj = FabricSkeleton(NumMachines)
            
            obj = obj@network_elements.Fabric(NumMachines);
            
            obj.machinesPorts = struct('ingress',cell(1,NumMachines),'egress',cell(1,NumMachines));
        end
        
        function setIngress(obj,id_machine,cap) % set ingress link of a given machine
            obj.machinesPorts(id_machine).ingress.id = id_machine;
            obj.machinesPorts(id_machine).ingress.linkCapacity = cap;
        end
        
        function setEgress(obj,id_machine,cap) % set egress link of a given machine
            obj.machinesPorts(id_machine).egress.id = id_machine + obj.numFabricPorts/2;
            obj.machinesPorts(id_machine).egress.linkCapacity = cap;
        end
    end
end

