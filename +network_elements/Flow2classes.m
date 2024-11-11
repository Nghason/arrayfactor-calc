classdef Flow2classes < network_elements.Flow

    % LAST MODIFIED: 01/01/2021
    
    methods
        function obj = Flow2classes(id,avgFlowVolume, coflowID, standardDivVolume, source, dest)
            obj = obj@network_elements.Flow(id, coflowID);
            obj.volume = abs(normrnd(avgFlowVolume,standardDivVolume)); % generate flow size randomly following a normal dist 
            
            obj.source = source;
            obj.destination = dest;
            
            
%             obj.source = fabricMachinesPorts(randi([1 numFabricMachinesPorts],1,1)).ingress; % assign a source randomly to each flow
%             
%             p = randi([1 numFabricMachinesPorts],1,1);
%             while (fabricMachinesPorts(p).ingress.id == obj.source.id)
%                 p = randi([1 numFabricMachinesPorts],1,1);
%             end
%             obj.destination = fabricMachinesPorts(p).egress; % assign a destination randomly to each flow (the destination should not be the same machine as the source)
            
            
            obj.links = [obj.source.id obj.destination.id];
            obj.remainingVolume = obj.volume;
            obj.volume_initial = obj.volume;

        end
        

    end
end

