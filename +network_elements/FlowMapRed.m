classdef FlowMapRed < network_elements.Flow

    % LAST MODIFIED: 26/01/2021
    
    methods
        function obj = FlowMapRed(id,avgFlowVolume, coflowID, source, dest)
            obj = obj@network_elements.Flow(id, coflowID);
            
            obj.volume = exprnd(avgFlowVolume);
            
            obj.source = source;
            obj.destination = dest;
            
            obj.links = [obj.source.id obj.destination.id];
            obj.remainingVolume = obj.volume;
            obj.volume_initial = obj.volume;

        end
        

    end
end


