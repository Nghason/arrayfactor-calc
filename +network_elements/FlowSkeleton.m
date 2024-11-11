classdef FlowSkeleton < network_elements.Flow

    % LAST MODIFIED: 17/12/2020
    
    methods
        function obj = FlowSkeleton(flowID,flowVolume,source,destination,coflowID)
            
            obj = obj@network_elements.Flow(flowID, coflowID);
            
            obj.volume = flowVolume;
            obj.source = source;
            obj.destination = destination;
            obj.links = [obj.source.id obj.destination.id];
            obj.remainingVolume = obj.volume;
            obj.volume_initial = obj.volume;
        end
    end
end

