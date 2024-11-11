function [config] = fabric_generation(config)

    % AUTHOR: Afaf + Cedric
    % LAST MODIFIED: 10/03/2021
    
    switch config.architecture_type
        
%        case {'twoTypeCoflows_architecture', 'twoTypeCoflows_architectureVol'}
        case 'twoTypeCoflows_architectureVol'
            
        config.NumMachines = randi([config.minNumMachines config.maxNumMachines]); % Random # of machines
        config.fabric = network_elements.Fabric2Classes(config.NumMachines, config.linkCapacitiesAvailable); % generating the fabric
        
        case {'csv_architecture', 'skeleton_architecture'}
            
        config.fabric = network_elements.FabricSkeleton(config.NumMachines);
        
        case {'mapRed_architecture', 'mapRed_fullFlows_architecture'}
            
            % Generating the fabric given the input parameters.
            % For mapRed architectures, the fabric object is not
            % diferent from the towTypeCoflows_architecture
            config.fabric = network_elements.Fabric2Classes(config.NumMachines, config.linkCapacitiesAvailable);
    end

end

