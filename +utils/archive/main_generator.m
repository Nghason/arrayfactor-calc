function  [architecture_config, coflowStruct, coflowStruct2] = main_generator(architecture_config)

    % AUTHOR: Afaf
    % LAST MODIFIED: 10/12/2020


    architecture_config = utils.fabric_generation(architecture_config);
    architecture_config = utils.coflows_generation(architecture_config);
   
    
    
    %% Preparing outputs
    
    if ~strcmp(architecture_config.architecture_type, 'skeleton_architecture')
    
        coflowStruct = struct('n_flows',{}, 'f_vol',{}, 'indicator',{});

        for ii = 1:length(architecture_config.coflows)
            coflowStruct(end+1).n_flows = architecture_config.coflows(ii).numFlows;
            for jj = 1:architecture_config.coflows(ii).numFlows
                coflowStruct(end).f_vol = [coflowStruct(end).f_vol architecture_config.coflows(ii).flows(jj).volume];
            end
            coflowStruct(end).indicator = architecture_config.coflows(ii).indicator;
        end

        % Similar to C??dric structure
        coflowStruct2 = struct('n_coflows',{}, 'n_flows',{}, 'f_vol',{}, 'indicator',{}, 'n_links',{}, 'links_capacities',{});


        coflowStruct2(end+1).n_coflows = length(coflowStruct);
        coflowStruct2(end).n_flows = [coflowStruct.n_flows];
        for ii = 1:length(architecture_config.coflows)
            coflowStruct2(end).f_vol{ii} = [coflowStruct(ii).f_vol];
            coflowStruct2(end).indicator{ii} = [coflowStruct(ii).indicator];
        end

        coflowStruct2(end).n_links = architecture_config.fabric.numFabricPorts;

        for jj =1:length(architecture_config.fabric.machinesPorts)
            coflowStruct2(end).links_capacities = [coflowStruct2(end).links_capacities architecture_config.fabric.machinesPorts(jj).ingress.linkCapacity];
        end

        for jj =1:length(architecture_config.fabric.machinesPorts)
            coflowStruct2(end).links_capacities = [coflowStruct2(end).links_capacities architecture_config.fabric.machinesPorts(jj).egress.linkCapacity];
        end
        
        architecture_config.coflowStruct = coflowStruct;
        architecture_config.coflowStruct2 = coflowStruct2;
    end
end

