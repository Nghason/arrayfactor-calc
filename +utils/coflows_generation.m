function config = coflows_generation(config)

% AUTHOR: Afaf
% LAST MODIFIED: 10/03/2021

config.coflows = [];

switch config.architecture_type

    case 'twoTypeCoflows_architectureVol'
        config.NumCoflows = randi([config.minNumCoflows config.maxNumCoflows]); % random number of coflows
        % generating coflows
        for ii = 1:config.NumCoflows
            if (ii == 1) % insure that there is at least 1 coflow of class 1
                flag_type = 0;
                config.coflows = [config.coflows network_elements.Coflow2classesVol(ii, config,flag_type)];
            elseif (ii == 2) % insure that there is at least 1 coflow of class 2
                flag_type = 1;
                config.coflows = [config.coflows network_elements.Coflow2classesVol(ii, config,flag_type)];
            else % decide of the coflow class randomly
                flag_type = 2;
                config.coflows = [config.coflows network_elements.Coflow2classesVol(ii, config,flag_type)];
            end
        end
        
    case 'csv_architecture'
        filename = config.filename;
        if (config.format == 0)
            config = utils.csvToCoflows(filename,config);
        else
            config = utils.H_csvToCoflows(filename,config);
        end
        
    case 'mapRed_architecture'
        config.NumCoflows = randi([config.minNumCoflows config.maxNumCoflows]); % random number of coflows
        % generating coflows
        for ii = 1:config.NumCoflows
            config.coflows = [config.coflows network_elements.CoflowMapRed(ii, config)];
        end

    case 'mapRed_fullFlows_architecture'
        config.NumCoflows = randi([config.minNumCoflows config.maxNumCoflows]); % random number of coflows
        % generating coflows
        for ii = 1:config.NumCoflows
            config.coflows = [config.coflows network_elements.CoflowMapRedFF(ii, config)];
        end
end

end


