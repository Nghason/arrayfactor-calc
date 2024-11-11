function coflows = fromCSVToCoflows(filename, numMachines,varargin)

format = 0;
matfile_name = 'initial_coflows.mat';
switch nargin
    case 3
        format = varargin{1};
    case 4
        format = varargin{1};
        matfile_name = varargin{2};
end
% if nargin >= 3
%    H = varargin{1};
% end

architecture_config.architecture_type = 'csv_architecture';
architecture_config.NumMachines = numMachines;
architecture_config.filename = filename;
architecture_config.format = format;

architecture_config = utils.fabric_generation(architecture_config);
architecture_config = utils.coflows_generation(architecture_config);

%% adapt to greedy allocation online 
coflows = architecture_config.coflows;
%%

% if ~strcmp(architecture_config.architecture_type, 'skeleton_architecture')
%     
%     coflowStruct = struct('n_flows',{}, 'f_vol',{}, 'indicator',{});
%     
%     for ii = 1:length(architecture_config.coflows)
%         coflowStruct(end+1).n_flows = architecture_config.coflows(ii).numFlows;
%         for jj = 1:architecture_config.coflows(ii).numFlows
%             coflowStruct(end).f_vol = [coflowStruct(end).f_vol architecture_config.coflows(ii).flows(jj).volume];
%         end
%         coflowStruct(end).indicator = architecture_config.coflows(ii).indicator;
%     end
%     
%     % Similar to Cedric structure
%     coflowStruct2 = struct('n_coflows',{}, 'n_flows',{}, 'f_vol',{}, 'indicator',{}, 'n_links',{}, 'links_capacities',{});
%     
%     
%     coflowStruct2(end+1).n_coflows = length(coflowStruct);
%     coflowStruct2(end).n_flows = [coflowStruct.n_flows];
%     for ii = 1:length(architecture_config.coflows)
%         coflowStruct2(end).f_vol{ii} = [coflowStruct(ii).f_vol];
%         coflowStruct2(end).indicator{ii} = [coflowStruct(ii).indicator];
%     end
%     
%     coflowStruct2(end).n_links = architecture_config.fabric.numFabricPorts;
%     
%     for jj =1:length(architecture_config.fabric.machinesPorts)
%         coflowStruct2(end).links_capacities = [coflowStruct2(end).links_capacities architecture_config.fabric.machinesPorts(jj).ingress.linkCapacity];
%     end
%     
%     for jj =1:length(architecture_config.fabric.machinesPorts)
%         coflowStruct2(end).links_capacities = [coflowStruct2(end).links_capacities architecture_config.fabric.machinesPorts(jj).egress.linkCapacity];
%     end
%     
% %     architecture_config.coflowStruct = coflowStruct;
% %     architecture_config.coflowStruct2 = coflowStruct2;
% end
% 
% coflows = architecture_config.coflows;
% fabric = architecture_config.fabric;
% 
% save generatedCoflowsFabric coflows fabric coflowStruct coflowStruct2

save(matfile_name, '-struct', 'architecture_config', 'fabric', 'coflows', '-v7.3');

end