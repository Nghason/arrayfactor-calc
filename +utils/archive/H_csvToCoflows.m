function architecture_config =  H_csvToCoflows(filename,architecture_config)
% Generates an instance of coflows from a H CSV file
% Returns the array of generated coflow objects
% Parameters:
% - filename: the name of the csv file
% - n_machines: number of machine in the fabric
% Cedric Richier, LIA, 2020
% (c) LIA, 2020

T = readtable(filename);

%% For now, generating a n_machines x n_machines fabric with link capacities set to one
n_machines = architecture_config.NumMachines;

for m = 1:n_machines
    architecture_config.fabric.setIngress(m,1);
    architecture_config.fabric.setEgress(m,1);
end

coflows_ids = unique(T.coflowID)+1;
n_coflows = length(coflows_ids);

for k = 1:n_coflows
    c = network_elements.CoflowSkeleton(coflows_ids(k),architecture_config);
    architecture_config.coflows = [architecture_config.coflows c];
end

tmp_coflows = architecture_config.coflows;

id_coflow = T.coflowID(1)+1;
idx = [tmp_coflows.id] == id_coflow;
c = tmp_coflows(idx);
id_flow = c.numFlows+1;

for n_line = 1:size(T,1)
    if T.coflowID(n_line) ~= id_coflow-1
        id_coflow = T.coflowID(n_line)+1;
        idx = [tmp_coflows.id] == id_coflow;
        c = tmp_coflows(idx);
        id_flow = c.numFlows+1;
    end
    c.addFlow(id_flow,T.bandwidth(n_line),architecture_config.fabric.machinesPorts(T.src(n_line)).ingress,...
        architecture_config.fabric.machinesPorts(T.dst(n_line)).egress,n_machines);
    id_flow = id_flow+1;
end

for c = architecture_config.coflows
    c.update(n_machines);
end

% Some printings
%displayCoflows(coflows);

% outputs for other scripts
coflowStruct = struct('n_flows',{}, 'f_vol',{}, 'indicator',{});

for ii = 1:length(architecture_config.coflows)
    coflowStruct(end+1).n_flows = architecture_config.coflows(ii).numFlows;
    for jj = 1:architecture_config.coflows(ii).numFlows
        coflowStruct(end).f_vol = [coflowStruct(end).f_vol architecture_config.coflows(ii).flows(jj).volume];
    end
    coflowStruct(end).indicator = architecture_config.coflows(ii).indicator;
end

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

% Save outputs in a .mat file
%fprintf('* Saving outputs in ''generatedCoflowsFabric.mat'' \n')
%save(out, 'coflows', 'fabric', 'coflowStruct', 'coflowStruct2');

end