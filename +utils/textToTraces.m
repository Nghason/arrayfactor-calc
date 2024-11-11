function [fabric, coflows] = textToTraces(input_text_file)


%% Build cblock from formatted text
% c_blocks = build_cblocks_from_txt(input_file);

fid = fopen(input_text_file, 'r');
header = fscanf(fid,'%d',[1,2]);
n_coflows = header(2);
c_blocks.c_data = cell(n_coflows,1);
c_blocks.flows_block = cell(n_coflows,1);
for k = 1:n_coflows
    c_data = fscanf(fid,'%d %d %f %f %d',[5,1])';
    n_flows = c_data(end);
    c_blocks.c_data{k} = c_data;
    c_blocks.flows_block{k} = fscanf(fid,'%d %f %d %d %d',[5,n_flows])';
end

n_links = fscanf(fid,'%d',[1,1])';
n_machines = n_links/2;
portCapacity = zeros(1, n_links);
for k = 1:n_links
    fabric_data{k} = fscanf(fid,'%d %d',[2,1])';
    portCapacity(k) = fabric_data{k}(2);
end
        
fclose(fid);


%% From cblock build csv file

csv_tmp_file = 'csv_tmp.csv';
% [~,c_arrivals,c_deadlines] = fromCBlocksToCSV(csv_tmp_file,c_blocks);

time_unit = 1; % 1e-3 in simulator of Cedric

c_data = cell2mat(c_blocks.c_data);

n_coflows = size(c_data,1);
tot_n_flows = sum(c_data(:,end));

flowID    = zeros(tot_n_flows,1);
coflowID  = zeros(tot_n_flows,1);
src       = zeros(tot_n_flows,1);
dst       = zeros(tot_n_flows,1);
bandwidth = zeros(tot_n_flows,1);
weight    = zeros(tot_n_flows,1);

T = table(flowID,coflowID,src,dst,bandwidth,weight);

c_arrivals = zeros(n_coflows,1);
c_deadlines = -ones(n_coflows,1);

f_id_base = 0;

% Get coflow data
for k = 1:n_coflows
    
    % coflow data
    c_data = c_blocks.c_data{k};
    c_id = cast(c_data(1), 'int32');
    c_w  = 1;
    % arrival and deadline are expressed in number of time_units
    % c_arrivals(k) = ceil(c_data(3)/time_unit);
    % c_deadlines(k) = ceil(c_data(4)/time_unit);
    c_arrivals(k) = c_data(3);
    c_deadlines(k) = c_data(4);
    
    % flows block
    flows = c_blocks.flows_block{k};
    n_flows = cast(c_data(5), 'int32');
    for j = 1:n_flows
        f_id              = f_id_base + j;
        T.flowID(f_id)    = f_id;
        T.coflowID(f_id)  = c_id;  
        T.src(f_id)       = cast(flows(j,4),'int32');          
        T.dst(f_id)       = cast(flows(j,5),'int32');
        T.bandwidth(f_id) = flows(j,2);
        T.weight(f_id)    = c_w;
    end
    f_id_base = f_id_base + n_flows;
end


% Save table T to a temporary csv file
writetable(T, csv_tmp_file);


%% From csv file build fabric and coflows 
% utils.fromCSVToCoflows(csv_tmp_file,n_machines,0,output_file);

% Fabric setup
link_caps = [1, 1, 1];
config.NumMachines = n_machines;
config.minNumMachines = n_machines; 
config.maxNumMachines = n_machines;
config.linkCapacitiesAvailable = link_caps; % possible link capacities (should not exceed 1x3 vector)

% Coflow setup
config.architecture_type = 'csv_architecture';
config.filename = csv_tmp_file;
config.format = 0; % 0 by default
        
% Construct fabric and coflows objects
config = utils.fabric_generation(config);
config = utils.coflows_generation(config);

delete(csv_tmp_file);

% architecture_config = csvToCoflows(filename, architecture_config);

% Get main variables
fabric      = config.fabric;
coflows     = config.coflows;
n_links     = fabric.numFabricPorts;     % nb of fabric ports (ingress+egress)    
n_coflows   = length(coflows);           % nb of coflows


% load(output_file,'fabric','coflows');

%% Modify additional parameters
k = 0;
for c = coflows
    k = k+1;
    c.arrival = c_arrivals(k);
    c.deadline = c_deadlines(k);
    
    %%% ADD c.addParam fields
end


for k = 1:n_machines
   fabric.machinesPorts(k).ingress.linkCapacity = portCapacity(k);
   fabric.machinesPorts(k).egress.linkCapacity = portCapacity(n_machines+k);
end

% save(output_file, 'fabric', 'coflows');



end