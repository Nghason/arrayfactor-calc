%{
command:
nohup /usr/local/matlab/bin/matlab -nodisplay -nosplash -nodesktop -r \
"run('/home/qtluu/Downloads/workspace/coflow_2022_02_07/wrapperOlivier.m');exit;" > \
/home/qtluu/Downloads/workspace/coflow_2022_02_07/+ResourceAllocationAlgos/+Trung/ifip/avignon.out &

nohup /usr/local/matlab/bin/matlab -nodisplay -nosplash -nodesktop -r \
"run('/home/qtluu/Downloads/workspace/coflow_2020_02_07/wrapperOlivier.m');exit;" > \
/home/qtluu/Downloads/workspace/coflow_2020_02_07/avignon.out &
%}


% clc
clear all

% Display current time instant
fprintf(strcat(datestr(datetime(now,'ConvertFrom','datenum')), "\n"));


%% Select methods to run

% mod_names = {'maxCdsOptim'; 'maxCdsRelax'; 'cs_mha';
%              'dcoflow_min_sum_negative'; 'dcoflow_min_sum_congested';
%              'wDcoflow_min_link'; 'wDcoflow_min_sum_negative'; 'wDcoflow_min_sum_congested';
%              'sincronia'; 'varys'};     

mod_names = {'maxCdsOptim'; 'maxCdsRelax'; 'cs_mha'; 
             'dcoflow_min_link'; 'dcoflow_min_sum_negative'; 'dcoflow_min_sum_congested';
             'sincronia'; 'varys'};   
         
% mod_names = {'cs_mha'; 'heu_v2_min_sum_negative'; 'heu_v2_min_sum_congested';
%              'sincronia'; 'varys'};  

% mod_names = {'dcoflow_min_sum_negative'; 'wDcoflow_min_sum_negative'}; 
% mod_names = {'dcoflow_min_sum_negative'};



%% Global options

base_path = '';

inputs.fixRandomness = true; % fix the random seed 
inputs.seed_id = 1;          % random seed id

inputs.exportNetworkAsText = true;     % export network instance as txt
inputs.buildFinalTable = true;          % export the summarized result table
inputs.T_final_dir = 'results/';        % directory to save the summarized table

% Necessary files for facebook trace generator
inputs.csv_out = strcat(base_path, '+utils/fb_extract.csv');                                                                      
inputs.fb_raw = strcat(base_path, '+utils/fb_coflows_raw.mat'); 

% Select LP solver
inputs.solver = 'gurobi'; % options: {'cplexmilp', 'intlinprog', 'gurobi'}

% Solver paths
% addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio1210\cplex\matlab\x64_win64'); % cplex path
% addpath('/home/qtluu/Downloads/workspace/install/gurobi950/linux64/matlab/');
addpath('C:\gurobi950\win64\matlab\'); % gurobi path



%% Input parameters
trace_type  = 'random';       % options: {'random','facebook'}

%======================================
n_iters = 100;                   % nb of iterations
n_machines = 15;                 % nb of machines
n_coflows  = 60;                 % nb of coflows

inputs.flow_limit = [1, n_machines]; % range of nb of flow per coflow (for facebook traces only)
inputs.Lambda = 0;              % poisson arrival rates (for online simulation only)

inputs.n_priority_classes = 2;
inputs.priority_probas = [0.3, 0.7]; 
inputs.priority_weights = [2, 1]; 
%======================================


%======================================
% Modify deadline range [a,b] of coflow
inputs.deadline_min = 1; % a
inputs.deadline_max = 2; % b
inputs.ClassOneProba = 0.6; % prob of selecting Class-1 coflows (default: 0.6)
%======================================


fprintf('* using generator: %s\n', trace_type);
fprintf('* using LP solver: %s\n', inputs.solver);
fprintf("* using deadline limit: [%d, %d] * CCT0 \n", inputs.deadline_min, inputs.deadline_max);
fprintf("* using seed: %d\n", inputs.seed_id);
fprintf("* using config [%d, %d],", n_machines, n_coflows);
fprintf(" flow_limit = [%d, %d]\n", inputs.flow_limit(1), inputs.flow_limit(2));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if inputs.fixRandomness
    rng(inputs.seed_id); % disp(""FIX THE RANDOMNESS")
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Initialization
  
% inputs.zn_type = 'N';                  % either 'N' or 'S' (LP and LP Approx)
% inputs.zn_type = 'S';
inputs.options = [];
if 	strcmp(inputs.solver, 'intlinprog') 
    inputs.options = optimoptions('intlinprog');
    inputs.options.Display = 'off';
    % inputs.options.MaxTime = 10;      % max compute time         
    % inputs.options.MaxNodes = 1e6;  % max compute nodes
end

n_methods = length(mod_names);        
inputs.n_methods = n_methods;
inputs.mod_names = mod_names;


for iter = 1:n_iters 
    fprintf('%d-', iter);
    if mod(iter, 50) == 0; fprintf('\n'); end
    
    inputs.iter = iter;

	% generation of colfows and fabric
    inputs.scenario = 'offline'; 
    network = utils.generateTraces(trace_type, inputs, n_machines, n_coflows); % '0' indicates no Lambda used in offline
    
    fabric = network.fabric;
    coflows = network.coflows;
 
    n_links    = fabric.numFabricPorts;     % nb of fabric ports (ingress+egress)    
    n_flows    = [coflows.numFlows];        % nb of flows of each coflow
    n_flows_all = sum(n_flows);             % total nb of flows
    deadlines = [coflows.deadline];         % deadline of all coflows

    %======================================================
    % OUTPUT OLIVIER FORMAT: SAVE SIMULATION INSTANCE
    %======================================================
    
    if inputs.exportNetworkAsText 
            
        config_name = strcat('off_', trace_type, "_m", num2str(n_machines),"_c", num2str(n_coflows));
        olivier_path = strcat("instances/", config_name);
        dir_path = strcat(base_path, olivier_path);
        fname = strcat(dir_path, "/", config_name, "_in", num2str(iter, "%03d"), ".txt");
        % create folder
        if ~exist(dir_path, 'dir'); mkdir(dir_path); end
        fid = fopen(fname,'w');
        
        fprintf(fid, "0\n"); % nb of slots (not important)
        fprintf(fid, "%d\n\n", n_coflows);

        for c = coflows
            fprintf(fid, "%d %d %.3f %.5f %d\n", c.id, c.weight, c.arrival, c.deadline, c.numFlows);
            %fprintf("coflow %d\n", c.id);
            for f = c.flows
               fprintf(fid, "%d%d %.5f 2 %d %d\n", c.id, f.id, f.volume, f.links(1), f.links(2)); 
            end
            fprintf(fid, "\n");
        end

        portCapacity = [[fabric.machinesPorts.ingress] [fabric.machinesPorts.egress]];
        portCapacity = [portCapacity.linkCapacity];

        fprintf(fid, "%d\n", n_links);
        for i = 1:n_links
           fprintf(fid, "%d %d\n", i, portCapacity(i)); 
        end

        fclose(fid);
    
    end % end exportNetworkAsText
end % end for iter








