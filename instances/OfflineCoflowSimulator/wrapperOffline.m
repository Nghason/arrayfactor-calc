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


%% Select methods to run
% No weight: dcoflow, dcoflow_mh
% Weighted: dcoflow_dp, wdcoflow


% mod_names = {'dcoflow_min_all'; 'dcoflow_min_negative'; 'dcoflow_min_congested';...
%              'w_dcoflow_min_all'; 'w_dcoflow_min_negative'; 'w_dcoflow_min_congested';...
%              'dcoflow_dp'; 'dcoflow_mh'; 'dcoflow_br'; 'sincronia'; 'cs_mha'; 'cs_dp';
%              'varys'; 'CdsOptim'; 'CdsRelax'; 'w_CdsOptim'; 'w_CdsRelax'};     
          
% mod_names = {'dcoflow_min_all';'sincronia';'cs_mha';'lpovls';'varys';'CdsOptim';'CdsRelax'}; 
          
% mod_names = {'cs_mha'; 'dcoflow_min_negative'; 'dcoflow_min_congested';
%              'sincronia'; 'varys'};  
% mod_names = {'sincronia'; 'lpovls''CdsOptim'; 'CdsRelax'}; 
% mod_names = {'sincronia'; 'lpovls'}; 
% mod_names = {'dcoflow_min_negative'; 'cs_mha'}; 
% mod_names = {'dcoflow_dp'};
% mod_names = {'dcoflow_mh'};
mod_names = {'cs_mha'};
% mod_names = {'cs_dp'};
% mod_names = {'dcoflow_min_negative';'dcoflow_mh';'dcoflow_br'};


fprintf(strcat(datestr(datetime(now,'ConvertFrom','datenum')), "\n")); % current time instant



%% GLOBAL PARAMETERS

global params;

params.fix_randomness = true; % fix the random seed 
params.seed_id = 1;           % random seed id

params.debug_mode = false;                   % export detailed simulation info to log.txt
params.export_detailed_results = true;     % export results for comparison
params.export_network_as_text = false;       % export network instance as txt
params.build_recap_table      = true;       % export the summarized result table
params.export_recap_table     = true;       % export the summarized result table

% Paths
params.base_path = '';
params.debug_path = 'logs/log.txt';
params.instance_path = 'instances/weight/';    % path to save instances
params.table_path = 'results/';   	% directory to save the summarized table

% Necessary files for facebook trace generator 
params.csv_out = strcat(params.base_path, '+utils/fb_extract.csv');                                                                      
params.fb_raw  = strcat(params.base_path, '+utils/fb_coflows_raw.mat'); 

% Select LP solver
params.solver = 'intlinprog'; % options: {'cplexmilp', 'intlinprog', 'gurobi'}
params.export_lp = true;  % export LP model to file 

% Solver paths
% addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio1210\cplex\matlab\x64_win64'); % cplex path
addpath('C:\Program Files\IBM\ILOG\CPLEX_Enterprise_Server1210\CPLEX_Studio\cplex\matlab\x64_win64'); % cplex path
% addpath('/home/qtluu/Downloads/workspace/install/gurobi950/linux64/matlab/');
addpath('C:\gurobi950\win64\matlab\'); % gurobi path



%% INPUT PARAMETERS
trace_type  = 'random';     % options: {'random','facebook'}
n_iters     = 1;          % nb of iterations
n_machines  = 10;            % nb of machines
n_coflows   = 10;           % nb of coflows

inputs.FLOW_RANGE = [1, n_machines];    % range of nb of flow per coflow (for facebook traces only)
inputs.ARRIVAL_RATE = 0;                % poisson arrival rates (for online simulation only)

% Priority classes and weights
inputs.PRIO_PROBAS      = [0.7, 0.3]; 
inputs.PRIO_WEIGHTS     = [1, 1]; 
inputs.PRIO_CLASS_NUM   = length(inputs.PRIO_PROBAS);

if length(unique(inputs.PRIO_WEIGHTS)) == 1
    params.instance_path = 'instances/noweight/';
end

% Modify deadline range [a,b] of coflow
inputs.DEADLINE_RANGE  = [1, 2]; % [a,b]

% Probability of selecting Class-1 coflows (default: 0.6)
inputs.CLASS_ONE_PROBA = 0.6; 




%% ADDITIONAL PARAMETERS 
% runPredictCct: predict accept/reject based on CCT in bottleneck)
% runAllocation: run allocation algorithm after having order

for mm = 1:length(mod_names)
    mod_name = mod_names{mm};
    switch mod_name
        case {'dcoflow_min_all','dcoflow_min_negative','dcoflow_min_congested',...
              'w_dcoflow_min_all','w_dcoflow_min_negative','w_dcoflow_min_congested',...
              'dcoflow_dp','dcoflow_mh','dcoflow_br'}
            inputs.mod_option.(mod_name).runPredictCct = true;
            inputs.mod_option.(mod_name).runAllocation = true;
        case {'sincronia','cs_mha','cs_dp','lpovls'}
            inputs.mod_option.(mod_name).runPredictCct = false;
            inputs.mod_option.(mod_name).runAllocation = true;
        case {'varys','CdsOptim','CdsRelax','w_CdsOptim','w_CdsRelax'} 
            % varys and maxCds don't need Backtrack and Alloc
            inputs.mod_option.(mod_name).runPredictCct = false;
            inputs.mod_option.(mod_name).runAllocation = false;
    end
end



%------------------------------------------------------------------------------------------------------
if params.debug_mode
    fid = fopen(params.debug_path, 'w');
    fprintf(fid, strcat(datestr(datetime(now,'ConvertFrom','datenum')), "\n")); % current time instant
    fprintf(fid, '* use generator: %s\n',                   trace_type);
    fprintf(fid, '* use LP solver: %s\n',                   params.solver);
    fprintf(fid, '* use deadline range: [%d, %d]*CCT0 \n',  inputs.DEADLINE_RANGE(1), inputs.DEADLINE_RANGE(2));
    fprintf(fid, '* use seed: %d\n',                        params.seed_id);
    fprintf(fid, '* use config [%d, %d],',                  n_machines, n_coflows);
    fprintf(fid, ' with flow_range = [%d, %d]\n',           inputs.FLOW_RANGE(1), inputs.FLOW_RANGE(2));
    %fclose(fid); 
    fclose('all');
end
%------------------------------------------------------------------------------------------------------



% Call main file
main_offline





%%
% fid = fopen('logs/compare.txt', 'w');
% for iter = 1:n_iters
%     fprintf(fid, '%d\t%.4f\n', iter, Final.variant.(mod_names{1}).wrac(iter)); 
% end

% Compare two files
% visdiff('logs/log_old.txt', 'logs/log.txt')






