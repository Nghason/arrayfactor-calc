
% 40 instances, [M,N] = [50,8000]

for lambda_i = [5] % lambda: Poisson arrival rate
    

%% Display & randomness options
OptimIn.showCoflowInfo = false; % display coflow information
OptimIn.showEpochInfo = true;   % display information in each epoch
OptimIn.showSummary = true;     % display information in each epoch
OptimIn.fixRandomness = true;   % fix the random seed 

OptimIn.setZeroArrivalTime = false; % set true to set all arrival times to zero (offline mode)
OptimIn.exportOlivierInstance = true; % set true to export instances with Olivier's format
OptimIn.setSameMachine = false; % set true to fix ingress and egress port of each flow is one single machine


%% Take input parameters

base_path = '+ResourceAllocationAlgos/+Trung/';
in_fb_raw = strcat(base_path, '+online/fb_coflows_raw.mat');
in_csv_out = strcat(base_path, '+online/fb_extract.csv');

OptimIn.csv_out = in_csv_out;
OptimIn.fb_raw = in_fb_raw;

%==============================================
% coflow generation method: {'random','facebook','from_instance'}
% generate_method  = 'random'; 
generate_method  = 'facebook';

scenario = 'offline';
% scenario = 'online';

fprintf(strcat("* using ", generate_method, " generator with ", scenario, " scenario\n"));
%==============================================


% for Trung only
% [in_n_machines, in_n_arrivals, in_max_n_flows, in_lambda] = deal(20, 1000, 30, 5); % bug dans preemption
[n_machines, n_arrivals, lambda] = deal(100, 400, lambda_i); 

flow_limit = [1, n_machines];
in_slot_duration = 0.1; 
in_seed_id = 1; 
in_n_iters = 100; 

% set arrival time of all coflows to zero in offline scenario
if strcmp(scenario, 'offline')
    n_coflows = n_arrivals;
    OptimIn.setZeroArrival = true; 
end

fprintf("* using config [%d, %d], flow_limit = [%d, %d], lambda = %.2f\n", n_machines, n_arrivals, flow_limit(1), flow_limit(2), lambda);

%==============================================
% Modify deadline range of coflow
OptimIn.deadline_min = 1;
OptimIn.deadline_max = 2;
OptimIn.ClassOneProba = 0.6; % probability of selecting coflows of Class 1
%==============================================

fprintf("* using deadline limit: [%d, %d] * CCT0 \n", OptimIn.deadline_min, OptimIn.deadline_max);


OptimIn.seed_id = in_seed_id;
if OptimIn.fixRandomness
    rng(OptimIn.seed_id); % FIX THE RANDOMNESS
    fprintf(strcat("* using seed ", num2str(in_seed_id), "\n"));
end



% batch distribution
in_batch_dist = ''; % one coflow per arrival
% in_batch_dist = 'Uniform';

OptimIn.batch_dist = in_batch_dist;
switch in_batch_dist
    % Note: in_batch_dist = '' means using no distribution, take one coflow per arrival
    case ''
        fprintf('* using no batch distribution\n');
    case 'Poisson'
        POISSON_LAMBDA = 5;
        OptimIn.pd = makedist('Poisson', 'lambda', POISSON_LAMBDA);
        fprintf('* using poisson batch with lambda = %.2f\n', POISSON_LAMBDA);
    case 'Binomial'
        [BINOM_N, BINOM_P] = deal(5, 1);
        OptimIn.pd = makedist('Binomial', 'N', BINOM_N, 'p', BINOM_P);
        fprintf('* using binomial batch with B(N,p) = B(%.2f, %.2f)\n', BINOM_N, BINOM_P);
    case 'Uniform'
        [UNI_MIN, UNI_MAX] = deal(5, 15);
        %OptimIn.coflows_per_arrival = randi([UNI_MIN, UNI_MAX], 100000, 1);
        fprintf('* using uniform batch between [min, max] = [%.2f, %.2f]\n', UNI_MIN, UNI_MAX);
        
        tmp = randi([UNI_MIN, UNI_MAX], 100000, 1);
        % Trick to get exactly sum(n_arrivals) = n_coflows 
        B = cumsum(tmp);
        n_coflows_target = n_arrivals;
        
        if B(1) >= n_coflows_target
            OptimIn.n_coflows_per_arrival = n_coflows_target;
        else
            C = B(B <= n_coflows_target);
            if C(end) == n_coflows_target
                OptimIn.n_coflows_per_arrival = tmp(1 : length(C));
            else
                OptimIn.n_coflows_per_arrival = tmp(1 : length(C)+1);
                OptimIn.n_coflows_per_arrival(end) = n_coflows_target - C(end);
            end
        end
        clearvars tmp B C;
        % Get the true n_arrivals
        n_arrivals = length(OptimIn.n_coflows_per_arrival);
        
end


%==============================================
% GENERATOR
%==============================================
for iter = 1:in_n_iters
    
    %--------------------------------------------------------------------------
    % Run generator
    %--------------------------------------------------------------------------
    if strcmp(generate_method, 'random') && strcmp(scenario, 'offline')
        [fabric, coflows] = ResourceAllocationAlgos.Trung.offline.genArchitecture(OptimIn, n_machines, n_coflows);
        
    elseif strcmp(generate_method, 'random') && strcmp(scenario, 'online')
        outputs = ResourceAllocationAlgos.Trung.online.onGenArchitecture(OptimIn, n_machines, n_arrivals, lambda);
        [fabric, coflows, n_coflows_per_arrival] = deal(outputs.fabric, outputs.coflows, outputs.n_coflows_per_arrival);
        
        
    elseif strcmp(generate_method, 'facebook') && strcmp(scenario, 'offline')
        [fabric, coflows] = ResourceAllocationAlgos.Trung.offline.offFacebookData(OptimIn, n_machines, n_coflows, flow_limit);
        
    elseif strcmp(generate_method, 'facebook') && strcmp(scenario, 'online')
        outputs = ResourceAllocationAlgos.Trung.online.onFacebookData(OptimIn, n_machines, n_arrivals, flow_limit, lambda);
        [fabric, coflows, n_coflows_per_arrival] = deal(outputs.fabric, outputs.coflows, outputs.n_coflows_per_arrival);
    end

    portCapacity = [[fabric.machinesPorts.ingress] [fabric.machinesPorts.egress]];
    portCapacity = [portCapacity.linkCapacity];
    n_links = n_machines*2;
    
    %--------------------------------------------------------------------------
    % Modify coflows to adapt to single machine scenario
    %--------------------------------------------------------------------------
    if OptimIn.setSameMachine
       for c = coflows
           tmp = zeros(n_links, c.numFlows);
           for f = c.flows
               f.destination.id = n_machines + f.source.id;
               f.destination.linkCapacity = portCapacity(f.destination.id);
               f.links = [f.source.id, f.destination.id];
               tmp(f.links, f.id) = 1;
           end
           
            % Update indicator of used links of the whole coflow
            c.indicator = tmp;
            c.addParam.indicator = (sum(c.indicator, 2) > 0);
            c.addParam.used_links = find(c.addParam.indicator == 1);
       end
    end
    
    
    arrival_times = [coflows.arrival];
    n_links    = fabric.numFabricPorts;	% nb of fabric ports (ingress+egress)    
    n_machines = n_links/2;        		% nb of machines on the fabric
    n_coflows  = length(coflows);       % nb of coflows
    n_flows    = [coflows.numFlows];    % nb of flows of each coflow
    n_flows_all = sum(n_flows);         % total nb of flows
    deadlines = [coflows.deadline];     % deadline of all coflows
    
    %--------------------------------------------------------------------------
    % Export to Olivier instance
    %--------------------------------------------------------------------------
    if OptimIn.exportOlivierInstance
        
        if strcmp(scenario, 'offline')
            config_name = strcat(generate_method, "_m", num2str(n_machines),"_c", num2str(n_coflows));
            olivier_path = strcat("ifip/offline/", config_name); 
            
        elseif strcmp(scenario, 'online')
            if strcmp(OptimIn.batch_dist, '') 
                % config_name = strcat("m", num2str(n_machines),"_c", num2str(n_coflows),"_l", num2str(lambda));
                % dir_path = strcat("+ResourceAllocationAlgos/+Trung/ifip/online/no_batch/", config_name); 
                config_name = strcat("m", num2str(n_machines),"_c", num2str(n_coflows));
                olivier_path = strcat("ifip/parallel/", config_name);    
            else
                if lambda > 1
                    config_name = strcat("m", num2str(n_machines),"_c", num2str(n_coflows),"_l", num2str(lambda));
                else
                    config_name = strcat("m", num2str(n_machines),"_c", num2str(n_coflows),"_l", num2str(10*lambda));
                end
                olivier_path = strcat("ifip/online/batch/", config_name);
            end
        end
        
        dir_path = strcat(base_path, olivier_path);
        fname = strcat(dir_path, "/", config_name, "_in", num2str(iter, "%03d"), ".txt");

        % create folder
        if ~exist(dir_path, 'dir'); mkdir(dir_path); end
        fid = fopen(fname,'w');

        % fprintf(fileID, "%d\n", epoch);
        fprintf(fid, "0\n"); % nb of slots (not important)
        fprintf(fid, "%d\n\n", n_coflows);

        for c = coflows
            fprintf(fid, "%d 1 %.3f %.5f %d\n", c.id, c.arrival, c.deadline, c.numFlows);
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
    end

 
end


end % end loop of lambda
