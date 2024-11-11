

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if params.fix_randomness
    rng(params.seed_id); % disp(""FIX THE RANDOMNESS")
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Initialization
% inputs.zn_type = 'N';  % either 'N' or 'S' (LP and LP Approx)
inputs.options = [];
if 	strcmp(params.solver, 'intlinprog') 
    inputs.options = optimoptions('intlinprog');
    inputs.options.Display = 'off';
    % inputs.options.MaxTime = 10;    % max compute time         
    % inputs.options.MaxNodes = 1e6;  % max compute nodes
end

n_methods = length(mod_names);        
inputs.n_methods = n_methods;
inputs.mod_names = mod_names;




t1 = cputime; % start

for iter = 1:n_iters 
    fprintf('%d-', iter);
    if mod(iter, 50) == 0; fprintf('\n'); end
    if n_iters < 50 && iter == n_iters; fprintf('\n'); end
    
    inputs.iter = iter;

	% generation of colfows and fabric
    inputs.scenario_onoff = 'offline'; 
    network = utils.generateTraces(trace_type, inputs, n_machines, n_coflows); % '0' indicates no Lambda used in offline
    % saved_networks{iter} = network;
    fabric = network.fabric;
    coflows = network.coflows;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % if false
    %     WEIGHTS = [1,1,10,1,1,3,1,1];
    %     DEADLINES = [6,8,9,11,20,25,28,35];
    %     VOLS = [4,1,6,3,6,8,7,10];
    %     for c = coflows
    %         c.weight = WEIGHTS(c.id);
    %         c.deadline = DEADLINES(c.id);
    %         c.flows(1).volume = VOLS(c.id);
    %         c.flows(1).volume_initial = VOLS(c.id);
    %     end
    %     clearvars WEIGHTS DEADLINES VOLS
    % end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    n_links    = fabric.numFabricPorts;     % nb of fabric ports (ingress+egress)    
    n_flows    = [coflows.numFlows];        % nb of flows of each coflow
    n_flows_all = sum(n_flows);             % total nb of flows
    deadlines = [coflows.deadline];         % deadline of all coflows
        
    %======================================================
    % OUTPUT OLIVIER FORMAT: SAVE SIMULATION INSTANCE
    %======================================================
    if params.export_network_as_text 
            
        config_name = strcat('off_', trace_type, "_m", num2str(n_machines),"_c", num2str(n_coflows));
        dir_path = strcat(params.base_path, params.instance_path, config_name);
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

        %fclose(fid);
        fclose('all');
    
    end % end EXPORT_NETWORK_AS_TEXT
    
    %================================================
    % RUN SCHEDULING ALGORITHMS
    %================================================
    mod_names = inputs.mod_names;

    for mm = 1:length(mod_names)

        mod_name = mod_names{mm}; 
        inputs.mod_name = mod_name;
        
        % fprintf('%s\n', mod_name);
        t_mm_1 = cputime; % start time

        %=================================================================================
        % Run algorithms
        %=================================================================================
        switch mod_name
            
            % Optimization from [Tseng2019]: CDS-LP (binary zn), CDS-LPA (continuous zn), direct allocation 
            case {'CdsOptim','CdsRelax','w_CdsOptim','w_CdsRelax'} 
                temp = algorithms.offline.optim_cds(inputs, fabric, coflows, mod_name); 
                temp.method = mod_name;
                out{mm} = temp;
                clearvars tseng_options;
                
            % CS-MHA [Luo2016]: Ordering
            case {'cs_mha', 'cs_dp'}    
                out{mm} = algorithms.offline.csmha(mod_name, fabric, coflows);
                out{mm}.ini_order = out{mm}.order;
                out{mm}.method = mod_name; % save method name
        
            % DCoflow: all variants
            case {'dcoflow_min_all','dcoflow_min_negative','dcoflow_min_congested',...
                  'w_dcoflow_min_all','w_dcoflow_min_negative','w_dcoflow_min_congested',...
                  'dcoflow_dp','dcoflow_mh'}  
                out{mm} = algorithms.offline.dcoflow(mod_name, fabric, coflows);
                 
            % Sincronia [Agarwal2018]: Ordering
            case 'sincronia'
                out{mm}.ini_order = algorithms.offline.sincronia(fabric, coflows);
                out{mm}.order = out{mm}.ini_order;
                out{mm}.method = mod_name; % save method name
                
            % LP-OV-LS [Shafiee2018]: Ordering
            case 'lpovls'
                out{mm}.ini_order = algorithms.offline.lpovls(fabric, coflows);
                out{mm}.order = out{mm}.ini_order;
                out{mm}.method = mod_name; % save method name

            % Varys [Chowdhury2014]: Ordering + Allocation
            case 'varys'   
                temp = algorithms.offline.varys(fabric, coflows);    
                temp.zn  = (temp.cct <= deadlines); % temp.cct contains CCTs of all coflows
                temp.ccts = temp.cct(temp.zn); 
                temp.nac = sum(temp.zn);
                temp.rejects = find(temp.zn < 0.1); % for safety
                temp.accepts = find(temp.zn > 0.9); % for safety
                temp.method = 'varys';
                out{mm} = temp;
                
        end

        % Elapsed time for scheduling
        if inputs.mod_option.(mod_name).runAllocation 
            out{mm}.t_sch = cputime - t_mm_1; 
        else
            out{mm}.t_sch = 0;
        end

        %=================================================================================
        % GREEDY ALLOCATION
        %=================================================================================
        if inputs.mod_option.(mod_name).runAllocation
            greedy_out = algorithms.offline.greedy_allocation(fabric, coflows, out{mm}.order);
            greedy_out_ccts = greedy_out.ccts(out{mm}.ini_order);
            tmp = greedy_out_ccts;
            tmp(tmp == 0) = inf;
            
            out{mm}.zn                = (tmp <= deadlines(out{mm}.ini_order)); % attention: ccts contains only cct of accepted coflows
            out{mm}.accepts           = out{mm}.ini_order(out{mm}.zn); 
            out{mm}.rejects           = out{mm}.ini_order(~out{mm}.zn);
            out{mm}.ccts_after_greedy = greedy_out_ccts;
            out{mm}.ccts              = greedy_out_ccts(out{mm}.zn); % greedy_out.ccts
            out{mm}.nac               = sum(out{mm}.zn);

        end
        
        %=================================================================================
        % Update status of each coflow
        %=================================================================================
        for ii = out{mm}.accepts
            coflows(ii).addParam.status = "accepted";
        end
        for ii = out{mm}.rejects
            coflows(ii).addParam.status = "rejected";
        end
        
        %=================================================================================
        % Normalized CCT of accepted coflows (divided by CCT0)
        %=================================================================================     
		tmp = zeros(1, n_coflows);
        for ii = 1:out{mm}.nac % out{mm}.accepts
            ci = out{mm}.accepts(ii);
            c = coflows(ci);
            % cct and normalized cct of coflow c
            c.addParam.cct 		= out{mm}.ccts(ii);
            c.addParam.cct_norm = round(c.addParam.cct/c.addParam.CCT0, 4);
            tmp(ii) = c.addParam.cct_norm;  
        end
        
        out{mm}.cct_norm_vec = tmp; % vector of normalized ccts
		
		% average cct and normalized cct of accepted coflows
        out{mm}.cct = sum(out{mm}.ccts)/out{mm}.nac;
        out{mm}.cct_norm = sum(out{mm}.cct_norm_vec)/out{mm}.nac; 
        
        % elapsed time
        t_mm_2 = cputime;       % end time
        t_mm = t_mm_2 - t_mm_1; % computing time of method 'mm'
        out{mm}.time = t_mm;

    end % end loop of mm

    %=================================================================================
    % SUMMARIZED RESULTS
    %=================================================================================
    for mm = 1:n_methods
        t_mm_1 = cputime; %clock; % start 
        mod_name = mod_names{mm};
        Final.variant.(mod_name).out{iter} = out{mm};
        
        % n_accepts and r_accepts
        Final.variant.(mod_name).nac(iter) = out{mm}.nac;
        Final.variant.(mod_name).rac(iter) = out{mm}.nac/n_coflows;
        
        % weighted acceptance rate
        wrac = 0;
        for ii = 1:out{mm}.nac % out{mm}.accepts
            ci = out{mm}.accepts(ii);
            c = coflows(ci);
			% cct and normalized cct of coflow c
			wrac = wrac + c.weight;
        end
        wrac = wrac/sum([coflows.weight]);
        Final.variant.(mod_name).wrac(iter) = wrac;
        clearvars wrac
        
        % average cct and average normalized cct (per accepted coflow)
        Final.variant.(mod_name).cct(iter) 		= out{mm}.cct;
        Final.variant.(mod_name).cct_norm(iter) = out{mm}.cct_norm;
        
        % prediction error of other methods than sincronia/maxCdsOptim/maxCdsRelax
        if inputs.mod_option.(mod_name).runPredictCct  
            Final.variant.(mod_name).pred_err(iter) = (out{mm}.pred_nac - out{mm}.nac)/(out{mm}.pred_nac); 
        end

        % time
        Final.variant.(mod_name).t_sch(iter) = out{mm}.t_sch; % time consumed for scheduling
        Final.variant.(mod_name).time(iter) = out{mm}.time;   % computing time of mm   
        
        %======================================================
        % OUTPUT OLIVIER FORMAT: SAVE RESULT OF EACH INSTANCE
        %======================================================
        if params.export_detailed_results         
            mm_final = Final.variant.(mod_name);
            mm_out = Final.variant.(mod_name).out{iter};
            
            instance_name = strcat('off_', trace_type, "_m", num2str(n_machines), "_c", num2str(n_coflows), "_in", num2str(iter, '%03d'));
            if iter == 1; fid = fopen('trung.txt', 'w'); else; fid = fopen('trung.txt', 'a+'); end
            %fid = fopen('details.txt', 'a+');
            
            switch mod_name
                case {'dcoflow_min_all','dcoflow_min_negative','dcoflow_min_congested',...
                      'w_dcoflow_min_all','w_dcoflow_min_negative','w_dcoflow_min_congested',...
                      'dcoflow_dp','dcoflow_mh'}
                    % Printing
                    fprintf(fid, strcat('Processing file TEST/', instance_name, '.txt\n'));
                    fprintf(fid, 'Nombre de coflows acceptes: %d\n', mm_final.nac(iter));
                    fprintf(fid, strcat('sigma=[',  strjoin(string(mm_out.order), ', '), ']\n'));
                    fprintf(fid, 'Weighted CAR: %.3f\n', mm_final.wrac(iter));
                    
                    % fprintf(fid, strcat('sigma*=[', strjoin(string(tmp.ini_order_star), ', '), ']\n'));
                    % fprintf(fid, strcat('sigma=[',  strjoin(string(tmp.order), ', '), ']\n'));
                    % 
                    % tmp2 = {};
                    % for iii = 1:length(tmp.ccts_after_greedy) % ccts of all coflows in order (also rejected ones after greedy)
                    %     tmp2{iii} = num2str(round(tmp.ccts_after_greedy(iii), 2), '%.2f'); 
                    % end
                    % fprintf(fid, strcat('ccts=[', strjoin(tmp2, ', '), ']\n'));
                case {'sincronia', 'cs_mha'}
                    fprintf(fid, strcat('sigma=[',  strjoin(string(mm_out.order), ', '), ']\n'));
                    tmp2 = {};
                    for iii = 1:length(mm_out.ccts_after_greedy) % ccts of all coflows in order (also rejected ones after greedy)
                        tmp2{iii} = num2str(round(mm_out.ccts_after_greedy(iii), 2), '%.2f'); 
                    end
                    fprintf(fid, strcat('ccts=[', strjoin(tmp2, ', '), ']\n'));
                case {'varys','CdsOptim','CdsRelax','w_CdsOptim','w_CdsRelax'}
                    fprintf(fid, strcat('accepts=[',  strjoin(string(sort(mm_out.accepts)), ', '), ']\n'));
                    fprintf(fid, '%d\n', mm_out.nac); 
            end
            
            %fprintf(fid, '\n');
            fclose('all');

        end  % end olivier file
        
    end % end loop of mm
    
end % end loop of iter

t2 = cputime;      % finish
t_total = t2 - t1; % computing time


%% Post-processing

% Mean and std_dev of nac, rac, cct, and time
Final.avg_nac = zeros(n_methods, 2); % n_accept
Final.avg_rac = zeros(n_methods, 2); % r_accept
Final.avg_wrac = zeros(n_methods, 2); % weighted r_accept

Final.avg_cct = zeros(n_methods, 2); % cct
Final.avg_cct_norm = zeros(n_methods, 2); % cct
Final.avg_t_sch = zeros(n_methods, 2); % scheduling time
Final.avg_time = zeros(n_methods, 2); % time
   
% Comparison of rac and cct vs reference
Final.comp_rac = zeros(n_methods, n_iters);
Final.comp_cct = zeros(n_methods, n_iters);
Final.avg_comp_rac = zeros(n_methods, 1);
Final.avg_comp_cct = zeros(n_methods, 1);

% Percentile
Final.Pr_rac      = zeros(n_methods, 5); 
Final.Pr_cct      = zeros(n_methods, 5); 
Final.Pr_comp_rac = zeros(n_methods, 5); 
Final.Pr_comp_cct = zeros(n_methods, 5); 

% Prediction error of other methods than sincronia/maxCdsOptim/maxCdsRelax
Final.pred_err = zeros(n_methods, 1);

% Reference method
if ismember('maxCdsOptim', inputs.mod_names)
	reference = 'maxCdsOptim'; 
else
    reference = mod_names{end};
end

for mm = 1:n_methods
    mod_name = mod_names{mm};
    
    %%% Average n_accept, r_accept, cct
    Final.avg_nac(mm,:)  = [mean(Final.variant.(mod_name).nac), std(Final.variant.(mod_name).nac)];
    Final.avg_rac(mm,:)  = [mean(Final.variant.(mod_name).rac), std(Final.variant.(mod_name).rac)]; 
    Final.avg_wrac(mm,:) = [mean(Final.variant.(mod_name).wrac), std(Final.variant.(mod_name).wrac)]; 
    
    tmp = Final.variant.(mod_name).cct;
    tmp(isnan(tmp)) = mean(tmp(~isnan(tmp)));
    
    Final.variant.(mod_name).cct = tmp; % replace NaN by average values
    Final.avg_cct(mm,:) 	 = [mean(Final.variant.(mod_name).cct), std(Final.variant.(mod_name).cct)];
    Final.avg_cct_norm(mm,:) = [mean(Final.variant.(mod_name).cct_norm), std(Final.variant.(mod_name).cct_norm)]; 

    %---------------------------------------------------------------------------------------
    %%% Compare with reference method
    Final.comp_rac(mm,:) = (Final.variant.(mod_name).rac - Final.variant.(reference).rac) ...
                          ./Final.variant.(reference).rac;
    
    Final.comp_cct(mm,:) = (Final.variant.(mod_name).cct - Final.variant.(reference).cct) ...
                          ./Final.variant.(reference).cct;
						  
	Final.comp_cct_norm(mm,:) = (Final.variant.(reference).cct_norm - Final.variant.(mod_name).cct_norm) ...
								./Final.variant.(reference).cct_norm;
    %---------------------------------------------------------------------------------------
    
    %%% Average rapport of rac and cct
    % Final.avg_comp_rac(mm,:) = [mean(Final.comp_rac(mm,:)), std(Final.comp_rac(mm,:))];
    % Final.avg_comp_cct(mm,:) = [mean(Final.comp_cct(mm,:)), std(Final.comp_cct(mm,:))];
    Final.avg_comp_rac(mm) = sum(Final.comp_rac(mm,:))/n_iters;
    Final.avg_comp_cct(mm) = sum(Final.comp_cct(mm,:))/n_iters;
    
    %%% Percentiles
    Final.Pr_rac(mm, :) = [ prctile(Final.variant.(mod_name).rac, 01), ...
                            prctile(Final.variant.(mod_name).rac, 10), ...
                            prctile(Final.variant.(mod_name).rac, 50), ...
                            prctile(Final.variant.(mod_name).rac, 90), ...
                            prctile(Final.variant.(mod_name).rac, 99) ];
                                
    Final.Pr_cct(mm, :) = [ prctile(Final.variant.(mod_name).cct, 01), ...
                            prctile(Final.variant.(mod_name).cct, 10), ...
                            prctile(Final.variant.(mod_name).cct, 50), ...
                            prctile(Final.variant.(mod_name).cct, 90), ...
                            prctile(Final.variant.(mod_name).cct, 99) ];
							
	Final.Pr_cct_norm(mm, :) = [prctile(Final.variant.(mod_name).cct_norm, 01), ...
								prctile(Final.variant.(mod_name).cct_norm, 10), ...
								prctile(Final.variant.(mod_name).cct_norm, 50), ...
								prctile(Final.variant.(mod_name).cct_norm, 90), ...
								prctile(Final.variant.(mod_name).cct_norm, 99) ];

    Final.Pr_comp_rac(mm, :) = [prctile(Final.comp_rac(mm,:), 01), ...
                                prctile(Final.comp_rac(mm,:), 10), ...
                                prctile(Final.comp_rac(mm,:), 50), ...
                                prctile(Final.comp_rac(mm,:), 90), ...
                                prctile(Final.comp_rac(mm,:), 99) ];
                                              
    Final.Pr_comp_cct(mm, :) = [prctile(Final.comp_cct(mm,:), 01), ...
                                prctile(Final.comp_cct(mm,:), 10), ...
                                prctile(Final.comp_cct(mm,:), 50), ...
                                prctile(Final.comp_cct(mm,:), 90), ...
                                prctile(Final.comp_cct(mm,:), 99) ];   

	Final.Pr_comp_cct_norm(mm, :) = [prctile(Final.comp_cct_norm(mm,:), 01), ...
									prctile(Final.comp_cct_norm(mm,:), 10), ...
									prctile(Final.comp_cct_norm(mm,:), 50), ...
									prctile(Final.comp_cct_norm(mm,:), 90), ...
									prctile(Final.comp_cct_norm(mm,:), 99) ];							
                            
   %%% Prediction error of other methods than sincronia/maxCdsOptim/maxCdsRelax
    if inputs.mod_option.(mod_name).runPredictCct
        Final.pred_err(mm) = mean(Final.variant.(mod_name).pred_err);
    end
        
   %%% Computing time
   Final.avg_t_sch(mm,:) = [mean(Final.variant.(mod_name).t_sch), std(Final.variant.(mod_name).t_sch)];
   Final.avg_time(mm,:) = [mean(Final.variant.(mod_name).time), std(Final.variant.(mod_name).time)];
   
end

fprintf('\n[%d,%d]\n', n_machines, n_coflows);
fprintf('t_total = %.5fs (%.5fh), n_iters = %d\n', t_total, t_total/3600, n_iters);


%% SUMMARIZED TABLE

if params.build_recap_table 
    Final.mod_names = mod_names;
    preci = 4; % rounding results up to a certain number
    config_name = strcat('m',num2str(n_machines), '_c', num2str(n_coflows));
    if strcmp(trace_type, 'random')
        T_final_name = strcat('off_', num2str(n_iters), '_T_final_', trace_type, '_', config_name, '.txt');
    elseif strcmp(trace_type, 'facebook')
        T_final_name = strcat('off_', num2str(n_iters), '_T_final_', trace_type, '_', config_name, ...
                              '_f', num2str(inputs.flow_limit(2)), '.txt');
    end
    % Save T_final to text file
    T_final = utils.buildFinalTable(Final, preci);
    
    if params.export_recap_table
        writetable(T_final, strcat(params.base_path, params.table_path, T_final_name), 'Delimiter', ' ');
    end
    
    disp(T_final);
    
    % Debug
    if params.debug_mode
        fid = fopen(params.debug_path, 'a+');
        fprintf(fid, '\n[%d,%d]\n', n_machines, n_coflows);
        fprintf(fid, 't_total = %.5fs (%.5fh), n_iters = %d\n', t_total, t_total/3600, n_iters);
        fprintf(fid, 'method\trac_avg\twrac_avg\tcct_norm_avg\n');
        for i = 1:length(mod_names)
           fprintf(fid, '%s\t%.2f\t%.2f\t%.2f\n', T_final.method{i}, T_final.rac_avg(i), T_final.wrac_avg(i), T_final.cct_norm_avg(i)); 
        end
        %fclose(fid);
        fclose('all');
    end

end









