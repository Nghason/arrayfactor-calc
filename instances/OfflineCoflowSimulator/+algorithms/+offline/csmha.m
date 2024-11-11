% (c) 2022 Quang-Trung Luu (LAAS-CNRS)
% Implementation of CS-MHA algorithm
% Reference: S. Luo et al., “Decentralized deadline-aware coflow scheduling for datacenter networks,” in Proc. IEEE ICC, 2016, pp. 1–6.
% Inputs:
% - fabric: the Fabric ojbect representing the network
% - coflows: an array of Coflow objects within the fabric
% Outputs:
% - order: scheduling order of coflows


function outputs = csmha(mod_name, fabric, coflows)


global params;


%% Initializations:

n_links = fabric.numFabricPorts; % nb of links
n_coflows = length(coflows);     % nb of coflows
S = [coflows.id];                % unscheduled coflows IDs:

% Capacity of each link
portCapacity = [[fabric.machinesPorts.ingress] [fabric.machinesPorts.egress]];
portCapacity = [portCapacity.linkCapacity];

% Matrix of processing time of each coflow on each link:
D = zeros(n_links, n_coflows);
for c = coflows
    D(:, c.id) = c.indicator*[c.flows.volume]' ./ portCapacity';
end

% Initialization
S_final = []; % set of admitted coflows
E_final = [];  % set of rejected coflows

if params.debug_mode
    fid = fopen(params.debug_path, 'a+');
    fprintf(fid, '\nEXECUTE CS-MHA ALGORITHM\n\n');
    fprintf(fid, 'Processing time matrix D: [n_links * n_coflow]\n');
    dlmwrite(params.debug_path, D, 'delimiter', '\t', '-append', 'precision', '%.5f')
    fprintf(fid, '\n');
    fclose('all');
end
    
for ell = 1:n_links
    
    switch mod_name
        case 'cs_mha' % use Moore-Hodgson algorithm on link ell
            [S_ell, E_ell] = algorithms.offline.moore_hodgson(fabric, coflows, ell);
            
        case 'cs_dp' % use DP on link ell
            % Coflows that use ell
            C_ell = D(ell,:) > 0; % set of coflows using ell
            [S_ell, E_ell] = algorithms.offline.DP(fabric, coflows(C_ell), ell);
    end
    
    % Set of rejected coflows on all links
    E_final = union(E_final, E_ell);
end

% Set of accepted coflows in all links
S_final = setdiff(S, E_final);

if params.debug_mode
    fid = fopen(params.debug_path, 'a+');
    fprintf(fid, '\nPreliminary results:\n');
    fprintf(fid, strcat('\tE_final = union(all E_ell) = [', strjoin(string(E_final),','),']\n'));
    fprintf(fid, strcat('\tS_final = S\\E_final = [', strjoin(string(S_final),','),']\n\n'));
    fclose('all');
end



%% Sort S_final w.r.t deadlines
EDD_table = utils.sort_EDD_order(coflows); % sort coflows wrt deadline (EDD order)
S_final = (EDD_table.EDD_id(ismember([EDD_table.EDD_id], S_final)))';

%----------------------------------------------------------------------------------------------------
if params.debug_mode
    fid = fopen(params.debug_path, 'a+');
    fprintf(fid, 'Sort S_final in increasing order of deadline (EDD):\n');
    fprintf(fid, strcat('\tDeadlines of coflows [', strjoin(string(S_final),','),'] are '));
    fprintf(fid, strcat(' [', strjoin(string( EDD_table.deadline(ismember(EDD_table.EDD_id, S_final)) ),','),']'));
    fprintf(fid, strcat(' => S_final = [', strjoin(string(S_ell),','),']\n'));
    fprintf(fid, '\nSort E_final in increasing order of max_{all ell} vol_{k,ell} / (T_k + epsilon)\n');
    fclose('all');
end
%----------------------------------------------------------------------------------------------------


%% Sort E_final w.r.t max(all ell) p_{ell,i}/T_i

if ~isempty(E_final)
    sort_E_table = table;
    
    if size(E_final, 1) == 1
        sort_E_table.id = E_final';
    else
        sort_E_table.id = E_final;
    end

    tmp = zeros(length(E_final), 1);
    for i = 1:length(E_final)
        c = coflows([coflows.id] == E_final(i));
        c_vols = D(:,c.id);
        
        switch mod_name
            case 'cs_mha' 
                c_criteria = max(c_vols)/(c.deadline + 1e-6);
                %----------------------------------------------------------------------------------------------------
                if params.debug_mode
                    fid = fopen(params.debug_path, 'a+');
                    fprintf(fid, '\tc%d:', c.id);
                    fprintf(fid, strcat(' volume on link [',  strjoin(string(1:n_links),','),  '] are'));
                    fprintf(fid, strcat(' [',  strjoin(string(D(:,c.id)),','),  ']\n'));
                    fprintf(fid, '\t\tCriteria for c%d = max(c_vols)/(c.deadline + 1e-6)\n', c.id);
                    fprintf(fid, strcat('\t\t\t\t\t\t= max(',  strjoin(string(D(:,c.id)),','),  ')'));
                    fprintf(fid, '/(%.5f + 1e-6)\n', c.deadline);
                    fprintf(fid, '\t\t\t\t\t\t= %.5f/(%.5f + 1e-6)\n', max(c_vols), c.deadline);
                    fprintf(fid, '\t\t\t\t\t\t= %.5f\n', c_criteria);
                    fclose('all');
                end
                %----------------------------------------------------------------------------------------------------
                
            case 'cs_dp'
                c_criteria = c.weight * max(c_vols)/(c.deadline + 1e-6);
                %----------------------------------------------------------------------------------------------------
                if params.debug_mode
                    fid = fopen(params.debug_path, 'a+');
                    fprintf(fid, '\tc%d:', c.id);
                    fprintf(fid, strcat(' volume on link [',  strjoin(string(1:n_links),','),  '] are'));
                    fprintf(fid, strcat(' [',  strjoin(string(D(:,c.id)),','),  ']\n'));
                    fprintf(fid, '\t\tCriteria for c%d = c.weight * max(c_vols)/(c.deadline + 1e-6)\n', c.id);
                    fprintf(fid, '\t\t\t\t\t\t= %d *', c.weight);
                    fprintf(fid, strcat(' max(',  strjoin(string(D(:,c.id)),','),  ')'));
                    fprintf(fid, '/(%.5f + 1e-6)\n', c.deadline);
                    fprintf(fid, '\t\t\t\t\t\t= %d * %.5f/(%.5f + 1e-6)\n', c.weight, max(c_vols), c.deadline);
                    fprintf(fid, '\t\t\t\t\t\t= %.5f\n', c_criteria);
                    fclose('all');
                end
        end
    
        
        tmp(i) = c_criteria;
        
        
    end
    sort_E_table.criteria = tmp;
    
    % Short tmp wrt criteria
    sort_E_table = sortrows(sort_E_table, 'criteria'); % sort the table by 'criteria'
    E_final = (sort_E_table.id)';
end
clearvars c_vols c_criteria tmp;

% Finally
order = [S_final, E_final];
outputs.order = order;

if params.debug_mode
    fid = fopen(params.debug_path, 'a+');
    fprintf(fid, '\nFinally:\n');
    fprintf(fid, strcat('\tS_final = [', strjoin(string(S_final),','),'], '));
    fprintf(fid, strcat(' E_final = [', strjoin(string(E_final),','),'], '));
    fprintf(fid, strcat(' sigma = [',   strjoin(string(order),','),']\n'));
    fclose('all');
end
        




end

