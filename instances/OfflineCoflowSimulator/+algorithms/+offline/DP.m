% (c) 2022 Quang-Trung Luu (LAAS-CNRS)
% Dynamic programming to solve the problem 1||\sum w_j U_j 
% 1||\sum w_j U_j : Graham notation indicating 1 machine and objective \sum w_j U_j 
% Inputs:
% - fabric: the Fabric ojbect representing the network
% - Sb: set of coflows that uses the bottleneck link
% - b: id of the bottleneck link
% Outputs:
% - accepts: list of accepted coflows
% - rejects: list of rejected coflows


function [accepts, rejects] = DP(fabric, Sb, b)

global params;

n_coflows = length(Sb);             % nb of coflows in Sb
n_links = fabric.numFabricPorts;    % nb of links
Sb_id = [Sb.id];                    % id of coflows in Sb
W = sum([Sb.weight]);               % sum of coflow weights

% Capacity of each link
portCapacity = [[fabric.machinesPorts.ingress] [fabric.machinesPorts.egress]];
portCapacity = [portCapacity.linkCapacity];

% Demand of each coflow on each link:
D = zeros(n_links, n_coflows);
for cid = 1:n_coflows
    c = Sb(cid);
    D(:, cid) = c.indicator*[c.flows.volume]';
end


if params.debug_mode
    fid = fopen(params.debug_path, 'a+');
    fprintf(fid, 'Dynamic programming: W = %d\n', W);
    fprintf(fid, '\nProcessing time matrix D of coflow set Sb: [n_links * n_coflow]\n');
    dlmwrite(params.debug_path, D, 'delimiter', '\t', '-append', 'precision', '%.5f')
    fprintf(fid, '\n');
    fclose('all');
end


% Sort coflows wrt EDD
EDD_table = utils.sort_EDD_order(Sb);

% Initialization
P             = zeros(1, W+1);
P_prev        = inf*ones(1, W+1);
selected      = cell(1, W+1);
selected_prev = cell(1, W+1);
  
sum_w = 0;
iter = 1;
P_prev(1) = 0;
    

%% DP iterations

for i = 1:n_coflows

    cid = EDD_table.EDD_id(i);          % get real id of coflow in Sb with EDD order
    cid_norm =  EDD_table.norm_id(i);   % norm_id of coflow cid
    c = Sb([Sb.id] == cid);             % take coflow of id = cid in Sb
    w_k = c.weight;
    
    % Processing time of coflow cid on bottleneck
    p_k_on_b = D(b, cid_norm) / portCapacity(b); 
    sum_w = sum_w + w_k;
    
    % Dynamic programming equation
    for w = 0:W
        if w > sum_w
            P(w+1) = inf;
            selected{w+1} = [];    
        elseif (w >= w_k) && (p_k_on_b + P_prev(w-w_k+1) <= c.deadline)
            cct = p_k_on_b + P_prev(w - w_k + 1);
            if (P_prev(w+1) <= cct) 
                selected{w+1} = selected_prev{w+1};
                P(w+1) = P_prev(w+1);
            else 
                P(w+1) = cct;
                selected{w+1} = selected_prev{w-w_k+1};
                selected{w+1} = [selected{w+1}, c.id]; %EDD_id(i)]; %cid];
            end
        else
            P(w+1) = P_prev(w+1);
            selected{w+1} = selected_prev{w+1};
        end     
    end
    
    % Update P_prev and selected_prev
    for w = 0:W
        P_prev(w+1) = P(w+1);
        selected_prev{w+1} = selected{w+1};
    end
    
    % Debug
    if params.debug_mode
        fid = fopen(params.debug_path, 'a+');
        fprintf(fid, 'Iteration %d of dynamic programming (coflow %d), ', iter, cid);
        fprintf(fid, ' Sb = '); fprintf(fid, strcat('[',strjoin(string(Sb_id),','),']\n')); 
        fprintf(fid, '\tw=');  for w = 0:W;  fprintf(fid, '%d\t', w);         end; fprintf(fid, '\n');
        fprintf(fid, '\tP=');  for w = 0:W;  fprintf(fid, '%.2f\t', P(w+1));  end; fprintf(fid, '\n');
        fprintf(fid, '\tS=');  for w = 0:W;  fprintf(fid, strcat('{',strjoin(string(selected{w+1}),','),'}\t')); end
        fprintf(fid, '\n');
        % fclose(fid);
        fclose('all');
    end
    
    iter = iter + 1;

end % end loop of iterations 


%% Find the largest finite value of W that is feasible

for i = W+1:-1:1
    if P(i) < inf
        accepts = selected{i};
        rejects = Sb_id(~ismember(Sb_id, accepts)); 
        
        % Remove rejected coflows from Sb
        if ~isempty(accepts)
            for j = accepts
                Sb([Sb.id] == j) = [];
                Sb_id(Sb_id == j) = [];
            end
        end
        
        % Display
        if params.debug_mode
            fid = fopen(params.debug_path, 'a+');
            fprintf(fid, 'Max weight = %d\n', i-1);
            fprintf(fid, strcat('Accepted coflows: [',strjoin(string(accepts),','),']\n'));
            fprintf(fid, strcat('Rejected coflows: [',strjoin(string(rejects),','),']\n'));
            % fclose(fid);
            fclose('all');
        end
        
        break;
    end
end

    
% outputs.accepts = accepts;
% outputs.rejects = rejects;


end


