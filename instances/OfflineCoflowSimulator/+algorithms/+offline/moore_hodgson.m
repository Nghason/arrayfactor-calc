% (c) 2022 Quang-Trung Luu (LAAS-CNRS)
% Implementation of Moore-Hodgson algorithm
% Reference:    J.M.Moore, “An n job, one machine sequencing algorithm for minimizing the number of late jobs,” 
%               Management Science, 15(1):102–109, 1968
% Inputs:
% - fabric: the Fabric ojbect representing the network
% - coflows: an array of Coflow objects within the fabric
% - ell: id of the considered port
% Outputs:
% - S_ell: set of accepted coflows on link ell
% - E_ell: set of rejected coflows on link ell


function [S_ell, E_ell] = moore_hodgson(fabric, coflows, ell)

global params;

%% Initializations:

n_links = fabric.numFabricPorts; % nb of links
n_coflows = length(coflows);     % nb of coflows

% Capacity of each link
portCapacity = [[fabric.machinesPorts.ingress] [fabric.machinesPorts.egress]];
portCapacity = [portCapacity.linkCapacity];

% Sort coflows wrt deadline (EDD order)
EDD_table = utils.sort_EDD_order(coflows);

% Matrix of processing time of each coflow on each link:
D = zeros(n_links, n_coflows);
for cid = 1:n_coflows
    c = coflows(cid);
    D(:, cid) = c.indicator*[c.flows.volume]' ./ portCapacity';
end


%% Initialization
    
% Processing time of sorted coflows on ell
D_ell = D(ell, [EDD_table.norm_id]); 

% Keep only sorted coflows having positive load on link ell
tmp = D_ell > 0;  % set of sorted coflows on ell

D_ell = D_ell(tmp);                 % processing times
S_ell = (EDD_table.norm_id(tmp))';  % real id in coflows
T_ell = (EDD_table.deadline(tmp))'; % deadlines

E_ell = [];

while true

    % Completion times (cumulative)
    D_ell_cum = cumsum(D_ell);  % completion times 
    dl_violation = D_ell_cum > T_ell;

    if sum(dl_violation) == 0
        break;
    end
    k = find(dl_violation == 1, 1); % first violated coflow

    % take the coflow having largest volume in 1:k
    % x* <- argmax(p_i, i <= k)
    S_until_k = S_ell(1:k);
    D_until_k = D_ell(1:k);
    x_star = S_until_k(D_until_k == max(D_until_k));

    x_star = x_star(1); % take the first in x_star (update 07/10/2021)

    % Remove xstar from S and D
    xstar_id = (S_ell == x_star);
    S_ell(xstar_id) = [];
    D_ell(xstar_id) = [];
    T_ell(xstar_id) = [];

    E_ell(end+1) = x_star;
end

% Return true id of coflows
S_ell = (EDD_table.EDD_id(ismember(EDD_table.norm_id, S_ell)))';
E_ell = (EDD_table.EDD_id(ismember(EDD_table.norm_id, E_ell)))';

if params.debug_mode
    fid = fopen(params.debug_path, 'a+');
    fprintf(fid, '\nRun Moore-Hodgson on link %d', ell);
    fprintf(fid, strcat(' with coflows [', strjoin(string([coflows.id]),','),']\n'));
    fprintf(fid, strcat('\taccepted coflows: S_{', num2str(ell), '} = [', strjoin(string(S_ell),','),']\n'));
    fprintf(fid, strcat('\trejected coflows: E_{', num2str(ell), '} = [', strjoin(string(E_ell),','),']\n'));
    fclose('all');
end



end

