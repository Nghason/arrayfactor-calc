% Sort coflows wrt EDD (Earliest Due Date) order
% (c) 2022 Quang-Trung Luu

function tab = sort_EDD_order(coflows)

global params;

n_coflows = length(coflows);
tab.norm_id = (1:n_coflows)';   % normalized id: from 1 to n_coflows
tab.EDD_id = [coflows.id]';     % real id of coflows
tab.deadline = [coflows.deadline]';   % due date (deadline) of coflows

% Short tmp wrt criteria
tab = struct2table(tab);         % convert the struct array to a table
tab = sortrows(tab, 'deadline'); % sort the table by 'criteria'

% if params.debug_mode
%     fid = fopen(params.debug_path, 'a+');   
%     fprintf(fid, strcat('EDD order: [',  strjoin(string(tab.EDD_id), ', '), ']\n'));
% end


end



