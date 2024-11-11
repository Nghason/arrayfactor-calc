function displayCoflows(coflows)
% Displays information about coflows structure
% Parameters:
% - coflows: an array of coflow objets
% Cedric Richier, LIA
% (c) LIA, 2020


for c = coflows
    fprintf('Coflow %d: deadline=%.4f\n', c.id, c.deadline);
    for f = c.flows
        fprintf('\t flow %d: vol = %.4f, IN=%d , OUT = %d\n\n',...
            f.id,f.volume,f.source.id,f.destination.id);
    end
end

%%

% for c = coflows
%     tbl = table;
%     for f = c.flows
%        tab = table(c.id, c.arrival, c.deadline, c.weight, c.numFlows, strcat(num2str(c.id), num2str(f.id)), f.volume, f.links(1), f.links(2));
%        tbl = [tbl; tab];
%     end
%     %tbl{c.id+1,:} = missing; 
%     tbl.Properties.VariableNames = {'cid','c_arrival','c_deadline','c_weight','n_flows','fid','f_vol','f_in','f_out'};
%     disp(tbl);
% end

% end