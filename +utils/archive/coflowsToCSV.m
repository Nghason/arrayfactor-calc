function T = coflowsToCSV(n_machines,coflows,filename,varargin)
% Generates a csv file from an array of Coflow objects
% Parameters:
% - n_machines: total number of machines in the fabric
% - coflows: An array of Coflow objects
% - filename: name of the output csv file (must contain .csv extension)
% Cedric Richier, LIA
% (c) LIA, 2020

H = 0;
if nargin == 4
    H = varargin{1};
end

tot_n_flows = sum([coflows.numFlows]);


flowID    = zeros(tot_n_flows,1);
coflowID  = zeros(tot_n_flows,1);
src       = zeros(tot_n_flows,1);
dst       = zeros(tot_n_flows,1);
bandwidth = zeros(tot_n_flows,1);
weight    = zeros(tot_n_flows,1);

T = table(flowID,coflowID,src,dst,bandwidth,weight);

f_id_base = 0;
j = 1;

for c = coflows
    c_id = c.id;
    c_w  = c.weight;
    for f = c.flows
        T.src(j)       = f.source.id;
        if H
            T.flowID(j)    = f_id_base + f.id - 1;
            T.coflowID(j)  = c_id - 1;            
            T.dst(j)       = f.destination.id - n_machines;
        else
            T.flowID(j)    = f_id_base + f.id;
            T.coflowID(j)  = c_id;            
            T.dst(j)       = f.destination.id;
        end
        T.bandwidth(j) = f.volume;
        T.weight(j)    = c_w;
        j = j+1;
    end
    f_id_base = f_id_base + c.numFlows;
end

writetable(T,filename);

end