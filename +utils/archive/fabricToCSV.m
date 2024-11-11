function F = fabricToCSV(fabric,filename)
% Generates a csv file that depict a fabric from a fabric object
% Returns the table that correspond to the csv format
% Parameters:
% - fabric: the fabric object
% - filename: the name of the csv file

n_links = fabric.numFabricPorts;
n_machines = n_links/2;
net_node = 0;

linkID = zeros(n_links,1);
src     = zeros(n_links,1);
dst     = zeros(n_links,1);
cap     = zeros(n_links,1);

F = table(linkID,src,dst,cap);
l_id = 0;
for p = [fabric.machinesPorts]
    F.linkID(l_id+1) = l_id;
    F.src(l_id+1) = p.ingress.id;
    F.dst(l_id+1) = net_node;
    F.cap(l_id+1) = p.ingress.linkCapacity;
    F.linkID(l_id+2) = l_id+1;
    F.src(l_id+2) = net_node;
    F.dst(l_id+2) = p.egress.id-n_machines;
    F.cap(l_id+2) = p.egress.linkCapacity;
    l_id = l_id+2;
end

writetable(F,filename);

end