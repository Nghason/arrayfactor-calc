% (c) 2022 Cedric Richier (LIA-Avignon)
% Evaluate CCTs of coflows from a given order
% Inputs: fabric, coflows, order
% Outpus: ordered CCTs 


function cct = evalCCTfromOrder(fabric, coflows, order)

%% Initializations:
% Number of colfows
n_coflows = length(coflows);

% Number of links
n_links = fabric.numFabricPorts;

% Link capacities
tmp = [[fabric.machinesPorts.ingress] [fabric.machinesPorts.egress]];
B = [tmp.linkCapacity];

cct = zeros(1,n_coflows);

D = zeros(n_links, n_coflows);
for cc = coflows(order)
    D(:,cc.id) = cc.indicator*[cc.getFlowsVolume]';
end

D = D./repmat(B',1,n_coflows);

for c = coflows(flip(order))
    %v_c = zeros(n_links,1);
    if c.numFlows > 1
        v_c = max(c.indicator, [], 2); % max(c.indicator,2);
    else
        v_c = c.indicator;
    end
    cct(c.id) = max(sum(D,2).*v_c);
    
    D(:,c.id) = 0;
end

% reorder vector cct w.r.t the order (19/11/2021)
cct = cct(order);

end