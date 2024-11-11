function [CCTs] = computeCCTs(coflows)

% LAST MODIFIED: 19/11/2020


%% Initalizations
n_coflows = length(coflows);
CCTs = zeros(n_coflows,1);

for c = coflows
    % !!!!  change it according to the version of pricing used !!! %
    CCTs(c.id) = max([c.flows.fct_pricing_2]);
end

end

