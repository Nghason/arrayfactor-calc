function [coflows] = orderingCoflows(coflows, n_coflows,order)
%function [coflows] = orderingCoflows(coflows, n_coflows,Order) % 1PARIS

% LAST MODIFIED: 28/10/2020

%coflowsOrder = randperm(n_coflows);
% ID of coflows in random order

% coflows = coflows(coflowsOrder);
coflows = coflows(order); % 1PARIS
% ordering coflows in the generated order

%coflowsPriority = num2cell(1:length(coflowsOrder));
coflowsPriority = num2cell(1:length(order)); % 1PARIS


[coflows.priority] = coflowsPriority{:};
% assigning priority to each coflow


end

