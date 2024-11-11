% Computes the value of CCT0 for a given coflow, that is its completion time when being alone in the fabric (isolation)
% Returns the value of CCT0 (notation 'gamma' in Chowdhury's paper related to Varys system)
% Parameters:
% - c: a coflow object
% - c_rvol: array of remaining volumes of all flows in coflow c
% - Rem: array of remaining bandwidth of each link in the fabric
% Cedric Richier, LIA
% (c) LIA, 2020

function gamma = computeCCT0(c,c_rvols,Rem)

% Get the indicator matrix to know sources and destinations of flows
indic = c.indicator;

% Compute the aggregated volume on each link
sums = indic*c_rvols';

% Compute gamma (defined in Chowdhury's paper related to Varys system)
gamma = max(sums./Rem');

end