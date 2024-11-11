function [flows] = updateParamsFlows(flows,delta,t)

% LAST MODIFIED: 13/01/2021


%% Eventually, managing some precission level
precision = 5;
[flows.ad_rate] = utils.num2csl(round([flows.ad_rate],precision));
[flows.d_rate_old] = utils.num2csl([flows.d_rate]);


%% Updating the completion times according to the duration
time_to_completion = [flows.remainingVolume]./[flows.ad_rate]; %!!!!!!!!!!!!! division round

% !!!!  change it according to the version of pricing used !!! %
[flows.fct_pricing_2] = utils.num2csl([flows.fct_pricing_2] + min(time_to_completion,delta));


%% Computing remaining volumes
[flows.remainingVolume] = utils.num2csl(max([flows.remainingVolume] - delta * [flows.ad_rate],0));
[flows.remainingVolume] = utils.num2csl(round([flows.remainingVolume],precision)); % managing precision level on volumes


%% Deactivate flows that finished their sessions
[flows.state_f] = utils.num2csl(int8([flows.remainingVolume] > 0));

if ~isempty(flows([flows.state_f] == 0))
    [flows([flows.state_f] == 0).d_rate] = utils.num2csl(zeros(length([flows([flows.state_f] == 0).d_rate]),1));
    [flows([flows.state_f] == 0).ad_rate] = utils.num2csl(zeros(length([flows([flows.state_f] == 0).ad_rate]),1));
    [flows([flows.state_f] == 0).d_rate_old] = utils.num2csl(zeros(length([flows([flows.state_f] == 0).d_rate]),1));
end

if ~isempty(flows([flows.state_f] == 0 & [flows.departure] == -1))
    [flows([flows.state_f] == 0).departure] = utils.num2csl(t*ones(length([flows([flows.state_f] == 0).departure]),1));
end

end

