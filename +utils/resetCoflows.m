function resetCoflows(coflows)
%
% Cedric Richier, LIA
% (c) LIA, 2021

for c = coflows
    c.state_c = 1;
    %c.priority = 0;
    c.prices = zeros(1,length(c.prices));
    c.maxPrices = c.prices;
    c.prices_prev = c.prices;
    c.stability = 0;
    c.current_w = 0;
    c.max_diff_prices = 0;
    c.ts_counter = 0;
    c.remaning_vol = 0;
    c.departure = -1;
    for f = c.flows
        f.remainingVolume = f.volume;
        f.d_rate = 1;
        f.fct = 0;
        f.state_f = 1;
    end
end

end