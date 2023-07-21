function [ip] = compute_ip(input)
    ip = computeAnalyticIP(input);
    if isnan(ip(1,1)) % if any part of the simulation resulted in numerical instabilties, just skip
        ip = NaN(num_trial,length(input.Frequency));
    end
    clear computeAnalyticIP
end
