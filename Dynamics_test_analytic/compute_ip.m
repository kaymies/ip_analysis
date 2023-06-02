function [ip,torques] = compute_ip(input,num_trial)
    ip = zeros(num_trial,length(input.Frequency));
%     torques = zeros(num_trial,3);
    for i = 1:num_trial
        ip(i,:) = computeAnalyticIP(input);
        if isnan(ip(i,1)) % if any part of the simulation resulted in numerical instabilties, just skip
            ip = NaN(num_trial,length(input.Frequency));
            break
        end

%         i
%         [~,warnID] = lastwarn;
%         if strcmp(warnID,'MATLAB:illConditionedMatrix')
%             ip(:) = NaN;
%             break
%         end
    end
    clear computeAnalyticIP
end
