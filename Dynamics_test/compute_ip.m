function [ip] = compute_ip(input,num_trial)
    ip = zeros(num_trial,length(input.Frequency));
    for i = 1:num_trial
        ip(i,:) = main_ip(input);
%         i
%         [~,warnID] = lastwarn;
%         if strcmp(warnID,'MATLAB:illConditionedMatrix')
%             ip(:) = NaN;
%             break
%         end
    end
    clear main_ip
end
