freq = load('freq.mat');
simulation = load('simulation.mat');
human = load('human.mat');
freq = freq.freq;
simulation = simulation.IP_1e6_03_1;
human = human.human;
%% 
g = fittype('a-b*exp(-x/tau)');
[f0, gof] = fit(freq,human,g,'StartPoint',[[ones(size(freq)), -exp(-freq)]\human; 1])
xx = linspace(0,8,50);
plot(freq,human,'o',xx,f0(xx),'r-');
%% Neglect Outliers
[numsim, numf] = size(simulation);
meansim = mean(simulation);
stdsim = std(simulation);
lowBound = meansim - 3*stdsim;
upBound = meansim + 3*stdsim;
for i = 1:numsim
    for j = 1:numf
        if simulation(i,j) < lowBound(j)
            simulation(i,j) = NaN;
        end
        if simulation(i,j) > upBound(j)
            simulation(i,j) = NaN;
        end 
    end
end        
%% RMSE
% rmse = sqrt(mean((simulation - human).^2,2));
% nanmean(rmse)
% nanstd(rmse)
% figure;
% histogram(rmse)
%% RMSE Noise
% rmse = sqrt(mean((simulation(:,24:38) - human(24:38)).^2,2));
% nanmean(rmse)
% nanstd(rmse)
% figure;
% histogram(rmse)
%% Difference
% diff = 1./simulation - 1./human;
% meansimdiff = nanmean(diff);
% nanmean(meansimdiff)
% nanstd(meansimdiff)
% figure;
% histogram(meansimdiff)
%% Difference 1.2-2.6 Hz
% diff = 1./simulation - 1./human;
% meansimdiff = nanmean(diff(:,5:11));
% nanmean(meansimdiff)
% nanstd(meansimdiff)
% figure;
% histogram(meansimdiff)
%% Difference Noise
diff = simulation(:,24:38) - human(24:38);
meansimdiff = nanmean(diff);
nanmean(meansimdiff)
nanstd(meansimdiff)
figure;
histogram(meansimdiff)
%% 95% Confidence Interval
% Determine if human data is within 95% Confidence Interval
% meanSim = mean(simulation);
% % medSim = median(simulation(1:10,:));
% % errorSim = 1.96*std(simulation)/sqrt(length(simulation));
% errorSim = std(simulation);
% lowBound = meanSim - errorSim;
% upBound = meanSim + errorSim;
% 
% inBound = zeros(1,length(human));
% for i = 1:length(human)
%     if human(i) > lowBound(i) && human(i) < upBound(i)
%         inBound(i) = 1;
%     else
%         inBound(i) = 0;
%     end
% end
% 
% percentInBound = 100*sum(inBound)/length(inBound)
% figure;
% hold on
% errorbar(freq,meanSim,errorSim,'o')
% plot(freq,human)
% legend('error','human')