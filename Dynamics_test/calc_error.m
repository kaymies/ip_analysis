clear all;
clc;
close all;
setpath
%% Load Data
% Human Data
subject_type = 'nonpar_7sub';
human_struct = load(sprintf('%s.mat',subject_type));
freq = human_struct.Frequency;
human_data = human_struct.IPDataAverage;
human_sd = human_struct.StandardDeviation;
% Best-Fit Simulation Data
controller = 'LQR_relative';
file = 'E4_1p9_3';
folder = fullfile('Data','Simulation',controller);
addpath(folder);
sim_struct = load(sprintf('%s.mat',file));
sim_data = sim_struct.DataWithoutOutliers;
generate_plot(human_struct,sim_struct)
%% Variance Accounted For (Rika's Method)
norm_difference = ((sim_data - human_data).^2)./human_data.^2
norm_variance = (human_sd./human_data).^2
vaf = 1 - sum(norm_difference)/sum(norm_variance)
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
%% RMSE
root_mean = sqrt(mean((meansim-human).^2));
%% Difference Beta (1.2-2.6 Hz)
% diffBeta = data(:,5:11) - human(5:11);
% meansimdiffBeta = mean(diffBeta, 'omitnan');
% meanBeta = mean(meansimdiffBeta, 'omitnan')
% CIBeta = 1.96*std(meansimdiffBeta, 'omitnan')/sqrt(7)
% figure;
% histogram(meansimdiffBeta)
%% RMSE 1.2-2.6 Hz
% rmse = sqrt(mean((data(:,5:11) - human(5:11)).^2,2));
% rmseBeta = mean(rmse, 'omitnan')
% rmseCIBeta = 1.96*std(rmse, 'omitnan')/sqrt(40)
% figure;
% histogram(rmse)
% %% Difference Noise (3-8 Hz)
% diff = data(:,14:38) - human(14:38); 
% % diff = data(:,24:38) - human(24:38); 
% meansimdiff = mean(diff, 'omitnan');
% meanSigma = mean(meansimdiff, 'omitnan')
% CISigma = 1.96* std(meansimdiff, 'omitnan')/sqrt(25)
% figure;
% histogram(meansimdiff)
%% RMSE 3-8 Hz
% rmse = sqrt(mean((data(:,14:38) - human(14:38)).^2,2));
% rmseSigma = mean(rmse, 'omitnan')
% rmseCISigma = 1.96*std(rmse, 'omitnan')/sqrt(40)
% figure;
% histogram(rmse)