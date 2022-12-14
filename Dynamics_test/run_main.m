close all;
clear all;
clc;
setpath
%% Set Model Parameters
input.TotalMass = 71.1; %71.1 Gruben
input.TotalHeight = 1.75; %1.75 Gruben    
input.gender = 'M';
input.plane = 'sgt';
input.model = 'DIP';
input.pose = 'pose_I';
input.FreqSampKin = 100;
input.trialDuration = 60;
input.CoordinateFrame = 'relative';
population = 'nonpar_7sub';
%% Set Controller Parameters
input.Controller.alpha = 1e-4;
input.Controller.beta = .2;
input.Controller.gamma = 1;
input.NoiseRatio = 0.2;
input.Controller.type = 'LQR_int';
%% Visualization
input.PostProc.AnimOn = 0;
input.PostProc.PlotOn = 0;
%% Run Simulation for n trials
numtrial = 40;
[f,IP_ratio] = main_ip(input);
data = zeros(numtrial,length(f));
torques = zeros(numtrial, 3);
data(1,:) = IP_ratio;
% torques(1,:) = torque_rms;
for i = 1: numtrial-1
    close all;
    [~,IP_ratio] = main_ip(input);
    data(i+1,:) = IP_ratio;
%     torques(i+1,:) = torque_rms;
    i
end
figure; boxplot(data,f+0.1)
xlabel('Frequency [Hz]')
ylabel('IP (Fraction of CoM)')
ylim([0 2.5])
% title('gamma = 1e6, beta = 2, t = 30 sec, freq = 100 Hz, 0.1x noise')
set(gca,'Fontsize',15);

% mean_torque_ratio = mean(torques)
% std_torque_ratio = std(torques)
%% Make sure to change variable names
% Save simulation results
% folder = '\..\Data';
filename = sprintf('simulation_%s_%s',population,input.Controller.type);
IP_1en4_02_02 = data;
save(filename,'IP_1en4_02_02','-append')
% rmpath([pwd,'\..\Data']);
%% Neglect Outliers
% data = frt_nobeam_sub;
% data = IP_1en4_02_02;
[numsim, numf] = size(data);
meansim = mean(data);
stdsim = std(data);
lowBound = meansim - 3*stdsim;
upBound = meansim + 3*stdsim;
for i = 1:numsim
    for j = 1:numf
        if data(i,j) < lowBound(j)
            data(i,j) = NaN;
        end
        if data(i,j) > upBound(j)
            data(i,j) = NaN;
        end 
    end
end
meansim = mean(data,'omitnan'); % mean of simulation without outliers
%% Plot human vs. simulation
% freq = load('freq.mat'); f = freq.freq;
human = load('human.mat');
human = eval(sprintf('human.%s',population))';

figure(2);
plot(f+0.1,human)
hold on
plot(f+0.1,meansim)
legend('human', 'simdata')
xlabel('Frequency [Hz]')
ylabel('IP (Fraction of CoM)')
xlim([0.5 5.3])
ylim([0 2.5])
% ave = mean(data);
% stdev = mean(data);
% figure;
% boxplot(data,f)
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