close all;
clear all;
clc;
setpath
%% Load Human Data
subject_type = 'marta_sgt_ground';
human_struct = load(sprintf('%s.mat',subject_type));
input.method = 'bpf'; %bpf or cpsd
% input.Frequency = eval(sprintf('human_struct.Frequency_%s',input.method));
input.Frequency = eval(sprintf('human_struct.Frequency'));
input.FrequencyWindow = 0.2;
% human_mean = eval(sprintf('human_struct.IPDataAverage_%s',input.method));
human_mean = eval(sprintf('human_struct.IPDataAverage'));
% human_data_sub = eval(sprintf('human_struct.IPAveSubject_%s',input.method));
% human_sd = std(human_data_sub);
%% Set Model Parameters
input.TotalMass = human_struct.MeanMass_kg;
input.TotalHeight = human_struct.MeanHeight_m;  
input.gender = 'M';
input.plane = human_struct.Plane;
input.model = 'DIP';
input.pose = human_struct.Pose;
input.FreqSampKin = 500;
input.trialDuration = 24;
input.CoordinateFrame = 'relative';
%% Set Controller Parameters
%-----% MAKE SURE TO CHANGE SIM FILENAME!!! %-----%
filename = 'test'; %alpha_beta_sigmar
folder_name = 'LQR_relative';
%-------------------------------------------------%
input.Controller.alpha = 1e6;
input.Controller.beta = 0.3;
input.NoiseRatio = 1.5;
input.Controller.gamma = 1;
input.Controller.kappa = 1;
input.Controller.eta = 1;
input.Controller.type = 'LQR';
%% Visualization
input.PostProc.AnimOn = 0;
input.PostProc.PlotOn = 0;
%% Run Simulation for n trials
numtrial = 40;
[data,torques,f] = compute_ip(input,numtrial);

% figure; boxplot(data,input.Frequency)
% xlabel('Frequency [Hz]')
% ylabel('IP (Fraction of CoM)')
% ylim([0 2.5])
% set(gca,'Fontsize',15);

mean_torque_ratio = mean(torques)
std_torque_ratio = std(torques)
%% Neglect Outliers
[data] = neglect_outliers(data);
meansim = mean(data,'omitnan'); % mean of simulation without outliers
stdsim = std(data,'omitnan');
%% Save File:
sim_struct = save_sim_file(input,filename,data,meansim,numtrial,folder_name,stdsim);
%% Plot human vs. simulation
generate_plot(input.Frequency,human_mean,[],meansim,stdsim)
% generate_plot(input.Frequency,human_mean,human_sd,meansim,stdsim)
%% Error
rmse = sqrt(mean((meansim-human_mean).^2))
%% Variance Accounted For (Rika's Method)
norm_difference = ((meansim - human_mean).^2);%./human_mean.^2;
% norm_variance = (human_sd).^2;
norm_variance = (human_mean - mean(human_mean)).^2;
vaf = 1 - sum(norm_difference)/sum(norm_variance)
% vaf = zeros(77,1);
% for i = 1:77
%     diff = (human_data_sub(:,i) - meansim(i)).^2;
%     var = (human_data_sub(:,i) - human_mean(i)).^2;
%     vaf(i) = 1 - sum(diff)/sum(var);
% end
% mean_vaf = mean(vaf)
% diffBeta = data(:,5:11) - human_struct.IPDataAverage(5:11);
% meansimdiffBeta = mean(diffBeta, 'omitnan');
% meanBeta = mean(meansimdiffBeta, 'omitnan')
% % CIBeta = 1.96*std(meansimdiffBeta, 'omitnan')/sqrt(7)
% diff = data(:,24:end) - human_struct.IPDataAverage(24:end);
% meansimdiffSigma = mean(diff,'omitnan');
% meanSigma = mean(meansimdiffSigma,'omitnan')
%% Plot IP Distribution
% for j = 1:length(input.Frequency)
%     figure(3)
%     subplot(2,5,j)
%     h = histogram(data(:,j),7);
%     hold on
%     pd = fitdist(data(:,j),'Normal');
%     x = linspace(h.BinEdges(1),h.BinEdges(end),100);
%     plot(x,normpdf(x,pd.mu,pd.sigma))
%     title(sprintf('Frequency: %0.5g - %0.5g Hz',input.Frequency(j),input.Frequency(j)+0.2))
%     figure(4)
%     subplot(2,5,j)
%     qqplot(data(:,j))
%     title(sprintf('Frequency: %0.5g - %0.5g Hz',input.Frequency(j),input.Frequency(j)+0.2))
% end