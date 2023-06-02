close all;
clear all;
clc;
setpath
%% Load Human Data
subject_type = 'marta_sgt_ground';
human_struct = load(sprintf('%s.mat',subject_type));
% f_i = 0.65; f_int = 0.5; f_end = 5.15;
% input.Frequency = f_i:f_int:f_end;
input.Frequency = human_struct.Frequency;
input.FrequencyWindow = 0.2;
%% Set Model Parameters
input.TotalMass = 68.7;
input.TotalHeight = 1.61;  
input.gender = 'M';
input.plane = 'sgt';
input.model = 'DIP';
input.pose = 'pose_T';
% input.FreqSampKin = 100;
% input.trialDuration = 60;
% input.CoordinateFrame = 'relative';
%% Set Controller Parameters
%-----% MAKE SURE TO CHANGE SIM FILENAME!!! %-----%
filename = 'test'; %alpha_beta_sigmar
folder_name = 'LQR_relative';
%-------------------------------------------------%
input.Controller.alpha = 1e6;
input.Controller.beta = 2.4;
input.NoiseRatio = 1.6;
input.Controller.gamma = 1;
input.Controller.kappa = 1;
input.Controller.eta = 1;
input.Controller.type = 'LQR';
input.simSensorNoiseRatio = 1;
input.simMotorNoise = 10;
input.simSensorNoise = 0;
input.Controller.param.delay =0; % 2023-05-18
%% Visualization
input.PostProc.AnimOn = 0;
input.PostProc.PlotOn = 0;
%% Run Simulation for n trials
numtrial = 1;
[data] = compute_ip(input,numtrial);

figure; boxplot(data,input.Frequency)
xlabel('Frequency [Hz]')
ylabel('IP (Fraction of CoM)')
ylim([0 2.5])
% title('gamma = 1e6, beta = 2, t = 30 sec, freq = 100 Hz, 0.1x noise')
set(gca,'Fontsize',15);

% mean_torque_ratio = mean(torques)
% std_torque_ratio = std(torques)
%% Neglect Outliers
[data] = neglect_outliers(data);
meansim = mean(data,'omitnan'); % mean of simulation without outliers
stdsim = std(data,'omitnan');
%% Save File:
sim_struct = save_sim_file(input,filename,data,meansim,numtrial,folder_name,stdsim);
%% Plot human vs. simulation
generate_plot(human_struct,sim_struct)
%% Error
rmse = sqrt(mean((meansim-human_struct.IPDataAverage).^2))

diffBeta = data(:,5:11) - human_struct.IPDataAverage(5:11);
meansimdiffBeta = mean(diffBeta, 'omitnan');
meanBeta = mean(meansimdiffBeta, 'omitnan')
% CIBeta = 1.96*std(meansimdiffBeta, 'omitnan')/sqrt(7)
diff = data(:,24:end) - human_struct.IPDataAverage(24:end);
meansimdiffSigma = mean(diff,'omitnan');
meanSigma = mean(meansimdiffSigma,'omitnan')
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