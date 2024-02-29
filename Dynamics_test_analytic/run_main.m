% close all;
% clear all;
% clc;
% setpath
%% Load Human Data
subject_type = 'duarte_old';
human_struct = load(sprintf('%s.mat',subject_type));
% f_i = 0.65; f_int = 0.5; f_end = 5.15;
% input.Frequency = f_i:f_int:f_end;
input.Frequency = human_struct.Frequency_cpsd;
%% Set Model Parameters
input.TotalMass = human_struct.MeanMass_kg;
input.TotalHeight = human_struct.MeanHeight_m;  
input.gender = 'M';
input.plane = human_struct.Plane;
input.model = 'DIP';
input.pose = human_struct.Pose;
%% Set Controller Parameters
%-----% MAKE SURE TO CHANGE SIM FILENAME!!! %-----%
filename = 'E6_2p3_1p7'; %alpha_beta_sigmar
folder_name = 'LQR_relative';
%-------------------------------------------------%
input.Controller.alpha = 1e6;
input.Controller.beta = 2.3;
input.NoiseRatio = 1.7;
input.Controller.gamma = 1;
input.Controller.kappa = 1;
input.Controller.eta = 1;
input.Controller.type = 'LQR';
input.simSensorNoiseRatio = 1;
input.simMotorNoise = 10;
input.simSensorNoise = 0;
input.Controller.param.delay =0; % 2023-05-18
%% Compute IP
data = compute_ip(input);
%% Plot human vs. simulation
figure();
plot(input.Frequency,human_struct.IPDataAverage_cpsd,'x')
hold on
plot(input.Frequency,data,'o')
% plot(sim_struct.Frequency,sim_struct.DataWithoutOutliers)
legend('human', 'simdata')
xlabel('Frequency [Hz]')
ylabel('IP (Fraction of CoM)')
ylim([0 3.5])
%% Error
rmse = sqrt(mean((data-human_struct.IPDataAverage_cpsd).^2))
%% Save File:
sim_struct = save_sim_file(input,filename,data,meansim,numtrial,folder_name,stdsim);
%% Variance Accounted For (Rika's Method)
human_data = human_struct.IPDataAverage;
human_sd = human_struct.StandardDeviation;
norm_difference = ((data - human_data).^2)./human_data.^2;
norm_variance = (std(human_struct.IPAveSubject)./human_data).^2;
vaf = 1 - sum(norm_difference)/sum(norm_variance)