% close all;
% clear all;
% clc;
% setpath
%% Load Human Data
subject_type = 'duarte_old';
human_struct = load(sprintf('%s.mat',subject_type));
% f_i = 0.65; f_int = 0.5; f_end = 5.15;
% input.Frequency = f_i:f_int:f_end;
input.Frequency = human_struct.Frequency;
%% Set Model Parameters
input.TotalMass = human_struct.MeanMass_kg;
input.TotalHeight = human_struct.MeanHeight_m;  
input.gender = 'M';
input.plane = human_struct.Plane;
input.model = 'DIP';
input.pose = human_struct.Pose;
%% Set Controller Parameters
%-----% MAKE SURE TO CHANGE SIM FILENAME!!! %-----%
filename = 'test'; %alpha_beta_sigmar
folder_name = 'LQR_relative';
%-------------------------------------------------%
input.Controller.alpha = 1e6;
input.Controller.beta = 0.2;
input.NoiseRatio = 0.4;
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
plot(human_struct.Frequency,human_struct.IPDataAverage,'x')
hold on
plot(human_struct.Frequency,data,'o')
% plot(sim_struct.Frequency,sim_struct.DataWithoutOutliers)
legend('human', 'simdata')
xlabel('Frequency [Hz]')
ylabel('IP (Fraction of CoM)')
ylim([0 3.5])
%% Error
rmse = sqrt(mean((data-human_struct.IPDataAverage).^2))
%% Variance Accounted For (Rika's Method)
human_data = human_struct.IPDataAverage;
human_sd = human_struct.StandardDeviation;
norm_difference = ((data - human_data).^2)./human_data.^2;
norm_variance = (std(human_struct.IPAveSubject)./human_data).^2;
vaf = 1 - sum(norm_difference)/sum(norm_variance)