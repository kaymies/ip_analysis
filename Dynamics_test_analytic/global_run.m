close all;
clear all;
clc;
setpath
%% To edit each time you run:
subject_type = 'duarte_old';
filename_rmse = 'rmse_duarte_old_finer_20240212';
%% Load Human Data
human_struct = load(sprintf('%s.mat',subject_type));
human = human_struct.IPDataAverage_cpsd;
% human = human_struct.IPDataAverage;
% f_i = 0.65; f_int = 0.5; f_end = 5.15;
% input.Frequency = f_i:f_int:f_end;
input.Frequency = human_struct.Frequency_cpsd;
% input.Frequency = human_struct.Frequency;
%% Set Model Parameters
input.TotalMass = human_struct.MeanMass_kg;
input.TotalHeight = human_struct.MeanHeight_m;  
input.gender = 'M';
input.plane = human_struct.Plane;
input.model = 'DIP';
input.pose = human_struct.Pose;
%% Set Controller Parameters
input.Controller.type = 'LQR';
input.simSensorNoiseRatio = 1;
input.simMotorNoise = 10;
input.simSensorNoise = 0;
input.Controller.param.delay =0; % 2023-05-18
%% Create empty arrays to log errors
% params.beta = [0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3];
% params.sigma_r = [0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3,5,10,15,20];
params.beta = 0.15:0.01:0.45;
% params.beta = 0.1:0.01:0.2;
params.sigma_r = [0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1];
% params.alpha = [1e-4,1e-3,1e-2,1,1e6]; %5
% params.beta = [0.3];
% params.sigma_r = [1];
% params.gamma = [0.1,1,10,35]; %4
% params.kappa = [1e-5,1e-4,1e-3,1e-2,1e-1,1,10,100,1e3,1e4,1e5]; %11
% params.eta = [1e-1,1,10,1e2]; %4

% rmse_nonpar7sub_lqr_6params = zeros(length(alpha),length(beta),length(sigma_r),...
%     length(gamma),length(kappa),length(eta));\
% load('C:\Users\kaymi\Dropbox (MIT)\Balance Project\Intersection Point Analysis\ip_analysis\Dynamics_test\Data\Error\rmse20230223_3.mat')
% rmse_nonpar7sub_lqr_6params = Error;
rmse = zeros(length(params.beta),length(params.sigma_r));
%% Run Simulation for n trials across different parameter values
count = 0;
tic
input.Controller.alpha = 1e6;
input.Controller.gamma = 1;
input.Controller.kappa = 1;
input.Controller.eta = 1;
for b = 1:length(params.beta)
    for s = 1:length(params.sigma_r)
        input.Controller.beta = params.beta(b);
        input.NoiseRatio = params.sigma_r(s);
        ip = compute_ip(input);
        rmse(b,s) = sqrt(mean((ip-human).^2));
%         figure(1);
%         plot(human_struct.Frequency,human_struct.IPDataAverage)
%         hold on
%         plot(human_struct.Frequency,ip)
%         legend('human', 'simdata')
        count = count+1
    end
end
toc
save_error_file(input,subject_type,params,filename_rmse,rmse,human_struct)