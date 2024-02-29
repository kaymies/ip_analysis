close all;
clear all;
clc;
setpath
%% To edit each time you run:
subject_type = 'duarte_young';
filename_rmse = 'rmse_duarte_young_finer_subs_20240212';
input.method = 'cpsd';
%% Load Human Data
human_struct = load(sprintf('%s.mat',subject_type));
human = eval(sprintf('human_struct.IPAveSubject_%s',input.method));
% f_i = 0.65; f_int = 0.5; f_end = 5.15;
% input.Frequency = f_i:f_int:f_end;
input.Frequency = eval(sprintf('human_struct.Frequency_%s',input.method));
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
params.beta = 0.15:0.01:0.45;
params.sigma_r = [0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1];
% params.beta = [0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3];
% params.sigma_r = [0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3,5,10,15,20];
% params.beta = [0.2,0.5];
% params.sigma_r = [0.4,2];
rmse = zeros(human_struct.NumSubjects,length(params.beta),length(params.sigma_r));
%% Run Simulation for n trials across different parameter values
count = 0;
point = 0;
tic
input.Controller.alpha = 1e6;
input.Controller.gamma = 1;
input.Controller.kappa = 1;
input.Controller.eta = 1;
for sub = 1:human_struct.NumSubjects
%     for trial = 1:3
%         point = point+1
        for b = 1:length(params.beta)
            for s = 1:length(params.sigma_r)
                input.Controller.beta = params.beta(b);
                input.NoiseRatio = params.sigma_r(s);
                ip = compute_ip(input);
                rmse(sub,b,s) = sqrt(mean((ip-human(sub,:)).^2));
%                 rmse(sub,b,s) = sqrt(mean((ip-reshape(human_struct.IP(sub,:),size(human_struct.Frequency))).^2));
    %             figure(1);
    %             plot(human_struct.Frequency,human_struct.IPAveSubject(sub,:))
    %             hold on
    %             plot(human_struct.Frequency,ip)
    %             legend('human', 'simdata')
                count = count+1
                %                         clear main_ip
    
            end
        end
%     end
end
toc
save_error_file(input,subject_type,params,filename_rmse,rmse,human_struct)