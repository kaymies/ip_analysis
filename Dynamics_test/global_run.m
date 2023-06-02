close all;
clear all;
clc;
setpath
%% Load Human Data
subject_type = 'marta_sgt_ground';
human_struct = load(sprintf('%s.mat',subject_type));
input.Frequency = human_struct.Frequency;
input.FrequencyWindow = 0.2;
%% Set Model Parameters
input.TotalMass = human_struct.MeanMass_kg; %71.1 Gruben, %76.9 M Marta, 57 F Marta
input.TotalHeight = human_struct.MeanHeight_m; %1.75 Gruben %1.72 M Marta, 1.61 F Marta   
input.gender = 'M';
input.plane = 'sgt';
input.model = 'DIP';
input.pose = 'pose_T';
input.FreqSampKin = 500;
input.trialDuration = 30;
input.CoordinateFrame = 'relative';
%% Set Controller Parameters
input.Controller.type = 'LQR';
%% Visualization
input.PostProc.AnimOn = 0;
input.PostProc.PlotOn = 0;
%% Create empty arrays to log errors
params.beta = [0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2];
params.sigma_r = [0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,3,5,10,15,20];
% params.alpha = [1e-4,1e-3,1e-2,1,1e6]; %5
% params.beta = [0.3]; %4
% params.sigma_r = [1.6]; %5
% params.gamma = [0.1,1,10,35]; %4
% params.kappa = [1e-5,1e-4,1e-3,1e-2,1e-1,1,10,100,1e3,1e4,1e5]; %11
% params.eta = [1e-1,1,10,1e2]; %4
% 
% alpha = params.alpha;
% beta = params.beta;
% sigma_r = params.sigma_r;
% gamma = params.gamma;
% kappa = params.kappa;
% eta = params.eta;
% num_alpha = 11;
% alpha = logspace(-4,6,num_alpha);
% beta = [0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2];
% sigma_r = [0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,3,5,10,15,20];
% rmse_nonpar7sub_lqr_6params = zeros(length(alpha),length(beta),length(sigma_r),...
%     length(gamma),length(kappa),length(eta));\
% load('C:\Users\kaymi\Dropbox (MIT)\Balance Project\Intersection Point Analysis\ip_analysis\Dynamics_test\Data\Error\rmse20230223_3.mat')
% rmse_nonpar7sub_lqr_6params = Error;
rmse = zeros(length(params.beta),length(params.sigma_r));
diff = zeros(length(params.beta),length(params.sigma_r));
diff_lowfreq = zeros(length(params.beta),length(params.sigma_r));
diff_highfreq = zeros(length(params.beta),length(params.sigma_r));
%% Run Simulation for n trials across different parameter values
num_trial = 40;
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

        [ip,~] = compute_ip(input,num_trial);
        ip_filtered = neglect_outliers(ip);
        meansim = mean(ip_filtered,'omitnan');
        rmse(b,s) = sqrt(mean((meansim-human_struct.IPDataAverage).^2));
        diff(b,s) = mean(mean(ip_filtered - human_struct.IPDataAverage,'omitnan'),'omitnan');
        diff_lowfreq(b,s) = mean(mean(ip_filtered(:,5:11) - human_struct.IPDataAverage(:,5:11),'omitnan'),'omitnan');
        diff_highfreq(b,s) = mean(mean(ip_filtered(:,16:38) - human_struct.IPDataAverage(:,16:38),'omitnan'),'omitnan');
%         figure(1);
%         plot(human_struct.Frequency,human_struct.IPDataAverage)
%         hold on
%         plot(human_struct.Frequency,meansim)
%         legend('human', 'simdata')
        count = count+1
        %                         clear main_ip

    end
end
toc
save_error_file(input,subject_type,num_trial,params,'rmse_martasgtground_20230601',rmse)
save_error_file(input,subject_type,num_trial,params,'diff_martasgtground_20230601',diff)
save_error_file(input,subject_type,num_trial,params,'diff_lowfreq_martasgtground_20230601',diff_lowfreq)
save_error_file(input,subject_type,num_trial,params,'diff_highfreq_martasgtground_20230601',diff_highfreq)