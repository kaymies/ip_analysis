close all;
clear all;
clc;
setpath
%% Load Human Data
subject_type = 'nonpar_7sub';
human_struct = load(sprintf('%s.mat',subject_type));
input.Frequency = human_struct.Frequency;
input.FrequencyWindow = 0.2;
%% Set Model Parameters
input.TotalMass = 87.5; %71.1 Gruben, %76.9 M Marta, 57 F Marta
input.TotalHeight = 1.67; %1.75 Gruben %1.72 M Marta, 1.61 F Marta   
input.gender = 'M';
input.plane = 'sgt';
input.model = 'DIP';
input.pose = 'pose_I';
input.FreqSampKin = 100;
input.trialDuration = 60;
input.CoordinateFrame = 'relative';
%% Set Controller Parameters
input.Controller.type = 'LQR';
%% Visualization
input.PostProc.AnimOn = 0;
input.PostProc.PlotOn = 0;
%% Create empty arrays to log errors
params.alpha = [1e-4,1e-3,1e-2,1,1e6]; %5
params.beta = [0.1,0.5,1,2]; %4
params.sigma_r = [0.01,0.5,1,2,5]; %5
params.gamma = [0.1,1,10,35]; %4
params.kappa = [1e-5,1e-4,1e-3,1e-2,1e-1,1,10,100,1e3,1e4,1e5]; %11
params.eta = [1e-1,1,10,1e2]; %4

alpha = params.alpha;
beta = params.beta;
sigma_r = params.sigma_r;
gamma = params.gamma;
kappa = params.kappa;
eta = params.eta;
% num_alpha = 11;
% num_alpha = 1;
% alpha = logspace(-4,6,num_alpha);
% alpha = 1e-4;
% beta = [0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2];
% beta = [5, 10];
% sigma_r = [0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,3,5,10,15,20];
% rmse_nonpar7sub_lqr_6params = zeros(length(alpha),length(beta),length(sigma_r),...
%     length(gamma),length(kappa),length(eta));\
load('C:\Users\kaymi\Dropbox (MIT)\Balance Project\Intersection Point Analysis\ip_analysis\Dynamics_test\Data\Error\rmse20230223_3.mat')
rmse_nonpar7sub_lqr_6params = Error;
%% Run Simulation for n trials across different parameter values
num_trial = 40;
count = 0;
tic
for a = 2%:length(alpha)
    for b = 3:length(beta)
        for s = 4:length(sigma_r)
            for g = 4:length(gamma)
                for k = 6:length(kappa)
                    for e = 1:length(eta)
                        input.Controller.alpha = alpha(a);
                        input.Controller.beta = beta(b);
                        input.NoiseRatio = sigma_r(s);
                        input.Controller.gamma = gamma(g);
                        input.Controller.kappa = kappa(k);
                        input.Controller.eta = eta(e);
                        
                        ip = compute_ip(input,num_trial);
                        ip_filtered = neglect_outliers(ip);
                        meansim = mean(ip_filtered,'omitnan');
                        rmse_nonpar7sub_lqr_6params(a,b,s,g,k,e) = sqrt(mean((meansim-human_struct.IPDataAverage).^2));
            %             figure(1);
            %             plot(f,human)
            %             hold on
            %             plot(f,meansim)
            %             legend('human', 'simdata')
                        count = count+1
%                         clear main_ip
                    end
                end
            end
        end
    end
end
toc
save_error_file(input,subject_type,num_trial,params,'rmse20230227_2',rmse_nonpar7sub_lqr_6params)