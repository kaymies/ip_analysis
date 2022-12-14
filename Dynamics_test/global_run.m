close all;
clear all;
clc;
%% Set Model Parameters
input.TotalMass = 71.1; %71.1 Gruben, %76.9 M Marta, 57 F Marta
input.TotalHeight = 1.75; %1.75 Gruben %1.72 M Marta, 1.61 F Marta   
input.gender = 'M';
input.plane = 'sgt';
input.model = 'DIP';
input.pose = 'pose_I';
input.FreqSampKin = 100;
input.trialDuration = 15;
input.CoordinateFrame = 'spatial';
population = 'nonpar_7sub';
%% Set Controller Parameters
input.Controller.gamma = 1;
input.Controller.type = 'LQR';
%% Visualization
input.PostProc.AnimOn = 0;
input.PostProc.PlotOn = 0;
%% Load Human Data
freq = load('freq.mat'); f = freq.freq;
human = load('human.mat');
human = eval(sprintf('human.%s',population))'; 
%% Create empty arrays to log errors
% alpha = [1e-4,1e6];
% beta = [0.3,2];
% sigma_r = [0.9,10];

% num_alpha = 11;
num_alpha = 1;
% alpha = logspace(-4,6,num_alpha);
alpha = 1e-4;
% beta = [0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2];
beta = [5, 10];
sigma_r = [0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,3,5,10,15,20];

% num_alpha = 1;
% alpha = 1e-4;
% beta = 0.2;
% sigma_r = 0.2;
root_mean_nonpar7sub_lqr_spatial2 = zeros(length(alpha), length(beta), length(sigma_r));
% diff = zeros(num_alpha, num_beta, num_sigma);
%% Run Simulation for n trials across different parameter values
num_trial = 40;
count = 1;
tic
for a = 1:length(alpha)
    for b = 1:length(beta)
        for s = 1:length(sigma_r)
            input.Controller.alpha = alpha(a);
            input.Controller.beta = beta(b);
            input.NoiseRatio = sigma_r(s);
            [~,ip] = main_ip(input);
            data = zeros(num_trial,length(ip));
            data(1,:) = ip;
            for i = 1:num_trial-1
                [~,ip] = main_ip(input);
                data(i+1,:) = ip;
%                 i
            end
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
%             diff(a,b,s) = mean(mean(data-human,2,'omitnan'));
%             root_mean(a,b,s) = mean(sqrt(mean((data - human).^2,2,'omitnan')));
            meansim = mean(data,'omitnan');
            root_mean_nonpar7sub_lqr_spatial2(a,b,s) = sqrt(mean((meansim-human).^2));
%             figure(1);
%             plot(f,human)
%             hold on
%             plot(f,meansim)
%             legend('human', 'simdata')
            count = count+1
            clear main_ip
        end
    end
end
toc
save('error','root_mean_nonpar7sub_lqr_spatial2','-append');