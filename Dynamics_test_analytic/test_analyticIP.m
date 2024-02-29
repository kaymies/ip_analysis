% Set figure properties
set(0, 'DefaultLineLineWidth', 2);
set(groot,'defaultAxesFontSize',16);
set(0,'defaultfigurecolor',[1 1 1]); % white figure background
set(groot,'defaultAxesBox','on');    % box on
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaultTextInterpreter','latex');
%% Load Human Data
subject_type = 'duarte_old';
human_struct = load(sprintf('%s.mat',subject_type));
input.method = 'cpsd';
input.Frequency = eval(sprintf('human_struct.Frequency_%s',input.method));
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
input.Controller.beta = 1;
input.NoiseRatio = 1;
input.Controller.gamma = 1;
input.Controller.kappa = 1;
input.Controller.eta = 1;
input.Controller.type = 'LQR';
input.simSensorNoiseRatio = 1;
input.simMotorNoise = 10;
input.simSensorNoise = 0;
input.Controller.param.delay =0; % 2023-05-18
% % Simulation LQR controller parameters
% gamma = 1;
% alpha = 1e6; 
% beta = 0.2;
% motorNoiseRatio = 0.2; % sigma_r
% controller_param = struct('Q',gamma*eye(4),'R',alpha*diag([beta,1/beta]));
% 
% controller_param.delay = 0; 
% 
% params =...
%     struct(...
%     'TotalMass', 68,...  % in kg
%     'TotalHeight', 1.7,...% in meters
%     'gender', 'M',...   % 'F','M'
%     'plane', 'sgt',...  % 'sgt','frt'
%     'model','DIP',...   % 'SIP','DIP'
%     'pose','pose_I',... %
%     'FreqSampKin', 100,...
%     'measurements', cell(1),...
%     'simMotorNoiseRatio',motorNoiseRatio,...
%     'simSensorNoiseRatio',1,...
%     'simMotorNoise', 10,...
%     'simSensorNoise', 0,...
%     'Controller',...
%         struct('type','LQR',...
%                'param',controller_param)...
%     );
%%
figure;
beta = 0.3; sigma = 0.6;
input.Controller.beta = beta;
input.NoiseRatio = sigma;
[IP_ratio_a] = computeAnalyticIP(input);
hold on;p_ankle = plot(input.Frequency,IP_ratio_a,'Linewidth',1.2);
input.TotalHeight = 1.68;  
input.TotalMass = 62.6; 
[IP_ratio_a] = computeAnalyticIP(input);
hold on;p_ankle = plot(input.Frequency,IP_ratio_a,'Linewidth',1.2);

beta = 0.3; sigma = 0.8;
input.Controller.beta = beta;
input.NoiseRatio = sigma;
[IP_ratio_a] = computeAnalyticIP(input);
p_ankle = plot(input.Frequency,IP_ratio_a,'Linewidth',1.2);
beta = 2.8; sigma = 0.8;
input.Controller.beta = beta;
input.NoiseRatio = sigma;
[IP_ratio_a] = computeAnalyticIP(input);
p_ankle = plot(input.Frequency,IP_ratio_a,'Linewidth',1.2);
% beta = 2.2; sigma = 2.1;
% input.Controller.beta = beta;
% input.NoiseRatio = sigma;
% [IP_ratio_a] = computeAnalyticIP(input);
% p_ankle = plot(input.Frequency,IP_ratio_a,'Linewidth',1.2);
plot(input.Frequency,human_struct.IPDataAverage)
legend('b=0.2,s=0.01','b=0.3,s=0.8','b=2.8,s=0.8','Younger')

for sigma = flip([0,0.3,0.5,0.7,1,2,10])
params.simMotorNoiseRatio = sigma;
[IP_ratio_a] = computeAnalyticIP(input);
hold on;p_ankle = plot(input.Frequency,IP_ratio_a,'Linewidth',1.2);
end
yline(0,'k--');
yline(1,'k--');
legend('10','2','1','0.7','0.5','0.3','0')
ylim([0,3])
xlim([0,8])
xlabel('Frequency (Hz)');
ylabel('$IP_z/CoM$');
title('LQR model: effect of $\sigma_r$');


params.simMotorNoiseRatio = motorNoiseRatio;
%%
figure;
for beta = flip(0.1:0.01:0.2)
    input.Controller.beta = beta;
    [IP_ratio_a] = computeAnalyticIP(input);
    hold on;p_ankle = plot(input.Frequency,IP_ratio_a,'Linewidth',4);
end
yline(0,'k--');
yline(1,'k--');
% legend('3','2','1','0.5','0.3','0.1','0.01')
ylim([0,3])
xlim([0,8])
xlabel('Frequency (Hz)');
ylabel('IP as a Fraction of CoM');
title('LQR model: effect of \beta');
mycolors = flip(parula(11));
ax = gca; 
ax.ColorOrder = mycolors;
improvePlot;
%%
beta = flip([0.1,0.3,0.5,0.7,1,2,3,10]);
sigma = flip([0,0.3,0.5,0.7,1,2,10]);
IPsameCoM = zeros(size(beta,2),size(sigma,2));
figure;
hold on;
for b = 1:length(beta)
    for s = 1:length(sigma)
    params.Controller.param.R = alpha*diag([beta(b),1/beta(b)]);
    params.simMotorNoiseRatio = sigma(s);
    [IP_ratio_a] = computeAnalyticIP(params);
    IPsameCoM(b,s) = input.Frequency(find((IP_ratio_a<1) == 1, 1));
    end
end
for i = 1:size(IPsameCoM,2)
    plot(beta,IPsameCoM(:,i),'-o');
end
legend('$\sigma$ = 10','$\sigma$ = 2','$\sigma$ = 1','$\sigma$ = 0.7','$\sigma$ = 0.5','$\sigma$ = 0.3','$\sigma$ = 0')
xlabel('$\beta$');
ylabel('$IP_z = CoM$ Frequency (Hz)');
title('Effect of $\sigma$ on IP Cross Over Frequency');
%%
beta = flip([0.1,0.3,0.5,0.7,1,2,3,10]);
sigma = flip([0,0.3,0.5,0.7,1,2,10]);
IPhighfreq = zeros(size(beta,2),size(sigma,2));
figure;
hold on;
for b = 1:length(beta)
    for s = 1:length(sigma)
    params.Controller.param.R = alpha*diag([beta(b),1/beta(b)]);
    params.simMotorNoiseRatio = sigma(s);
    [f_a,IP_ratio_a] = computeAnalyticIP(params);
    IPhighfreq(b,s) = IP_ratio_a(end);
    end
end
for i = 1:size(IPhighfreq,1)
    plot(sigma,IPhighfreq(i,:),'-o');
end
legend('$\beta$ = 10','$\beta$ = 3','$\beta$ = 2','$\beta$ = 1','$\beta$ = 0.7','$\beta$ = 0.5','$\beta$ = 0.3','$\beta$ = 0.1')
xlabel('$\sigma$');
ylabel('$IP_z$ Highest Frequency Value');
title('Effect of $\beta$ on IP High Frequency Asymptote');