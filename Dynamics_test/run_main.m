close all;
clear all;
clc;
setpath
%% Load Human Data
subject_type = 'old_dom';
human_struct = load(sprintf('%s.mat',subject_type));
% f_i = 0.65; f_int = 0.5; f_end = 5.15;
% input.Frequency = f_i:f_int:f_end;
input.Frequency = human_struct.Frequency;
input.FrequencyWindow = 0.2;
%% Set Model Parameters
input.TotalMass = 71.1; %71.1 Gruben
input.TotalHeight = 1.75; %1.75 Gruben    
input.gender = 'M';
input.plane = 'sgt';
input.model = 'DIP';
input.pose = 'pose_I';
input.FreqSampKin = 100;
input.trialDuration = 60;
input.CoordinateFrame = 'relative';
%% Set Controller Parameters
%-----% MAKE SURE TO CHANGE SIM FILENAME!!! %-----%
filename = 'E4_1p9_3'; %alpha_beta_sigmar
%-------------------------------------------------%
input.Controller.alpha = 1e4;
input.Controller.beta = 1.9;
input.NoiseRatio = 4;
input.Controller.gamma = 1;
input.Controller.type = 'LQR';
%% Visualization
input.PostProc.AnimOn = 0;
input.PostProc.PlotOn = 0;
%% Run Simulation for n trials
numtrial = 40;

data = zeros(numtrial,length(input.Frequency));
% torques = zeros(numtrial, 3);
% torques(1,:) = torque_rms;
for i = 1:numtrial
    close all;
    [IP_ratio] = main_ip(input);
    data(i,:) = IP_ratio;
%     torques(i,:) = torque_rms;
    i
end
figure; boxplot(data,input.Frequency)
xlabel('Frequency [Hz]')
ylabel('IP (Fraction of CoM)')
ylim([0 2.5])
% title('gamma = 1e6, beta = 2, t = 30 sec, freq = 100 Hz, 0.1x noise')
set(gca,'Fontsize',15);

% mean_torque_ratio = mean(torques)
% std_torque_ratio = std(torques)
%% Neglect Outliers
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
meansim = mean(data,'omitnan'); % mean of simulation without outliers
%% Save File:
sim_struct = save_sim_file(input,filename,data,meansim,numtrial);
%% Plot human vs. simulation
generate_plot(human_struct,sim_struct)