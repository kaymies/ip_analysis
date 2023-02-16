close all;
clear all;
clc;
setpath
%% Load Human Data
subject_type = 'nonpar_7sub';
human_struct = load(sprintf('%s.mat',subject_type));
% f_i = 0.65; f_int = 0.5; f_end = 5.15;
% input.Frequency = f_i:f_int:f_end;
input.Frequency = human_struct.Frequency;
input.FrequencyWindow = 0.2;
%% Set Model Parameters
input.TotalMass = 87.5; %71.1 Gruben
input.TotalHeight = 1.67; %1.75 Gruben    
input.gender = 'M';
input.plane = 'sgt';
input.model = 'DIP';
input.pose = 'pose_I';
input.FreqSampKin = 100;
input.trialDuration = 60;
input.CoordinateFrame = 'relative';
%% Set Controller Parameters
%-----% MAKE SURE TO CHANGE SIM FILENAME!!! %-----%
filename = 'En4_0p2_0p2'; %alpha_beta_sigmar
folder_name = 'LQR_relative_varyQ';
%-------------------------------------------------%
input.Controller.alpha = 1e-4;
input.Controller.beta = 0.2;
input.NoiseRatio = 0.2;
input.Controller.gamma = 5;
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
sim_struct = save_sim_file(input,filename,data,meansim,numtrial,folder_name);
%% Plot human vs. simulation
generate_plot(human_struct,sim_struct)
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