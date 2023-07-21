setpath;

human_struct = load("duarte_old.mat");
bestparams = load('subs_bestparams_duarte_old_20230705.mat');
numSubplots = human_struct.NumSubjects;
input.Frequency = human_struct.Frequency;
input.TotalMass = human_struct.MeanMass_kg;
input.TotalHeight = human_struct.MeanHeight_m;  
input.gender = 'M';
input.plane = human_struct.Plane;
input.model = 'DIP';
input.pose = human_struct.Pose;
input.Controller.alpha = 1e6;
input.Controller.gamma = 1;
input.Controller.kappa = 1;
input.Controller.eta = 1;
input.Controller.type = 'LQR';
input.simSensorNoiseRatio = 1;
input.simMotorNoise = 10;
input.simSensorNoise = 0;
input.Controller.param.delay =0;

% [betas,TFrm1] = rmoutliers(log(bestparams.BestBetas),'percentiles',[5 95]);
% [sigmas,TFrm2] = rmoutliers(log(bestparams.BestSigmas),'percentiles',[5 95]);
% TFrm = logical(TFrm1+TFrm2);
outliers = [1 7 8 9 16 20 24 26 28 35 36 37];

gridRows = 5;
gridCols = 8;
% Create a figure
figure;
% Iterate over each subplot
for i = 1:numSubplots
    % Select the current subplot
    subplot(gridRows, gridCols, i);
    % Plot the data in the current subplot
    input.Controller.beta = bestparams.BestBetas(i);
    input.NoiseRatio = bestparams.BestSigmas(i);
    if i == outliers
        plot(input.Frequency, compute_ip(input),'r-',LineWidth=2);
        hold on
        plot(human_struct.Frequency,human_struct.IPAveSubject(i,:),'b-',LineWidth=1.5)
    else
        plot(input.Frequency, compute_ip(input),'r--',LineWidth=2);
        hold on
        plot(human_struct.Frequency,human_struct.IPAveSubject(i,:),'b--',LineWidth=1.5)
    end
    % Add labels and title to the subplot
    xlabel('Frequency [Hz]')
    ylabel('zIP/zCoM')
    title(['Older Sub', num2str(i)]);
    ylim([0 4])
end
legend('Fitted (Analytic)','Raw (Duarte)')

% IP_reshape = zeros(38,3,29);
% for col = 1:38
%     % Extract every three rows from the original matrix
%     extractedRows = IP((col-1)*3+1:col*3, :);
%     
%     % Assign the extracted rows as a column in the reshaped matrix
%     IP_reshape(col, :, :) = extractedRows;
% end