setpath;

human_struct = load("duarte_young.mat");
bestparams = load('subs_bestparams_duarte_young_20230627.mat');
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
outliers = [6 9 19 46 61];
outliers_close = [6 46];
params_close = [0.3 0.3 ;
                0.4 0.6];

% gridRows = 3;
% gridCols = 4;
% % Create a figure
% figure;
% o_c = 1;
% % Iterate over each subplot
% for i = 1:12
%     % Select the current subplot
%     subplot(gridRows, gridCols, i);
%     % Plot the data in the current subplot
%     input.Controller.beta = bestparams.BestBetas(i);
%     input.NoiseRatio = bestparams.BestSigmas(i);
%     if ismember(i,outliers)
%         plot(input.Frequency, compute_ip(input),'r-',LineWidth=2);
%         hold on
%         plot(human_struct.Frequency,human_struct.IPAveSubject(outliers(i),:),'b-',LineWidth=1.5)
%         if ismember(i,outliers_close)
%             input.Controller.beta = params_close(1,o_c);
%             input.NoiseRatio = params_close(2,o_c);
%             o_c = o_c + 1;
%             plot(input.Frequency,compute_ip(input),'g-',LineWidth=2)
%             if i == 12
%                 legend('Fitted (Analytic)','Raw (Duarte)','Best Fit with Small \beta')
%             end
%         end
%     else
%         plot(input.Frequency, compute_ip(input),'r--',LineWidth=2);
%         hold on
%         plot(human_struct.Frequency,human_struct.IPAveSubject(i,:),'b--',LineWidth=1.5)
%     end
%     % Add labels and title to the subplot
%     xlabel('Frequency [Hz]')
%     ylabel('zIP/zCoM')
%     title(['Older Sub', num2str(i)]);
%     ylim([0 4])
% end

gridRows = 2;
gridCols = 3;
% Create a figure
figure;
o_c = 1;
% Iterate over each subplot
for i = 1:5
    % Select the current subplot
    subplot(gridRows, gridCols, i);
    % Plot the data in the current subplot
    input.Controller.beta = bestparams.BestBetas(outliers(i));
    input.NoiseRatio = bestparams.BestSigmas(outliers(i));
    plot(input.Frequency, compute_ip(input),'-','Color','#648FFF',LineWidth=2);
    hold on
    plot(human_struct.Frequency,human_struct.IPAveSubject(outliers(i),:),'-','Color','#000000',LineWidth=1.5)
    if ismember(outliers(i),outliers_close)
        input.Controller.beta = params_close(1,o_c);
        input.NoiseRatio = params_close(2,o_c);
        o_c = o_c + 1;
        plot(input.Frequency,compute_ip(input),'--','Color','#DC267F',LineWidth=2)
        if i == 4
            legend('Fitted (Analytic)','Raw (Duarte)','Best Fit with Small \beta')
        end
    end
    % Add labels and title to the subplot
    xlabel('Frequency [Hz]')
    ylabel('zIP/zCoM')
    title(['Younger Sub', num2str(outliers(i))]);
    ylim([0 4])
end
% legend('Fitted (Analytic)','Raw (Duarte)')

% IP_reshape = zeros(38,3,29);
% for col = 1:38
%     % Extract every three rows from the original matrix
%     extractedRows = IP((col-1)*3+1:col*3, :);
%     
%     % Assign the extracted rows as a column in the reshaped matrix
%     IP_reshape(col, :, :) = extractedRows;
% end