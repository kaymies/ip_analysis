function [] = generate_plot(human,sim)
%% Load Data
% Human Data
setpath;
human_struct = human;
% Best-Fit Simulation Data
sim_struct = sim;
%% Plot human vs. simulation
figure();
plot(human_struct.Frequency,human_struct.IPDataAverage,'x')
hold on
errorbar(sim_struct.Frequency,sim_struct.DataWithoutOutliers,sim_struct.StandardDeviation,'o')
% plot(sim_struct.Frequency,sim_struct.DataWithoutOutliers)
legend('human', 'simdata')
xlabel('Frequency [Hz]')
ylabel('IP (Fraction of CoM)')
xlim([0 6.3])
ylim([0 3.5])
end
%% Plot Small alpha vs Large alpha RMSE
% rmses = [0.188 0.065;
%         0.166 0.067;
%         0.104 0.067;
%         0.122 0.051]';
% x = 1:4;
% xl = ["Sagittal Ground" "Sagittal Beam" "Frontal Ground" "Frontal Beam"];
% bar(condition,rmses)
% legend('\alpha = 10^{-4}', '\alpha = 10^6')
% xlabel('Balance Condition')
% ylabel('RMSE')
% set(gca,'XTickLabel',xl)
% improvePlot;
%% Plot with standard deviation
% old = load('duarte_old.mat');
% young = load('duarte_young.mat');
% mean = [old.IPDataAverage; young.IPDataAverage];
% std = [old.StandardDeviation; young.StandardDeviation];
% high = zeros(size(mean));
% low = zeros(size(mean));
% for i = 1:height(mean)
%     high(i,:) = [mean(i,:) + std(i,:)];
%     low(i,:) = [mean(i,:) - std(i,:)];
% end
% freq = old.Frequency;
% x = [freq, fliplr(freq)];
% inBetween = [high(1,:), fliplr(low(1,:))];
% h = fill(x, inBetween, [255/255 0/255 0/255], 'LineStyle','none');
% set(h,'facealpha',.5)
% hold on;
% plot(freq, mean(1,:), 'color', [255/255 0/255 0/255], 'LineWidth', 2);
% inBetween2 = [high(2,:), fliplr(low(2,:))];
% h2 = fill(x, inBetween2, [0 255/255 0/255], 'LineStyle','none');
% set(h2,'facealpha',.5)
% plot(freq, mean(2,:), 'color', [0 255/255 0/255], 'LineWidth', 2);
% % plot(freq,young.IPDataAverage,'color', [0 0 255/255], 'LineWidth', 2);
% legend('','Old','','Young')
% xlabel('Frequency (Hz)')
% ylabel('IP (Fraction of CoM)')

%% Plot each human
% old = load('duarte_old.mat');
% young = load('duarte_young.mat');
% for i = 1:3:old.NumSubjects*3
%     ave = mean(old.IP(i:i+2,:),'omitnan');
%     plot(old.Frequency,ave)
%     hold on
%     xlabel('Frequency (Hz)')
%     ylabel('IP (Fraction of CoM)')
% end
% figure();
% for i = 1:3:young.NumSubjects*3
%     ave = mean(young.IP(i:i+2,:),'omitnan');
%     plot(young.Frequency,ave)
%     hold on
%     xlabel('Frequency (Hz)')
%     ylabel('IP (Fraction of CoM)')
% end
%% Plot Human vs. Sim (Multiple Conditions)
% setpath;
% figure;
% % % Sagittal Side-by-Side
% % human_struct = load('gruben2019.mat');
% % sim_struct = load('E6_0p3_0p9.mat');
% % subplot(2, 3, 1);
% % plot(human_struct.Frequency, human_struct.IPDataAverage,'x');
% % hold on
% % errorbar(sim_struct.Frequency,sim_struct.DataWithoutOutliers,sim_struct.StandardDeviation,'o')
% % xlabel('Frequency [Hz]')
% % ylabel('IP (Fraction of CoM)')
% % xlim([0 8])
% % ylim([0 2.5])
% % % title('Plot 1');
% 
% % Sagittal Ground
% human_struct = load('marta_sgt_ground.mat');
% sim_struct = load('E6_0p3_1p6.mat');
% freq = human_struct.Frequency;
% subplot(2, 2, 1);
% plot(freq, human_struct.IPDataAverage,'x');
% hold on
% errorbar(sim_struct.Frequency,sim_struct.DataWithoutOutliers,sim_struct.StandardDeviation,'o')
% xlabel('Frequency [Hz]')
% ylabel('IP (Fraction of CoM)')
% xlim([0 8])
% ylim([0 2.5])
% % title('Plot 2');
% 
% % Sagittal Beam
% human_struct = load('marta_sgt_beam.mat');
% sim_struct = load('E6_0p2_1p1.mat');
% subplot(2, 2, 2);
% plot(freq, human_struct.IPDataAverage,'x');
% hold on
% errorbar(sim_struct.Frequency,sim_struct.DataWithoutOutliers,sim_struct.StandardDeviation,'o')
% xlabel('Frequency [Hz]')
% ylabel('IP (Fraction of CoM)')
% xlim([0 8])
% ylim([0 2.5])
% % title('Plot 3');
% 
% % Frontal Ground
% human_struct = load('marta_frt_ground.mat');
% sim_struct = load('E6_0p2_0p1.mat');
% subplot(2, 2, 3);
% plot(freq, human_struct.IPDataAverage,'x');
% hold on
% errorbar(sim_struct.Frequency,sim_struct.DataWithoutOutliers,sim_struct.StandardDeviation,'o')
% xlabel('Frequency [Hz]')
% ylabel({'Intersection Point Height'; '(Normalized by CoM)'})
% xlim([0 8])
% ylim([0 2.5])
% % title('Plot 4');
% 
% % Frontal Beam
% subplot(2, 2, 4);
% human_struct = load('marta_frt_beam.mat');
% sim_struct = load('E6_1p2_0p5.mat');
% plot(freq, human_struct.IPDataAverage,'x');
% hold on
% errorbar(sim_struct.Frequency,sim_struct.DataWithoutOutliers,sim_struct.StandardDeviation,'o')
% 
% legend('Human Data', 'Simulated Data')
% xlabel('Frequency [Hz]')
% % ylabel('IP (Fraction of CoM)')
% xlim([0 8])
% ylim([0 2.5])
% improvePlot;
%% Plot Human (Multiple Conditions)
setpath;
figure;
% % Sagittal Side-by-Side
% human_struct = load('gruben2019.mat');
% sim_struct = load('E6_0p3_0p9.mat');
% subplot(2, 3, 1);
% plot(human_struct.Frequency, human_struct.IPDataAverage,'x');
% hold on
% errorbar(sim_struct.Frequency,sim_struct.DataWithoutOutliers,sim_struct.StandardDeviation,'o')
% xlabel('Frequency [Hz]')
% ylabel('IP (Fraction of CoM)')
% xlim([0 8])
% ylim([0 2.5])
% % title('Plot 1');

% Sagittal Ground
human_struct = load('marta_sgt_ground.mat');
freq = human_struct.Frequency;
subplot(2, 2, 1);
errorbar(freq, human_struct.IPDataAverage,human_struct.StandardDeviation,'x');
ylabel({'Intersection Point Height'; '(Normalized by CoM)'})
xlim([0 8])
ylim([0 2.5])

% Sagittal Beam
human_struct = load('marta_sgt_beam.mat');
subplot(2, 2, 2);
errorbar(freq, human_struct.IPDataAverage,human_struct.StandardDeviation,'x');
xlim([0 8])
ylim([0 2.5])

% Frontal Ground
human_struct = load('marta_frt_ground.mat');
subplot(2, 2, 3);
errorbar(freq, human_struct.IPDataAverage,human_struct.StandardDeviation,'x');
xlabel('Frequency [Hz]')
ylabel({'Intersection Point Height'; '(Normalized by CoM)'})
xlim([0 8])
ylim([0 2.5])

% Frontal Beam
subplot(2, 2, 4);
human_struct = load('marta_frt_beam.mat');
errorbar(freq, human_struct.IPDataAverage,human_struct.StandardDeviation,'x');
xlabel('Frequency [Hz]')
% ylabel('IP (Fraction of CoM)')
xlim([0 8])
ylim([0 2.5])
improvePlot;

%% Plot Global Param Search (Multiple Conditions)
minColorLimit = 0.051;                   % determine colorbar limits from data
maxColorLimit = 0.5;
crange = [minColorLimit,maxColorLimit];
fig=figure(1);
% subplot(2,3,1,'Parent',fig);
% filename = 'rmse_gruben2019_20230615';
% error = load(sprintf('%s.mat',filename));
% beta = error.Parameters.beta;
% sigma = error.Parameters.sigma_r;
% error = error.Error;
% heatmap(round(sigma,2),round(beta,2),abs(error),'Colormap',parula,'ColorLimits',crange)
% colorbar off

subplot(2,2,1,'Parent',fig);
filename = 'rmse_marta_sgt_ground_20231025';
error = load(sprintf('%s.mat',filename));
beta = error.Parameters.beta(3:end);
sigma = error.Parameters.sigma_r(3:end-4);
error = error.Error(3:end,3:end-4);
h = heatmap(round(sigma,2),flip(round(beta,2)),abs(flip(error)),'Colormap',parula,'ColorLimits',crange)
h.YDisplayLabels = {'3.0','','','','','2.5','','','','','2.0','','','','','1.5','','','','','1.0','','','','','0.5','','','','0.1'};
h.XDisplayLabels = {'0.1','','','','0.5','','','','','1.0','','','','','1.5','','','','','2.0','','','','','2.5','','','','','3.0'};
colorbar off

subplot(2,2,2,'Parent',fig);
filename = 'rmse_marta_sgt_beam_20231025';
error = load(sprintf('%s.mat',filename));
error = error.Error(3:end,3:end-4);
h = heatmap(round(sigma,2),flip(round(beta,2)),abs(flip(error)),'Colormap',parula,'ColorLimits',crange)
h.YDisplayLabels = {'3.0','','','','','2.5','','','','','2.0','','','','','1.5','','','','','1.0','','','','','0.5','','','','0.1'};
h.XDisplayLabels = {'0.1','','','','0.5','','','','','1.0','','','','','1.5','','','','','2.0','','','','','2.5','','','','','3.0'};
colorbar off

subplot(2,2,3,'Parent',fig);
filename = 'rmse_marta_frt_ground_20231025';
error = load(sprintf('%s.mat',filename));
error = error.Error(3:end,3:end-4);
h = heatmap(round(sigma,2),flip(round(beta,2)),abs(flip(error)),'Colormap',parula,'ColorLimits',crange)
h.YDisplayLabels = {'3.0','','','','','2.5','','','','','2.0','','','','','1.5','','','','','1.0','','','','','0.5','','','','0.1'};
h.XDisplayLabels = {'0.1','','','','0.5','','','','','1.0','','','','','1.5','','','','','2.0','','','','','2.5','','','','','3.0'};
colorbar off

subplot(2,2,4,'Parent',fig);
filename = 'rmse_marta_frt_beam_20231025';
error = load(sprintf('%s.mat',filename));
error = error.Error(3:end,3:end-4);
h = heatmap(round(sigma,2),flip(round(beta,2)),abs(flip(error)),'Colormap',parula,'ColorLimits',crange)
h.YDisplayLabels = {'3.0','','','','','2.5','','','','','2.0','','','','','1.5','','','','','1.0','','','','','0.5','','','','0.1'};
h.XDisplayLabels = {'0.1','','','','0.5','','','','','1.0','','','','','1.5','','','','','2.0','','','','','2.5','','','','','3.0'};
colorbar off

% for i=1:6
%     if i == 5
%     else
% 
%     sph{i} = subplot(2,3,i,'Parent',fig);
%     contourf(sph{i},x,y,i.*z(:,:,i))                % changes due to illustration purposes
%     caxis(sph{i},[minColorLimit,maxColorLimit]);    % set colorbar limits
% end
h = axes(fig,'visible','off'); 
h.XLabel.Visible = 'on';
h.YLabel.Visible = 'on';
% ylabel(h,'\beta','FontWeight','bold');
% xlabel(h,'\sigma_r','FontWeight','bold');
c = colorbar(h,'Position',[0.93 0.168 0.022 0.7]);  % attach colorbar to h
colormap(c,'Parula')
caxis(h,crange);             % set colorbar limits