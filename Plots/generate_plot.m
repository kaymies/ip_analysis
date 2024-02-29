%% Load Data
% Human Data
setpath;
% human_struct = human;
% Best-Fit Simulation Data
% sim_struct = sim;
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
%% Plot with standard deviation
setpath;
old = load('duarte_old.mat');
young = load('duarte_young.mat');
mean = [old.IPDataAverage_cpsd; young.IPDataAverage_cpsd];
stdev = [std(old.IPAveSubject_cpsd); std(young.IPAveSubject_cpsd)];
high = zeros(size(mean));
low = zeros(size(mean));
for i = 1:height(mean)
    high(i,:) = [mean(i,:) + stdev(i,:)];
    low(i,:) = [mean(i,:) - stdev(i,:)];
end
freq = old.Frequency_cpsd;
x = [freq, fliplr(freq)];
inBetween = [high(1,:), fliplr(low(1,:))];
h = fill(x, inBetween, [255/255 0/255 0/255], 'LineStyle','none');
set(h,'facealpha',.5)
hold on;
plot(freq, mean(1,:), 'color', [255/255 0/255 0/255], 'LineWidth', 2);
inBetween2 = [high(2,:), fliplr(low(2,:))];
h2 = fill(x, inBetween2, [0 255/255 0/255], 'LineStyle','none');
set(h2,'facealpha',.5)
plot(freq, mean(2,:), 'color', [0 255/255 0/255], 'LineWidth', 2);
% plot(freq,young.IPDataAverage,'color', [0 0 255/255], 'LineWidth', 2);
legend('','Older subjects (n = 38)','','Younger Subjects (n = 65)')
xlabel('Frequency (Hz)')
ylabel({'Intersection Point Height'; '(Normalized by CoM)'})
yline(1);
improvePlot;
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
subplot(2, 2, 2);
errorbar(freq, human_struct.IPDataAverage,human_struct.StandardDeviation,'x',Color,'k');
xlim([0 8])
ylim([-0.2 3.5])

% Sagittal Beam
human_struct = load('marta_sgt_beam.mat');
subplot(2, 2, 1);
errorbar(freq, human_struct.IPDataAverage,human_struct.StandardDeviation,'x');
ylabel({'Intersection Point Height'; '(Normalized by CoM)'})
xlim([0 8])
ylim([-0.2 3.5])

% Frontal Ground
human_struct = load('marta_frt_ground.mat');
subplot(2, 2, 4);
errorbar(freq, human_struct.IPDataAverage,human_struct.StandardDeviation,'x');
xlabel('Frequency [Hz]')
xlim([0 8])
ylim([-0.2 3.5])

% Frontal Beam
subplot(2, 2, 3);
human_struct = load('marta_frt_beam.mat');
errorbar(freq, human_struct.IPDataAverage,human_struct.StandardDeviation,'x');
xlabel('Frequency [Hz]')
ylabel({'Intersection Point Height'; '(Normalized by CoM)'})
xlim([0 8])
ylim([-0.2 3.5])
improvePlot;
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

% Condition 1
human_struct = load('duarte_old.mat');
sim_struct = load('E6_0p17_0p01.mat');
freq = human_struct.Frequency_cpsd;
subplot(1, 2, 1);
% plot(freq,human_struct.IPDataAverage,'x')
errorbar(freq, human_struct.IPDataAverage_cpsd,std(human_struct.IPAveSubject_cpsd),'x');
hold on
errorbar(sim_struct.Frequency,sim_struct.DataWithoutOutliers,sim_struct.StandardDeviation,'o')
ylabel({'Intersection Point Height'; '(Normalized by CoM)'})
xlabel('Frequency [Hz]')
xlim([0 8])
ylim([0 3.5])

% Condition 2
human_struct = load('duarte_young.mat');
sim_struct = load('E6_0p25_0p5.mat');
subplot(1, 2, 2);
% plot(freq,human_struct.IPDataAverage,'x')
errorbar(freq, human_struct.IPDataAverage_cpsd,std(human_struct.IPAveSubject_cpsd),'x');
hold on
errorbar(sim_struct.Frequency,sim_struct.DataWithoutOutliers,sim_struct.StandardDeviation,'o')
xlabel('Frequency [Hz]')
xlim([0 8])
ylim([0 3.5])

% % Condition 3
% human_struct = load('marta_frt_beam.mat');
% sim_struct = load('E6_1p1_0p4.mat');
% subplot(2, 2, 3);
% plot(freq, human_struct.IPDataAverage,'x');
% hold on
% errorbar(sim_struct.Frequency,sim_struct.DataWithoutOutliers,sim_struct.StandardDeviation,'o')
% ylabel({'Intersection Point Height'; '(Normalized by CoM)'})
% xlabel('Frequency [Hz]')
% xlim([0 8])
% ylim([0 2.5])
% 
% % Condition 4
% subplot(2, 2, 4);
% human_struct = load('marta_frt_ground.mat');
% sim_struct = load('E6_0p2_0p01.mat');
% plot(freq, human_struct.IPDataAverage,'x');
% hold on
% errorbar(sim_struct.Frequency,sim_struct.DataWithoutOutliers,sim_struct.StandardDeviation,'o')
% xlabel('Frequency [Hz]')
% xlim([0 8])
% ylim([0 2.5])

legend('Human Data', 'Simulated Data')
improvePlot;
%% Plot Global Param Search (Multiple Conditions)
setpath;
minColorLimit = 0.1279;                   % determine colorbar limits from data
maxColorLimit = 0.28;
crange = [minColorLimit,maxColorLimit];
fig=figure();
% subplot(2,3,1,'Parent',fig);
% filename = 'rmse_gruben2019_20230615';
% error = load(sprintf('%s.mat',filename));
% beta = error.Parameters.beta;
% sigma = error.Parameters.sigma_r;
% error = error.Error;
% heatmap(round(sigma,2),round(beta,2),abs(error),'Colormap',parula,'ColorLimits',crange)
% colorbar off

subplot(2,2,1,'Parent',fig);
filename = 'rmse_duarte_old_20240212';
error = load(sprintf('%s.mat',filename));
beta = error.Parameters.beta(3:end);
sigma = error.Parameters.sigma_r(3:end-4);
error = error.Error(3:end,3:end-4);
h = heatmap(round(sigma,2),flip(round(beta,2)),abs(flip(error)),'Colormap',parula,'ColorLimits',crange)
h.YDisplayLabels = {'3.0','','','','','2.5','','','','','2.0','','','','','1.5','','','','','1.0','','','','','0.5','','','','0.1'};
h.XDisplayLabels = {'0.1','','','','0.5','','','','','1.0','','','','','1.5','','','','','2.0','','','','','2.5','','','','','3.0'};
% h.YDisplayLabels = flip({'0.15','','','','','0.20','','','','','0.25','','','','','0.30','','','','','0.35','','','','','0.40','','','','','0.45'});
% h.XDisplayLabels = {'0.01','0.05','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0'};
colorbar off

subplot(2,2,3,'Parent',fig);
filename = 'rmse_duarte_old_finer_20240212';
error = load(sprintf('%s.mat',filename));
beta = error.Parameters.beta;
sigma = error.Parameters.sigma_r;
error = error.Error;
h = heatmap(round(sigma,2),flip(round(beta,2)),abs(flip(error)),'Colormap',parula,'ColorLimits',crange)
h.YDisplayLabels = {'0.15','','','','','0.20','','','','','0.25','','','','','0.30','','','','','0.35','','','','','0.40','','','','','0.45'};
h.XDisplayLabels = {'0.01','0.05','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0'};
colorbar off

subplot(2,2,2,'Parent',fig);
filename = 'rmse_duarte_young_20240212';
error = load(sprintf('%s.mat',filename));
beta = error.Parameters.beta(3:end);
sigma = error.Parameters.sigma_r(3:end-4);
error = error.Error(3:end,3:end-4);
h = heatmap(round(sigma,2),flip(round(beta,2)),abs(flip(error)),'Colormap',parula,'ColorLimits',crange)
h.YDisplayLabels = {'3.0','','','','','2.5','','','','','2.0','','','','','1.5','','','','','1.0','','','','','0.5','','','','0.1'};
h.XDisplayLabels = {'0.1','','','','0.5','','','','','1.0','','','','','1.5','','','','','2.0','','','','','2.5','','','','','3.0'};
colorbar off

subplot(2,2,4,'Parent',fig);
filename = 'rmse_duarte_young_finer_20240212';
error = load(sprintf('%s.mat',filename));
beta = error.Parameters.beta;
sigma = error.Parameters.sigma_r;
error = error.Error;
h = heatmap(round(sigma,2),flip(round(beta,2)),abs(flip(error)),'Colormap',parula,'ColorLimits',crange)
h.YDisplayLabels = {'0.15','','','','','0.20','','','','','0.25','','','','','0.30','','','','','0.35','','','','','0.40','','','','','0.45'};
h.XDisplayLabels = {'0.01','0.05','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0'};
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
improvePlot;
%% Plot Small alpha vs Large alpha RMSE
rmses = [0.153 0.068;
         0.128 0.051;
         0.175 0.072;
         0.088 0.068]';
x = 1:4;
xl = ["Beam Sagittal" "Beam Frontal" "Ground Sagittal" "Ground Frontal"];
figure();
bar(x,rmses)
legend('\alpha = 10^{-4}', '\alpha = 10^6')
xlabel('Balance Condition')
ylabel('Best-Fit RMSE')
set(gca,'XTickLabel',xl)
improvePlot;
%% T-test
setpath;
old = load('bestparams_duarte_old_finer_subs_20240212.mat');
young = load('bestparams_duarte_young_finer_subs_20240212.mat');
Y = [old.BestBetas];
figure();
violin(Y)
%% Plot Best-Fit Parameter Relation to Other Information
% Best-fit parameters:
old = load('bestparams_duarte_old_finer_subs_20240212.mat');
young = load('bestparams_duarte_young_finer_subs_20240212.mat');
old_betas = old.BestBetas;
old_sigmas = old.BestSigmas;
young_betas = young.BestBetas;
young_sigmas = young.BestSigmas;

% Other demographic information:
old = load('duarte_old.mat');
young = load('duarte_young.mat');
old_demo = old.DataInfo.FES_T;
young_demo = young.DataInfo.FES_T;

% Plot
figure();
sz = 50;
subplot(1, 2, 1);
hold on
scatter(old_demo,old_betas,sz,'filled')
scatter(young_demo,young_betas,sz,'filled')
ylabel('\beta')
%%----- Change! ------%%
xlabel('FES')

subplot(1, 2, 2);
hold on
scatter(old_demo,oldsigmas,sz,'filled')
scatter(young_demo,young_sigmas,sz,'filled')
ylabel('\sigma_r')
%%----- Change! ------%%
xlabel('FES')
legend('Older Subjects','Younger Subjects')
improvePlot;