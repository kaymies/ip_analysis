function [] = generate_plot(frequency,human_mean,human_sd,sim_mean,sim_sd)
%% Load Data
% Human Data
setpath;
%% Plot human vs. simulation
figure();
errorbar(frequency,human_mean,human_sd,'x')
% plot(frequency,human_mean,'x')
hold on
errorbar(frequency,sim_mean,sim_sd,'o')
% plot(sim_struct.Frequency,sim_struct.DataWithoutOutliers)
legend('Human Data', 'Simulated Data')
xlabel('Frequency (Hz)')
ylabel('Intersection Point Height (Normalized by CoM)')
xlim([0 8])
ylim([0 3.5])
improvePlot;
end
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
% 
% figure;
% % Sagittal Side-by-Side
% % human_struct = load('gruben2019.mat');
% % sim_struct = load('E6_0p3_0p9.mat');
% % subplot(2, 2, 1);
% % plot(human_struct.Frequency, human_struct.IPDataAverage,'x');
% % hold on
% % errorbar(sim_struct.Frequency,sim_struct.DataWithoutOutliers,sim_struct.StandardDeviation,'o')
% % xlabel('Frequency [Hz]')
% % ylabel({'Intersection Point Height';'(Normalized by CoM)'})
% % xlim([0 8])
% % ylim([0 2.5])
% % % title('Plot 1');
% 
% % Sagittal Ground
% human_struct = load('marta_sgt_ground.mat');
% sim_struct = load('E6_0p3_1p6.mat');
% sim_struct_smalla = load('En4_1p3_0p6');
% subplot(2, 2, 1);
% % errorbar(sim_struct_smalla.Frequency,sim_struct_smalla.DataWithoutOutliers,sim_struct_smalla.StandardDeviation,'o')
% hold on
% plot(human_struct.Frequency, human_struct.IPDataAverage,'x');
% errorbar(sim_struct.Frequency,sim_struct.DataWithoutOutliers,sim_struct.StandardDeviation,'o')
% xlabel('Frequency [Hz]')
% ylabel({'Intersection Point Height';'(Normalized by CoM)'})
% xlim([0 8])
% ylim([0 2.5])
% % title('Plot 2');
% 
% % Sagittal Beam
% human_struct = load('marta_sgt_beam.mat');
% sim_struct = load('E6_0p2_1p1.mat');
% sim_struct_smalla = load('En4_1p3_0p6');
% subplot(2, 2, 2);
% % errorbar(sim_struct_smalla.Frequency,sim_struct_smalla.DataWithoutOutliers,sim_struct_smalla.StandardDeviation,'o')
% hold on
% plot(human_struct.Frequency, human_struct.IPDataAverage,'x');
% errorbar(sim_struct.Frequency,sim_struct.DataWithoutOutliers,sim_struct.StandardDeviation,'o')
% xlabel('Frequency [Hz]')
% ylabel({'Intersection Point Height';'(Normalized by CoM)'})
% xlim([0 8])
% ylim([0 2.5])
% % title('Plot 3');
% 
% % Frontal Ground
% human_struct = load('marta_frt_ground.mat');
% sim_struct = load('E6_0p2_0p01.mat');
% sim_struct_smalla = load('En4_2p3_0p4');
% subplot(2, 2, 3);
% % errorbar(sim_struct_smalla.Frequency,sim_struct_smalla.DataWithoutOutliers,sim_struct_smalla.StandardDeviation,'o')
% hold on
% plot(human_struct.Frequency, human_struct.IPDataAverage,'x');
% errorbar(sim_struct.Frequency,sim_struct.DataWithoutOutliers,sim_struct.StandardDeviation,'o')
% xlabel('Frequency [Hz]')
% ylabel({'Intersection Point Height';'(Normalized by CoM)'})
% xlim([0 8])
% ylim([0 2.5])
% % title('Plot 4');
% 
% % Frontal Beam
% subplot(2, 2, 4);
% human_struct = load('marta_frt_beam.mat');
% sim_struct = load('E6_1p2_0p5.mat');
% sim_struct_smalla = load('En4_3_0p1');
% % errorbar(sim_struct_smalla.Frequency,sim_struct_smalla.DataWithoutOutliers,sim_struct_smalla.StandardDeviation,'o')
% hold on
% plot(human_struct.Frequency, human_struct.IPDataAverage,'x');
% errorbar(sim_struct.Frequency,sim_struct.DataWithoutOutliers,sim_struct.StandardDeviation,'o')
% legend('Human Data', 'Best Fit Simulation Data')
% xlabel('Frequency [Hz]')
% ylabel({'Intersection Point Height';'(Normalized by CoM)'})
% xlim([0 8])
% ylim([0 2.5])
% 
% improvePlot;