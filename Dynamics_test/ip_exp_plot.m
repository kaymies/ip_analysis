close all;
clear all;
clc;
setpath;
%%
d1 = load(sprintf('%s.mat','paretic_7sub'));
d2 = load(sprintf('%s.mat','nonpar_7sub'));
d3 = load(sprintf('%s.mat','gruben2019'));
% d4 = load(sprintf('%s.mat','duarte_young_fitted'));
%% Plot
figure();
% plot(d1.Frequency,d1.IPDataAverage)
% hold on
% plot(d2.Frequency,d2.IPDataAverage)
% plot(d3.Frequency(2:24),d3.IPDataAverage(2:24))
% plot(d3.Frequency,d3.IPDataAverage)
% plot(d4.Frequency,d4.IPDataAverage)
% legend('Duarte Old','Duarte Young','Duarte Old (partial)')%,'Duarte Old (Fitted)','Duarte Young (Fitted)')

std_error1 = d1.StandardDeviation;
std_error2 = d2.StandardDeviation;
% std_error3 = std(mean(d1.IPAveSubject(ans,:)))/sqrt(38);
mean1 = d1.IPDataAverage;
mean2 = d2.IPDataAverage;
% mean3 = mean(d1.IPAveSubject(ans,:));
freq = d1.Frequency;
up1 = mean1 + std_error1;
low1 = mean1 - std_error1;
up2 = mean2 + std_error2;
low2 = mean2 - std_error2;
% up3 = mean3 + std_error3;
% low3 = mean3 - std_error3;
x2 = [freq, fliplr(freq)];
inBetween = [up1, fliplr(low1)];
h = fill(x2, inBetween,[0 255/255 0/255],'LineStyle','none');
set(h,'facealpha',.5)
hold on;
plot(freq, mean1,'LineWidth', 2);
x_par = [freq, fliplr(freq)];
inBetween2 = [up2, fliplr(low2)];
h2 = fill(x_par, inBetween2, [0 255/255 0/255], 'LineStyle','none');
set(h2,'facealpha',.5)
hold on;
plot(freq, mean2, 'LineWidth', 2);
plot(d3.Frequency(2:24),d3.IPDataAverage(2:24),'LineWidth',2)
% x_3 = [freq, fliplr(freq)];
% inBetween3 = [up3, fliplr(low3)];
% h3 = fill(x_3, inBetween3, [0 255/255 255/255], 'LineStyle','none');
% set(h3,'facealpha',.5)
% hold on;
% plot(freq, mean3, 'color', [0 255/255 255/255], 'LineWidth', 2);
% plot(freq,young.IPDataAverage,'color', [0 0 255/255], 'LineWidth', 2);
legend('','Paretic Limb','','Non-Paretic Limb','Younger Subjects','','','')
xlabel('Frequency (Hz)')
ylabel('IP (Normalized by CoM)')