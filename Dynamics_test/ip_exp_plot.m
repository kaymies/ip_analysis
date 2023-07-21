close all;
clear all;
clc;
setpath;
%%
d1 = load(sprintf('%s.mat','duarte_old'));
d2 = load(sprintf('%s.mat','duarte_young'));
% d3 = load(sprintf('%s.mat','duarte_old_fitted'));
% d4 = load(sprintf('%s.mat','duarte_young_fitted'));
%% Plot
figure();
plot(d1.Frequency,d1.IPDataAverage)
hold on
plot(d2.Frequency,d2.IPDataAverage)
plot(d1.Frequency,mean(d1.IPAveSubject(ans,:)))
% plot(d3.Frequency,d3.IPDataAverage)
% plot(d4.Frequency,d4.IPDataAverage)
legend('Duarte Old','Duarte Young','Duarte Old (partial)')%,'Duarte Old (Fitted)','Duarte Young (Fitted)')

std_error1 = std(d1.IPAveSubject)/sqrt(d1.NumSubjects);
std_error2 = std(d2.IPAveSubject)/sqrt(d2.NumSubjects);
std_error3 = std(mean(d1.IPAveSubject(ans,:)))/sqrt(38);
mean1 = d1.IPDataAverage;
mean2 = d2.IPDataAverage;
mean3 = mean(d1.IPAveSubject(ans,:));
freq = d1.Frequency;
up1 = mean1 + std_error1;
low1 = mean1 - std_error1;
up2 = mean2 + std_error2;
low2 = mean2 - std_error2;
up3 = mean3 + std_error3;
low3 = mean3 - std_error3;
x2 = [freq, fliplr(freq)];
inBetween = [up1, fliplr(low1)];
h = fill(x2, inBetween, [255/255 0/255 0/255], 'LineStyle','none');
set(h,'facealpha',.5)
hold on;
plot(freq, mean1, 'color', [255/255 0/255 0/255], 'LineWidth', 2);
x_par = [freq, fliplr(freq)];
inBetween2 = [up2, fliplr(low2)];
h2 = fill(x_par, inBetween2, [0 255/255 0/255], 'LineStyle','none');
set(h2,'facealpha',.5)
hold on;
plot(freq, mean2, 'color', [0 255/255 0/255], 'LineWidth', 2);
x_3 = [freq, fliplr(freq)];
inBetween3 = [up3, fliplr(low3)];
h3 = fill(x_3, inBetween3, [0 255/255 255/255], 'LineStyle','none');
set(h3,'facealpha',.5)
hold on;
plot(freq, mean3, 'color', [0 255/255 255/255], 'LineWidth', 2);
% plot(freq,young.IPDataAverage,'color', [0 0 255/255], 'LineWidth', 2);
legend('','','','Duarte Older','','Duarte Younger','','Duarte Older (partial)')
xlabel('Frequency (Hz)')
ylabel('IP (Fraction of CoM)')