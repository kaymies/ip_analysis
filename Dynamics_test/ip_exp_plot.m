close all;
clear all;
clc;
setpath;
%%
d1 = load(sprintf('%s.mat','duarte_old'));
d2 = load(sprintf('%s.mat','duarte_young'));
d3 = load(sprintf('%s.mat','gruben2019'));
d4 = load(sprintf('%s.mat','young_ave'));
%% Plot
figure();
plot(d1.Frequency,d1.IPDataAverage)
hold on
plot(d2.Frequency,d2.IPDataAverage)
plot(d3.Frequency,d3.IPDataAverage)
plot(d4.Frequency,d4.IPDataAverage)
legend('Duarte Old','Duarte Young','Gruben2019','Previous Young')

std_nonpar = nonpar.StandardDeviation;
std_paretic = paretic.StandardDeviation;
up_nonpar = mean_nonpar + std_nonpar;
low_nonpar = mean_nonpar - std_nonpar;
up_paretic = mean_paretic + std_paretic;
low_paretic = mean_paretic - std_paretic;
x2 = [freq, fliplr(freq)];
inBetween = [up_nonpar, fliplr(low_nonpar)];
h = fill(x2, inBetween, [255/255 0/255 0/255], 'LineStyle','none');
set(h,'facealpha',.5)
hold on;
plot(freq, mean_nonpar, 'color', [255/255 0/255 0/255], 'LineWidth', 2);
x_par = [freq, fliplr(freq)];
inBetween2 = [up_paretic, fliplr(low_paretic)];
h2 = fill(x_par, inBetween2, [0 255/255 0/255], 'LineStyle','none');
set(h2,'facealpha',.5)
hold on;
plot(freq, mean_paretic, 'color', [0 255/255 0/255], 'LineWidth', 2);
plot(freq,young.IPDataAverage,'color', [0 0 255/255], 'LineWidth', 2);
legend('','Non-Paretic','','Paretic','Young')
xlabel('Frequency (Hz)')
ylabel('IP (Fraction of CoM)')