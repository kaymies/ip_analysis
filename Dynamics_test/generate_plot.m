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
xlim([0.5 5.3])
ylim([0 2.5])
end
%% Plot with standard deviation
old = load('duarte_old.mat');
young = load('duarte_young.mat');
mean = [old.IPDataAverage; young.IPDataAverage];
std = [old.StandardDeviation; young.StandardDeviation];
high = zeros(size(mean));
low = zeros(size(mean));
for i = 1:height(mean)
    high(i,:) = [mean(i,:) + std(i,:)];
    low(i,:) = [mean(i,:) - std(i,:)];
end
freq = old.Frequency;
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
legend('','Old','','Young')
xlabel('Frequency (Hz)')
ylabel('IP (Fraction of CoM)')

%% Plot each human
old = load('duarte_old.mat');
young = load('duarte_young.mat');
for i = 1:3:old.NumSubjects*3
    ave = mean(old.IP(i:i+2,:),'omitnan');
    plot(old.Frequency,ave)
    hold on
    xlabel('Frequency (Hz)')
    ylabel('IP (Fraction of CoM)')
end
figure();
for i = 1:3:young.NumSubjects*3
    ave = mean(young.IP(i:i+2,:),'omitnan');
    plot(young.Frequency,ave)
    hold on
    xlabel('Frequency (Hz)')
    ylabel('IP (Fraction of CoM)')
end