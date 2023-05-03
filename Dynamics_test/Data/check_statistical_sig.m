N_bins = 20;
% Older subjects
old = load('duarte_old.mat');
for i = 1:size(old.IP,2)
    plot_histogram(N_bins,IP(:,i))
end
% Younger subjects
young = load('duarte_young.mat');
for i = 1:size(young.IP,2)
    plot_histogram(N_bins,IP(:,i))
end
for i = 1:size(old.IP,2)
    ttest_plot(old.IP(:,i), young.IP(:,i), 0.95, {sprintf('Frequency: %02d',i)})
end