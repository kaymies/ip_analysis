close all;
clear all;
clc;
setpath;
%% Raw Data
d1 = load(sprintf('%s.mat','duarte_old'));
d2 = load(sprintf('%s.mat','duarte_young'));
% By Frequency
figure;
N_bins = round(sqrt(height(d1.IPAveSubject)+height(d2.IPAveSubject)));
for i = 1:length(d1.Frequency)
    subplot(5,8,i)
    ip = rmoutliers(real(log([d1.IPAveSubject(:,i);d2.IPAveSubject(:,i)])));
    plot_histogram([],N_bins,ip)
end
%%
% By Subject Group
figure;
N_bins1 = round(sqrt(height(d1.IPAveSubject)+length(d1.Frequency)));
N_bins2 = round(sqrt(height(d2.IPAveSubject)+length(d2.Frequency)));
subplot(2,1,1)
ip = rmoutliers(real(log(d1.IPAveSubject(:))));
plot_histogram([],N_bins1,ip)
xlim([-3 1])
subtitle('Older Subjects')
subplot(2,1,2)
ip = rmoutliers(log(d2.IPAveSubject(:)));
plot_histogram([],N_bins2,ip)
xlim([-3 1])
subtitle('Younger Subjects')
%% Run ANOVA
old = d1;
young = d2;
% zIP = zeros(length(old.Frequency)*(height(old.IPAveSubject)+height(young.IPAveSubject)),1);
zIP = [];
frequency = [];
subs = [];
% frequency = repmat(old.Frequency,1,height(old.IPAveSubject)+height(young.IPAveSubject))';
% subs = vertcat(repmat('Older  ',height(old.IPAveSubject)*length(old.Frequency),1),repmat('Younger',height(young.IPAveSubject)*length(young.Frequency),1));
for f = 1:length(old.Frequency)
    [zIP_new,TFrm] = rmoutliers(real(log(old.IPAveSubject(:,f))));
    frequency_new = repmat(old.Frequency(f),sum(~TFrm),1);
    subs_new = repmat('Older  ',sum(~TFrm),1);
    zIP = [zIP;zIP_new];
    frequency = [frequency;frequency_new];
    subs = [subs;subs_new];
%     i = i + length(old.Frequency);
end
%%
for f = 1:length(old.Frequency)
    [zIP_new,TFrm] = rmoutliers(real(log(young.IPAveSubject(:,f))));
    frequency_new = repmat(young.Frequency(f),sum(~TFrm),1);
    subs_new = repmat('Younger',sum(~TFrm),1);
    zIP = [zIP;zIP_new];
    frequency = [frequency;frequency_new];
    subs = [subs;subs_new];
%     i = i + length(old.Frequency);
end
%%
[p,tbl,stats,terms] = anovan(zIP,{subs,frequency},'model',2,'varnames',{'Age','Frequency'})
% On only a subset of the data f = 0.5 Hz, 1.1 Hz, 2.1 Hz, 4.1 Hz
zIP = [];
frequency = [];
subs = [];
f_idx = [1,4,9,19];
for f = 1:length(f_idx)
    [zIP_new,TFrm] = rmoutliers(real(log(old.IPAveSubject(:,f_idx(f)))));
    frequency_new = repmat(old.Frequency(f_idx(f)),sum(~TFrm),1);
    subs_new = repmat('Older  ',sum(~TFrm),1);
    zIP = [zIP;zIP_new];
    frequency = [frequency;frequency_new];
    subs = [subs;subs_new];
%     i = i + length(old.Frequency);
end

for f = 1:length(f_idx)
    [zIP_new,TFrm] = rmoutliers(real(log(young.IPAveSubject(:,f_idx(f)))));
    frequency_new = repmat(young.Frequency(f_idx(f)),sum(~TFrm),1);
    subs_new = repmat('Younger',sum(~TFrm),1);
    zIP = [zIP;zIP_new];
    frequency = [frequency;frequency_new];
    subs = [subs;subs_new];
%     i = i + length(old.Frequency);
end
[p,tbl,stats,terms] = anovan(zIP,{subs,frequency},'model',2,'varnames',{'Age','Frequency'})

% T test
figure;
[c,~,~,gnames] = multcompare(stats,'Dimension',[1 2],'CriticalValueType','bonferroni');
alpha_corrected = 0.05/4;
for f = 1:length(f_idx)
    ttest_plot(rmoutliers(real(log(old.IPAveSubject(:,f_idx(f))))),...
        rmoutliers(real(log(young.IPAveSubject(:,f_idx(f))))),...
        1-alpha_corrected, sprintf('log(zIP) at %0.5g Hz',old.Frequency(f_idx(f))))
end
    % legend('older','younger')
%% Best Fit Parameters By Subject
% Older subjects
figure;
old = load('subs_bestparams_duarte_old_20230705.mat');
old_betas = rmoutliers(old.BestBetas,'percentiles',[5 95]);
old_sigmas = rmoutliers(old.BestSigmas,'percentiles',[5 95]);
subplot(2,2,1)
% plot_histogram([0.1, 0.2, 0.3, 0.4],[],old_betas(old_betas < 1.5))
plot_histogram([],round(sqrt(length(old_betas))),old_betas)
xlim([-0.2 2.9])
title('Older Subjects \beta')
subplot(2,2,2)
% plot_histogram(log(0:0.1:0.9),[],log(old_sigmas(old_sigmas < 1.5)))
plot_histogram([],round(sqrt(length(old_sigmas))),old_sigmas)
xlim([-0.2 2.9])
title('Older Subjects \sigma_r')
% Younger subjects
young = load('subs_bestparams_duarte_young_20230627.mat');
young_betas = rmoutliers(young.BestBetas,'percentiles',[5 95]);
young_sigmas = rmoutliers(young.BestSigmas,'percentiles',[5 95]);
subplot(2,2,3)
plot_histogram([],round(sqrt(length(young_betas))),young_betas)
xlim([-0.2 2.4])
% plot_histogram(round(sqrt(length(young_betas(young_betas < 1.5)))),log(young_betas(young_betas < 1.5)))
title('Younger Subjects by Trial \beta')
subplot(2,2,4)
plot_histogram([],round(sqrt(length(young_sigmas))),young_sigmas)
xlim([-0.2 2.4])
% plot_histogram(round(sqrt(length(young_sigmas(young_sigmas < 1.5)))),log(young_sigmas(young_sigmas < 1.5)))
title('Younger Subjects by Trial \sigma_r')
sgtitle('Best-Fit Histogram by Subject')
%% Boxplot
[old_betas,TFrm1] = rmoutliers(log(old.BestBetas),'percentiles',[5 95]);
[old_sigmas,TFrm2] = rmoutliers(log(old.BestSigmas),'percentiles',[5 95]);
[young_betas,TFrm3] = rmoutliers(log(young.BestBetas),'percentiles',[5 95]);
[young_sigmas,TFrm4] = rmoutliers(log(young.BestSigmas),'percentiles',[5 95]);
figure;
subplot(1,2,1)
label = [repmat({'Older'},length(old_betas),1); repmat({'Younger'},length(young_betas),1)];
boxplot([old_betas; young_betas],label)
subtitle('\beta')
subplot(1,2,2)
label = [repmat({'Older'},length(old_sigmas),1); repmat({'Younger'},length(young_sigmas),1)];
boxplot([old_sigmas; young_sigmas],label)
subtitle('\sigma_r')
%% By Trial
% Older subjects by trial
figure;
old = load('trials_bestparams_duarte_old_20230705.mat.mat');
old_betas = old.BestBetas;
old_sigmas = old.BestSigmas;
subplot(2,2,1)
plot_histogram(round(sqrt(length(old_betas(old_betas < 1.5)))),log(old_betas(old_betas < 1.5)))
title('Older Subjects \beta')
subplot(2,2,2)
plot_histogram(round(sqrt(length(old_sigmas(old_sigmas < 1.5)))),log(old_sigmas(old_sigmas < 1.5)))
title('Older Subjects \sigma_r')
% Younger subjects by trial
young = load('trials_bestparams_duarte_young_20230706.mat.mat');
young_betas = young.BestBetas;
young_sigmas = young.BestSigmas;
subplot(2,2,3)
plot_histogram(round(sqrt(length(young_betas(young_betas < 1.5)))),log(young_betas(young_betas < 1.5)))
title('Younger Subjects by Trial \beta')
subplot(2,2,4)
plot_histogram(round(sqrt(length(young_sigmas(young_sigmas < 1.5)))),log(young_sigmas(young_sigmas < 1.5)))
title('Younger Subjects by Trial \sigma_r')
sgtitle('Best-Fit Histogram by Trial')
%%
old = load('subs_bestparams_duarte_old_20230705.mat');
young = load('subs_bestparams_duarte_young_20230627.mat');
dold = load(sprintf('%s.mat','duarte_old'));
dyoung = load(sprintf('%s.mat','duarte_young'));

group1 = dold.AgebySub < 70;
sigmas = old.BestSigmas;
betas = old.BestBetas;
figure;
hold on
histogram(betas(group1),linspace(0,2.94,7),EdgeColor='blue',FaceColor='blue')
histogram(betas(~group1),linspace(0,2.94,7),EdgeColor='red',FaceColor='red')
legend('Younger than 70','Older than 70')
xlabel('\beta')
improvePlot;
%%
younger = Age < 25;
older = Age >= 25;
figure;
hold on
histogram(young_sigmas(younger),N_bins,EdgeColor='blue',FaceColor='blue')
histogram(young_sigmas(older),N_bins,EdgeColor='red',FaceColor='red')
%% Kreg's Idea
% [betas,TFrm1] = rmoutliers(old.BestBetas,'percentiles',[5 95]);
% [sigmas,TFrm2] = rmoutliers(old.BestSigmas,'percentiles',[5 95]);
% [betas,TFrm3] = rmoutliers(young.BestBetas,'percentiles',[5 95]);
% [sigmas,TFrm4] = rmoutliers(young.BestSigmas,'percentiles',[5 95]);

% TFrm_old = logical(TFrm1+TFrm2);
% TFrm_old = ~TFrm_old;
x = old.BestBetas;%(TFrm_old);
y = old.BestSigmas;%(TFrm_old);

% Count the number of points on each coordinate
%uniqueCoords gives us all of the unique coordinates of x and y
% pointCounts gives us the index of uniqueCoords that each row of [x,y] is
% corresponds to
[uniqueCoords, ~, pointCounts] = unique([y, x], 'rows');
%countPerCoord gives us how many times each unique coordinate shows up in
%[x,y]
countPerCoord = accumarray(pointCounts, 1);

% Create the scatter plot
figure;
scatter(uniqueCoords(:,1),uniqueCoords(:,2),countPerCoord*50,'filled');
xlabel('\sigma_r');
ylabel('\beta');
title('\beta and \sigma_r Distribution');
hold on

% TFrm_young = logical(TFrm3+TFrm4);
% TFrm_young = ~TFrm_young;
x = young.BestBetas;%(TFrm_young);
y = young.BestSigmas;%(TFrm_young);

% Count the number of points on each coordinate
[uniqueCoords, ~, pointCounts] = unique([y, x], 'rows');
countPerCoord = accumarray(pointCounts, 1);

% Create the scatter plot
scatter(uniqueCoords(:,1),uniqueCoords(:,2),countPerCoord*50,'filled');
legend('Older','Younger')
%% T test
old = load('trials_bestparams_duarte_old_20230705.mat');
old_betas = old.BestBetas; old_sigmas = old.BestSigmas;
young = load('trials_bestparams_duarte_young_20230706.mat');
young_betas = young.BestBetas; young_sigmas = young.BestSigmas;
ttest_plot(old_betas(old_betas <= 1.5), young_betas(young_betas <= 1.5), 0.95, '\beta by Trials')
% legend('older','younger')
%%
zIP = [];
for f = 1:length(old.Frequency)
    zIP_freq = rmoutliers(real(old.IPAveSubject(:,f)));
    zIP_new = mean(zIP_freq);
    zIP = [zIP zIP_new];
%     i = i + length(old.Frequency);
end