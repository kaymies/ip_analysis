% Example of ANOVA analysis
% ANOVA is short for "ANalysis Of VAriations"
% ANOVA "answers the question":

% ANOVA tests the hypothesis that the samples in the columns are 
% drawn from populations with the same mean against the alternative 
% hypothesis that the population means are not all the same. 
% If you have a critical p-value p_critical = 0.05 (by 2.671 convention)
% When ANOVA returns a p-value > 0.05, you can conclude that
% there is NO statistically significant difference between categories
% When ANOVA returns a p-value < 0.05, you can conclude that
% there is a statistically signfiicant difference between categories
clear all;
setpath;

old = load('duarte_old.mat');
young = load('duarte_young.mat');
zIP = zeros(length(old.Frequency)*(height(old.IPAveSubject)+height(young.IPAveSubject)),1);
frequency = repmat(old.Frequency,1,height(old.IPAveSubject)+height(young.IPAveSubject))';
subs = vertcat(repmat('Older  ',height(old.IPAveSubject)*length(old.Frequency),1),repmat('Younger',height(young.IPAveSubject)*length(young.Frequency),1));
i = 1;
for o = 1:height(old.IPAveSubject)
    zIP(i:i+length(old.Frequency)-1) = old.IPAveSubject(o,:)';
    i = i + length(old.Frequency);
end

for y = 1:height(young.IPAveSubject)
    zIP(i:i+length(young.Frequency)-1) = young.IPAveSubject(y,:)';
    i = i + length(young.Frequency);
end
%%
[p,tbl,stats,terms] = anovan(zIP,{subs,frequency},'model',2,'varnames',{'Age','Frequency'})
figure;
[c,~,~,gnames] = multcompare(stats,'Estimate','column');
%%
% 
% x = 1:6;   % 6 Measurements
% n_cat = 5; % 5 categories
% 
% % Array of category of every datapoint
% cat = repmat([1:n_cat], length(x), 1);
% 
% 
% %% Random data
% rng(216); % Seed random number for repeatable randomness
% y = randn(size(cat));
% [p_value, tbl, stats] = anova1(y, 1:n_cat, 'off');
% 
% % See MATLAB documentation for multcompare for details on how to use and
% % read the multcompare plot.
% figure;
% [c,~,~,gnames] = multcompare(stats);
% 
% figure;
% hold on;
% for idx = 1:n_cat
%     h = scatter(cat(:, idx), y(:, idx), 'SizeData', 100);
%     h.MarkerFaceColor = h.MarkerEdgeColor;
%     h.MarkerFaceAlpha = 0.5;
%     h.MarkerEdgeAlpha = 0.5;
% end
% xlabel('Category');
% ylabel('Y Arb');
% title({'Random data', sprintf('p=%.1e', round(p_value, 2, 'significant'))});
% xlim([0, n_cat + 1]);
% improvePlot;
% 
% %% Add trend to same random data
% % Note that the data looks like it has a linear trend with respect to
% % category number, but category order does not matter to ANOVA
% rng(216); % Seed random number for repeatable randomness
% y = randn(size(cat)) + 0.5*cat;
% [p_value, tbl, stats] = anova1(y, 1:n_cat, 'off');
% 
% % See MATLAB documentation for multcompare for details on how to use and
% % read the multcompare plot.
% figure;
% [c,~,~,gnames] = multcompare(stats);
% 
% figure;
% hold on;
% for idx = 1:n_cat
%     h = scatter(cat(:, idx), y(:, idx), 'SizeData', 100);
%     h.MarkerFaceColor = h.MarkerEdgeColor;
%     h.MarkerFaceAlpha = 0.5;
%     h.MarkerEdgeAlpha = 0.5;
% end
% xlabel('Category');
% ylabel('Y Arb');
% title({'Random data with trend', sprintf('p=%.1e', round(p_value, 2, 'significant'))});
% xlim([0, n_cat + 1]);
% improvePlot;
% 
% 
% %% Same random data as trend example, different category order
% % This demonstrates how category order does not matter to ANOVA
% rng(216); % Seed random number for repeatable randomness
% y = randn(size(cat)) + 0.5*cat;
% shuffle_cat = repmat([3,4,2,5,1], length(x), 1);
% [p_value, tbl, stats] = anova1(y, 1:n_cat, 'off');
% 
% % See MATLAB documentation for multcompare for details on how to use and
% % read the multcompare plot.
% figure;
% [c,~,~,gnames] = multcompare(stats);
% 
% figure;
% hold on;
% for idx = 1:n_cat
%     h = scatter(shuffle_cat(:, idx), y(:, idx), 'SizeData', 100);
%     h.MarkerFaceColor = h.MarkerEdgeColor;
%     h.MarkerFaceAlpha = 0.5;
%     h.MarkerEdgeAlpha = 0.5;
% end
% xlabel('Category');
% ylabel('Y Arb');
% title({'Random data with trend', sprintf('p=%.1e', round(p_value, 2, 'significant'))});
% xlim([0, n_cat + 1]);
% improvePlot;
% 
% %% Add random data, one column is obviously different
% rng(216); % Seed random number for repeatable randomness
% y = randn(size(cat));
% y(:, 4) = y(:, 4) + 2;
% [p_value, tbl, stats] = anova1(y, 1:n_cat, 'off');
% 
% % See MATLAB documentation for multcompare for details on how to use and
% % read the multcompare plot.
% figure;
% [c,~,~,gnames] = multcompare(stats);
% 
% figure;
% hold on;
% for idx = 1:n_cat
%     h = scatter(cat(:, idx), y(:, idx), 'SizeData', 100);
%     h.MarkerFaceColor = h.MarkerEdgeColor;
%     h.MarkerFaceAlpha = 0.5;
%     h.MarkerEdgeAlpha = 0.5;
% end
% xlabel('Category');
% ylabel('Y Arb');
% title({'Random data with a distinct category', sprintf('p=%.1e', round(p_value, 2, 'significant'))});
% xlim([0, n_cat + 1]);
% improvePlot;
