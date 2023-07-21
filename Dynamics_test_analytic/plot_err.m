close all;
clear all;
clc;
setpath;
%% Load Error
filename = 'rmse_duarte_old37_20230720';
error = load(sprintf('%s.mat',filename));
beta = error.Parameters.beta;
sigma = error.Parameters.sigma_r;
% numsubs = 38;
% numsubs = error.NumSubjects*3;
error = error.Error;
% rmse_avb = root_mean_nonpar7sub_lqr_spatial2; %each 2D array is alpha v beta, for given sigma
% rmse_avs = permute(rmse_avb,[1 3 2]); %each 2D array is alpha v sigma, for given beta
% rmse_bvs = permute(rmse_avb,[2 3 1]); %each 2D array is beta v sigma, for given alpha
%% For Outlier Elders:
% outliers = [1 7 8 9 16 20 24 26 28 35 36 37];
% figure;
% for o = 1:length(outliers)
%     filename = sprintf('rmse_duarte_old%d_20230720',outliers(o));
%     error = load(sprintf('%s.mat',filename));
%     beta = error.Parameters.beta;
%     sigma = error.Parameters.sigma_r;
%     error = error.Error;
%     crange = [0.1425 0.7553];
%     subplot(3,4,o)
%     heatmap(round(sigma,2),round(flip(beta),2),abs(flip(error)),'Colormap',parula,'ColorLimits',crange) % first 2 entries is column, row
%     xlabel('\sigma_r');
%     ylabel('\beta')
%     title(sprintf('Old Sub %d',outliers(o)))
% end
%% For Multiple Subjects:
% best_betas = zeros(numsubs,1);
% best_sigmas = zeros(numsubs,1);
% for i = 1:numsubs
%     error_sub = squeeze(error(i,:,:));
%     [v_min,loc_min] = min(abs(error_sub(:)));
%     [b_min,s_min] = ind2sub(size(error_sub),loc_min);
%     best_betas(i) = beta(b_min); best_sigmas(i) = sigma(s_min);
% end
% params.BestBetas = best_betas;
% params.BestSigmas = best_sigmas;
% filename = 'trials_bestparams_duarte_young_20230706.mat';
% folder = fullfile('Data','BestParams'); addpath(folder);
% file = sprintf('%s.mat',filename);
% 
% % Save the struct to a file
% save(fullfile(folder,file), '-struct', 'params');
%% Find Minimum Error
[v_min,loc_min] = min(abs(error(:)))
[b_min,s_min] = ind2sub(size(error),loc_min);
% [v_max,loc_max] = max(error(:));
[v_max,loc_max] = max(error(find(abs(error(:))<5)))
[b_max,s_max] = ind2sub(size(error),loc_max);
% [v_min,loc_min] = min(rmse_avb(:));
% [a_min,b_min,s_min] = ind2sub(size(error),loc_min);
% [v_max,loc_max] = max(rmse_avb(:));
% [i_max,j_max,k_max] = ind2sub(size(error),loc_max);
%% Parameters
best_params = [beta(b_min), sigma(s_min)]%, gamma(g_min),...
%     kappa(k_min), eta(e_min)] %minimum error parameters
worst_params = [beta(b_max), sigma(s_max)]%, gamma(g_max),...
%     kappa(k_max), eta(e_max)] %maximum error parameters
%% Heatmap RMSE
crange = [v_min v_max];
figure();
heatmap(round(sigma,2),round(flip(beta),2),abs(flip(error)),'Colormap',parula,'ColorLimits',crange) % first 2 entries is column, row
xlabel('\sigma_r');
ylabel('\beta')
% %% Heatmap Difference
% crange = [0 max(abs(v_min),abs(v_max))];
% figure();
% heatmap(round(sigma,2),round(beta,2),abs(error),'Colormap',parula,'ColorLimits',crange) % first 2 entries is column, row
% xlabel('\sigma_r');
% ylabel('\beta')
% 
% %%
% % For given alpha:
% for a = 1:length(alpha)
%     figure();
%     heatmap(round(sigma,2),round(beta,2),rmse_bvs(:,:,a),'Colormap',parula,'ColorLimits',crange) % first 2 entries is column, row
%     xlabel('\sigma_r');
%     ylabel('\beta')
% end
% % For given beta:
% for b = 1:length(beta)
%     figure();
%     heatmap(round(sigma,2),alpha,rmse_avs(:,:,b),'Colormap',jet,'ColorLimits',crange) % first 2 entries is column, row
%     xlabel('\sigma_r');
%     ylabel('\alpha')
% end
% % For given sigma_r:
% for s = 1:length(sigma)
%     figure();
%     heatmap(round(beta,2),alpha,rmse_avb(:,:,s),'Colormap',jet,'ColorLimits',crange) % first 2 entries is column, row
%     xlabel('\beta');
%     ylabel('\alpha')
% end
% %% RMSE Sensitivity
% % root_mean_nonpar = root_mean_nonpar7sub_lqr_spatial;
% [~,loc_min] = min(root_mean(:));
% [a_min_par,b_min_par,s_min_par] = ind2sub(size(root_mean),loc_min);
% [~,loc_min] = min(root_mean_nonpar(:));
% [a_min_nonpar,b_min_nonpar,s_min_nonpar] = ind2sub(size(root_mean_nonpar),loc_min);
% [~,loc_min] = min(root_mean_young(:));
% [a_min_young,b_min_young,s_min_young] = ind2sub(size(root_mean_young),loc_min);
% % [~,loc_min] = min(root_mean_old(:));
% % [a_min_old,b_min_old,s_min_old] = ind2sub(size(root_mean_old),loc_min);
% 
% alpha_2 = logspace(-4,6,11)';
% beta_2 = [0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2]';
% sigma_2 = [0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,3,5,10]';
% %%
% figure();
% t = tiledlayout(2,3);
% % Alpha sweep:
% ax1 = nexttile;
% hold on
% plot(alpha_2,root_mean_nonpar(:,b_min_nonpar,s_min_nonpar),'--.', 'MarkerSize',20,'DisplayName','Old Non-Paretic')
% plot(alpha_2,root_mean_young(:,b_min_young,s_min_young),'--.', 'MarkerSize',20,'DisplayName','Young')
% plot(alpha_2,root_mean(:,b_min_par,s_min_par),'--.', 'MarkerSize',20,'DisplayName','Paretic')
% plot(alpha,error(:,b_min,s_min,g_min,k_min,e_min),'--.', 'MarkerSize',20,'DisplayName','New Non-Paretic')
% set(gca,'Xscale','log')
% xlabel('\alpha');
% ylabel('RMSE Error')
% 
% % Beta sweep:
% ax2 = nexttile;hold on
% plot(beta_2,root_mean_nonpar(a_min_nonpar,:,s_min_nonpar),'--.', 'MarkerSize',20,'DisplayName','Old Non-Paretic')
% plot(beta_2,root_mean_young(a_min_young,:,s_min_young),'--.', 'MarkerSize',20,'DisplayName','Young')
% plot(beta_2,root_mean(a_min_par,:,s_min_par),'--.', 'MarkerSize',20,'DisplayName','Paretic')
% plot(beta,error(a_min,:,s_min,g_min,k_min,e_min),'--.', 'MarkerSize',20,'DisplayName','New Non-Paretic')
% set(gca,'Xscale','log')
% xlabel('\beta');
% ylabel('RMSE Error')
% % Sigma sweep:
% ax3 = nexttile;
% hold on
% plot(sigma_2,squeeze(root_mean_nonpar(a_min_nonpar,b_min_nonpar,:)),'--.', 'MarkerSize',20,'DisplayName','Old Non-Paretic')
% plot(sigma_2,squeeze(root_mean_young(a_min_young,b_min_young,:)),'--.', 'MarkerSize',20,'DisplayName','Young')
% plot(sigma_2,squeeze(root_mean(a_min_par,b_min_par,:)),'--.', 'MarkerSize',20,'DisplayName','Paretic')
% plot(sigma_r,squeeze(error(a_min,b_min,:,g_min,k_min,e_min)),'--.', 'MarkerSize',20,'DisplayName','New Non-Paretic')
% set(gca,'Xscale','log')
% xlabel('\sigma_r');
% ylabel('RMSE Error')
% legend('show')
% % Gamma sweep:
% ax4 = nexttile;
% plot(gamma,squeeze(error(a_min,b_min,s_min,:,k_min,e_min)),'--.', 'MarkerSize',20)
% set(gca,'Xscale','log')
% xlabel('\gamma');
% ylabel('RMSE Error')
% % Kappa sweep:
% ax5 = nexttile;
% plot(kappa,squeeze(error(a_min,b_min,s_min,g_min,:,e_min)),'--.', 'MarkerSize',20)
% set(gca,'Xscale','log')
% xlabel('\kappa');
% ylabel('RMSE Error')
% % Eta sweep:
% ax6 = nexttile;
% plot(eta,squeeze(error(a_min,b_min,s_min,g_min,k_min,:)),'--.', 'MarkerSize',20)
% set(gca,'Xscale','log')
% xlabel('\eta');
% ylabel('RMSE Error')
% linkaxes([ax1 ax2 ax3 ax4 ax5 ax6],'y')
% title(t,'Best Fit Params: \alpha = 0.01, \beta = 0.5, \sigma_r = 1, \gamma = 35, \kappa = 0.1, \eta = 100')
% 
% %%
% % For best beta and sigma:
% figure();
% plot(alpha,root_mean(:,b_min_par,s_min_par),'--.', 'MarkerSize',20)
% hold on
% plot(alpha,root_mean_nonpar(:,b_min_nonpar,s_min_nonpar),'--.', 'MarkerSize',20)
% plot(alpha,root_mean_young(:,b_min_young,s_min_young),'--.', 'MarkerSize',20)
% plot(alpha,root_mean_old(:,b_min_old,s_min_old),'--.', 'MarkerSize',20)
% xlabel('\alpha power');
% ylabel('RMSE Error')
% 
% % For best alpha and sigma:
% figure();
% plot(beta,root_mean(a_min_par,:,s_min_par),'--.', 'MarkerSize',20)
% hold on
% plot(beta,root_mean_nonpar(a_min_nonpar,:,s_min_nonpar),'--.', 'MarkerSize',20)
% plot(beta,root_mean_young(a_min_young,:,s_min_young),'--.', 'MarkerSize',20)
% plot(beta,root_mean_old(a_min_old,:,s_min_old),'--.', 'MarkerSize',20)
% xlabel('\beta');
% ylabel('RMSE Error')
% legend('Paretic','Non-Paretic','Young','Old')
% % For best alpha and beta:
% figure();
% error = root_mean(a_min_par,b_min_par,:);
% plot(sigma,error(:),'--.', 'MarkerSize',20)
% hold on
% error = root_mean_nonpar(a_min_nonpar,b_min_nonpar,:);
% plot(sigma,error(:),'--.', 'MarkerSize',20)
% error = root_mean_young(a_min_young,b_min_young,:);
% plot(sigma,error(:),'--.', 'MarkerSize',20)
% error = root_mean_old(a_min_old,b_min_old,:);
% plot(sigma,error(:),'--.', 'MarkerSize',20)
% xlabel('\sigma_r');
% ylabel('RMSE Error')
% legend('Paretic','Non-Paretic','Young','Old')
% %% 3D Scatter
% % close all
% % [A,B,S] = meshgrid(beta,alpha,sigma(1:23));
% [G,K,E] = meshgrid(kappa,gamma,eta);
% % [A,G,E] = meshgrid(gamma,alpha,eta);
% % rmse_bva = permute(root_mean2,[2 1 3]); %each 2D array is beta v alpha, for given sigma
% % A = A(:);
% % B = B(:);
% % S = S(:);
% G = G(:);
% K = K(:);
% E = E(:);
% % vals = error(:,:,1:23);
% % vals = vals(:);
% % vals = rmse_fake(:);
% figure();
% % scatter3(A,B,S,[],vals,'filled'), grid on
% sample = 1; % downsample 
% % scatter3(downsample(A,sample),downsample(B,sample),downsample(S,sample),[],downsample(vals,sample),'filled'), grid on
% B = reshape(error(a_min,b_min,s_min,:,:,:),[length(gamma),length(kappa),length(eta)]);
% % B = reshape(error(:,b_min,s_min,:,k_min,:),[length(alpha),length(gamma),length(eta)]);
% % scatter3(downsample(G,sample),downsample(A,sample),downsample(E,sample),[],downsample(B(:),sample),'filled')
% s = scatter3(downsample(K,sample),downsample(G,sample),downsample(E,sample),[],downsample(B(:),sample),'filled')
% s.SizeData = 100;
% set(gca,'Xscale','log')
% set(gca,'ZScale','log')
% set(gca,'YScale','log')
% grid on
% colormap parula
% caxis([min(B(:)) max(B(:))])
% colorbar
% xlabel('Gamma')
% ylabel('Kappa')
% zlabel('Eta')