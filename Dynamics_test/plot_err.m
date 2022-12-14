close all;
clear all;
clc;
%% RMSE
load error
rmse_avb = root_mean_nonpar7sub_lqr_spatial2; %each 2D array is alpha v beta, for given sigma
rmse_avs = permute(rmse_avb,[1 3 2]); %each 2D array is alpha v sigma, for given beta
rmse_bvs = permute(rmse_avb,[2 3 1]); %each 2D array is beta v sigma, for given alpha
%% Find Minimum RMSE
[v_min,loc_min] = min(rmse_avb(:));
[a_min,b_min,s_min] = ind2sub(size(rmse_avb),loc_min);
[v_max,loc_max] = max(rmse_avb(:));
[i_max,j_max,k_max] = ind2sub(size(rmse_avb),loc_max);
%% Parameters
% alpha = linspace(-4,6,11)'; %values of alpha
% beta = [0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2]'; %values of beta
% sigma = [0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,3,5,10]'; %values of sigma
alpha = 1e-4;
beta = [5, 10];
sigma = [0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,3,5,10,15,20];

best_alpha = alpha(a_min),best_beta = beta(b_min),best_sigma = sigma(s_min) %minimum error parameters
alpha(i_max),beta(j_max),sigma(k_max) %maximum error parameters
%% Heatmap RMSE
% % crange = [v_min v_max];
% crange = [0.025 1.35];
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
%% RMSE Sensitivity
root_mean_nonpar = root_mean_nonpar7sub_lqr_spatial;
[~,loc_min] = min(root_mean(:));
[a_min_par,b_min_par,s_min_par] = ind2sub(size(root_mean),loc_min);
[~,loc_min] = min(root_mean_nonpar(:));
[a_min_nonpar,b_min_nonpar,s_min_nonpar] = ind2sub(size(root_mean_nonpar),loc_min);
[~,loc_min] = min(root_mean_young(:));
[a_min_young,b_min_young,s_min_young] = ind2sub(size(root_mean_young),loc_min);
[~,loc_min] = min(root_mean_old(:));
[a_min_old,b_min_old,s_min_old] = ind2sub(size(root_mean_old),loc_min);

% For best beta and sigma:
figure();
plot(alpha,root_mean(:,b_min_par,s_min_par),'--.', 'MarkerSize',20)
hold on
plot(alpha,root_mean_nonpar(:,b_min_nonpar,s_min_nonpar),'--.', 'MarkerSize',20)
plot(alpha,root_mean_young(:,b_min_young,s_min_young),'--.', 'MarkerSize',20)
plot(alpha,root_mean_old(:,b_min_old,s_min_old),'--.', 'MarkerSize',20)
xlabel('\alpha power');
ylabel('RMSE Error')
legend('Paretic','Non-Paretic','Young','Old')
% For best alpha and sigma:
figure();
plot(beta,root_mean(a_min_par,:,s_min_par),'--.', 'MarkerSize',20)
hold on
plot(beta,root_mean_nonpar(a_min_nonpar,:,s_min_nonpar),'--.', 'MarkerSize',20)
plot(beta,root_mean_young(a_min_young,:,s_min_young),'--.', 'MarkerSize',20)
plot(beta,root_mean_old(a_min_old,:,s_min_old),'--.', 'MarkerSize',20)
xlabel('\beta');
ylabel('RMSE Error')
legend('Paretic','Non-Paretic','Young','Old')
% For best alpha and beta:
figure();
error = root_mean(a_min_par,b_min_par,:);
plot(sigma,error(:),'--.', 'MarkerSize',20)
hold on
error = root_mean_nonpar(a_min_nonpar,b_min_nonpar,:);
plot(sigma,error(:),'--.', 'MarkerSize',20)
error = root_mean_young(a_min_young,b_min_young,:);
plot(sigma,error(:),'--.', 'MarkerSize',20)
error = root_mean_old(a_min_old,b_min_old,:);
plot(sigma,error(:),'--.', 'MarkerSize',20)
xlabel('\sigma_r');
ylabel('RMSE Error')
legend('Paretic','Non-Paretic','Young','Old')
%% 3D Scatter
% close all
[A,B,S] = meshgrid(beta,alpha,sigma(1:23));
% rmse_bva = permute(root_mean2,[2 1 3]); %each 2D array is beta v alpha, for given sigma
A = A(:);
B = B(:);
S = S(:);
vals = rmse_avb(:,:,1:23);
vals = vals(:);
% vals = rmse_fake(:);
figure();
% scatter3(A,B,S,[],vals,'filled'), grid on
sample = 1; % downsample 
scatter3(downsample(A,sample),downsample(B,sample),downsample(S,sample),[],downsample(vals,sample),'filled'), grid on
colormap parula
caxis([0.025 1.35])
colorbar
xlabel('Beta')
ylabel('Alpha')
zlabel('Sigma')