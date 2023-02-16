clear all;
clc;
close all;
%% Load Human Data
subject_type = 'paretic';
paretic_struct = load(sprintf('%s.mat',subject_type));
subject_type = 'nonpar';
nonpar_struct = load(sprintf('%s.mat',subject_type));
%% Plot Force AP
figure(1);
plot(paretic_struct.ForceAP(:,:))
hold on
plot(nonpar_struct.ForceAP(:,:),'--')
xlabel('Time (\times 10^{-2} s)')
ylabel('AP Force (N)')
legend('Subject 1','Subject 2','Subject 3','Subject 4','Subject 5',...
    'Subject 7','Subject 8','Subject 9','Subject 10','Subject 1','Subject 2',...
    'Subject 3','Subject 4','Subject 5','Subject 7','Subject 8','Subject 9',...
    'Subject 10')
mean_forceAP = zeros(9,2);
sd_forceAP = zeros(9,2);
for i = 1:9
    mean_forceAP(i,1) = mean(paretic_struct.ForceAP(:,i));
    mean_forceAP(i,2) = mean(nonpar_struct.ForceAP(:,i));
    sd_forceAP(i,1) = std(paretic_struct.ForceAP(:,i));
    sd_forceAP(i,2) = std(nonpar_struct.ForceAP(:,i));
end
figure(2);
b = bar([1,2,3,4,5,7,8,9,10],mean_forceAP, 'grouped');
hold on
[ngroups, nbars] = size(mean_forceAP);
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    % Calculate center of each bar
    x(1:5) = (1:5) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    x(6:9) = 1 + (6:9) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, mean_forceAP(:,i), sd_forceAP(:,i), 'k', 'linestyle', 'none');
end
hold off
xlabel('Subject Number')
ylabel('AP Force (N)')
legend('Paretic','Non-Paretic')
%% Plot Force Vertical
figure(3);
plot(paretic_struct.ForceVert(:,:),'LineWidth',2)
hold on
plot(nonpar_struct.ForceVert(:,:),'--','LineWidth',1)
xlabel('Time (\times 10^{-2} s)')
ylabel('Vertical Force (N)')
legend('Subject 1','Subject 2','Subject 3','Subject 4','Subject 5',...
    'Subject 7','Subject 8','Subject 9','Subject 10','Subject 1','Subject 2',...
    'Subject 3','Subject 4','Subject 5','Subject 7','Subject 8','Subject 9',...
    'Subject 10')
mean_forceVert = zeros(9,2);
sd_forceVert = zeros(9,2);
for i = 1:9
    mean_forceVert(i,1) = mean(paretic_struct.ForceVert(:,i));
    mean_forceVert(i,2) = mean(nonpar_struct.ForceVert(:,i));
    sd_forceVert(i,1) = std(paretic_struct.ForceVert(:,i));
    sd_forceVert(i,2) = std(nonpar_struct.ForceVert(:,i));
end
figure(4);
b = bar([1,2,3,4,5,7,8,9,10],mean_forceVert, 'grouped');
hold on
[ngroups, nbars] = size(mean_forceVert);
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    % Calculate center of each bar
    x(1:5) = (1:5) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    x(6:9) = 1 + (6:9) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, mean_forceVert(:,i), sd_forceVert(:,i), 'k', 'linestyle', 'none');
end
hold off
xlabel('Subject Number')
ylabel('AP Force (N)')
legend('Paretic','Non-Paretic')
%% Plot CoP AP
figure(5);
plot(paretic_struct.CoPAP(:,:),'LineWidth',2)
hold on
plot(nonpar_struct.CoPAP(:,:),'--','LineWidth',1)
xlabel('Time (\times 10^{-2} s)')
ylabel('AP CoP (m)')
legend('Subject 1','Subject 2','Subject 3','Subject 4','Subject 5',...
    'Subject 7','Subject 8','Subject 9','Subject 10','Subject 1','Subject 2',...
    'Subject 3','Subject 4','Subject 5','Subject 7','Subject 8','Subject 9',...
    'Subject 10')
mean_CoPAP = zeros(9,2);
sd_CoPAP = zeros(9,2);
for i = 1:9
    mean_CoPAP(i,1) = mean(paretic_struct.CoPAP(:,i));
    mean_CoPAP(i,2) = mean(nonpar_struct.CoPAP(:,i));
    sd_CoPAP(i,1) = std(paretic_struct.CoPAP(:,i));
    sd_CoPAP(i,2) = std(nonpar_struct.CoPAP(:,i));
end
figure(6);
b = bar([1,2,3,4,5,7,8,9,10],mean_CoPAP, 'grouped');
hold on
[ngroups, nbars] = size(mean_CoPAP);
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    % Calculate center of each bar
    x(1:5) = (1:5) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    x(6:9) = 1 + (6:9) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, mean_CoPAP(:,i), sd_CoPAP(:,i), 'k', 'linestyle', 'none');
end
hold off
xlabel('Subject Number')
ylabel('AP CoP (m)')
legend('Paretic','Non-Paretic')
% Compare Mean Paretic vs. Mean Non-Paretic
all_mean_ForceAP = mean(paretic_struct.ForceAP(:))