close all;
clear all;
clc;
setpath
%% Load Human Data
subject_type = 'duarte_old';
human_struct = load(sprintf('%s.mat',subject_type));
input.method = 'cpsd';
input.Frequency = eval(sprintf('human_struct.Frequency_%s',input.method));
input.FrequencyWindow = 0.2;
human_mean = eval(sprintf('human_struct.IPDataAverage_%s',input.method));
human_data_sub = eval(sprintf('human_struct.IPAveSubject_%s',input.method));
human_sd = std(human_data_sub);
%% Set Model Parameters
input.TotalMass = human_struct.MeanMass_kg;
input.TotalHeight = human_struct.MeanHeight_m;  
input.gender = 'M';
input.plane = human_struct.Plane;
input.model = 'DIP';
input.pose = human_struct.Pose;
input.FreqSampKin = 1000;
input.trialDuration = 60;
input.CoordinateFrame = 'relative';
%% Set Controller Parameters
%-----% MAKE SURE TO CHANGE SIM FILENAME!!! %-----%
filename = 'test'; %alpha_beta_sigmar
folder_name = 'LQR_relative';
%-------------------------------------------------%
% input.Controller.alpha = 1e6;
% input.Controller.beta = 0.19;
input.NoiseRatio = 1;
input.Controller.gamma = 1;
input.Controller.kappa = 1;
input.Controller.eta = 1;
input.Controller.type = 'LQR';
%% Visualization
input.PostProc.AnimOn = 0;
input.PostProc.PlotOn = 0;
%% Average Gain
%%%%%% Make sure to change population %%%%%%%
input.Controller.alpha = 10^6;
input.Controller.beta = 0.25;
[~, ~, model] = test_model(input);
% Gain
Ks = model.K;
% Eigenstructure
[stiff_V,stiff_eig] = eig(model.K(:,1:2));
[damp_V,damp_eig] = eig(model.K(:,3:4));
% Symmetric + Antisymmetric
stiffness_sym = (Ks(:,1:2)+Ks(:,1:2)')/2;
stiffness_asym = (Ks(:,1:2)-Ks(:,1:2)')/2;
damping_sym = (Ks(:,3:4)+Ks(:,3:4)')/2;
damping_asym = (Ks(:,3:4)-Ks(:,3:4)')/2;
% Ratio
K_ratio_sqrt = sqrt(abs(det(stiffness_asym)))/sqrt(abs(det(stiffness_sym)));
K_ratio = abs(det(stiffness_asym))/abs(det(stiffness_sym));
B_ratio_sqrt = sqrt(abs(det(damping_asym)))/sqrt(abs(det(damping_sym)));
B_ratio = abs(det(damping_asym))/abs(det(damping_sym));
%% Save File
filename = 'eigen_bestparams_avg';
eigen = load(sprintf('%s.mat',filename));
folder = fullfile('Data','Eigenstructure');
% addpath(folder);
file = sprintf('%s.mat',filename);
folder = fullfile(folder,file);
%Save data as struct
% eigen.D_ratio = D_ratio;
% eigen.V_angle = V_angle;
% eigen.Values = Ds;
% eigen.Vectors = Vs;
eigen.Young.Gain = Ks;
eigen.Young.StiffnessEig = stiff_eig;
eigen.Young.DampingEig = damp_eig;
eigen.Young.StiffnessV = stiff_V;
eigen.Young.DampingV = damp_V;
eigen.Young.StiffnessSym = stiffness_sym;
eigen.Young.StiffnessAsym = stiffness_asym;
eigen.Young.DampingSym = damping_sym;
eigen.Young.DampingAsym = damping_asym;
eigen.Young.KRatioSqrt = K_ratio_sqrt;
eigen.Young.KRatio = K_ratio;
eigen.Young.BRatioSqrt = B_ratio_sqrt;
eigen.Young.BRatio = B_ratio;
save(folder, '-struct', 'eigen');
%% Plot
figure();
hold on;
drawStiffnessEllipseSym(Young.Gain,1,[],[],[],[],[]);
drawStiffnessEllipseAsym(Young.Gain,1,[],[],[],[],[]);
drawStiffnessEllipseSym(Old.Gain,1,[],[],[],[],[]);
drawStiffnessEllipseAsym(Old.Gain,1,[],[],[],[],[]);
improvePlot;
legend('Younger','','Older','')
%% Per-Subject Gain
%%%%%% Make sure to change population %%%%%%%
bestparams = load('bestparams_duarte_old_finer_subs_20240212.mat');
input.Controller.alpha = 10^6;
beta = bestparams.BestBetas;
numsubs = length(beta);
stiffness_sym = zeros(2,2,numsubs);
stiffness_asym = zeros(2,2,numsubs);
damping_sym = zeros(2,2,numsubs);
damping_asym = zeros(2,2,numsubs);
stiff_V = zeros(2,2,numsubs);
stiff_eig = zeros(2,2,numsubs);
damp_V = zeros(2,2,numsubs);
damp_eig = zeros(2,2,numsubs);
K_ratio_sqrt = zeros(numsubs,1);
K_ratio = zeros(numsubs,1);
B_ratio_sqrt = zeros(numsubs,1);
B_ratio = zeros(numsubs,1);
Ks = zeros(2,4,length(beta));
for b = 1:numsubs
    input.Controller.beta = beta(b);
    clear test_model
    [~, ~, model] = test_model(input);
    % Gain
    Ks(:,:,b) = model.K;
    % Eigenstructure
    [stiff_V(:,:,b),stiff_eig(:,:,b)] = eig(model.K(:,1:2));
    [damp_V(:,:,b),damp_eig(:,:,b)] = eig(model.K(:,3:4));
    % Symmetric + Antisymmetric
    stiffness_sym(:,:,b) = (model.K(:,1:2)+model.K(:,1:2)')/2;
    stiffness_asym(:,:,b) = (model.K(:,1:2)-model.K(:,1:2)')/2;
    damping_sym(:,:,b) = (model.K(:,3:4)+model.K(:,3:4)')/2;
    damping_asym(:,:,b) = (model.K(:,3:4)-model.K(:,3:4)')/2;
    % Ratio
    K_ratio_sqrt(b,:) = sqrt(abs(det(stiffness_asym(:,:,b))))/sqrt(abs(det(stiffness_sym(:,:,b))));
    K_ratio(b,:) = abs(det(stiffness_asym(:,:,b)))/abs(det(stiffness_sym(:,:,b)));
    B_ratio_sqrt(b,:) = sqrt(abs(det(damping_asym(:,:,b))))/sqrt(abs(det(damping_sym(:,:,b))));
    B_ratio(b,:) = abs(det(damping_asym(:,:,b)))/abs(det(damping_sym(:,:,b)));
end
%% Save File
filename = 'eigen_bestparams';
eigen = load(sprintf('%s.mat',filename));
folder = fullfile('Data','Eigenstructure');
% addpath(folder);
file = sprintf('%s.mat',filename);
folder = fullfile(folder,file);
%Save data as struct
% eigen.D_ratio = D_ratio;
% eigen.V_angle = V_angle;
% eigen.Values = Ds;
% eigen.Vectors = Vs;
eigen.Young.Gain = Ks;
eigen.Young.StiffnessEig = stiff_eig;
eigen.Young.DampingEig = damp_eig;
eigen.Young.StiffnessV = stiff_V;
eigen.Young.DampingV = damp_V;
eigen.Young.StiffnessSym = stiffness_sym;
eigen.Young.StiffnessAsym = stiffness_asym;
eigen.Young.DampingSym = damping_sym;
eigen.Young.DampingAsym = damping_asym;
eigen.Young.KRatioSqrt = K_ratio_sqrt;
eigen.Young.KRatio = K_ratio;
eigen.Young.BRatioSqrt = B_ratio_sqrt;
eigen.Young.BRatio = B_ratio;
save(folder, '-struct', 'eigen');
%% Plot
eigen = load('eigen_bestparams.mat');
old_Kratio = eigen.Old.KRatioSqrt;
old_Bratio = eigen.Old.BRatioSqrt;
young_Kratio = eigen.Young.KRatioSqrt;
young_Bratio = eigen.Young.BRatioSqrt;
for i = 1:length(old_Kratio)
    if old_Kratio(i) > 3
        old_Kratio(i) = NaN;
    end
end
old_Bratio = neglect_outliers(old_Bratio);
young_Kratio = neglect_outliers(young_Kratio);
young_Bratio = neglect_outliers(young_Bratio);
[mean_old_kratio,unc_old_kratio,std_old_kratio,mean_young_kratio,unc_young_kratio,std_young_kratio] = ttest_plot(old_Kratio,young_Kratio,0.95,'KRatioSqrt')
[mean_old_bratio,unc_old_bratio,std_old_bratio,mean_young_bratio,unc_young_bratio,std_young_bratio] = ttest_plot(old_Bratio,young_Bratio,0.95,'BRatioSqrt')
%%----------------------------------------------------------------------%%
%% Below is draft code
for b = 1:length(beta)
    figure(1);
    hold on;
    drawStiffnessEllipseSym(squeeze(Ks(:,:,b)),1,[],[],[],[],length(beta));
    improvePlot;
%     subplot(1,2,1)
    colorbar
    cmap = colormap(parula(length(beta))) ; %Create Colormap
    set(gca,'CLim',[0.01 3]);
    cbh = colorbar ; %Create Colorbar
    cbh.TickLabels = beta;
    cbh.FontSize = 10;
    cbh.Ticks = linspace(0,3,length(beta));
%     ylim([-2500,2500])
%     xlim([-2500,2500])
%     subplot(1,2,2)
%     ylim([-1000,1000])
%     xlim([-1000,1000])

%     
%     figure(1);
%     drawStiffnessEllipseAsym(squeeze(Ks(:,:,b)),1,[],[],[],[],2);
%     improvePlot;
%     subplot(1,2,1)
% %     ylim([-1500,1500])
% %     xlim([-1500,1500])
%     subplot(1,2,2)
% %     ylim([-300,300])
% %     xlim([-300,300])
end
filename = 'eigen_bestparams_avg';
eigen = load(sprintf('%s.mat',filename));
folder = fullfile('Data','Eigenstructure');
% addpath(folder);
file = sprintf('%s.mat',filename);
folder = fullfile(folder,file);
%Save data as struct
% eigen.D_ratio = D_ratio;
% eigen.V_angle = V_angle;
% eigen.Values = Ds;
% eigen.Vectors = Vs;
eigen.D_stiff = D_stiff;
eigen.D_damp = D_damp;
eigen.V_stiff = V_stiff;
eigen.V_damp = V_damp;
% save(folder, '-struct', 'eigen');
%% Subject Compute Gains of Best Param Beta
%%%%%% Make sure to change population %%%%%%%
alpha = 10^6;
bestparams = load('bestparams_duarte_old_finer_subs_20231018.mat');
beta = bestparams.BestBetas;
Ks = zeros(2,4,length(beta));
for b = 1:length(beta)
    input.Controller.alpha = alpha;
    input.Controller.beta = beta(b);
    clear test_model
    [~, ~, model] = test_model(input);
    Ks(:,:,b) = model.K;
end

% filename = 'eigen_bestparams';
% eigen = load(sprintf('%s.mat',filename));
% folder = fullfile('Data','Eigenstructure');
% % addpath(folder);
% file = sprintf('%s.mat',filename);
% folder = fullfile(folder,file);
% %Save data as struct
% % eigen.D_ratio = D_ratio;
% % eigen.V_angle = V_angle;
% % eigen.Values = Ds;
% % eigen.Vectors = Vs;
% eigen.GainsOld = Ks;
% % eigen.GainsYoung = Ks;
% save(folder, '-struct', 'eigen');
%% Subject Compute Principal Eigenvector/value for Stiffness/Damping
gains = GainsOld;
numsubs = 38;
stiffness_sym = zeros(2,2,numsubs);
stiffness_asym = zeros(2,2,numsubs);
damping_sym = zeros(2,2,numsubs);
damping_asym = zeros(2,2,numsubs);
% V_stiff_asym = zeros(2,2,numsubs);
% D_stiff = zeros(2,2,numsubs);
% D_damp = zeros(2,2,numsubs);
for i = 1:numsubs
%     stiffness_sym(:,:,i) = (gains(:,1:2,i)+gains(:,1:2,i)')/2;
%     stiffness_asym(:,:,i) = (gains(:,1:2,i)-gains(:,1:2,i)')/2;
%     damping_sym(:,:,i) = (gains(:,3:4,i)+gains(:,3:4,i)')/2;
%     damping_asym(:,:,i) = (gains(:,3:4,i)-gains(:,3:4,i)')/2;

    [~,D_stiff(:,:,i)] = eig(gains(:,1:2,i));
    [~,D_damp(:,:,i)] = eig(gains(:,3:4,i));
    [V_stiff_sym,D_stiff_sym] = eig((gains(:,1:2,i)+gains(:,1:2,i)')/2);
    stiffness_sym(i,1) = abs(D_stiff_sym(1,1)*D_stiff_sym(2,2));
    [d,ind] = sort(abs(diag(D_stiff_sym)));
    Vs_stiff_sym = V_stiff_sym(:,ind);
    V_princ = Vs_stiff_sym(:,2);
    if V_princ(1) * V_princ(2) > 0
        V_princ = abs(V_princ);
        stiffness_sym(i,2) = real(acosd(max(min(real(dot(V_princ,[1,0])),1),-1)));
    else
        if V_princ(1) < 0
            V_princ = -V_princ;
        end
        stiffness_sym(i,2) = -real(acosd(max(min(real(dot(V_princ,[1,0])),1),-1)));
    end
    [V_stiff_asym(:,:,i),D_stiff_asym] = eig((gains(:,1:2,i)-gains(:,1:2,i)')/2);
    stiffness_asym(i) = abs(D_stiff_asym(1,1)*D_stiff_asym(2,2));
    [V_damp_sym,D_damp_sym] = eig((gains(:,3:4,i)+gains(:,3:4,i)')/2);
    damping_sym(i) = abs(D_damp_sym(1,1)*D_damp_sym(2,2));
    [V_damp_asym,D_damp_asym] = eig((gains(:,3:4,i)-gains(:,3:4,i)')/2);
    damping_asym(i) = abs(D_damp_asym(1,1)*D_damp_asym(2,2));
end
filename = 'eigen_bestparams';
eigen = load(sprintf('%s.mat',filename));
folder = fullfile('Data','Eigenstructure');
% addpath(folder);
file = sprintf('%s.mat',filename);
folder = fullfile(folder,file);
%Save data as struct
% eigen.GainsOld = GainsOld;
% eigen.GainsYoung = GainsYoung;
eigen.StiffnessSymEigOld = stiffness_sym;
eigen.StiffnessAsymEigOld = stiffness_asym;
eigen.DampingSymEigOld = damping_sym;
eigen.DampingAsymEigOld = damping_asym;
save(folder, '-struct', 'eigen');
%% Run Simulation for n trials
alpha = 10^6;
beta = [0.19,0.25];
% alpha = 10.^(-4:1:6);
% beta = [0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3];
% D_ratio = zeros(length(alpha),length(beta),2);
% V_angle = zeros(length(alpha),length(beta),2);
Ds = zeros(4,length(beta),length(alpha));
Vs = zeros(4,4,length(beta),length(alpha));
Ks = zeros(2,4,length(beta),length(alpha));
count = 1;
for a = 1:length(alpha)
    for b = 1:length(beta)
        input.Controller.alpha = alpha(a);
        input.Controller.beta = beta(b);
        clear test_model
        [~, ~, model] = test_model(input);
        Ks(:,:,b,a) = model.K
        A_cl = model.A_ol-model.B*model.K;
        
%         figure(1);
%         drawStiffnessEllipseSym(model.K,1,[],[],[],[],length(alpha));
%         figure(2);
%         drawStiffnessEllipseAsym(model.K,2,[],[],[],[],length(alpha));
%         colorbar
%         cmap = colormap(parula(32)) ; %Create Colormap
%         set(gca,'CLim',[0.01 3]);
%         cbh = colorbar ; %Create Colorbar
%         cbh.TickLabels = [0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3];    %Replace the labels of these 8 ticks with the numbers 1 to 8
%         cbh.FontSize = 10;
%         cbh.Ticks = linspace(0,3,32); %Create 8 ticks from zero to 1
% 
        if count == 1
            [V_op,D_op] = eig(model.A_ol);
        end
        [V_cl,D_cl] = eig(A_cl)
        [D_cl,I] = sort(diag(D_cl))
        V_cl = V_cl(:, I)
        Ds(:,b,a) = D_cl;
        Vs(:,:,b,a) = V_cl;
%         
% %         figure();
% %         sys_open = zpk([],diag(D_op)./(2*pi)',1);
% %         sys_closed = zpk([],diag(D_cl)./(2*pi)',1);
% %         h = pzplot(sys_open,'b',sys_closed,'r');
% %         xlabel('Real Axis')
% %         ylabel('Imaginary Axis')
% %         legend('Open Loop','Closed Loop')
% %         h.AxesGrid.XUnits = 'Hz';
% %         h.AxesGrid.YUnits = 'Hz';
% %         grid on
%         

%         D_ratio(a,b,:) = [D_cl(1)/D_cl(2);D_cl(3)/D_cl(4)];
%         V_angle(a,b,:) = [real(acosd(max(min(real(dot(V_cl(:,1),V_cl(:,2))),1),-1)));
%                           real(acosd(max(min(real(dot(V_cl(:,3),V_cl(:,4))),1),-1)))];
        alpha(a)
        beta(b)
    end
end
% filename = 'eigen_sweep';
% eigen = load(sprintf('%s.mat',filename));
% folder = fullfile('Data','Eigenstructure');
% % addpath(folder);
% file = sprintf('%s.mat',filename);
% folder = fullfile(folder,file);
% %Save data as struct
% % eigen.D_ratio = D_ratio;
% % eigen.V_angle = V_angle;
% eigen.Values = Ds;
% eigen.Vectors = Vs;
% eigen.Gains = Ks;
% save(folder, '-struct', 'eigen');
%% Plot Results
figure();
subplot(1,2,1)
semilogx(alpha,squeeze(Ds(1,12,:)))
hold on
semilogx(alpha,squeeze(Ds(2,12,:)))
xlabel('\alpha')
ylabel('Eigenvalue')
title('First Pole Pair')
subplot(1,2,2)
semilogx(alpha,squeeze(Ds(3,12,:)))
hold on
semilogx(alpha,squeeze(Ds(4,12,:)))
xlabel('\alpha')
ylabel('Eigenvalue')
title('Second Pole Pair')
improvePlot;
%%
V_0 = squeeze(Vs(:,:,12,1));
V_anglediff = zeros(4,length(alpha));
for a = 1:length(alpha)
    for i = 1:4
        if sign(V_0(1,i)) == sign(Vs(1,i,12,a))
            V = Vs(:,i,12,a);
        else
            V = -Vs(:,i,12,a);
        end
        V_anglediff(i,a) = real(acosd(max(min(real(dot(V_0(:,i),V)),1),-1)));
    end
end
figure();
subplot(1,2,1)
semilogx(alpha,V_anglediff(1,:))
hold on
semilogx(alpha,V_anglediff(2,:))
xlabel('\alpha')
ylabel('Eigenvector''s Angle Difference from \alpha = 10^{-4} Eigenvector (deg)')
title('First Pole Pair') 
% subplot(1,3,2) 
% plot(alpha(10:end),V_anglediff(1,10:end)) 
% hold on 
% plot(alpha(10:end),V_anglediff(2,10:end)) 
% xlabel('\alpha') 
% ylabel('Eigenvector''s Angle Difference from \alpha = 10^{-4} Eigenvector (deg)')
% title('First Pole Pair (Zoomed In)') 
subplot(1,2,2) 
semilogx(alpha,V_anglediff(3,:)) 
hold on 
semilogx(alpha,V_anglediff(4,:)) 
xlabel('\alpha') 
ylabel('Eigenvector''s Angle Difference from \alpha = 10^{-4} Eigenvector (deg)')
title('Second Pole Pair')
improvePlot;
%%
close all;
for a = 1:length(alpha)
    figure(1);
    drawStiffnessEllipseSym(squeeze(Ks(:,:,1,a)),1,[],[],[],[],length(alpha))
    colorbar
    cmap = colormap(parula(length(alpha))) ; %Create Colormap
    set(gca,'CLim',[1e-4 1e6]);
    cbh = colorbar ; %Create Colorbar
    cbh.TickLabels = 10.^(-4:1:6);
    cbh.FontSize = 10;
    cbh.Ticks = linspace(1e-4,1e6,11); %Create 8 ticks from zero to 1
    improvePlot;
%     subplot(1,2,1)
%     ylim([-1000,1000])
%     xlim([-1000,1000])
%     subplot(1,2,2)
%     ylim([-400,400])
%     xlim([-400,400])
%     colorbar
%     cmap = colormap(parula(length(alpha))) ; %Create Colormap
%     set(gca,'CLim',[0.01 3]);
%     cbh = colorbar ; %Create Colorbar
%     cbh.TickLabels = 10.^(-4:1:6);
% %     cbh.TickLabels = [0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3];    %Replace the labels of these 8 ticks with the numbers 1 to 8
%     cbh.FontSize = 10;
%     cbh.Ticks = linspace(0,3,length(alpha)); %Create 8 ticks from zero to 1
%     improvePlot;
%     
    figure(2);
    drawStiffnessEllipseAsym(Ks(:,:,12,a),2,[],[],[],[],length(alpha));
    colorbar
    cmap = colormap(parula(length(alpha))) ; %Create Colormap
    set(gca,'CLim',[1e-4 1e6]);
    cbh = colorbar ; %Create Colorbar
    cbh.TickLabels = 10.^(-4:1:6);
    cbh.FontSize = 10;
    cbh.Ticks = linspace(1e-4,1e6,11); %Create 8 ticks from zero to 1
    improvePlot;
    subplot(1,2,1)
    ylim([-550,550])
    xlim([-550,550])
    subplot(1,2,2)
    ylim([-175,175])
    xlim([-175,175])
end
%%
Ks = Gains;
for a = 1:length(beta)
    figure(7);
    drawStiffnessEllipseSym(squeeze(Ks(:,:,a,11)),7,[],[],[],[],length(beta))
    colorbar
    cmap = colormap(parula(length(beta))) ; %Create Colormap
    set(gca,'CLim',[0.01 3]);
    cbh = colorbar ; %Create Colorbar
    cbh.TickLabels = [0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3];    %Replace the labels of these 8 ticks with the numbers 1 to 8
    cbh.FontSize = 10;
    cbh.Ticks = linspace(0,3,length(beta)); %Create 8 ticks from zero to 1
    improvePlot;
    subplot(1,2,1)
    ylim([-2500,2500])
    xlim([-2500,2500])
    subplot(1,2,2)
    ylim([-1000,1000])
    xlim([-1000,1000])

%     
    figure(8);
    drawStiffnessEllipseAsym(Ks(:,:,a,11),8,[],[],[],[],length(beta));
    colorbar
    cmap = colormap(parula(length(beta))) ; %Create Colormap
    set(gca,'CLim',[0.01 3]);
    cbh = colorbar ; %Create Colorbar
    cbh.TickLabels = [0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3];    %Replace the labels of these 8 ticks with the numbers 1 to 8
    cbh.FontSize = 10;
    cbh.Ticks = linspace(0,3,length(beta)); %Create 8 ticks from zero to 1
    improvePlot;
    subplot(1,2,1)
    ylim([-1500,1500])
    xlim([-1500,1500])
    subplot(1,2,2)
    ylim([-300,300])
    xlim([-300,300])
end
%% Compute K_ratio
StiffnessAsym = StiffnessAsymYoung;
StiffnessSym = StiffnessSymYoung;
DampingAsym = DampingAsymYoung;
DampingSym = DampingSymYoung;
numsubs = 65;
for i = 1:numsubs
    K_ratio_sqrt(i,1) = sqrt(abs(det(StiffnessAsym(:,:,i))))/sqrt(abs(det(StiffnessSym(:,:,i))))*100;
    K_ratio(i,1) = abs(det(StiffnessAsym(:,:,i)))/abs(det(StiffnessSym(:,:,i)))*100;
    B_ratio_sqrt(i,1) = sqrt(abs(det(DampingAsym(:,:,i))))/sqrt(abs(det(DampingSym(:,:,i))))*100;
    B_ratio(i,1) = abs(det(DampingAsym(:,:,i)))/abs(det(DampingSym(:,:,i)))*100;
end