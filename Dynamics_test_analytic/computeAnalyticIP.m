function [IP_ratio] = computeAnalyticIP(input)
% Compute intersection-point curve over frequency analytically using the
% cross-spectral-density (CPSD) matrix between COP and GRF angle
%
% Inputs:
% params = struct containing parameters for IP analysis
%     'TotalMass':   in kg
%     'TotalHeight': in meters
%     'gender': 'F','M'
%     'plane' : 'sgt','frt'
%     'model' : 'SIP','DIP'
%     'pose'  : 'T-pose','Arms_adducted_original'
%     'Controller.type' : 'LQR','ICFV'
%     'Controller.param': struct, including:
%                         - Q, R: matrices for LQR
%                         - delay
%     'simMotorNoiseRatio' = ankle-to-hip motor noise ratio (sigma_r = 0 -> hip noise only)
%     'simSensorNoiseRatio' = ankle-to-hip forward-path noise ratio
%     'simMotorNoise' = motor noise level
%     'simSensorNoise' = sensor noise level
%
% Rika Sugimoto Dimitrova
% 2021-05-28
%
% Modified:
% 2021-06-01 - corrected calculation errors; now the results are the same
%              as Jongwoo's analytic method in ../JW - Analytic IP/main2D.m
% 2023-03-28 - updated to have sigma_r defined the same as Kaymie's
% 2023-04-04 - new analytic IP calculation based on cross spectrum
% 2023-04-05 - added functionality to compute analytic IP for the ICFV
%              controller
%            - updated the function arguments and outputs
% 2023-05-18 - added effect of delay
% 2023-06-01 - cleaning up code
%            - feedback gains computed within code (no need to pass as arg)

%% Load data
% Input parameters
if nargin ~=0
    TotalMass = input.TotalMass; 
    TotalHeight = input.TotalHeight;    
    gender = input.gender;
    plane = input.plane;
    model_type = input.model;
    pose = input.pose;
    f = input.Frequency;
else
    TotalMass = 68.5; 
    TotalHeight = 1.66;    
    gender = 'M';
    plane = 'sgt';
    model_type = 'DIP'; % 'SIP','DIP'
    pose = 'pose_T'; % 'pose_T' 'pose_I'
    f_i = 0.5; f_int = 0.2; f_end = 7.9;
	f = f_i:f_int:f_end;
end
%% Define Controller
if nargin ~=0
    struct_Controller = input.Controller;
    alpha = struct_Controller.alpha;
    beta = struct_Controller.beta;
    gamma = struct_Controller.gamma;
    kappa = struct_Controller.kappa;
    eta = struct_Controller.eta;
    sigma_r = input.NoiseRatio;
    sigma_f = input.simSensorNoiseRatio;
    motorNoiseLvL = input.simMotorNoise;
    fpNoiseLvL = input.simSensorNoise;
    delay = input.Controller.param.delay; % 2023-05-18
    controller_type = struct_Controller.type;

else
    alpha = 1e6;
    beta = 0.3;
    gamma = 1;
    kappa = 1;
    eta = 1;
    sigma_r = 0.9;
    sigma_f = 1;
    motorNoiseLvL = 10;
    fpNoiseLvL = 0;
    delay = 0;
    controller_type = 'LQR';
end
switch model_type
    case 'SIP'
        param.R = alpha*beta;
    case 'DIP'
        param.R = alpha*diag([beta, 1/beta]);
end
param.Q = gamma*[kappa 0 0 0;
                 0 1/kappa 0 0;
                 0 0 eta 0;
                 0 0 0 1/eta];
struct_Controller.param = param;


addpath([pwd,'/AutoDerived/',pose,'/',model_type,'/',gender,'_',plane]);

% Frequency over which to compute IP curve
w = 2*pi*f;

%% Generate model dynamics equations, A & B
% linearize dynamics
model_param = [TotalMass; TotalHeight];
switch model_type
    case 'SIP'
        n_q = 1;
    case 'DIP'
        n_q = 2;
end
q_eq = zeros(n_q,1); Dq_eq = zeros(n_q,1);
[M_eq,~,~,~] = auto_Dynamics(q_eq, Dq_eq, model_param);
[J_CoM_eq,J_G_eq] = auto_Jacobian(q_eq, Dq_eq, model_param);
[DJ_CoM_eq] = auto_DJacobian_CoM(q_eq, Dq_eq, model_param);

A_lin = [zeros(n_q), eye(n_q); 
    -inv(M_eq)*J_G_eq, zeros(n_q)];
B_lin = [zeros(n_q); inv(M_eq)];

[m_1,c_1,~,L_1,m_2,c_2,~,~] = auto_LumpParam(model_param);
% c_1 = c_1(2);c_2 = c_2(2);L_1 = L_1(2);
L_COM = (m_1*c_1 + m_2*(L_1+c_2))/(m_1+m_2);

% C_lin and D_lin based on Jongwoo's compact calcs:
C_1 = [0 0 0 0];
D_1 = [1 0];

J_CoM_x = J_CoM_eq(1,:);
dJ_CoM_x = DJ_CoM_eq(1,:);
J_2 = -(m_1+m_2)*[dJ_CoM_x J_CoM_x];
C_2 = J_2*A_lin;
D_2 = J_2*B_lin;
C_lin = [C_1; C_2]; D_lin = [D_1; D_2];

K_x = lqr(A_lin,B_lin,struct_Controller.param.Q,struct_Controller.param.R); % 2023-06-01

%%
s = tf('s');
A_CL = A_lin - B_lin * K_x * exp(-s*delay); % 2023-05-18: testing with delay
C_CL = C_lin - D_lin * K_x * exp(-s*delay);
B_CL = [B_lin, -B_lin*K_x];
D_CL = [D_lin, -D_lin*K_x];

if strcmp(controller_type,'LQR')
    H = C_CL * inv(s*eye(2*n_q) - A_CL) * B_lin + D_lin;
    
    Gww = [sigma_r^2 0; 0 1];

    % New computation method, introduced to work with delays (2023-05-18)
    [H_mag,H_phase,~] = bode(H,w);
    [HT_mag,HT_phase,~] = bode(Gww*(H.'),w);
    HT_frf = HT_mag.*exp(1j*pi/180*HT_phase);
    HC_frf = H_mag.*exp(-1j*pi/180*H_phase);

    IP_ratio = f*0;
    for i = 1:length(w)
        Gcpsd_frf(:,:,i) = HC_frf(:,:,i)*HT_frf(:,:,i); % added for with-delay computations 2023-05-19
        [V,D] = eig(real(Gcpsd_frf(:,:,i)));
        [~,i_max] = max(diag(abs(D)));
        IP_ratio(i) = V(1,i_max)/V(2,i_max) / L_COM;
    end

%%
elseif strcmp(controller_type,'ICFV')
    H = C_CL * inv(s*eye(2*n_q) - A_CL) * B_CL + D_CL;
    
    w1 = motorNoiseLvL*sigma_r;
    w2 = motorNoiseLvL;
    n1 = fpNoiseLvL*sigma_f;
    n2 = fpNoiseLvL;
    n3 = fpNoiseLvL*sigma_f;
    n4 = fpNoiseLvL;
    Gww_sqrt = [w1 0 0 0 0 0;
           0 w2 0 0 0 0;
           0 0 n1/s 0 0 0;
           0 0 0 n2/s 0 0;
           0 0 0 0 n3 0;
           0 0 0 0 0 n4];
    Gww = Gww_sqrt'*Gww_sqrt;
    
    % New computation method, introduced to work with delays (2023-05-18)
    [H_mag,H_phase,~] = bode(H,w);
    [HT_mag,HT_phase,~] = bode(Gww*(H.'),w);
    HT_frf = HT_mag.*exp(1j*pi/180*HT_phase);
    HC_frf = H_mag.*exp(-1j*pi/180*H_phase);

    IP_ratio = f*0;
    for i = 1:length(w)
        Gcpsd_frf(:,:,i) = HC_frf(:,:,i)*HT_frf(:,:,i); % added for with-delay computations 2023-05-19
        [V,D] = eig(real(Gcpsd_frf(:,:,i)));
        [~,i_max] = max(diag(abs(D)));
        IP_ratio(i) = V(1,i_max)/V(2,i_max) / L_COM;
    end
end

%% Clean up
rmpath([pwd,'/AutoDerived/',pose,'/',model_type,'/',gender,'_',plane]);

end