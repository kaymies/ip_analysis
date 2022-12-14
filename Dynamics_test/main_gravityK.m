function [x,Dx,A_cl,K,A_ol,B_ol,f_IP,IP_ratio,output_dynamics] = ...
    main_gravityK(input,K_ref,ifig,mkrcolor)
% Based on main_KS.m from Kaymie's project
%
% Rika Sugimoto Dimitrova
% 2020-12-18
% Modified:
% 2021-01-05 - modified code to also enable simulation of single inverted pendulum
% 2021-01-12 - modified to include simDuration as input parameter
% 2021-02-05 - added pose to input, to access dynamics files in
%              Dynamics_gen_JW
% 2021-02-09 - added K_ref as input, to use a pre-computed gain if desired
% 2021-04-20 - experimenting with noise
% 2021-05-27 - cleaned up IP calc & plot code;
%            - adding f_IP and IP_ratio to output.
% 2021-06-01 - adding output_dynamics to output.
% 2021-06-21 - adding COM_z and COM_x to output_dynamics.

%% Model Parameters
if nargin ~=0
    TotalMass = input.TotalMass; 
    TotalHeight = input.TotalHeight;    
    gender = input.gender;
    plane = input.plane;
    model_type = input.model;
    pose = input.pose; % added 2021-02-05
    Fs_Hz = input.FreqSampKin;
    motorNoiseLvL = input.simMotorNoise;
    sensorNoiseLvL = input.simSensorNoise;
    t_0 = 0; t_f = input.trialDuration;
    sigma_r = input.simMotorNoiseRatio;
    sensorNoiseRatio = input.simSensorNoiseRatio;
    if isfield(input,'noise_type')
        noise_type = input.noise_type;
    else
        noise_type = 'w';
    end
else
    TotalMass = 71.13; 
    TotalHeight = 1.75;    
    gender = 'M';
    plane = 'sgt';
    model_type = 'DIP'; % 'SIP','DIP'
    pose = 'T-pose';
    Fs_Hz = 100; % Sampling frequency
    motorNoiseLvL = 0.01; % changed default from 10 to 0.01 (2021-05-27)
    sensorNoiseLvL = 0;
    t_0 = 0; t_f = 30;
    sigma_r = 1.6;
    noise_type = 'w';
end
model_param = [TotalMass; TotalHeight];
addpath([pwd,'\Dynamics_gen_JW\',pose,'\',model_type,'\',gender,'_',plane]);

% lumped model parameters for the table
% switch model_type
%     case 'SIP'
%         [m_1,c_1,j_1,L_1] = auto_LumpParam(model_param);
%         c_1 = c_1(2); L_1 = L_1(2);
%         [m_1,c_1,j_1,L_1]
%     case 'DIP'
%         [m_1,c_1,j_1,L_1,m_2,c_2,j_2,L_2] = auto_LumpParam(model_param);
%         c_1 = c_1(2); c_2 = c_2(2); L_1 = L_1(2); L_2 = L_2(2);
%         [m_1,c_1,j_1,L_1,m_2,c_2,j_2,L_2]
% end

%% Unwrap Input
if nargin ~=0 % if there is an input into the function  
    struct_Controller = input.Controller;
    isAnimOn = input.PostProc.AnimOn;
    isPlotOn = input.PostProc.PlotOn;    
else 
    % LQR
    alpha = 1; gamma = 1e6; beta = 0.3;
    param.Q = alpha*eye(4); param.R = gamma*diag([beta, 1/beta]);
    
    struct_Controller.type = 'LQR'; % designate controller type used here
    struct_Controller.param = param;
        
    isAnimOn = 1; % change to 0 if you don't want to see animation
    isPlotOn = 1; % change to 0 if you don't want to see plots   
end
%% Determine Controller - defined below
switch struct_Controller.type 
    case 'JIC'
        FuncHandle_Controller = @JIC;      
    case 'CIC'
        FuncHandle_Controller = @CIC;        
    case 'LQR'
        FuncHandle_Controller = @LQR;
end
%% Simulation Setup
% initial condition
% q: joint ('relative, not spatial') angles
switch model_type
    case 'SIP'
        % q=[q_ankle]
        q_0 = [0];
        Dq_0 = [0];
    case 'DIP'
        % q=[q_ankle, q_hip]
        q_0 = [0;0];
        %q_0 = [2;-8]*pi/180; %2021-05-24
        Dq_0 = [0;0];
end
n_q = length(q_0); % number of joints

% x: state = [q;Dq]
x_0 = [q_0; Dq_0]; % initial state
n_x = length(x_0); % number of states

%% Run simulation
% Here I decided to use semi-implicit Euler integrator due to its ease of
% handling noise and calculating output and Dx simultaneously. 
% later on, when you add time delay, semi-implicit Euler is easier to
% handle than ode45. (you cannot use ode45 but should use dde23 for
% instance) - Jongwoo
dt = 1/Fs_Hz;    % timestep
t = t_0:dt:t_f; % time
N = length(t);  % sample number

% memory allocation 
x = zeros(n_x, N);    % 4xN matrix - trajectory 
Dx = zeros(n_x, N);
y = cell(1, N);       % 1x Nt structures - output 
ankle = zeros(1,N); % ankle torque
ankle_noise = zeros(1,N); % ankle actuator noise

% Uncorrelated Gaussian noise
sensorNoise_white = sensorNoiseLvL*randn([n_x,N]);
motorNoise_white = motorNoiseLvL*randn([n_q,N]);
sensorNoise(1:n_q,:) = sensorNoiseRatio*sensorNoise_white(1:n_q,:); % Scale joint positions
sensorNoise(n_q+1:n_x,:) = sensorNoise_white(n_q+1:n_x,:);          % Leave joint velocities unscaled
motorNoise = motorNoise_white;

switch noise_type
    case 'b'
        % Integrate white noise % 2021-05-05 % 2022-01-17 - changed to q only (no dq, no motor noise)
        sensorNoise(1:n_q,:) = sensorNoiseRatio*sensorNoise_white(1:n_q,:);
        sensorNoise((n_q+1):end,:) = sensorNoise_white((n_q+1):end,:);
        motorNoise = motorNoise_white;
        for i = 2:length(t)
        sensorNoise(1:n_q,i) = sensorNoise(1:n_q,i) + sensorNoise(1:n_q,i-1);
        sensorNoise((n_q+1):end,i) = sensorNoise((n_q+1):end,i) + sensorNoise((n_q+1):end,i-1);
        motorNoise(:,i) = motorNoise(:,i) + motorNoise(:,i-1);
        end
    case 'p'
        pinkNoise = dsp.ColoredNoise(1,N,n_q);
        sensorNoise = sensorNoiseLvL*[pinkNoise()';pinkNoise()'];
        motorNoise = motorNoiseLvL*pinkNoise()';
end

% % low-passed white noise
% lpf_fpass_Hz = 1/80*2*pi;% Peterka's: 80s time constant %0.1;  
% sensorNoise = lowpass(sensorNoise_white',lpf_fpass_Hz,Fs_Hz);
% motorNoise = lowpass(motorNoise_white',lpf_fpass_Hz,Fs_Hz);
% sensorNoise = sensorNoise';
% motorNoise = motorNoise';

% % 1/f noise
% pinkNoise = dsp.ColoredNoise(1,N,n_q);
% sensorNoise = sensorNoiseLvL*[pinkNoise()';pinkNoise()'];
% motorNoise = motorNoiseLvL*pinkNoise()';

% % 1/f^2 noise
% brownianNoise = dsp.ColoredNoise(2,N,n_q);
% sensorNoise = sensorNoiseLvL*[brownianNoise()';brownianNoise()'];
% motorNoise = motorNoiseLvL*brownianNoise()';

% run simulation
x_i = x_0; % initialize state
if nargin < 2
        K_ref = [];
end
[u_i,A_cl,K] = FuncHandle_Controller(...
    t(1), x_i, struct_Controller.param, model_param, [], model_type, K_ref);
[A_ol,B_ol] = getLinearDynamics(model_param, model_type);

for i = 1:N
    t_i = t(i);
    % calculate Dx and output y at time t; also output ankle torque input
    % and ankle noise for graphing purposes later
    u_i = FuncHandle_Controller(...
        t_i, x_i, struct_Controller.param, model_param, sensorNoise(:,i), model_type);
    [Dx_i, y_i, ankle_i, ankle_noise_i] = func_Dynamics(...
        t_i, x_i, u_i, model_param, motorNoise(:,i), model_type, sigma_r); % these functions are all below
    % log state and output
    x(:, i) = x_i;
    Dx(:, i) = Dx_i;
    y{i} = y_i;
    ankle(i) = ankle_i;
    ankle_noise(i) = ankle_noise_i;
    
    % integrating
    % unwrapping states
    q_t = x_i(1:n_q, 1);
    Dq_t = x_i(n_q+1:end, 1);
    DDq_t = Dx_i(n_q+1:end, 1);
    % semi-implicit Euler integration
    Dq_next = Dq_t + dt*DDq_t;
    q_post = q_t + dt*Dq_next;
    x_i = [q_post; Dq_next];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% post process
if isAnimOn
    addpath([pwd,'\..\KS - MATLAB Balance No Cane Template\Dynamics_test']);
    b_saveanim= 0; 
    title_ = 'test';
    t_anim = t(1:100:end); 
    x_anim = x(:, 1:100:end);
    postproc_animation(t_anim, x_anim, model_param, model_type, b_saveanim, title_ ); 
    rmpath([pwd,'\..\KS - MATLAB Balance No Cane Template\Dynamics_test']);
end

%-- COMPUTE IP --%
COM_z = zeros(1,N); % vertical CoM
COM_x = zeros(1,N); % horizontal CoM % 2021-06-21
MLCOP = zeros(1,N); % CoP
MLFx = zeros(1,N); % horizontal ground reactin force
MLFz = zeros(1,N); % vertical ground reaction force
% renaming component of output y from simulation above
for i = 1:N
    MLCOP(:,i) = y{i}.COP;
    COM_z(:,i) = y{i}.COM_z;
    COM_x(:,i) = y{i}.COM_x; % 2021-06-21
    MLFx(:,i) = y{i}.Fx;
    MLFz(:,i) = y{i}.Fz;
end
MLCOP = MLCOP - mean(MLCOP); % detrend CoP data for PCA later

% Obtain IP at different frequencies (inspired by Marta's code)
f_i = 0.4; f_int = 0.2; f_range = 7.4; % starting freq, interval between tested frequencies, frequency range
f = f_i:f_int:f_i+f_range; % frequencies are defined so that the center point is at 0.5 to 7.9 Hz
IP = zeros(1,length(f)); % IP
theta = -MLFx./MLFz;
theta = theta - mean(theta);
cop_ham = MLCOP'.*hamming(length(MLCOP)); % apply Hamming window to CoP data
theta_ham = theta'.*hamming(length(theta)); % apply Hamming window to the force angle
for n = 1:length(f) % for all frequencies designated above:
    [B,A] = butter(2,[f(n) f(n)+f_int]/(1/dt/2)); % define the Butterworth filter
    COP_f = filtfilt(B,A,cop_ham); % apply filter to CoP data
    theta_f = filtfilt(B,A,theta_ham); % apply filter to force angle data
    [coeff,~,~,~,~,~] = pca([COP_f theta_f]); % PCA
    IP(n) = coeff(1,1)/coeff(2,1); % compute IP and store
end
IP_ratio = IP/mean(COM_z);
mean(COM_z);

f_IP = f + 0.1;
output_dynamics.COP = MLCOP;
output_dynamics.Fx = MLFx;
output_dynamics.Fz = MLFz;
output_dynamics.thetaF = theta;
output_dynamics.COM = [COM_x; COM_z]; % 2021-06-21

if isPlotOn
    % Plot a scatter plot of IP frequency response
    figure(ifig);
    if nargin > 3 && ~isempty(mkrcolor)
        scatter(f_IP, IP_ratio,10,mkrcolor,'filled') % 0.87 m is the height of the CoM and it stays roughly constant throughout the entire simulation
    else
        scatter(f_IP, IP_ratio,10,'k','filled') % 0.87 m is the height of the CoM and it stays roughly constant throughout the entire simulation    
    end
    ylabel('IP/CoM')
    xlabel('Frequency [Hz]')
end

%% sub function definitions
% Controllers 
% (1) Joint stiffness %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function u_JIC = JIC(t, x, param)
        % unwrap states
        q = x(1:2,:);
        Dq = x(3:4,:);
        % controller        
%         err = q-q_des; Derr = Dq-Dq_des;
        err = q; Derr = Dq;
        u_JIC = -param.Kp*err + -param.Kd*Derr;        
    end
% (2) Cartesian stiffness %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function u_CIC = CIC(t, x, param)
        % unwrap states
        q = x(1:2,:);
        Dq = x(3:4,:);      
        % Jacobian
        [~,~,c] = auto_COMPosition(q);  % COM position
        [J_c,~,~] = auto_Jacobian(q);   % COM Jacobian        
        Dc = J_c*Dq;
        Dc_0 = [0;0];        
        % controller                
        err = c-param.c_0; Derr = Dc - Dc_0;
        u_CIC = J_c'*(-param.Kp*err + -param.Kd*Derr);  
    end
% (3) LQR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [u_LQR, A_cl, K] = LQR(t, x, param, model_param, sensorNoise, model_type, K_ref)
        if nargin < 5 || isempty(sensorNoise)
            sensorNoise = zeros(size(x));
        end
        if nargin < 6 || isempty(model_type)
            model_type = 'DIP';
        end
        Q = param.Q; R = param.R;
        persistent K_lqr; % calculate K_lqr only once for all time
        %persistent K_pp; % pole placement gain (either use K_lqr OR K_pp)
        if isempty(K_lqr) && (nargin < 7 || isempty(K_ref)) % finding K_lqr for u = -Kx
            [A_lin,B_lin] = getLinearDynamics(model_param, model_type);
            %[V_op,D_op] = eig(A_lin); % linearized system open loop poles
            %D_op_hz = diag(D_op)./(2*pi);
            % design lqr control
            K_lqr = lqr(A_lin,B_lin,Q,R);
            %B_s = 0.5*(K_lqr(:,3:4) + (K_lqr(:,3:4))');
            %B_a = 0.5*(K_lqr(:,3:4) - (K_lqr(:,3:4))');
            %B_s_norm = norm(B_s);
            %B_a_norm = norm(B_a);
            % compute closed-loop poles
            %A_cl = A_lin-B_lin*K_lqr;
            %[V_cl,D_cl] = eig(A_cl);
            %D_cl_hz = diag(D_cl)./(2*pi);
        elseif isempty(K_lqr) && nargin > 6 && ~isempty(K_ref)
            K_lqr = K_ref;
        end
        K = K_lqr;
        [A_lin,B_lin,Kg] = getLinearDynamics(model_param, model_type);

        %if ( abs(x(1)) < 1*pi/180 ) && ( abs(x(2)) < 1*pi/180 )
        %   K(1:2,1:2) = -Kg;            
        %end
        %K(:,3:4) = K(:,3:4)/2;

        A_cl = A_lin-B_lin*K;
        
        q_eq = [0; 0]*pi/180; Dq_eq = [0;0]; %2021-05-24
        x_eq = [q_eq; Dq_eq];
        u_LQR = -K*((x-x_eq)+sensorNoise);
        %u_LQR = -K_lqr*(x+[1;1;0.01;0.001].*sensorNoiseLvL.*randn(size(x)));
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plant Dynamics
    function [Dx, output_struct, ankle, ankle_noise] = func_Dynamics(t, x, u_ctrl, param, motorNoise, model_type, sigma_r)
        
        if nargin < 5 || isempty(motorNoise)
            motorNoise = 10*randn(length(q),1);
        end
        if nargin < 6 || isempty(model_type)
            model_type = 'DIP';
        end
        if nargin < 7 || isempty(sigma_r)
            sigma_r = 1.6; % ankle:hip noise ratio; edit as necessary (2021-02-05)
        end
        switch model_type
            case 'SIP'
                q = x(1,:);
                Dq = x(2,:);
            case 'DIP'
                % unwrap states
                q = x(1:2,:);
                Dq = x(3:4,:);
        end
        % eom Mq'' + Cq' + G = u_ctrl + u_noise
        [M, C, G, ~] = auto_Dynamics(q, Dq, param);
        %G = G*0; % 2022-07-09 Simulate no gravity
        
        % input noise, multiplier can be changed
        u_noise = motorNoise; % noise size can be determined by you
        u_noise(1) = sigma_r*u_noise(1);
        ankle_noise = u_noise(1);
        u = u_ctrl + u_noise;
        ankle = u(1);
        
        % solve accl 
        DDq = M\(-C*Dq -G + u);
        Dx = [Dq; DDq];
        
        % outputs
        if nargout > 1           
            % vertical force
            [Fx,Fz] = auto_GroundReactionForce(q,Dq,DDq, param);
            
            output_struct.Fz = Fz;
            output_struct.Fx = Fx;
            
            % theta_F
            output_struct.theta_F = Fx/Fz;
            
            % COP
            output_struct.COP = u(1)/Fz;
            
            % COMz
            p_O_CoM = auto_COMinformation(q,Dq,param);
            output_struct.COM_z = p_O_CoM(2);
            output_struct.COM_x = p_O_CoM(1); % 2021-06-21
            
            % Kaymie point (intersection point)
            output_struct.KP = -Fz/Fx*output_struct.COP;            
        end 
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linearized Dynamics Matrices (Rika 2020-12-18)
    function [A_lin,B_lin,Kg] = getLinearDynamics(model_param,model_type)
        switch model_type
            case 'SIP'
                q_eq = [0]; Dq_eq = [0];
            case 'DIP'
                % equilibrium position
                q_eq = [0; 0]; Dq_eq = [0;0];
                q_eq = [0; 0]*pi/180; Dq_eq = [0;0]; %2021-05-24
        end
        % linearize dynamics
        [M_eq, ~, ~, ~] = auto_Dynamics(q_eq, Dq_eq, model_param);
        [~,Jac_G_eq] = auto_Jacobian(q_eq, Dq_eq, model_param);
        %Jac_G_eq = Jac_G_eq*0; % 2022-07-09 Simulate no gravity
        switch model_type
            case 'SIP'
                A_lin = [0, 1; 
                    -inv(M_eq)*Jac_G_eq, 0];
                B_lin = [0; inv(M_eq)];
            case 'DIP'
                A_lin = [zeros(2),  eye(2); 
                    -inv(M_eq)*Jac_G_eq, zeros(2)];
                B_lin = [zeros(2,2); inv(M_eq)];
        end
        Kg = Jac_G_eq;
    end

rmpath([pwd,'\Dynamics_gen_JW\',pose,'\',model_type,'\',gender,'_',plane]);

end