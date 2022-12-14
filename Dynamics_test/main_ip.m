function [f, IP_ratio] = main_ip(input)
%% Model Parameters and Set Up
if nargin ~=0
    TotalMass = input.TotalMass; 
    TotalHeight = input.TotalHeight;    
    gender = input.gender;
    plane = input.plane;
    model_type = input.model;
    pose = input.pose;
    Fs_Hz = input.FreqSampKin; % sampling frequency
    t_f = input.trialDuration;
    coord = input.CoordinateFrame;
else
    TotalMass = 71.1; 
    TotalHeight = 1.75;    
    gender = 'M';
    plane = 'sgt';
    model_type = 'DIP'; % 'SIP','DIP'
    pose = 'pose_I'; % 'pose_T' 'pose_I'
    Fs_Hz = 100;
    t_f = 15;
    coord = 'spatial'; % 'relative','spatial'
end
model_param = [TotalMass; TotalHeight];
addpath([pwd,'\dynamics_models\',pose,'\',model_type,'\',gender,'_',plane]);

% lumped model parameters
% switch model_type
%     case 'SIP'
%         [m_1,c_1,j_1,L_1] = auto_LumpParam(model_param);
%         c_1 = c_1(2); L_1 = L_1(2);
%         [m_1,c_1,j_1,L_1];
%     case 'DIP'
%         [m_1,c_1,j_1,L_1,m_2,c_2,j_2,L_2] = auto_LumpParam(model_param);
%         c_1 = c_1(2); c_2 = c_2(2); L_1 = L_1(2); L_2 = L_2(2);
%         [m_1,c_1,j_1,L_1,m_2,c_2,j_2,L_2];
% end
%% Visualization
if nargin ~=0 % if there is an input into the function  
    isAnimOn = input.PostProc.AnimOn;
    isPlotOn = input.PostProc.PlotOn;
else 
    isAnimOn = 1; % change to 0 if you don't want to see animation
    isPlotOn = 1; % change to 0 if you don't want to see plots   
end
%% Simulation Setup
% initial condition
% q: joint ('relative, not spatial') angles q=[q_ankle, q_hip]
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

% simulation time
t_0 = 0; % initial time
%% Define Controller
if nargin ~=0
    struct_Controller = input.Controller;
    alpha = struct_Controller.alpha;
    beta = struct_Controller.beta;
    gamma = struct_Controller.gamma;
    noise_ratio = input.NoiseRatio;
else
    alpha = 1e-4;
    beta = 0.2;
    gamma = 1;
    noise_ratio = 0.2;
    struct_Controller.type = 'LQR_int'; % designate controller type used here
end
switch model_type
    case 'SIP'
        param.R = alpha*beta;
    case 'DIP'
        param.R = alpha*diag([beta, 1/beta]);
end
param.Q = gamma*eye(2*n_q);
struct_Controller.param = param;
%% Determine Controller - defined below
switch struct_Controller.type
    case 'JIC'
        FuncHandle_Controller = @JIC;
    case 'CIC'
        FuncHandle_Controller = @CIC;
    case 'LQR'
        FuncHandle_Controller = @LQR;
    case 'LQR_int'
        FuncHandle_Controller = @LQR_int;
end
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
y = cell(1, N);       % 1x Nt structures - output 
% torque = zeros(2,N); % 2xN matrix - ankle, hip torque

% run simulation
x_i = x_0; % initialize state
T = [1 0;
     1 1]; % coordinate transformatino matrix from relative to spatial
T_state = [1 0 0 0;
           1 1 0 0;
           0 0 1 0;
           0 0 1 1]; 
       
for i = 1:N
    t_i = t(i);
    % calculate Dx and output y at time t; also output ankle torque input
    % and ankle noise for graphing purposes later
    torque_i = FuncHandle_Controller(t_i, x_i, struct_Controller.param, model_param, model_type, T, T_state);
    [Dx_i, y_i] = func_Dynamics(t_i, x_i, torque_i, model_param, model_type, noise_ratio); % these functions are all below
    % log state and output
    switch coord
        case 'spatial'
%             torque(:,i) = T.'\torque_i; % spatial torque
            x(:, i) = T_state*x_i; % spatial states
        case 'relative'
%             torque(:,i) = torque_i;
            x(:, i) = x_i;
    end
    y{i} = y_i;
    
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
% outputs
    % Angular momentum, CoM, CoP, IP, ground reaction forces    
    COM_z = zeros(1,N); % vertical CoM
    MLCOP = zeros(1,N); % CoP
    MLFx = zeros(1,N); % horizontal ground reactin force
    MLFz = zeros(1,N); % vertical ground reaction force
    % renaming component of output y from simulation above
    for i = 1:N
        MLCOP(:,i) = y{i}.COP;
        COM_z(:,i) = y{i}.COM_z;
        MLFx(:,i) = y{i}.Fx;
        MLFz(:,i) = y{i}.Fz;
    end

    % Obtain IP at different frequencies (inspired by Marta's code)
    f_i = 0.55; f_int = 0.5; f_end = 5.05; % starting freq, interval between tested frequencies, frequency range
%     f_i = 0.5; f_int = 0.2; f_range = 7.4;
    f = f_i:f_int:f_end; % frequencies
    IP = zeros(1,length(f)); % IP
    theta = -MLFx./MLFz;
    theta = theta - mean(theta); % detrend
    MLCOP = MLCOP - mean(MLCOP); % detrend
    cop_ham = MLCOP'.*hamming(length(MLCOP)); % apply Hamming window to CoP data
    theta_ham = theta'.*hamming(length(theta)); % apply Hamming window to the force angle
    for n = 1:length(f) % for all frequencies designated above:
        [B,A] = butter(2,[f(n) f(n)+.2]/(1/dt/2)); % define the Butterworth filter
        COP_f = filtfilt(B,A,cop_ham); % apply filter to CoP data
        theta_f = filtfilt(B,A,theta_ham); % apply filter to force angle data
        [coeff,~,~,~,~,~] = pca([COP_f theta_f]); % PCA
        IP(n) = coeff(1,1)/coeff(2,1); % compute IP and store
    end
    IP_ratio = IP/mean(COM_z);
%     mean(COM_z);
    
%     torque_rms = [rms(torque(1,:)), rms(torque(2,:)), rms(torque(1,:))/rms(torque(2,:))];

% animation
if isAnimOn    
    b_saveanim= 0; 
    title_ = 'test';
    t_anim = t(1:100:end); 
    x_anim = x(:, 1:100:end);
    postproc_animation(t_anim, x_anim, model_param, model_type, b_saveanim, title_ ); 
end

% plot
if isPlotOn
    % Plot joint angles wrt time
    figure(1);
        switch model_type
            case 'SIP'
                plot(t, x(1,:),'r');
                ylabel('Ankle Joint Angle [rad]')
                xlabel('Time [s]')
            case 'DIP'
                plot(t, x(1,:),'r');
                hold on
                plot(t,x(2,:),'b');
                ylabel('Joint Angles [rad]')
                xlabel('Time [s]')
                legend('q1','q2')
        end
    % Plot a scatter plot of IP frequency response
    figure(2);  
        scatter(f+0.1, IP_ratio,30,'k','filled')      
        ylabel('IP/CoM')
        xlabel('Frequency [Hz]')
    
    % Example linear fit of theta_f vs. CoP relationship for a given
    % frequency band
%     f_start = 7; f_end = f_start + f_int; % which frequency band are you interested in?
%     [B,A] = butter(2,[f_start f_end]/(1/dt/2)); % define the Butterworth filter
%     COP_ex = filtfilt(B,A,cop_ham); % apply filter to CoP data
%     theta_ex = filtfilt(B,A,theta_ham); % apply filter to force angle data
%     [coeff,~,latent,~,explained,~] = pca([COP_ex theta_ex]); % PCA
%     slope_ex = coeff(2,1)/coeff(1,1); % compute IP and store
%     figure(4);
%         plot(COP_ex, theta_ex, '.', 'Color', uint8([0 0 255]))
%         ylabel('\theta_f [rad]')
%         xlabel('CoP [m]')
end
%% sub function definitions
% Controllers 
% (1) Joint stiffness %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     function u_JIC = JIC(t, x, param)
%     end
% (2) Cartesian stiffness %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     function u_CIC = CIC(t, x, param)
%     end
% (3) LQR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function u_LQR = LQR(t, x, param, model_param, model_type, T, T_state)
        Q = param.Q;
        switch coord
            case 'spatial'
                R = T*param.R*T.';
            case 'relative'
                R = param.R;
        end
%         persistent K_rika;
        persistent K_lqr; % calculate K_lqr only once for all time
%         persistent K_pp; % pole placement gain (either use K_lqr OR K_pp)
        if isempty(K_lqr) % finding K_lqr for u = -Kx
            [A_lin,B_lin,~] = getLinearDynamics(model_param,model_type);
            switch coord
                case 'spatial'
                    A_lin = T_state*A_lin/T_state;
                    B_lin = T_state*B_lin*T.';
            end
%             design lqr control
            K_lqr = lqr(A_lin,B_lin,Q,R);
            % compute stiffness and damping matrices of gain matrix
%             K_s = 0.5*(K_lqr(:,1:2) + (K_lqr(:,1:2))') % symmetric stiffness matrix
%             K_a = 0.5*(K_lqr(:,1:2) - (K_lqr(:,1:2))') % antisymmetric stiffness matrix
%             K_s_norm = norm(K_s)
%             K_a_norm = norm(K_a)
%             B_s = 0.5*(K_lqr(:,3:4) + (K_lqr(:,3:4))') % symmetric damping matrix
%             B_a = 0.5*(K_lqr(:,3:4) - (K_lqr(:,3:4))') % antisymmetric damping matrix
%             B_s_norm = norm(B_s)
%             B_a_norm = norm(B_a)
%             K_lqr = [4*K_lqr(1:2,1:2), 2*K_lqr(1:2,3:end)]; % nx stiffness/damping matrix
%             % compute open-loop poles
%             [V_op,D_op] = eig(A_lin); % linearized system open loop eigenstructure
%             D_op_hz = diag(D_op)./(2*pi); % open loop eigenvalues in Hz
%             % compute closed-loop poles
%             [V_cl,D_cl] = eig(A_lin-B_lin*K_lqr);
%             D_cl_hz = diag(D_cl)./(2*pi)
%             % plot open/closed loop poles
%             figure; 
%             plot(real(diag(D_op))./(2*pi),imag(diag(D_op))./(2*pi),'x','MarkerSize',10,'LineWidth',3,'MarkerEdgeColor','b')
%             hold on
%             plot(real(diag(D_cl))./(2*pi),imag(diag(D_cl))./(2*pi),'+','MarkerSize',10,'LineWidth',3,'MarkerEdgeColor','r')
%             legend('Open Loop Poles','Closed Loop Poles')
%             hold off
%             xlabel('Real [Hz]')
%             ylabel('Imaginary [Hz]')

%             % Pole Placement
%             zeta = 0.5; % damping ratio
%             p = [zeta*D_op(3,3)+D_op(3,3)*sqrt(zeta^2-1) zeta*D_op(3,3)-D_op(3,3)*sqrt(zeta^2-1) zeta*D_op(4,4)+D_op(4,4)*sqrt(zeta^2-1) zeta*D_op(4,4)-D_op(4,4)*sqrt(zeta^2-1)] % pole locations
%             K_pp = place(A_lin,B_lin,p); % pole placement gain
%             [V_pp,D_pp] = eig(A_lin - B_lin*K_pp); % pole placement closed loop poles
            
            % Rika's Estimated K-matrix:
%             K_rika = [1725.83125910413,328.235139259213,568.837767710976,96.4166202704406;363.924144334841,243.860833486195,149.495001245118,44.5900110519538]
        end
%         K = K_rika;
%         K = K_pp;
        K = K_lqr;
%         K
        switch coord
            case 'spatial'
                x = T_state*x;
                u_LQR = T.'*(-K*x); % back to relative torques
            case 'relative'
                u_LQR = -K*x;
        end
    end
% (4) LQR_Intermittent %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function u_LQR_int = LQR_int(t, x, param, model_param, model_type, T, T_state)
        Q = param.Q;
        switch coord
            case 'spatial'
                R = T*param.R*T.';
            case 'relative'
                R = param.R;
        end
        persistent K_lqr; % calculate K_lqr only once for all time
        if isempty(K_lqr) % finding K_lqr for u = -Kx
            [A_lin,B_lin,~] = getLinearDynamics(model_param,model_type);
            switch coord
                case 'spatial'
                    A_lin = T_state*A_lin/T_state;
                    B_lin = T_state*B_lin*T.';
            end
%             design lqr control
            K_lqr = lqr(A_lin,B_lin,Q,R);
        end
        K = K_lqr;
    %--Rika/Federico Intermittent K--%
        [~,~,Kg] = getLinearDynamics(model_param, model_type);
        if (abs(x(1)) < 1*pi/180 ) && (abs(x(2)) < 1*pi/180 )
          K(1:2,1:2) = -Kg;            
        end
        switch coord
            case 'spatial'
                x = T_state*x;
                u_LQR_int = T.'*(-K*x); % back to relative torques
            case 'relative'
                u_LQR_int = -K*x;
        end
%         K(:,3:4) = K(:,3:4)/2;
%         u_LQR_int = -K*x;
    end    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plant Dynamics
    function [Dx, output_struct] = func_Dynamics(t, x, u_ctrl, param, model_type, sigma_r)
       switch model_type
            case 'SIP'
                q = x(1,:);
                Dq = x(2,:);
                q_n = length(q);
            case 'DIP'
                q = x(1:2,:);
                Dq = x(3:4,:);
                q_n = length(q);
       end
        % eom Mq'' + Cq' + G = u_ctrl + u_noise
        [M, C, G, ~] = auto_Dynamics(q, Dq, param);
        
        % input noise, multiplier can be changed
        u_noise = 10*randn(q_n,1); % noise size can be determined by you
%         sigma_r = 0.5; % change this to change noise ratio
        u_noise(1) = sigma_r*u_noise(1);
        u = u_ctrl + u_noise;
        
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
%                 q_eq = [0; 0]*pi/180; Dq_eq = [0;0]; %2021-05-24
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
rmpath([pwd,'\dynamics_models\',pose,'\',model_type,'\',gender,'_',plane]);
end