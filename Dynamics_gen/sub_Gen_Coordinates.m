% Joints
% syms q_A q_K q_H real
% syms q_RS q_LS q_RE q_LE q_RW q_LW real
%% Joint Labelling  
% here you determine the degree-of-freedom of the model

% eg double inverted pendulum model with two arms
% q = [q_A; q_H; q_RS; q_LS];
% Dq = [Dq_A; Dq_H; Dq_RS; Dq_LS];
% nq = size(q,1);
% q_K = 0; q_RE = 0; q_LE = 0; q_RW = 0; q_LW = 0;
% Dq_K = 0; Dq_RE = 0; Dq_LE = 0; Dq_RW = 0; Dq_LW = 0;
% q_NE = 0; Dq_NE = 0; % neck

% double inverted pendulum model
% home: superman
% relative angles
q = [q_A; q_H];
Dq = [Dq_A; Dq_H];
nq = size(q,1);
q_K = 0; q_RS = pi; q_LS = -pi; q_RE = 0; q_LE = 0; q_RW = 0; q_LW = 0;
Dq_K = 0; Dq_RS = 0; Dq_LS = 0; Dq_RE = 0; Dq_LE = 0; Dq_RW = 0; Dq_LW = 0;
q_NE = 0; Dq_NE = 0; % neck

% single inverted pendulum model (Rika 2021-01-05)
% q = [q_A];
% Dq = [Dq_A];
% q_H = 0; dq_H = 0;
% nq = size(q,1);

% T-pose (arms out to the side)
% %q_RS = pi/2; q_LS = -pi/2; % Rika 2020-12-10: Marta's data has arms in abducted position
q_RS = -pi/2; q_LS = pi/2; % Rika 2020-12-13: Corrected angles
%% Actuated joint : assuming torque-actuated model
q_act = q;
%% Gen.Positions, Gen.Velocities, Gen.Accelerations
% double inverted pendulum model
q_list = {'q_A', 'q(1)'
            'q_H', 'q(2)'};
Dq_list = {'Dq_A', 'Dq(1)'
            'Dq_H', 'Dq(2)'};  
DDq = [DDq_A; DDq_H];
% DDq = [DDq_A; DDq_A + DDq_H];  

DDq_list = {'DDq_A', 'DDq(1)'
            'DDq_H', 'DDq(2)'};               

% single inverted pendulum model (Rika 2021-01-05)
% q_list = {'q_A', 'q(1)'};
% Dq_list = {'Dq_A', 'Dq(1)'};        
% DDq = [DDq_A];        
% DDq_list = {'DDq_A', 'DDq(1)'};         
        