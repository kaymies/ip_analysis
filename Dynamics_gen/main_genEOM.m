% Generic balance model
% In sub_Gen_Coordinates and sub_Gen_ForwardKin_eqns
% One should determine the degree of freedom of the model
% And in sub_Gen_MechEnergy

clc; close all;
clear; 

%% Parameters
gender = 'M';
plane = 'sgt';
syms totalHeight totalMass positive
% parameters
param_list = {'totalMass', 'param(1)'
            'totalHeight', 'param(2)'};       
        
% totalHeight = 1.73
% totalMass = 70
sub_Gen_Parameters;
sub_Gen_Variables;
%% Coordinates/States
sub_Gen_Coordinates; % fully_actuated / underactuated <- The most important part
%% Kinematics
% output: COM position, joint position for visualization purpose. 
sub_Gen_ForwardKin_eqns;


%% Geometric Constraints
% sub_Gen_LoopConstraint; % closed kinematic chain constraint
%% Inverse kinematics for dependent joints
% Express dependent joints as functions of independent joints
% sub_Gen_DependentJointValues;
%% Kinetic Energy and Potential Energy to construct Lagrangian
sub_Gen_MechEnergy
%% M*ddq+C*dq+G = Torq
% Solve for Lagrangian Mechanics
[M, C, G, B] = LagrangianDynamics(KE_tot, PE_tot, q, Dq, q_act); 

%%
[m_tot, p_O_CoM] = sub_Get_Com([2*Shank.m, 2*Thigh.m, Trunk.m, Head.m, Upperarm.m, Forearm.m, Hand.m, Upperarm.m, Forearm.m, Hand.m],...
    [SE2_O_Shankc.p, SE2_O_Thighc.p, SE2_O_Trunkc.p, SE2_O_Headc.p, SE2_O_LUarmc.p, SE2_O_LFarmc.p, SE2_O_LHandc.p, SE2_O_RUarmc.p, SE2_O_RFarmc.p, SE2_O_RHandc.p]);

J_CoM = jacobian(p_O_CoM,q);    J_CoM = simplify(J_CoM, 'Steps', 3);
v_O_CoM = simplify(jacobian(p_O_CoM,q)*Dq,3);
DJ_CoM = jacobian(v_O_CoM, q);
%%
%% Write Functions
% Joint Positions
write_fcn('auto_JointPosition.m', {'q', 'param'}, ...
    [q_list; param_list],...
    {sym(SE2_O_CF.p),         'p_O_CF'    
     SE2_O_CH.p,         'p_O_CH'
     SE2_O_CK.p,         'p_O_CK'     
     SE2_O_NE.p,         'p_O_NE'
     SE2_O_HE.p,         'p_O_HE'
     SE2_O_RS.p,         'p_O_RS'
     SE2_O_RE.p,         'p_O_RE'
     SE2_O_RW.p,         'p_O_RW'     
     SE2_O_RH.p,         'p_O_RH'
     SE2_O_LS.p,         'p_O_LS'
     SE2_O_LE.p,         'p_O_LE'
     SE2_O_LW.p,         'p_O_LW'
     SE2_O_LH.p,         'p_O_LH'     
     });
 
 
% write_fcn('auto_COMPosition.m', {'q', 'param'}, ...
%     [q_list; param_list],...
%     {eval(p_O_lbc),       'p_O_lbc'
%      eval(p_O_ubc),       'p_O_ubc'
%      eval(p_O_armrc),       'p_O_armrc'
%      eval(p_O_armlc),       'p_O_armlc'
%      });

write_fcn('auto_COMinformation.m',...
    {'q','Dq','param'},[q_list; Dq_list; param_list],...
    {p_O_CoM,  'p_O_CoM'
     v_O_CoM,  'v_O_CoM'     
     sym(m_tot), 'm_tot'});

write_fcn('auto_MechEnergy.m',...
    {'q','Dq','param'},[q_list; Dq_list; param_list],...
    {KE_tot,'KE_tot'; PE_tot,'PE_tot'});


write_fcn('auto_Dynamics.m',...
    {'q','Dq','param'},[q_list; Dq_list; param_list],...
    {M,'M'; C,'C'; G,'G'; B,'B'});


J_G = jacobian(G,q);
write_fcn('auto_Jacobian.m', {'q','Dq','param'},[q_list; Dq_list; param_list],...
    {J_CoM, 'J_CoM'
    J_G, 'J_G'});


write_fcn('auto_DJacobian_CoM.m', {'q','Dq','param'},[q_list; Dq_list; param_list],...
    {DJ_CoM, 'DJ_CoM'});


a_O_CoM = J_CoM*DDq + DJ_CoM*Dq;
Fz = m_tot*(a_O_CoM(2)+g);
Fx = m_tot*(a_O_CoM(1));
write_fcn('auto_GroundReactionForce.m', ...    
    {'q','Dq','DDq','param'},[q_list; Dq_list; DDq_list; param_list],...
    {Fx,'Fx'; Fz,'Fz'});
%% Lumped Parameter Values
q_A = 0; q_H = 0; q_RE = 0; q_RS = 0; q_RW = 0; q_LE = 0; q_LS = 0; q_LW =0;
[m_1, c_1, j_1] = sub_Get_Com([Thigh.m, Shank.m], eval([SE2_O_Thighc.p, SE2_O_Shankc.p]), [Thigh.j, Shank.j]);
m_1 = 2*m_1; j_1 = 2*j_1; % two legs
% eval([m_1, c_1, j_1]);
L_1 = eval(SE2_O_CH.p);
[m_2, c_2, j_2] = sub_Get_Com(...
    [Trunk.m, Head.m, ...
    Upperarm.m, Forearm.m, Hand.m, ...
    Upperarm.m, Forearm.m, Hand.m],...
    eval([SE2_O_Trunkc.p, SE2_O_Headc.p, ...
    SE2_O_LUarmc.p, SE2_O_LFarmc.p, SE2_O_LHandc.p, ...
    SE2_O_RUarmc.p, SE2_O_RFarmc.p, SE2_O_RHandc.p]-SE2_O_CH.p), ...
    [Trunk.j, Head.j, ...
    Upperarm.j, Forearm.j, Hand.j, ...
    Upperarm.j, Forearm.j, Hand.j]);
% eval([m_2, c_2, j_2]);
L_2 = eval(SE2_O_HE.p-SE2_O_CH.p);


write_fcn('auto_LumpParam.m',...
    {'param'},param_list,...
    {m_1,'m_1'; c_1(2),'c_1'; j_1,'j_1'; L_1(2),'L_1'; m_2,'m_2'; c_2(2),'c_2'; j_2,'j_2'; L_2(2),'L_2'});

% single-inverted pendulum
% [m_sip, c_sip, j_sip] = sub_Get_Com([m_1, m_2], eval([c_1, c_2]), [j_1, j_2]);
% L_sip = eval(SE2_O_HE.p);
% write_fcn('auto_LumpParam.m',...
%     {'param'},param_list,...
%     {m_sip,'m_sip'; c_sip,'c_sip'; j_sip,'j_sip'; L_sip,'L_sip'});

% [m_1,c_1,j_1,L_1,m_2,c_2,j_2,L_2] = auto_LumpParam([70;1.70])