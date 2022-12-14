%%%%%%% q_a = [q_A; q_H]; %%%%%%% active joints
%%%%%%% q_p = [q_S; q_W; q_C]; %%%%%%%% passive joints
%%%%%%% q_g = [q_A; q_H]; %%%%%%% independent joints
%%%%%%% q_s = [q_S; q_W; q_C]; %%%%%%%% dependent joints
%%%%%%% q_j = [q_g; q_s]; %%%%%%%% all joints
%% q_s = q_s(q_g);
e1 = p_O_Cr - p_O_Sr;
e2 = -p_O_Sr;
e3 = -p_O_H;


kappa1 = acos( (L_arm^2+norm(e1)^2-L_C^2)/(2*L_arm*norm(e1)) ); kappa1 = simplify(kappa1);
kappa2 = atan2(e2(2),e2(1)) - atan2(e1(2), e1(1)); kappa2 = simplify(kappa2);
kappa3 = acos( (norm(e2)^2+L_sh^2-norm(e3)^2)/(2*norm(e2)*L_sh) ); kappa3 = simplify(kappa3);

q_Sr = pi/2 - (kappa1+kappa2+kappa3);

zeta = acos((L_C^2+L_arm^2-norm(e1)^2)/(2*L_C*L_arm));

q_Wr = pi/2 - zeta;
% q_Cr = q_A+q_H+q_Sr+q_Wr;

% left side
e1 = p_O_Cl - p_O_Sl;
e2 = -p_O_Sl;
e3 = -p_O_H;

kappa1 = acos((L_arm^2+norm(e1)^2-L_C^2)/(2*L_arm*norm(e1)) ); kappa1 = simplify(kappa1);
kappa2 = atan2(e1(2),e1(1)) - atan2(e2(2), e2(1)); kappa2 = simplify(kappa2);
kappa3 = acos( (norm(e2)^2+L_sh^2-norm(e3)^2)/(2*norm(e2)*L_sh) ); kappa3 = simplify(kappa3);
% q_Sl = (kappa1+kappa2+kappa3)-pi/2;
q_Sl = - pi/2 + (kappa1+kappa2-kappa3);

zeta = acos((L_C^2+L_arm^2-norm(e1)^2)/(2*L_C*L_arm));
q_Wl = zeta - pi/2;
% q_Cl = q_A+q_H+q_Sl+q_Wl;
%% Express Dq_v as functions of Dq_u
% Dq_v = PHI*Dq_u  shall be done numerically

%% Write
write_fcn('auto_GetDependentJoint.m',{'q_g'},q_g_list, {(eval(q_s)), 'q_s'});