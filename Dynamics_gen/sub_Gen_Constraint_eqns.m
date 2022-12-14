% function sub_Gen_Constraint_eqns
%% Get Forward Jacobian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% param.x_c = 2.5; % constraint
syms x_c real
param_list = {'x_c','param.x_c'};


Con = p_O_EE(1) - x_c;
J_Con = jacobian(Con, q);
dJ_Con = jacobian(J_Con*Dq, q);

simplify(J_Con, 'Steps', 10);
simplify(dJ_Con, 'Steps', 10);

