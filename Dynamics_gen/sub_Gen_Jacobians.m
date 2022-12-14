% function sub_Gen_Jacobians
%% Get Forward Jacobian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% com 
Jac_CoM = jacobian(p_O_CoM, r)+jacobian(p_O_CoM, s)*H;
dJac_CoM = jacobian(Jac_CoM*Dr,r)+jacobian(Jac_CoM*Dr,s)*H;
% simplify(Jac_CoM, 'Steps', 2);
% simplify(dJac_CoM, 'Steps', 2);

% EE
Jac_HE = jacobian(p_O_HE, r) + jacobian(p_O_HE, s)*H;
dJac_HE = jacobian(Jac_HE*Dr,r) + jacobian(Jac_HE*Dr,s)*H;

% Cane
Jac_CaneR = jacobian(p_O_RC, q); 
Jac_CaneL = jacobian(p_O_LC, q); 

Jac_CaneRr = jacobian(p_O_RC, r) + jacobian(p_O_RC, s)*H; 
Jac_CaneLr = jacobian(p_O_LC, r) + jacobian(p_O_LC, s)*H; 

%%
Gr = Lambda'*G;
Jac_Gr = jacobian(Gr,r)+jacobian(Gr,s)*H;
% Jac_Gr = simplify(Jac_Gr,'Steps',1);

Jac_G = jacobian(G,q); 
% Jac_G = simplify(Jac_G, 'Steps', 1);
%% 
