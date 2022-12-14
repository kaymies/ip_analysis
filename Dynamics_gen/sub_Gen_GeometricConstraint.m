%%%%%%% q_a = [q_A; q_H]; %%%%%%% active joints
%%%%%%% q_p = [q_S; q_W; q_C]; %%%%%%%% passive joints
%%%%%%% q_g = [q_A; q_H]; %%%%%%% independent joints
%%%%%%% q_s = [q_S; q_W; q_C]; %%%%%%%% dependent joints
%%%%%%% q_j = [q_g; q_s]; %%%%%%%% all joints
%% Geometrical Constraint Equations
% g(q) = 0;
g(1:2) = T_O_W1(1:2,3) - T_O_W2(1:2,3); % position
g(3) = 