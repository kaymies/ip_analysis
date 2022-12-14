%% Compute Kinetic Energy
% KE = KE_translation + KE_rotation;
KE_Shank = KE_trans(Shank.m, SE2_O_Shankc.p, q, Dq) +  1/2*Shank.j*Dq_A^2;
KE_Thigh = KE_trans(Thigh.m, SE2_O_Thighc.p, q, Dq) +  1/2*Thigh.j*(Dq_A+Dq_K)^2;
KE_Trunk = KE_trans(Trunk.m, SE2_O_Trunkc.p, q, Dq) +  1/2*Trunk.j*(Dq_A+Dq_K+Dq_H)^2;
KE_Head  = KE_trans(Head.m, SE2_O_Headc.p, q, Dq) +  1/2*Head.j*(Dq_A+Dq_K+Dq_H+Dq_NE)^2;
KE_LUarm = KE_trans(Upperarm.m, SE2_O_LUarmc.p, q, Dq) +  1/2*Upperarm.j*(Dq_A+Dq_K+Dq_H+Dq_LS)^2;
KE_LFarm = KE_trans(Forearm.m, SE2_O_LFarmc.p, q, Dq) +  1/2*Forearm.j*(Dq_A+Dq_K+Dq_H+Dq_LS+Dq_LE)^2;
KE_LHand = KE_trans(Forearm.m, SE2_O_LHandc.p, q, Dq) +  1/2*Forearm.j*(Dq_A+Dq_K+Dq_H+Dq_LS+Dq_LE+Dq_LW)^2;

KE_RUarm = KE_trans(Upperarm.m, SE2_O_RUarmc.p, q, Dq) +  1/2*Upperarm.j*(Dq_A+Dq_K+Dq_H+Dq_RS)^2;
KE_RFarm = KE_trans(Forearm.m, SE2_O_RFarmc.p, q, Dq) +  1/2*Forearm.j*(Dq_A+Dq_K+Dq_H+Dq_RS+Dq_RE)^2;
KE_RHand = KE_trans(Forearm.m, SE2_O_RHandc.p, q, Dq) +  1/2*Forearm.j*(Dq_A+Dq_K+Dq_H+Dq_RS+Dq_RE+Dq_RW)^2;

%% Compute Total Kinetic Energy
KE_tot = KE_Shank*2+KE_Thigh*2+KE_Trunk+KE_Head+KE_LUarm+KE_LFarm+KE_LHand + KE_RUarm+KE_RFarm+KE_RHand;
KE_tot = simplify(KE_tot, 'Steps', 3);  % iterate 3 times
%% Compute Potential Energy
PE_Shank = g*Shank.m*SE2_O_Shankc.p(2);
PE_Thigh = g*Thigh.m*SE2_O_Thighc.p(2);
PE_Trunk = g*Trunk.m*SE2_O_Trunkc.p(2);
PE_Head = g*Head.m*SE2_O_Headc.p(2);

PE_LUarm = g*Upperarm.m*SE2_O_LUarmc.p(2);
PE_LFarm = g*Forearm.m*SE2_O_LFarmc.p(2);
PE_LHand = g*Hand.m*SE2_O_LHandc.p(2);

PE_RUarm = g*Upperarm.m*SE2_O_RUarmc.p(2);
PE_RFarm = g*Forearm.m*SE2_O_RFarmc.p(2);
PE_RHand = g*Hand.m*SE2_O_RHandc.p(2);

PE_tot = PE_Shank*2+PE_Thigh*2+PE_Trunk+PE_Head+PE_LUarm+PE_LFarm+PE_LHand + PE_RUarm+PE_RFarm+PE_RHand;
PE_tot = simplify(PE_tot, 'Steps', 3);     

%%





