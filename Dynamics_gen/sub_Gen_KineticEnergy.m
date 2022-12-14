%% Compute Kinetic Energy
% KE = KE_translation + KE_rotation;
KE_lb = KE_trans(m_lb, p_O_Ac, q, Dq) + 1/2*(j_lb)*Dq_A^2;
KE_ub = KE_trans(m_ub, p_O_Hc, q, Dq) + 1/2*(j_ub)*(Dq_A+Dq_H)^2;
KE_arm = KE_trans(m_arm, p_O_Sc, q, Dq) + 1/2*(j_arm)*(Dq_A+Dq_H+Dq_S)^2;
KE_cane = KE_trans(m_C, p_O_Cc, q, Dq) + 1/2*(j_arm)*(Dq_A+Dq_H+Dq_S)^2;

%% Compute Total Kinetic Energy
KE_tot = KE_lb+KE_ub+KE_arm+KE_cane;
KE_tot = simplify(KE_tot, 'Steps', 3);  % iterate 3 times
    
    

