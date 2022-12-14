%% Total mass, Total com
[m_tot p_O_CoM]  = sub_Get_Com([m_1 m_2], [p_O_J1c p_O_J2c]);
p_O_CoM = simplify(p_O_CoM,'Steps',3);

%% Compute Total Potential Energy
PE_tot = m_tot*g*p_O_CoM(2);
PE_tot = simplify(PE_tot, 'Steps', 3); % iterate  3 times

