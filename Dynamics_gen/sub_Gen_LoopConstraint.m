%% Loop Constraint Equations
% initialize
c_L1 = sym(zeros(2,1));
c_L2 = sym(zeros(2,1));
% cLoop(q) = 0;
c_L1 = p_O_LC - dvector.p_O_LC;
c_L2 = p_O_RC - dvector.p_O_RC;

% stack
c_L = [c_L1;
    c_L2];
c_L = simplify(c_L, 'Steps', 3);

nc = size(c_L,1);
%% Constraint Jacobian
J_L = jacobian(c_L, q);
J_L = simplify(J_L, 'Steps', 3);

% decompose constraint jacobian
J_Lr = J_L*jacobian(q, r);   J_Lr = simplify(J_Lr, 'Steps', 3);
J_Ls = J_L*jacobian(q, s);   J_Ls = simplify(J_Ls, 'Steps', 3);
J_sr = -J_Ls\J_Lr; J_sr = simplify(J_sr, 'Steps', 3);
%% induced jacobians
Lambda = [eye(nr); -J_Ls\J_Lr];     Lambda = simplify(Lambda, 'Steps', 3);
%% Time-derivative of constraint jacobians to compute acceleration
Dc_L = jacobian(c_L, q)*Dq;
DJ_L = jacobian(Dc_L, q);              DJ_L = simplify(DJ_L, 'Steps', 3);

Omega = [zeros(nr, nq); -J_Ls\DJ_L];
Omega = simplify(Omega, 'Steps', 3);

%% to use write_fcn, symbolize numeric matrices
if isnumeric(J_L)
    J_L = sym(J_L);
    J_L = simplify(J_L, 'Steps', 3);
end

if isnumeric(DJ_L)
    DJ_L = sym(DJ_L);
    DJ_L = simplify(DJ_L, 'Steps',3);
end

if isnumeric(J_Lr)
    J_Lr = sym(J_Lr);
    J_Lr = simplify(J_Lr, 'Steps', 3);
end
if isnumeric(J_Ls)
    J_Ls = sym(J_Ls);
    J_Ls = simplify(J_Ls, 'Steps', 3);
end
if isnumeric(Lambda)
    Lambda = sym(Lambda);
    Lambda = simplify(Lambda, 'Steps', 3);
end  
if isnumeric(Omega)
    Omega = sym(Omega);
    Omega = simplify(Omega, 'Steps', 3);
end
