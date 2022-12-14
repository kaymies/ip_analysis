% Function to output the dynamics matrices. Uses the lagrangian method.
% Inputs:
%   K: Kinetic Energy scalar
%   U: Potential Energy scalar
%   q: Generalized coordinates
%   dq: Time-derivative of the generalized coordinates
%   q_rel: Relative angles of the system (N-1 dimension)
% Outputs:
%   D: D(q) Inertia matrix
%   C: C(q,dq) Coriollis matrix
%   G: G(q) Gravity matrix
%   B: B(q) Input Matrix?
function [D, C, G, B] = LagrangianDynamics(K, U, q, dq, q_act)

D = simplify( jacobian(jacobian(K,dq), dq) ) ;
C = sym(zeros(length(q), length(q)));
for k=1:length(q)
    for j=1:length(q)
%         C(k,j) = sym(0) ;
        for i=1:length(q)
            C(k,j) = C(k,j) + 1/2 * ( diff(D(k,j),q(i)) + diff(D(k,i),q(j)) - diff(D(i,j),q(k)) ) * dq(i);
        end
    end
end
C = simplify(C, 'Steps', 10) ;
G = simplify( jacobian(U,q) , 'Steps', 10)' ;
B = jacobian(q_act, q)' ;