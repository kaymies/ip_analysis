% Function to find the kinetic energy of a link
% Inputs:
%   m: Mass of the link
%   x: Position of the link - used to compute velocity
%   q: Generalized coordinates
%   dq: Time-derivative of the generalized coordinates
% Outputs:
%   K: Kinetic Energy of the Link
function K = KE_trans(m, x, q, dq)

v = jacobian(x,q)*dq ; % velocity
K = 1/2 * m * (v' * v) ; % kinetic energy