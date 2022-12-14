function [J_CoM,J_G] = auto_Jacobian(q,Dq,param)

J_CoM = zeros([2,1]) ;
  J_CoM(1,1)=-(467907820151*param(2)*cos(q(1)))/846648300000;
  J_CoM(2,1)=-(467907820151*param(2)*sin(q(1)))/846648300000;

 J_G = zeros([1,1]) ;
  J_G(1,1)=-(459017571568131*param(1)*param(2)*cos(q(1)))/87050000000000;

 