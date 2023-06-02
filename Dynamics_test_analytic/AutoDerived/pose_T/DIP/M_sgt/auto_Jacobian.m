function [J_CoM,J_G] = auto_Jacobian(q,Dq,param)

J_CoM = zeros([2,2]) ;
  J_CoM(1,1)=-(param(2)*(114352244432*cos(q(1) + q(2)) + 368292075455*cos(q(1))))/846648300000;
  J_CoM(1,2)=-(7147015277*param(2)*cos(q(1) + q(2)))/52915518750;
  J_CoM(2,1)=-(param(2)*(114352244432*sin(q(1) + q(2)) + 368292075455*sin(q(1))))/846648300000;
  J_CoM(2,2)=-(7147015277*param(2)*sin(q(1) + q(2)))/52915518750;

 J_G = zeros([2,2]) ;
  J_G(1,1)=-(981*param(1)*param(2)*(114352244432*cos(q(1) + q(2)) + 368292075455*cos(q(1))))/87050000000000;
  J_G(1,2)=-(7011221986737*param(1)*param(2)*cos(q(1) + q(2)))/5440625000000;
  J_G(2,1)=-(7011221986737*param(1)*param(2)*cos(q(1) + q(2)))/5440625000000;
  J_G(2,2)=-(7011221986737*param(1)*param(2)*cos(q(1) + q(2)))/5440625000000;

 