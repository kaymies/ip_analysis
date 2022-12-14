function [DJ_CoM] = auto_DJacobian_CoM(q,Dq,param)

DJ_CoM = zeros([2,1]) ;
  DJ_CoM(1,1)=(467907820151*Dq(1)*param(2)*sin(q(1)))/846648300000;
  DJ_CoM(2,1)=-(467907820151*Dq(1)*param(2)*cos(q(1)))/846648300000;

 