function [Fx,Fz] = auto_GroundReactionForce(q,Dq,DDq,param)

Fx = zeros([1,1]) ;
  Fx(1,1)=(4863*param(1)*((467907820151*Dq(1)^2*param(2)*sin(q(1)))/846648300000 - (467907820151*...
         DDq(1)*param(2)*cos(q(1)))/846648300000))/5000;

 Fz = zeros([1,1]) ;
  Fz(1,1)=-(4863*param(1)*((467907820151*Dq(1)^2*param(2)*cos(q(1)))/846648300000 + (467907820151*...
         DDq(1)*param(2)*sin(q(1)))/846648300000 - 981/100))/5000;

 