function [J_CoM,J_G] = auto_Jacobian(q,Dq,param)

J_CoM = zeros([2,2]) ;
  J_CoM(1,1)=- (2291718100814621124653788003875045731221609*param(2)*cos(q(1) +...
          q(2)))/17139280081590577448771387661911851008000000 - (21657464588*param(2)*cos(q(1)))/52814484375;
  J_CoM(1,2)=-(2291718100814621124653788003875045731221609*param(2)*cos(q(1) + q(2)))/17139280081590577448771387661911851008000000;
  J_CoM(2,1)=- (2291718100814621124653788003875045731221609*param(2)*sin(q(1) +...
          q(2)))/17139280081590577448771387661911851008000000 - (21657464588*param(2)*sin(q(1)))/52814484375;
  J_CoM(2,2)=-(2291718100814621124653788003875045731221609*param(2)*sin(q(1) + q(2)))/17139280081590577448771387661911851008000000;

 J_G = zeros([2,2]) ;
  J_G(1,1)=- (5311493190207*param(1)*param(2)*cos(q(1)))/1355468750000 -...
          (2248175456899143323285366031801419862328398429*param(1)*param(2)*cos(q(1) + q(2)))/1759499033116782409277424049061888000000000000;
  J_G(1,2)=-(2248175456899143323285366031801419862328398429*param(1)*param(2)*cos(q(1) + q(2)))/1759499033116782409277424049061888000000000000;
  J_G(2,1)=-(2248175456899143323285366031801419862328398429*param(1)*param(2)*cos(q(1) + q(2)))/1759499033116782409277424049061888000000000000;
  J_G(2,2)=-(2248175456899143323285366031801419862328398429*param(1)*param(2)*cos(q(1) + q(2)))/1759499033116782409277424049061888000000000000;

 