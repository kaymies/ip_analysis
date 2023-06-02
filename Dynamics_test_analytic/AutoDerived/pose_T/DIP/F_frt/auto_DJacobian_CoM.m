function [DJ_CoM] = auto_DJacobian_CoM(q,Dq,param)

DJ_CoM = zeros([2,2]) ;
  DJ_CoM(1,1)=(param(2)*(2291718100814621124653788003875045731221609*Dq(1)*sin(q(1) + q(2)) +...
          2291718100814621124653788003875045731221609*Dq(2)*sin(q(1) + q(2)) + 7028249084006354683098952670509263673622528*Dq(1)*sin(q(1))))/17139280081590577448771387661911851008000000;
  DJ_CoM(1,2)=(param(2)*(2291718100814621124653788003875045731221609*Dq(1)*sin(q(1) + q(2)) +...
          2291718100814621124653788003875045731221609*Dq(2)*sin(q(1) + q(2))))/17139280081590577448771387661911851008000000;
  DJ_CoM(2,1)=-(param(2)*(2291718100814621124653788003875045731221609*Dq(1)*cos(q(1) + q(2)) +...
          2291718100814621124653788003875045731221609*Dq(2)*cos(q(1) + q(2)) + 7028249084006354683098952670509263673622528*Dq(1)*cos(q(1))))/17139280081590577448771387661911851008000000;
  DJ_CoM(2,2)=-(param(2)*(2291718100814621124653788003875045731221609*Dq(1)*cos(q(1) + q(2)) +...
          2291718100814621124653788003875045731221609*Dq(2)*cos(q(1) + q(2))))/17139280081590577448771387661911851008000000;

 