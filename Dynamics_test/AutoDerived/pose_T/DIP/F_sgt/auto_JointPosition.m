function [p_O_CF,p_O_CH,p_O_CK,p_O_NE,p_O_HE,p_O_RS,p_O_RE,p_O_RW,p_O_RH,p_O_LS,p_O_LE,p_O_LW,p_O_LH] = auto_JointPosition(q,param)

p_O_CF = zeros([2,1]) ;
  p_O_CF(1,1)=0;
  p_O_CF(2,1)=0;

 p_O_CH = zeros([2,1]) ;
  p_O_CH(1,1)=-(8071*param(2)*sin(q(1)))/17350;
  p_O_CH(2,1)=(8071*param(2)*cos(q(1)))/17350;

 p_O_CK = zeros([2,1]) ;
  p_O_CK(1,1)=-(2193*param(2)*sin(q(1)))/8675;
  p_O_CK(2,1)=(2193*param(2)*cos(q(1)))/8675;

 p_O_NE = zeros([2,1]) ;
  p_O_NE(1,1)=-(param(2)*(6148*sin(q(1) + q(2)) + 8071*sin(q(1))))/17350;
  p_O_NE(2,1)=(param(2)*(6148*cos(q(1) + q(2)) + 8071*cos(q(1))))/17350;

 p_O_HE = zeros([2,1]) ;
  p_O_HE(1,1)=-(param(2)*(8585*sin(q(1) + q(2)) + 8071*sin(q(1))))/17350;
  p_O_HE(2,1)=(param(2)*(8585*cos(q(1) + q(2)) + 8071*cos(q(1))))/17350;

 p_O_RS = zeros([2,1]) ;
  p_O_RS(1,1)=-(param(2)*(4979*sin(q(1) + q(2)) + 8071*sin(q(1))))/17350;
  p_O_RS(2,1)=(param(2)*(4979*cos(q(1) + q(2)) + 8071*cos(q(1))))/17350;

 p_O_RE = zeros([2,1]) ;
  p_O_RE(1,1)=-(param(2)*(4979*sin(q(1) + q(2)) + 8071*sin(q(1))))/17350;
  p_O_RE(2,1)=(param(2)*(4979*cos(q(1) + q(2)) + 8071*cos(q(1))))/17350;

 p_O_RW = zeros([2,1]) ;
  p_O_RW(1,1)=-(param(2)*(4979*sin(q(1) + q(2)) + 8071*sin(q(1))))/17350;
  p_O_RW(2,1)=(param(2)*(4979*cos(q(1) + q(2)) + 8071*cos(q(1))))/17350;

 p_O_RH = zeros([2,1]) ;
  p_O_RH(1,1)=-(param(2)*(4979*sin(q(1) + q(2)) + 8071*sin(q(1))))/17350;
  p_O_RH(2,1)=(param(2)*(4979*cos(q(1) + q(2)) + 8071*cos(q(1))))/17350;

 p_O_LS = zeros([2,1]) ;
  p_O_LS(1,1)=-(param(2)*(4979*sin(q(1) + q(2)) + 8071*sin(q(1))))/17350;
  p_O_LS(2,1)=(param(2)*(4979*cos(q(1) + q(2)) + 8071*cos(q(1))))/17350;

 p_O_LE = zeros([2,1]) ;
  p_O_LE(1,1)=-(param(2)*(4979*sin(q(1) + q(2)) + 8071*sin(q(1))))/17350;
  p_O_LE(2,1)=(param(2)*(4979*cos(q(1) + q(2)) + 8071*cos(q(1))))/17350;

 p_O_LW = zeros([2,1]) ;
  p_O_LW(1,1)=-(param(2)*(4979*sin(q(1) + q(2)) + 8071*sin(q(1))))/17350;
  p_O_LW(2,1)=(param(2)*(4979*cos(q(1) + q(2)) + 8071*cos(q(1))))/17350;

 p_O_LH = zeros([2,1]) ;
  p_O_LH(1,1)=-(param(2)*(4979*sin(q(1) + q(2)) + 8071*sin(q(1))))/17350;
  p_O_LH(2,1)=(param(2)*(4979*cos(q(1) + q(2)) + 8071*cos(q(1))))/17350;

 