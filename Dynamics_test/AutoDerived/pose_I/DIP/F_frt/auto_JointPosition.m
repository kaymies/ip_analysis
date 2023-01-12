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
  p_O_RS(1,1)=-(param(2)*(19916*sin(q(1) + q(2)) - 7981*cos(q(1) + q(2)) + 32284*sin(q(1))))/69400;
  p_O_RS(2,1)=(param(2)*(19916*cos(q(1) + q(2)) + 7981*sin(q(1) + q(2)) + 32284*cos(q(1))))/69400;

 p_O_RE = zeros([2,1]) ;
  p_O_RE(1,1)=(80936955523371977160460348597670087*param(2)*cos(q(1) +...
          q(2)))/703799613246712963710969619624755200 - (1114*param(2)*sin(q(1) + q(2)))/8675 - (8071*param(2)*sin(q(1)))/17350;
  p_O_RE(2,1)=(1114*param(2)*cos(q(1) + q(2)))/8675 + (80936955523371977160460348597670087*param(2)*...
         sin(q(1) + q(2)))/703799613246712963710969619624755200 + (8071*param(2)*cos(q(1)))/17350;

 p_O_RW = zeros([2,1]) ;
  p_O_RW(1,1)=(40468477761685982015338505870408257*param(2)*cos(q(1) +...
          q(2)))/351899806623356481855484809812377600 + (83*param(2)*sin(q(1) + q(2)))/3470 - (8071*param(2)*sin(q(1)))/17350;
  p_O_RW(2,1)=(40468477761685982015338505870408257*param(2)*sin(q(1) +...
          q(2)))/351899806623356481855484809812377600 - (83*param(2)*cos(q(1) + q(2)))/3470 + (8071*param(2)*cos(q(1)))/17350;

 p_O_RH = zeros([2,1]) ;
  p_O_RH(1,1)=(40468477761685980077913041862018967*param(2)*cos(q(1) +...
          q(2)))/351899806623356481855484809812377600 + (239*param(2)*sin(q(1) + q(2)))/3470 - (8071*param(2)*sin(q(1)))/17350;
  p_O_RH(2,1)=(40468477761685980077913041862018967*param(2)*sin(q(1) +...
          q(2)))/351899806623356481855484809812377600 - (239*param(2)*cos(q(1) + q(2)))/3470 + (8071*param(2)*cos(q(1)))/17350;

 p_O_LS = zeros([2,1]) ;
  p_O_LS(1,1)=-(param(2)*(7981*cos(q(1) + q(2)) + 19916*sin(q(1) + q(2)) + 32284*sin(q(1))))/69400;
  p_O_LS(2,1)=(param(2)*(19916*cos(q(1) + q(2)) - 7981*sin(q(1) + q(2)) + 32284*cos(q(1))))/69400;

 p_O_LE = zeros([2,1]) ;
  p_O_LE(1,1)=- (80936955523371977160460348597670087*param(2)*cos(q(1) +...
          q(2)))/703799613246712963710969619624755200 - (1114*param(2)*sin(q(1) + q(2)))/8675 - (8071*param(2)*sin(q(1)))/17350;
  p_O_LE(2,1)=(1114*param(2)*cos(q(1) + q(2)))/8675 - (80936955523371977160460348597670087*param(2)*...
         sin(q(1) + q(2)))/703799613246712963710969619624755200 + (8071*param(2)*cos(q(1)))/17350;

 p_O_LW = zeros([2,1]) ;
  p_O_LW(1,1)=(83*param(2)*sin(q(1) + q(2)))/3470 - (40468477761685982015338505870408257*param(2)*...
         cos(q(1) + q(2)))/351899806623356481855484809812377600 - (8071*param(2)*sin(q(1)))/17350;
  p_O_LW(2,1)=(8071*param(2)*cos(q(1)))/17350 - (40468477761685982015338505870408257*param(2)*...
         sin(q(1) + q(2)))/351899806623356481855484809812377600 - (83*param(2)*cos(q(1) + q(2)))/3470;

 p_O_LH = zeros([2,1]) ;
  p_O_LH(1,1)=(239*param(2)*sin(q(1) + q(2)))/3470 - (40468477761685980077913041862018967*param(2)*...
         cos(q(1) + q(2)))/351899806623356481855484809812377600 - (8071*param(2)*sin(q(1)))/17350;
  p_O_LH(2,1)=(8071*param(2)*cos(q(1)))/17350 - (40468477761685980077913041862018967*param(2)*...
         sin(q(1) + q(2)))/351899806623356481855484809812377600 - (239*param(2)*cos(q(1) + q(2)))/3470;

 