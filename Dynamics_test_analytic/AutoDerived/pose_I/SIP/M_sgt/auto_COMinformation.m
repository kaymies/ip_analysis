function [p_O_CoM,v_O_CoM,m_tot] = auto_COMinformation(q,Dq,param)

p_O_CoM = zeros([2,1]) ;
  p_O_CoM(1,1)=-(5000*((61*param(1)*((1536770845809054406339*param(2)*...
         cos(q(1)))/35311675119957558208092164488953856000 + (379651*param(2)*sin(q(1)))/870500))/10000 - (61*param(1)*((1536770845809054406339*param(2)*...
         cos(q(1)))/35311675119957558208092164488953856000 - (379651*param(2)*sin(q(1)))/870500))/10000 - (271*param(1)*((20193591868813040730741*param(2)*...
         cos(q(1)))/1765583755997877910404608224447692800000 - (30385069*param(2)*sin(q(1)))/43525000))/10000 + (271*param(1)*((20193591868813040730741*param(2)*...
         cos(q(1)))/1765583755997877910404608224447692800000 + (30385069*param(2)*sin(q(1)))/43525000))/10000 - (81*param(1)*((100521298322732069556673*param(2)*...
         cos(q(1)))/3531167511995755820809216448895385600000 - (48665257*param(2)*sin(q(1)))/87050000))/5000 + (81*param(1)*((100521298322732069556673*param(2)*...
         cos(q(1)))/3531167511995755820809216448895385600000 + (48665257*param(2)*sin(q(1)))/87050000))/5000 + (414571119887*param(1)*param(2)*...
         sin(q(1)))/870500000000))/(4863*param(1));
  p_O_CoM(2,1)=(5000*((61*param(1)*((379651*param(2)*cos(q(1)))/870500 - (1536770845809054406339*...
         param(2)*sin(q(1)))/35311675119957558208092164488953856000))/10000 + (61*param(1)*((379651*param(2)*...
         cos(q(1)))/870500 + (1536770845809054406339*param(2)*sin(q(1)))/35311675119957558208092164488953856000))/10000 + (271*...
         param(1)*((30385069*param(2)*cos(q(1)))/43525000 - (20193591868813040730741*param(2)*...
         sin(q(1)))/1765583755997877910404608224447692800000))/10000 + (271*param(1)*((30385069*param(2)*cos(q(1)))/43525000 + (20193591868813040730741*param(2)*...
         sin(q(1)))/1765583755997877910404608224447692800000))/10000 + (81*param(1)*((48665257*param(2)*cos(q(1)))/87050000 - (100521298322732069556673*param(2)*...
         sin(q(1)))/3531167511995755820809216448895385600000))/5000 + (81*param(1)*((48665257*param(2)*cos(q(1)))/87050000 + (100521298322732069556673*param(2)*...
         sin(q(1)))/3531167511995755820809216448895385600000))/5000 + (414571119887*param(1)*param(2)*cos(q(1)))/870500000000))/(4863*param(1));

 v_O_CoM = zeros([2,1]) ;
  v_O_CoM(1,1)=-(467907820151*Dq(1)*param(2)*cos(q(1)))/846648300000;
  v_O_CoM(2,1)=-(467907820151*Dq(1)*param(2)*sin(q(1)))/846648300000;

 m_tot = zeros([1,1]) ;
  m_tot(1,1)=(4863*param(1))/5000;

 