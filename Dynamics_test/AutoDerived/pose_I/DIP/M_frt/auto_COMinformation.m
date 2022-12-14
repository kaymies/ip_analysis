function [p_O_CoM,v_O_CoM,m_tot] = auto_COMinformation(q,Dq,param)

p_O_CoM = zeros([2,1]) ;
  p_O_CoM(1,1)=-(5000*((81*param(1)*((5540257*param(2)*sin(q(1) + q(2)))/87050000 -...
          (388428426319533039767715486646422859327*param(2)*cos(q(1) + q(2)))/3531167511995755820809216448895385600000 + (1725*param(2)*...
         sin(q(1)))/3482))/5000 + (81*param(1)*((388428426319533039767715486646422859327*param(2)*cos(q(1) +...
          q(2)))/3531167511995755820809216448895385600000 + (5540257*param(2)*sin(q(1) + q(2)))/87050000 + (1725*param(2)*sin(q(1)))/3482))/5000 + (271*...
         param(1)*((8822569*param(2)*sin(q(1) + q(2)))/43525000 - (194214213159766549950915035876205477259*param(2)*...
         cos(q(1) + q(2)))/1765583755997877910404608224447692800000 + (1725*param(2)*sin(q(1)))/3482))/10000 + (271*...
         param(1)*((194214213159766549950915035876205477259*param(2)*cos(q(1) +...
          q(2)))/1765583755997877910404608224447692800000 + (8822569*param(2)*sin(q(1) + q(2)))/43525000 + (1725*param(2)*sin(q(1)))/3482))/10000 + (61*...
         param(1)*((3884284263195329866119292284730517821*param(2)*cos(q(1) +...
          q(2)))/35311675119957558208092164488953856000 - (51599*param(2)*sin(q(1) + q(2)))/870500 + (1725*param(2)*sin(q(1)))/3482))/10000 - (61*param(1)*...
         ((3884284263195329866119292284730517821*param(2)*cos(q(1) + q(2)))/35311675119957558208092164488953856000 + (51599*param(2)*sin(q(1) +...
          q(2)))/870500 - (1725*param(2)*sin(q(1)))/3482))/10000 + (6519*param(1)*param(2)*(4888741*sin(q(1) + q(2)) +...
          14375000*sin(q(1))))/435250000000 + (1041*param(1)*param(2)*(12078357*sin(q(1) + q(2)) + 14375000*...
         sin(q(1))))/435250000000 + (21666915091*param(1)*param(2)*sin(q(1)))/174100000000))/(4863*param(1));
  p_O_CoM(2,1)=(5000*((81*param(1)*((5540257*param(2)*cos(q(1) + q(2)))/87050000 -...
          (388428426319533039767715486646422859327*param(2)*sin(q(1) + q(2)))/3531167511995755820809216448895385600000 + (1725*param(2)*...
         cos(q(1)))/3482))/5000 + (81*param(1)*((5540257*param(2)*cos(q(1) + q(2)))/87050000 +...
          (388428426319533039767715486646422859327*param(2)*sin(q(1) + q(2)))/3531167511995755820809216448895385600000 + (1725*param(2)*...
         cos(q(1)))/3482))/5000 + (271*param(1)*((8822569*param(2)*cos(q(1) + q(2)))/43525000 -...
          (194214213159766549950915035876205477259*param(2)*sin(q(1) + q(2)))/1765583755997877910404608224447692800000 + (1725*param(2)*...
         cos(q(1)))/3482))/10000 + (271*param(1)*((8822569*param(2)*cos(q(1) + q(2)))/43525000 +...
          (194214213159766549950915035876205477259*param(2)*sin(q(1) + q(2)))/1765583755997877910404608224447692800000 + (1725*param(2)*...
         cos(q(1)))/3482))/10000 + (61*param(1)*((3884284263195329866119292284730517821*param(2)*sin(q(1) +...
          q(2)))/35311675119957558208092164488953856000 - (51599*param(2)*cos(q(1) + q(2)))/870500 + (1725*param(2)*cos(q(1)))/3482))/10000 - (61*param(1)*...
         ((51599*param(2)*cos(q(1) + q(2)))/870500 + (3884284263195329866119292284730517821*param(2)*sin(q(1) +...
          q(2)))/35311675119957558208092164488953856000 - (1725*param(2)*cos(q(1)))/3482))/10000 + (21666915091*param(1)*param(2)*cos(q(1)))/174100000000 +...
          (6519*param(1)*param(2)*(4888741*cos(q(1) + q(2)) + 14375000*cos(q(1))))/435250000000 + (1041*param(1)*...
         param(2)*(12078357*cos(q(1) + q(2)) + 14375000*cos(q(1))))/435250000000))/(4863*param(1));

 v_O_CoM = zeros([2,1]) ;
  v_O_CoM(1,1)=-(param(2)*(99615744696*Dq(1)*cos(q(1) + q(2)) + 99615744696*Dq(2)*cos(q(1) + q(2)) +...
          368292075455*Dq(1)*cos(q(1))))/846648300000;
  v_O_CoM(2,1)=-(param(2)*(99615744696*Dq(1)*sin(q(1) + q(2)) + 99615744696*Dq(2)*sin(q(1) + q(2)) +...
          368292075455*Dq(1)*sin(q(1))))/846648300000;

 m_tot = zeros([1,1]) ;
  m_tot(1,1)=(4863*param(1))/5000;

 