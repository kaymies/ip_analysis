function [m_1,c_1,j_1,L_1,m_2,c_2,j_2,L_2] = auto_LumpParam(param)

m_1 = zeros([1,1]) ;
  m_1(1,1)=(1849*param(1))/5000;

 c_1 = zeros([1,1]) ;
  c_1(1,1)=(21666915091*param(2))/64382180000;

 j_1 = zeros([1,1]) ;
  j_1(1,1)=(18052288961074124251*param(1)*param(2)^2)/2802234384500000000000;

 L_1 = zeros([1,1]) ;
  L_1(1,1)=(1725*param(2))/3482;

 m_2 = zeros([1,1]) ;
  m_2(1,1)=(1507*param(1))/2500;

 c_2 = zeros([1,1]) ;
  c_2(1,1)=(491139150059638097397*param(2))/2253729971988070400000;

 j_2 = zeros([1,1]) ;
  j_2(1,1)=(160134257266760952562750583344901651816917*param(1)*param(2)^2)/8426175823884121748261778227200000000000000;

 L_2 = zeros([1,1]) ;
  L_2(1,1)=(4231*param(2))/8705;

 