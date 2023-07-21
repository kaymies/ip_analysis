function [m_1,c_1,j_1,L_1,m_2,c_2,j_2,L_2] = auto_LumpParam(param)

m_1 = zeros([1,1]) ;
  m_1(1,1)=(1959*param(1))/5000;

 c_1 = zeros([1,1]) ;
  c_1(1,1)=(13941533551*param(2))/42485812500;

 j_1 = zeros([1,1]) ;
  j_1(1,1)=(121133730011100474317*param(1)*param(2)^2)/18428221171875000000000;

 L_1 = zeros([1,1]) ;
  L_1(1,1)=(8071*param(2))/17350;

 m_2 = zeros([1,1]) ;
  m_2(1,1)=(5823*param(1))/10000;

 c_2 = zeros([1,1]) ;
  c_2(1,1)=(11299042597*param(2))/50514525000;

 j_2 = zeros([1,1]) ;
  j_2(1,1)=(770639201632582127401*param(1)*param(2)^2)/87642700875000000000000;

 L_2 = zeros([1,1]) ;
  L_2(1,1)=(1717*param(2))/3470;

 