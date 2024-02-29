function [KE_tot,PE_tot] = auto_MechEnergy(q,Dq,param)

KE_tot = zeros([1,1]) ;
  KE_tot(1,1)=(param(1)*param(2)^2*...
         (7271883847371331515412139819865750138640461064351346634905814922506573535836577423*Dq(1)^2 + 999561669344098160147462680559017060529392670935261928654225475718295031708065423*...
         Dq(2)^2 + 3497503326155654882698478887140083015945147069525769871969888498412609863680000000*Dq(1)^2*...
         cos(q(2)) + 1999123338688196320294925361118034121058785341870523857308450951436590063416130846*Dq(1)*Dq(2) +...
          3497503326155654882698478887140083015945147069525769871969888498412609863680000000*Dq(1)*Dq(2)*cos(q(2))))/62345719988871481643264076810438956574919700129256211072232064863436800000000000000;

 PE_tot = zeros([1,1]) ;
  PE_tot(1,1)=(981*param(1)*param(2)*(99615744696*cos(q(1) + q(2)) + 368292075455*cos(q(1))))/87050000000000;

 