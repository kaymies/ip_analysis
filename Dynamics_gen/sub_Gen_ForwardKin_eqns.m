%% To derive SE(2)
%% SO(2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%active
% base > FootFrame

% Rika 2020-12-13: upper body config dependent on plane when arms are
%                  abducted sideways, in the frontal plane
switch plane
    case 'frt'
        SE2_O_CF = ClassSE2(pi/2+q_A, 0);
        SE2_CF_CK = ClassSE2(q_K, Shank.L);
        SE2_CK_CH = ClassSE2(q_H, Thigh.L);
        SE2_CH_NE = ClassSE2(q_NE, Trunk.L);
        SE2_NE_HE = ClassSE2(0, Head.L);
        SE2_CH_RS = ClassSE2(q_RS, l_sjc, -w_sjc);
        SE2_CH_LS = ClassSE2(q_LS, l_sjc, w_sjc);
        SE2_RS_RE = ClassSE2(q_RE, Upperarm.L);
        SE2_LS_LE = ClassSE2(q_LE, Upperarm.L);
        SE2_RE_RW = ClassSE2(q_RW, Forearm.L);
        SE2_LE_LW = ClassSE2(q_LW, Forearm.L);
        SE2_RW_RH = ClassSE2(0, Hand.L);
        SE2_LW_LH = ClassSE2(0, Hand.L);

        SE2_CF_Shankc = ClassSE2(0, Shank.c);
        SE2_CK_Thighc = ClassSE2(0, Thigh.c);
        SE2_CH_Trunkc = ClassSE2(0, Trunk.c);
        SE2_NE_Headc = ClassSE2(0, Head.c);
        SE2_RS_RUarmc = ClassSE2(0, Upperarm.c);
        SE2_RE_RFarmc = ClassSE2(0, Forearm.c);
        SE2_RW_RHandc = ClassSE2(0, Hand.c);
        SE2_LS_LUarmc = ClassSE2(0, Upperarm.c);
        SE2_LE_LFarmc = ClassSE2(0, Forearm.c);
        SE2_LW_LHandc = ClassSE2(0, Hand.c);

        %% Forward kinematics using SE(2) for Joint
        SE2_O_CF = SE2_O_CF;
        SE2_O_CK = SE2_O_CF*SE2_CF_CK;
        SE2_O_CH = SE2_O_CK*SE2_CK_CH;
        SE2_O_NE = SE2_O_CH*SE2_CH_NE;
        SE2_O_HE = SE2_O_NE*SE2_NE_HE;
        SE2_O_RS = SE2_O_CH*SE2_CH_RS;
        SE2_O_RE = SE2_O_RS*SE2_RS_RE;
        SE2_O_RW = SE2_O_RE*SE2_RE_RW;
        SE2_O_RH = SE2_O_RW*SE2_RW_RH;
        SE2_O_LS = SE2_O_CH*SE2_CH_LS;
        SE2_O_LE = SE2_O_LS*SE2_LS_LE;
        SE2_O_LW = SE2_O_LE*SE2_LE_LW;
        SE2_O_LH = SE2_O_LW*SE2_LW_LH;
        % COM link
        SE2_O_Shankc = SE2_O_CF*SE2_CF_Shankc;
        SE2_O_Thighc = SE2_O_CK*SE2_CK_Thighc;
        SE2_O_Trunkc = SE2_O_CH*SE2_CH_Trunkc;
        SE2_O_Headc =  SE2_O_NE*SE2_NE_Headc;
        SE2_O_RUarmc = SE2_O_RS*SE2_RS_RUarmc; 
        SE2_O_RFarmc = SE2_O_RE*SE2_RE_RFarmc; 
        SE2_O_RHandc = SE2_O_RW*SE2_RW_RHandc; 
        SE2_O_LUarmc = SE2_O_LS*SE2_LS_LUarmc; 
        SE2_O_LFarmc = SE2_O_LE*SE2_LE_LFarmc; 
        SE2_O_LHandc = SE2_O_LW*SE2_LW_LHandc; 
%% maybe we can shorten the above a little more.
    case 'sgt'
        SE2_O_CF = ClassSE2(pi/2+q_A, 0);
        SE2_CF_CK = ClassSE2(q_K, Shank.L);
        SE2_CK_CH = ClassSE2(q_H, Thigh.L);
        SE2_CH_NE = ClassSE2(q_NE, Trunk.L);
        SE2_NE_HE = ClassSE2(0, Head.L);
        SE2_CH_RS = ClassSE2(q_RS, l_sjc, -w_sjc);
        SE2_CH_LS = ClassSE2(q_LS, l_sjc, w_sjc);
        SE2_RS_RE = ClassSE2(q_RE, 0);
        SE2_LS_LE = ClassSE2(q_LE, 0);
        SE2_RE_RW = ClassSE2(q_RW, 0);
        SE2_LE_LW = ClassSE2(q_LW, 0);
        SE2_RW_RH = ClassSE2(0, 0);
        SE2_LW_LH = ClassSE2(0, 0);

        SE2_CF_Shankc = ClassSE2(0, Shank.c);
        SE2_CK_Thighc = ClassSE2(0, Thigh.c);
        SE2_CH_Trunkc = ClassSE2(0, Trunk.c);
        SE2_NE_Headc = ClassSE2(0, Head.c);

        %% Forward kinematics using SE(2) for Joint
        SE2_O_CF = SE2_O_CF;
        SE2_O_CK = SE2_O_CF*SE2_CF_CK;
        SE2_O_CH = SE2_O_CK*SE2_CK_CH;
        SE2_O_NE = SE2_O_CH*SE2_CH_NE;
        SE2_O_HE = SE2_O_NE*SE2_NE_HE;
        SE2_O_RS = SE2_O_CH*SE2_CH_RS;
        SE2_O_RE = SE2_O_RS;
        SE2_O_RW = SE2_O_RS;
        SE2_O_RH = SE2_O_RS;
        SE2_O_LS = SE2_O_CH*SE2_CH_LS;
        SE2_O_LE = SE2_O_LS;
        SE2_O_LW = SE2_O_LS;
        SE2_O_LH = SE2_O_LS;
        % COM link
        SE2_O_Shankc = SE2_O_CF*SE2_CF_Shankc;
        SE2_O_Thighc = SE2_O_CK*SE2_CK_Thighc;
        SE2_O_Trunkc = SE2_O_CH*SE2_CH_Trunkc;
        SE2_O_Headc =  SE2_O_NE*SE2_NE_Headc;
        SE2_O_RUarmc = SE2_O_RS; 
        SE2_O_RFarmc = SE2_O_RS; 
        SE2_O_RHandc = SE2_O_RS; 
        SE2_O_LUarmc = SE2_O_LS; 
        SE2_O_LFarmc = SE2_O_LS; 
        SE2_O_LHandc = SE2_O_LS;
end