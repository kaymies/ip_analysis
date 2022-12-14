%% Symbolic Variables
%% Joints
% % independent joints, r
% syms q_A q_H real
% syms Dq_A Dq_H real
% syms DDq_A DDq_H real
% % passive joints, s
% syms q_RS q_RW real
% syms Dq_RS Dq_RW real
% syms DDq_RS DDq_RW real
% 
% syms q_LS q_LW real
% syms Dq_LS Dq_LW real
% syms DDq_LS DDq_LW real

%% Joints
syms q_A q_K q_H real
syms q_RS q_LS q_RE q_LE q_RW q_LW real
syms q_NE real % neck

syms Dq_A Dq_K Dq_H real
syms Dq_RS Dq_LS Dq_RE Dq_LE Dq_RW Dq_LW real
syms Dq_NE real % neck

syms DDq_A DDq_K DDq_H real
syms DDq_RS DDq_LS DDq_RE DDq_LE DDq_RW DDq_LW real
syms DDq_NE real % neck