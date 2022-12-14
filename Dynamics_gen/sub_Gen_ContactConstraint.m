% sub_Gen_ContactConstraint
% p_O_C*: constant
% p_O_C*1: tree-branch kinematics

ContactConstraint = [...
    p_O_RC;
    p_O_LC];
ContactJacobian = jacobian(ContactConstraint,q);
Jac_Con = ContactJacobian;

dContactJacobian = jacobian(Jac_Con*Dq,q);
dJac_Con = dContactJacobian;




