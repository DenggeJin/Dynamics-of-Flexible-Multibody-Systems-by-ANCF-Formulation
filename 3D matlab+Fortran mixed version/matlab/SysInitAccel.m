function [ddq0_sys,lamda0] = SysInitAccel(M_sys,K0_sys,Qg0_sys,Cq,q0_sys,dof_sys)
%--------------------------------------------------------------------------


Mlamda11 = M_sys;
Mlamda12 = Cq';
Mlamda21 = Cq;
Mlamda22 = zeros(6,6);
Mlamda = [Mlamda11,Mlamda12;Mlamda21,Mlamda22];
F1q = Qg0_sys - K0_sys*q0_sys; 
F1lamda = zeros(6,1);
F1 = [F1q; F1lamda]; 
qlamda0 = Mlamda\F1;
ddq0_sys = qlamda0(1:dof_sys,1);
lamda0 = qlamda0(dof_sys+1:end,1);

end