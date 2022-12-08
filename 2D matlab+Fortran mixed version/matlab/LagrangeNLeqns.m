function reserr = LagrangeNLeqns(qlamda,q0_sys,dq0_sys,ddq0_sys,beta,gamma,timestep,elmnum_1,elmnum_2,Beq,M_sys,LK_sys,Qg_sys,Cq,CK_elm_1,CK_elm_2,dof_sys)
%--------------------------------------------------------------------------

qti_sys = qlamda(1:dof_sys);
lgi_sys = qlamda(dof_sys+1:end);

ddqti_sys = 1/(beta*timestep^2)*(qti_sys-q0_sys) - 1/(beta*timestep)*dq0_sys - (1/(2*beta)-1)*ddq0_sys;
dqti_sys = gamma/(beta*timestep)*(qti_sys-q0_sys) + (1-gamma/beta)*dq0_sys + (1-gamma/(2*beta))*timestep*ddq0_sys;

NLKi_sys = SysNLStiffMatrix(Beq,CK_elm_1,CK_elm_2,qti_sys,elmnum_1,elmnum_2,dof_sys);
Fi_sys = (LK_sys+NLKi_sys)*qti_sys - Qg_sys;

reserr1 = M_sys*ddqti_sys + Cq'*lgi_sys + Fi_sys;
reserr2 = Cq*ddqti_sys;
reserr = [reserr1; reserr2];

end