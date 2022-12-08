function Cq = SysConsJacoGen(elmnum_1,dof_sys)
%--------------------------------------------------------------------------

Cq = zeros(4,dof_sys);
Cq(1,1) = 1;
Cq(2,2) = 1;
Cq(3,6*elmnum_1+1)=1;
Cq(3,6*(elmnum_1+1)+1)=-1;
Cq(4,6*elmnum_1+2)=1;
Cq(4,6*(elmnum_1+1)+2)=-1;

end