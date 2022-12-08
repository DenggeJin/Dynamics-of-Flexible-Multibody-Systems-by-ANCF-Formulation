function Cq = SysConsJacoGen(elmnum_1,dof_sys)
%--------------------------------------------------------------------------


Cq = zeros(6,dof_sys);
Cq(1,1) = 1;
Cq(2,2) = 1;
Cq(3,3) = 1;
Cq(4,12*elmnum_1+1) = 1;
Cq(4,12*(elmnum_1+1)+1) = -1;
Cq(5,12*elmnum_1+2) = 1;
Cq(5,12*(elmnum_1+1)+2) = -1;
Cq(6,12*elmnum_1+3) = 1;
Cq(6,12*(elmnum_1+1)+3) = -1;

end