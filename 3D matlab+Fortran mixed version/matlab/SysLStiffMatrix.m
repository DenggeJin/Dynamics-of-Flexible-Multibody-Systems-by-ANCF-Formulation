function LKn_sys = SysLStiffMatrix(Beqn,E_1,E_2,nu_1,nu_2,elmlen_1,elmlen_2,elmwid_1,elmwid_2,elmhgt_1,elmhgt_2,elmnum_1,elmnum_2,dof_sys)
%--------------------------------------------------------------------------


LK_elm_1 = ElmLStiffGen(E_1,nu_1,elmlen_1,elmwid_1,elmhgt_1); %生成部件1的单元维度下的单元刚度矩阵线性部分
LK_elm_2 = ElmLStiffGen(E_2,nu_2,elmlen_2,elmwid_2,elmhgt_2); %生成部件1的单元维度下的单元刚度矩阵线性部分
LKn_sys = zeros(dof_sys,dof_sys); %初始化部件维度下的部件刚度矩阵线性部分
for k = 1:elmnum_1
    LKn_sys = LKn_sys + Beqn(:,:,k)'*LK_elm_1*Beqn(:,:,k); %将单元刚度矩阵线性部分组装到部件刚度矩阵线性部分中
end
for k = elmnum_1+1:elmnum_1+elmnum_2
    LKn_sys = LKn_sys + Beqn(:,:,k)'*LK_elm_2*Beqn(:,:,k); %将单元刚度矩阵线性部分组装到部件刚度矩阵线性部分中
end

end