function Mn_sys = SysMassMatrix(Beqn,rho_1,rho_2,elmlen_1,elmlen_2,elmwid_1,elmwid_2,elmhgt_1,elmhgt_2,elmnum_1,elmnum_2,dof_sys)


mass_elm_1 = ElmMassGen(rho_1,elmlen_1,elmwid_1,elmhgt_1); %生成部件1的单元维度下的单元质量矩阵
mass_elm_2 = ElmMassGen(rho_2,elmlen_2,elmwid_2,elmhgt_2); %生成部件2的单元维度下的单元质量矩阵
Mn_sys = zeros(dof_sys,dof_sys); %初始化部件维度下的部件质量矩阵
for k = 1:elmnum_1
    Mn_sys = Mn_sys + Beqn(:,:,k)'*mass_elm_1*Beqn(:,:,k); %将单元质量矩阵组装到部件质量矩阵中
end
for k = elmnum_1+1:elmnum_1+elmnum_2
    Mn_sys = Mn_sys + Beqn(:,:,k)'*mass_elm_2*Beqn(:,:,k); %将单元质量矩阵组装到部件质量矩阵中
end

end