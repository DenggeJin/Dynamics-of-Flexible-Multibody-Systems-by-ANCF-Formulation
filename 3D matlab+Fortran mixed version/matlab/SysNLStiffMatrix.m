function NLKn_sys = SysNLStiffMatrix(Beqn,CK_elm_1,CK_elm_2,qn_sys,elmnum_1,elmnum_2,dof_sys)
%--------------------------------------------------------------------------


NLKn_sys = zeros(dof_sys,dof_sys);
for k = 1:elmnum_1
    NLK_elm_1 = ElmNLStiffGen(CK_elm_1,Beqn(:,:,k)*qn_sys); %生成部件1上系统维度下的单元刚度矩阵非线性部分
    NLKn_sys = NLKn_sys + Beqn(:,:,k)'*NLK_elm_1*Beqn(:,:,k); %将单元刚度矩阵非线性部分组装到系统刚度矩阵非线性部分中
end
for k = elmnum_1+1:elmnum_1+elmnum_2
    NLK_elm_2 = ElmNLStiffGen(CK_elm_2,Beqn(:,:,k)*qn_sys); %生成部件2上系统维度下的单元刚度矩阵非线性部分
    NLKn_sys = NLKn_sys + Beqn(:,:,k)'*NLK_elm_2*Beqn(:,:,k); %将单元刚度矩阵非线性部分组装到系统刚度矩阵非线性部分中
end

end