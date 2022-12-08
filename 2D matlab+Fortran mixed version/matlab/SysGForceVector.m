function Qgn_sys = SysGForceVector(Beqn,rho_1,rho_2,grav,elmlen_1,elmlen_2,elmwid_1,elmwid_2,elmhgt_1,elmhgt_2,elmnum_1,elmnum_2,dof_sys)
%--------------------------------------------------------------------------



Qg_elm_1 = ElmGForceGen(rho_1,grav,elmlen_1,elmwid_1,elmhgt_1); %生成部件1的单元维度下的单元广义外力向量
Qg_elm_2 = ElmGForceGen(rho_2,grav,elmlen_2,elmwid_2,elmhgt_2); %生成部件2的单元维度下的单元广义外力向量
Qgn_sys = zeros(dof_sys,1); %初始化部件维度下的部件广义外力向量
for k = 1:elmnum_1
    Qgn_sys = Qgn_sys + Beqn(:,:,k)'*Qg_elm_1; %将单元广义外力向量组装到部件广义外力向量中
end
for k = elmnum_1+1:elmnum_1+elmnum_2
    Qgn_sys = Qgn_sys + Beqn(:,:,k)'*Qg_elm_2; %将单元广义外力向量组装到部件广义外力向量中
end

end