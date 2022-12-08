function Qgn_sys = SysGForceVector(Beqn,rho_1,rho_2,grav,elmlen_1,elmlen_2,elmwid_1,elmwid_2,elmhgt_1,elmhgt_2,elmnum_1,elmnum_2,dof_sys)
%--------------------------------------------------------------------------



Qg_elm_1 = ElmGForceGen(rho_1,grav,elmlen_1,elmwid_1,elmhgt_1); %���ɲ���1�ĵ�Ԫά���µĵ�Ԫ������������
Qg_elm_2 = ElmGForceGen(rho_2,grav,elmlen_2,elmwid_2,elmhgt_2); %���ɲ���2�ĵ�Ԫά���µĵ�Ԫ������������
Qgn_sys = zeros(dof_sys,1); %��ʼ������ά���µĲ���������������
for k = 1:elmnum_1
    Qgn_sys = Qgn_sys + Beqn(:,:,k)'*Qg_elm_1; %����Ԫ��������������װ��������������������
end
for k = elmnum_1+1:elmnum_1+elmnum_2
    Qgn_sys = Qgn_sys + Beqn(:,:,k)'*Qg_elm_2; %����Ԫ��������������װ��������������������
end

end