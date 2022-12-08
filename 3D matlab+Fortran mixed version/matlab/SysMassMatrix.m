function Mn_sys = SysMassMatrix(Beqn,rho_1,rho_2,elmlen_1,elmlen_2,elmwid_1,elmwid_2,elmhgt_1,elmhgt_2,elmnum_1,elmnum_2,dof_sys)


mass_elm_1 = ElmMassGen(rho_1,elmlen_1,elmwid_1,elmhgt_1); %���ɲ���1�ĵ�Ԫά���µĵ�Ԫ��������
mass_elm_2 = ElmMassGen(rho_2,elmlen_2,elmwid_2,elmhgt_2); %���ɲ���2�ĵ�Ԫά���µĵ�Ԫ��������
Mn_sys = zeros(dof_sys,dof_sys); %��ʼ������ά���µĲ�����������
for k = 1:elmnum_1
    Mn_sys = Mn_sys + Beqn(:,:,k)'*mass_elm_1*Beqn(:,:,k); %����Ԫ����������װ����������������
end
for k = elmnum_1+1:elmnum_1+elmnum_2
    Mn_sys = Mn_sys + Beqn(:,:,k)'*mass_elm_2*Beqn(:,:,k); %����Ԫ����������װ����������������
end

end