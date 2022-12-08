function LKn_sys = SysLStiffMatrix(Beqn,E_1,E_2,nu_1,nu_2,elmlen_1,elmlen_2,elmwid_1,elmwid_2,elmhgt_1,elmhgt_2,elmnum_1,elmnum_2,dof_sys)
%--------------------------------------------------------------------------


LK_elm_1 = ElmLStiffGen(E_1,nu_1,elmlen_1,elmwid_1,elmhgt_1); %���ɲ���1�ĵ�Ԫά���µĵ�Ԫ�նȾ������Բ���
LK_elm_2 = ElmLStiffGen(E_2,nu_2,elmlen_2,elmwid_2,elmhgt_2); %���ɲ���1�ĵ�Ԫά���µĵ�Ԫ�նȾ������Բ���
LKn_sys = zeros(dof_sys,dof_sys); %��ʼ������ά���µĲ����նȾ������Բ���
for k = 1:elmnum_1
    LKn_sys = LKn_sys + Beqn(:,:,k)'*LK_elm_1*Beqn(:,:,k); %����Ԫ�նȾ������Բ�����װ�������նȾ������Բ�����
end
for k = elmnum_1+1:elmnum_1+elmnum_2
    LKn_sys = LKn_sys + Beqn(:,:,k)'*LK_elm_2*Beqn(:,:,k); %����Ԫ�նȾ������Բ�����װ�������նȾ������Բ�����
end

end