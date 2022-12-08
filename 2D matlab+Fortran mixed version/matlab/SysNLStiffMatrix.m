function NLKn_sys = SysNLStiffMatrix(Beqn,CK_elm_1,CK_elm_2,qn_sys,elmnum_1,elmnum_2,dof_sys)
%--------------------------------------------------------------------------


NLKn_sys = zeros(dof_sys,dof_sys);
for k = 1:elmnum_1
    NLK_elm_1 = ElmNLStiffGen(CK_elm_1,Beqn(:,:,k)*qn_sys); %���ɲ���1��ϵͳά���µĵ�Ԫ�նȾ�������Բ���
    NLKn_sys = NLKn_sys + Beqn(:,:,k)'*NLK_elm_1*Beqn(:,:,k); %����Ԫ�նȾ�������Բ�����װ��ϵͳ�նȾ�������Բ�����
end
for k = elmnum_1+1:elmnum_1+elmnum_2
    NLK_elm_2 = ElmNLStiffGen(CK_elm_2,Beqn(:,:,k)*qn_sys); %���ɲ���2��ϵͳά���µĵ�Ԫ�նȾ�������Բ���
    NLKn_sys = NLKn_sys + Beqn(:,:,k)'*NLK_elm_2*Beqn(:,:,k); %����Ԫ�նȾ�������Բ�����װ��ϵͳ�նȾ�������Բ�����
end

end