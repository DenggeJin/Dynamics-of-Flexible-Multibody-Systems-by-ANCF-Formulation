function [q0_sys, dq0_sys] = SysInitCondGen(len_1,elmlen_1,elmnum_1,elmlen_2,elmnum_2,dof_sys)
%--------------------------------------------------------------------------


q0_sys = zeros(dof_sys,1); %��ʼ��ϵͳ��ʼʱ�̹�����������
dq0_sys = zeros(dof_sys,1); %��ʼ��ϵͳ��ʼʱ�̹����ٶ�����
for i = 1:elmnum_1+1
    q0_sys(6*(i-1)+1) = (i-1)*elmlen_1;
    q0_sys(6*(i-1)+3) = 1;
    q0_sys(6*(i-1)+6) = 1;
end
for i = 1:elmnum_2+1
    q0_sys(6*(i+elmnum_1)+1) = (i-1)*elmlen_2 + len_1;
    q0_sys(6*(i+elmnum_1)+3) = 1;
    q0_sys(6*(i+elmnum_1)+6) = 1;
end

end