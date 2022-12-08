
clc
clear
startdatestr = datestr(now,'HH:MM:SS.FFF  ');
disp([startdatestr,'��ʼ��������˫������:'])
%----------����������-----------------------------------------------------
disp([datestr(now,'HH:MM:SS.FFF  '),'��ʼ�����������......']);
timestep = 0.001; %ʱ�䲽�������������
totaltime = 1; %��ʱ�������������
time = 0:timestep:totaltime;
timenum = length(time); %������ʱ�����
grav = 9.81; %�������ٶȣ����������
%----------���벿��1�Ĳ���--------------------------------------------------
E_1 = 7.0e9; %����1�Ĳ��ϵ���ģ�������������
nu_1 = 0.22; %����1�Ĳ��ϲ��ɱȣ����������
rho_1 = 6500; %����1�Ĳ����ܶȣ����������
len_1 = 0.5; %����1�ļ��γ��ȳߴ磨���������
wid_1 = 0.03; %����1�ļ��ο�ȳߴ磨���������
hgt_1 = 0.03; %����1�ļ��θ߶ȳߴ磨���������
elmnum_1 = 2; %����1���ֵ�Ԫ���������������
%----------���벿��2�Ĳ���--------------------------------------------------
E_2 = 7.0E9; %����2�Ĳ��ϵ���ģ�������������
nu_2 = 0.22; %����2�Ĳ��ϲ��ɱȣ����������
rho_2 = 6500; %����2�Ĳ����ܶȣ����������
len_2 = 1.5; %����2�ļ��γ��ȳߴ磨���������
wid_2 = 0.03; %����2�ļ��ο�ȳߴ磨���������
hgt_2 = 0.03; %����2�ļ��θ߶ȳߴ磨���������
elmnum_2 = 6; %����2���ֵ�Ԫ���������������
disp([datestr(now,'HH:MM:SS.FFF  '),'������ȡ���.']);
%----------���㲿���ĵ�Ԫ�ߴ�------------------------------------------------
elmlen_1 = len_1/elmnum_1; %���㲿��1�ĵ�Ԫ����
elmwid_1 = wid_1; %���㲿��1�ĵ�Ԫ���
elmhgt_1 = hgt_1; %���㲿��1�ĵ�Ԫ�߶�
dof_comp1 = (elmnum_1 + 1)*12; %���㲿��1�����ɶ���
elmlen_2 = len_2/elmnum_2; %���㲿��2�ĵ�Ԫ����
elmwid_2 = wid_2; %���㲿��2�ĵ�Ԫ���
elmhgt_2 = hgt_2; %���㲿��2�ĵ�Ԫ�߶�
dof_comp2 = (elmnum_2 + 1)*12; %���㲿��2�����ɶ���
dof_sys = dof_comp1 + dof_comp2; %����ϵͳ�����ɶ���
%----------���㵥Ԫ�ĸնȾ�������Բ��ֵĲ������-----------------------------
disp([datestr(now,'HH:MM:SS.FFF  '),'��ʼ���㵥Ԫ�նȾ�������Բ��ֵĲ������......']);
CK_elm_1 = CKelmGen(E_1,nu_1,elmlen_1,elmwid_1,elmhgt_1); %���㲿��1�ĵ�Ԫά���µĵ�Ԫ�նȾ�������Բ��ֵĲ������
CK_elm_2 = CKelmGen(E_2,nu_2,elmlen_2,elmwid_2,elmhgt_2); %���㲿��2�ĵ�Ԫά���µĵ�Ԫ�նȾ�������Բ��ֵĲ������
PrintTimeStr = datestr(now,'HH:MM:SS.FFF  ');
disp([datestr(now,'HH:MM:SS.FFF  '),'�������������.']);
%----------����ϵͳ��������-------------------------------------------------
disp([datestr(now,'HH:MM:SS.FFF  '),'��ʼ����ϵͳ��������......']);
Beq = LagrangeBooleanMatrix(dof_sys,elmnum_1,elmnum_2); %����ϵͳ����任����
disp([datestr(now,'HH:MM:SS.FFF  '),'ϵͳ��������������.']);
%----------����ϵͳ��ʼ����-------------------------------------------------
disp([datestr(now,'HH:MM:SS.FFF  '),'��ʼ����ϵͳ��ʼ����......']);
[q0_sys, dq0_sys] = SysInitCondGen(len_1,elmlen_1,elmnum_1,elmlen_2,elmnum_2,dof_sys); %����ϵͳ��ʼ����
disp([datestr(now,'HH:MM:SS.FFF  '),'ϵͳ��ʼ�����������.']);
%----------����ϵͳ�Ķ�����������--------------------------------------------
disp([datestr(now,'HH:MM:SS.FFF  '),'��ʼ����ϵͳ�Ķ�����������......']);
M_sys = SysMassMatrix(Beq,rho_1,rho_2,elmlen_1,elmlen_2,elmwid_1,elmwid_2,elmhgt_1,elmhgt_2,elmnum_1,elmnum_2,dof_sys); %����ϵͳ�Ķ�����������
disp([datestr(now,'HH:MM:SS.FFF  '),'ϵͳ��������������.']);
%----------����ϵͳ�ĸնȾ������Բ���----------------------------------------
disp([datestr(now,'HH:MM:SS.FFF  '),'��ʼ����ϵͳ�ĸնȾ������Բ���......']);
LK_elm_1 = ElmLStiffGen(E_1,nu_1,elmlen_1,elmwid_1,elmhgt_1); %���ɲ���1�ĵ�Ԫά���µĵ�Ԫ�նȾ������Բ���
LK_elm_2 = ElmLStiffGen(E_2,nu_2,elmlen_2,elmwid_2,elmhgt_2); %���ɲ���1�ĵ�Ԫά���µĵ�Ԫ�նȾ������Բ���
LK_sys = SysLStiffMatrix(Beq,E_1,E_2,nu_1,nu_2,elmlen_1,elmlen_2,elmwid_1,elmwid_2,elmhgt_1,elmhgt_2,elmnum_1,elmnum_2,dof_sys); %����ϵͳ�ĸնȾ������Բ���
disp([datestr(now,'HH:MM:SS.FFF  '),'ϵͳ�ĸնȾ������Բ��ּ������.']);
%----------����ϵͳ�Ĺ�����������--------------------------------------------
disp([datestr(now,'HH:MM:SS.FFF  '),'��ʼ����ϵͳ�Ĺ�����������......']);
Qg_sys = SysGForceVector(Beq,rho_1,rho_2,grav,elmlen_1,elmlen_2,elmwid_1,elmwid_2,elmhgt_1,elmhgt_2,elmnum_1,elmnum_2,dof_sys); %����ϵͳ�Ĺ�����������
disp([datestr(now,'HH:MM:SS.FFF  '),'ϵͳ�Ĺ������������������.']);
%----------�����ʼʱ��ϵͳ�նȾ�������Բ���---------------------------------
disp([datestr(now,'HH:MM:SS.FFF  '),'��ʼ�����ʼʱ��ϵͳ�նȾ�������Բ���......']);
NLK0_sys = SysNLStiffMatrix(Beq,CK_elm_1,CK_elm_2,q0_sys,elmnum_1,elmnum_2,dof_sys); %�����ʼʱ��ϵͳ�նȾ�������Բ���
disp([datestr(now,'HH:MM:SS.FFF  '),'��ʼʱ��ϵͳ�նȾ�������Բ��ּ������.']);
%----------�����ʼʱ��ϵͳ�նȾ���------------------------------------------
disp([datestr(now,'HH:MM:SS.FFF  '),'��ʼ�����ʼʱ��ϵͳ�նȾ���......']);
K0_sys = LK_sys + NLK0_sys; %�����ʼʱ��ϵͳ�նȾ���
disp([datestr(now,'HH:MM:SS.FFF  '),'��ʼʱ��ϵͳ�նȾ���������.']);
%----------����ϵͳԼ�����̵��ſ˱�------------------------------------------
disp([datestr(now,'HH:MM:SS.FFF  '),'��ʼ����ϵͳԼ�����̵��ſ˱�......']);
Cq = SysConsJacoGen(elmnum_1,dof_sys); %����ϵͳԼ�����̵��ſ˱�
disp([datestr(now,'HH:MM:SS.FFF  '),'ϵͳԼ�������ſ˱��������.']);
%----------�����ʼʱ��ϵͳ������ٶ��Լ��������ճ���-------------------------
disp([datestr(now,'HH:MM:SS.FFF  '),'��ʼ�����ʼʱ��ϵͳ������ٶ��Լ��������ճ���......']);
[ddq0_sys,lamda0] = SysInitAccel(M_sys,K0_sys,Qg_sys,Cq,q0_sys,dof_sys); %�����ʼʱ��ϵͳ������ٶ��Լ��������ճ���
disp([datestr(now,'HH:MM:SS.FFF  '),'��ʼʱ��ϵͳ������ٶ��Լ��������ճ��Ӽ������.']);

save('data.mat')
%%

save2txt( 'M.txt', M_sys);
save2txt( 'LK.txt', LK_sys);
save2txt( 'q0.txt', q0_sys);
save2txt('lamda0.txt',lamda0);
save2txt( 'dq0.txt', dq0_sys);
save2txt( 'ddq0.txt', ddq0_sys);
save2txt('dof_sys.txt',dof_sys);
save2txt_2( 'CK_1.txt', CK_elm_1);
save2txt_2('CK_2.txt',CK_elm_2);
save2txt_1( 'Beq.txt', Beq);
save2txt( 'Cq.txt', Cq);
save2txt('Qg.txt',Qg_sys);
save2txt('LK_1.txt',LK_elm_1);
save2txt('LK_2.txt',LK_elm_2);


%%
disp([datestr(now,'HH:MM:SS.FFF  '),'��ʼ���벢��ͼ..................']);

q_sys=read(elmlen_1,elmhgt_1,elmwid_1,elmnum_1,elmnum_2);

disp([datestr(now,'HH:MM:SS.FFF  '),'���벢��ͼ����']);






