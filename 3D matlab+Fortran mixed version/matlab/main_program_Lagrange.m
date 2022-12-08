
clc
clear
startdatestr = datestr(now,'HH:MM:SS.FFF  ');
disp([startdatestr,'开始计算柔性双摆算例:'])
%----------输入求解参数-----------------------------------------------------
disp([datestr(now,'HH:MM:SS.FFF  '),'开始读入输入参数......']);
timestep = 0.001; %时间步长（输入参数）
totaltime = 1; %总时长（输入参数）
time = 0:timestep:totaltime;
timenum = length(time); %计算总时间点数
grav = 9.81; %重力加速度（输入参数）
%----------输入部件1的参数--------------------------------------------------
E_1 = 7.0e9; %部件1的材料弹性模量（输入参数）
nu_1 = 0.22; %部件1的材料泊松比（输入参数）
rho_1 = 6500; %部件1的材料密度（输入参数）
len_1 = 0.5; %部件1的几何长度尺寸（输入参数）
wid_1 = 0.03; %部件1的几何宽度尺寸（输入参数）
hgt_1 = 0.03; %部件1的几何高度尺寸（输入参数）
elmnum_1 = 2; %部件1划分单元数量（输入参数）
%----------输入部件2的参数--------------------------------------------------
E_2 = 7.0E9; %部件2的材料弹性模量（输入参数）
nu_2 = 0.22; %部件2的材料泊松比（输入参数）
rho_2 = 6500; %部件2的材料密度（输入参数）
len_2 = 1.5; %部件2的几何长度尺寸（输入参数）
wid_2 = 0.03; %部件2的几何宽度尺寸（输入参数）
hgt_2 = 0.03; %部件2的几何高度尺寸（输入参数）
elmnum_2 = 6; %部件2划分单元数量（输入参数）
disp([datestr(now,'HH:MM:SS.FFF  '),'参数读取完成.']);
%----------计算部件的单元尺寸------------------------------------------------
elmlen_1 = len_1/elmnum_1; %计算部件1的单元长度
elmwid_1 = wid_1; %计算部件1的单元宽度
elmhgt_1 = hgt_1; %计算部件1的单元高度
dof_comp1 = (elmnum_1 + 1)*12; %计算部件1的自由度数
elmlen_2 = len_2/elmnum_2; %计算部件2的单元长度
elmwid_2 = wid_2; %计算部件2的单元宽度
elmhgt_2 = hgt_2; %计算部件2的单元高度
dof_comp2 = (elmnum_2 + 1)*12; %计算部件2的自由度数
dof_sys = dof_comp1 + dof_comp2; %计算系统的自由度数
%----------计算单元的刚度矩阵非线性部分的不变矩阵-----------------------------
disp([datestr(now,'HH:MM:SS.FFF  '),'开始计算单元刚度矩阵非线性部分的不变矩阵......']);
CK_elm_1 = CKelmGen(E_1,nu_1,elmlen_1,elmwid_1,elmhgt_1); %计算部件1的单元维度下的单元刚度矩阵非线性部分的不变矩阵
CK_elm_2 = CKelmGen(E_2,nu_2,elmlen_2,elmwid_2,elmhgt_2); %计算部件2的单元维度下的单元刚度矩阵非线性部分的不变矩阵
PrintTimeStr = datestr(now,'HH:MM:SS.FFF  ');
disp([datestr(now,'HH:MM:SS.FFF  '),'不变矩阵计算完成.']);
%----------计算系统布尔矩阵-------------------------------------------------
disp([datestr(now,'HH:MM:SS.FFF  '),'开始计算系统布尔矩阵......']);
Beq = LagrangeBooleanMatrix(dof_sys,elmnum_1,elmnum_2); %计算系统坐标变换矩阵
disp([datestr(now,'HH:MM:SS.FFF  '),'系统布尔矩阵计算完成.']);
%----------计算系统初始条件-------------------------------------------------
disp([datestr(now,'HH:MM:SS.FFF  '),'开始计算系统初始条件......']);
[q0_sys, dq0_sys] = SysInitCondGen(len_1,elmlen_1,elmnum_1,elmlen_2,elmnum_2,dof_sys); %计算系统初始条件
disp([datestr(now,'HH:MM:SS.FFF  '),'系统初始条件计算完成.']);
%----------计算系统的定常质量矩阵--------------------------------------------
disp([datestr(now,'HH:MM:SS.FFF  '),'开始计算系统的定常质量矩阵......']);
M_sys = SysMassMatrix(Beq,rho_1,rho_2,elmlen_1,elmlen_2,elmwid_1,elmwid_2,elmhgt_1,elmhgt_2,elmnum_1,elmnum_2,dof_sys); %计算系统的定常质量矩阵
disp([datestr(now,'HH:MM:SS.FFF  '),'系统质量矩阵计算完成.']);
%----------计算系统的刚度矩阵线性部分----------------------------------------
disp([datestr(now,'HH:MM:SS.FFF  '),'开始计算系统的刚度矩阵线性部分......']);
LK_elm_1 = ElmLStiffGen(E_1,nu_1,elmlen_1,elmwid_1,elmhgt_1); %生成部件1的单元维度下的单元刚度矩阵线性部分
LK_elm_2 = ElmLStiffGen(E_2,nu_2,elmlen_2,elmwid_2,elmhgt_2); %生成部件1的单元维度下的单元刚度矩阵线性部分
LK_sys = SysLStiffMatrix(Beq,E_1,E_2,nu_1,nu_2,elmlen_1,elmlen_2,elmwid_1,elmwid_2,elmhgt_1,elmhgt_2,elmnum_1,elmnum_2,dof_sys); %计算系统的刚度矩阵线性部分
disp([datestr(now,'HH:MM:SS.FFF  '),'系统的刚度矩阵线性部分计算完成.']);
%----------计算系统的广义外力向量--------------------------------------------
disp([datestr(now,'HH:MM:SS.FFF  '),'开始计算系统的广义外力向量......']);
Qg_sys = SysGForceVector(Beq,rho_1,rho_2,grav,elmlen_1,elmlen_2,elmwid_1,elmwid_2,elmhgt_1,elmhgt_2,elmnum_1,elmnum_2,dof_sys); %计算系统的广义外力向量
disp([datestr(now,'HH:MM:SS.FFF  '),'系统的广义外力向量计算完成.']);
%----------计算初始时刻系统刚度矩阵非线性部分---------------------------------
disp([datestr(now,'HH:MM:SS.FFF  '),'开始计算初始时刻系统刚度矩阵非线性部分......']);
NLK0_sys = SysNLStiffMatrix(Beq,CK_elm_1,CK_elm_2,q0_sys,elmnum_1,elmnum_2,dof_sys); %计算初始时刻系统刚度矩阵非线性部分
disp([datestr(now,'HH:MM:SS.FFF  '),'初始时刻系统刚度矩阵非线性部分计算完成.']);
%----------计算初始时刻系统刚度矩阵------------------------------------------
disp([datestr(now,'HH:MM:SS.FFF  '),'开始计算初始时刻系统刚度矩阵......']);
K0_sys = LK_sys + NLK0_sys; %计算初始时刻系统刚度矩阵
disp([datestr(now,'HH:MM:SS.FFF  '),'初始时刻系统刚度矩阵计算完成.']);
%----------生成系统约束方程的雅克比------------------------------------------
disp([datestr(now,'HH:MM:SS.FFF  '),'开始生成系统约束方程的雅克比......']);
Cq = SysConsJacoGen(elmnum_1,dof_sys); %生成系统约束方程的雅克比
disp([datestr(now,'HH:MM:SS.FFF  '),'系统约束方程雅克比生成完成.']);
%----------计算初始时刻系统广义加速度以及拉格朗日乘子-------------------------
disp([datestr(now,'HH:MM:SS.FFF  '),'开始计算初始时刻系统广义加速度以及拉格朗日乘子......']);
[ddq0_sys,lamda0] = SysInitAccel(M_sys,K0_sys,Qg_sys,Cq,q0_sys,dof_sys); %计算初始时刻系统广义加速度以及拉格朗日乘子
disp([datestr(now,'HH:MM:SS.FFF  '),'初始时刻系统广义加速度以及拉格朗日乘子计算完成.']);

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
disp([datestr(now,'HH:MM:SS.FFF  '),'开始读入并绘图..................']);

q_sys=read(elmlen_1,elmhgt_1,elmwid_1,elmnum_1,elmnum_2);

disp([datestr(now,'HH:MM:SS.FFF  '),'读入并绘图结束']);






