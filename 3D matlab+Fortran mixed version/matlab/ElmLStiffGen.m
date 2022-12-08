function LK_elm = ElmLStiffGen(E,nu,elmlen,elmwid,elmhgt)


lamda = nu*E/((1+nu)*(1-2*nu)); %计算拉梅常数
mu = E/(2*(1+nu)); %计算拉梅常数
halfelmhgt = elmhgt/2;
halfelmwid = elmwid/2; %计算单元的半宽度
%计算单元的半高度
syms x y z
s1 = 1 - 3*(x/elmlen)^2 + 2*(x/elmlen)^3;
s2 = elmlen*(x/elmlen - 2*(x/elmlen)^2 + (x/elmlen)^3);
s3 = elmlen*(y/elmlen - (x*y)/(elmlen^2));
s4 = elmlen*(z/elmlen - (x*z)/(elmlen^2));
s5 = 3*(x/elmlen)^2 - 2*(x/elmlen)^3;
s6 = elmlen*(-(x/elmlen)^2 + (x/elmlen)^3);
s7 = elmlen*((x*y)/(elmlen^2));
s8 = elmlen*((x*z)/(elmlen^2));
I = eye(3);
S = [s1*I, s2*I, s3*I, s4*I, s5*I, s6*I, s7*I, s8*I]; %生成单元的形函数矩阵
Sx = diff(S,x); %单元形函数分别对x的偏导数矩阵
Sy = diff(S,y); %单元形函数分别对y的偏导数矩阵
Sz = diff(S,z); %单元形函数分别对y的偏导数矩阵
LKint_elm_xyz = Sx'*Sx + Sy'*Sy + Sz'*Sz; %生成被积函数
LKint_elm_yz = int(LKint_elm_xyz, x, 0, elmlen); %对变量x进行积分
LKint_elm_z = int(LKint_elm_yz, y, -halfelmhgt, halfelmhgt); %对变量y进行积分
LKint_elm = int(LKint_elm_z, z, -halfelmwid, halfelmwid);  %对变量z进行积分
LK_elm = -(lamda + mu) * eval(LKint_elm); %计算单元维度下的单元刚度矩阵线性部分

end