function M_elm = ElmMassGen(rho,elmlen,elmwid,elmhgt)


halfelmhgt = elmhgt/2; %计算单元的半高度
halfelmwid = elmwid/2; %计算单元的半宽度
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
Mint_elm_xyz = S'*S; %生成被积函数
Mint_elm_yz = int(Mint_elm_xyz, x, 0, elmlen); %对变量x进行积分
Mint_elm_z = int(Mint_elm_yz, y, -halfelmhgt, halfelmhgt); %对变量y进行积分
Mint_elm = int(Mint_elm_z, z, -halfelmwid, halfelmwid); %对变量z进行积分
M_elm = rho * eval(Mint_elm); %计算单元维度下的单元质量矩阵

end