function M_elm = ElmMassGen(rho,elmlen,elmwid,elmhgt)


halfelmhgt = elmhgt/2; %���㵥Ԫ�İ�߶�
halfelmwid = elmwid/2; %���㵥Ԫ�İ���
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
S = [s1*I, s2*I, s3*I, s4*I, s5*I, s6*I, s7*I, s8*I]; %���ɵ�Ԫ���κ�������
Mint_elm_xyz = S'*S; %���ɱ�������
Mint_elm_yz = int(Mint_elm_xyz, x, 0, elmlen); %�Ա���x���л���
Mint_elm_z = int(Mint_elm_yz, y, -halfelmhgt, halfelmhgt); %�Ա���y���л���
Mint_elm = int(Mint_elm_z, z, -halfelmwid, halfelmwid); %�Ա���z���л���
M_elm = rho * eval(Mint_elm); %���㵥Ԫά���µĵ�Ԫ��������

end