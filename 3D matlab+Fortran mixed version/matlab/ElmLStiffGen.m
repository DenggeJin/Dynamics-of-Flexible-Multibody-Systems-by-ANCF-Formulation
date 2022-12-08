function LK_elm = ElmLStiffGen(E,nu,elmlen,elmwid,elmhgt)


lamda = nu*E/((1+nu)*(1-2*nu)); %������÷����
mu = E/(2*(1+nu)); %������÷����
halfelmhgt = elmhgt/2;
halfelmwid = elmwid/2; %���㵥Ԫ�İ���
%���㵥Ԫ�İ�߶�
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
Sx = diff(S,x); %��Ԫ�κ����ֱ��x��ƫ��������
Sy = diff(S,y); %��Ԫ�κ����ֱ��y��ƫ��������
Sz = diff(S,z); %��Ԫ�κ����ֱ��y��ƫ��������
LKint_elm_xyz = Sx'*Sx + Sy'*Sy + Sz'*Sz; %���ɱ�������
LKint_elm_yz = int(LKint_elm_xyz, x, 0, elmlen); %�Ա���x���л���
LKint_elm_z = int(LKint_elm_yz, y, -halfelmhgt, halfelmhgt); %�Ա���y���л���
LKint_elm = int(LKint_elm_z, z, -halfelmwid, halfelmwid);  %�Ա���z���л���
LK_elm = -(lamda + mu) * eval(LKint_elm); %���㵥Ԫά���µĵ�Ԫ�նȾ������Բ���

end