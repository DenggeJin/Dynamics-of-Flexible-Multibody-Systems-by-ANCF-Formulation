function LK_elm = ElmLStiffGen(E,nu,elmlen,elmwid,elmhgt)


lamda = nu*E/((1+nu)*(1-2*nu)); %������÷����
mu = E/(2*(1+nu)); %������÷����
halfelmhgt = elmhgt/2; %���㵥Ԫ�İ�߶�
syms x y
s1 = 1 - 3*(x/elmlen)^2 + 2*(x/elmlen)^3;
s2 = elmlen*(x/elmlen - 2*(x/elmlen)^2 + (x/elmlen)^3);
s3 = elmlen*(y/elmlen - (x*y)/(elmlen^2));
s4 = 3*(x/elmlen)^2 - 2*(x/elmlen)^3;
s5 = elmlen*(-(x/elmlen)^2 + (x/elmlen)^3);
s6 = elmlen*((x*y)/(elmlen^2));
I = eye(2);
S = [s1*I, s2*I, s3*I, s4*I, s5*I, s6*I]; %���ɵ�Ԫ���κ�������
Sx = diff(S,x); %��Ԫ�κ����ֱ��x��ƫ��������
Sy = diff(S,y); %��Ԫ�κ����ֱ��y��ƫ��������
LKint_elm_xyz = Sx'*Sx + Sy'*Sy; %���ɱ�������
LKint_elm_yz = int(LKint_elm_xyz, x, 0, elmlen); %�Ա���x���л���
LKint_elm_z = int(LKint_elm_yz, y, -halfelmhgt, halfelmhgt); %�Ա���y���л���
LKint_elm = elmwid * LKint_elm_z; %�Ա���z���л���
LK_elm = -(lamda + mu) * eval(LKint_elm); %���㵥Ԫά���µĵ�Ԫ�նȾ������Բ���

end