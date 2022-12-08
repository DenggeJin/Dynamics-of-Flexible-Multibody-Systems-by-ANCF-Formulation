function CK_elm = CKelmGen(E,nu,elmlen,elmwid,elmhgt)


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
S11 = Sx'*Sx; 
S22 = Sy'*Sy;
S12 = Sx'*Sy;
S21 = Sy'*Sx;
CK1int = zeros(12,12,12,12); %��ʼ��CK1int
CK2int = zeros(12,12,12,12); %��ʼ��CK2int
CK3int = zeros(12,12,12,12); %��ʼ��CK3int
for i = 1:12
    for j = 1:12
        %----------����CK1int----------------------------------------------
        CK1int_elm_xyz = S11(i,:)'*S11(j,:) + S22(i,:)'*S22(j,:); %���ɱ�������
        CK1int_elm_yz = int(CK1int_elm_xyz, x, 0, elmlen); %�Ա���x���л���
        CK1int_elm_z = int(CK1int_elm_yz, y, -halfelmhgt, halfelmhgt); %�Ա���y���л���
        CK1int_elm = elmwid * CK1int_elm_z; %�Ա���z���л���
        CK1int(:,:,i,j) = eval(CK1int_elm);
        %----------����CK2int----------------------------------------------
        CK2int_elm_xyz = S11(i,:)'*S22(j,:) + S22(i,:)'*S11(j,:); %���ɱ�������
        CK2int_elm_yz = int(CK2int_elm_xyz, x, 0, elmlen); %�Ա���x���л���
        CK2int_elm_z = int(CK2int_elm_yz, y, -halfelmhgt, halfelmhgt); %�Ա���y���л���
        CK2int_elm = elmwid * CK2int_elm_z; %�Ա���z���л���
        CK2int(:,:,i,j) = eval(CK2int_elm);
        %----------����CK3int----------------------------------------------
        CK3int_elm_xyz = S12(i,:)'*S21(j,:) + S21(i,:)'*S12(j,:); %���ɱ�������
        CK3int_elm_yz = int(CK3int_elm_xyz, x, 0, elmlen); %�Ա���x���л���
        CK3int_elm_z = int(CK3int_elm_yz, y, -halfelmhgt, halfelmhgt); %�Ա���y���л���
        CK3int_elm = elmwid * CK3int_elm_z; %�Ա���z���л���
        CK3int(:,:,i,j) = eval(CK3int_elm);
    end
end
CK_elm = (lamda + 2*mu)/2*CK1int + lamda/2*CK2int + mu*CK3int; %���㵥Ԫά���µĵ�Ԫ�նȾ�������Բ��ֶ�Ӧ�Ĳ������

end