function CK_elm = CKelmGen(E,nu,elmlen,elmwid,elmhgt)


lamda = nu*E/((1+nu)*(1-2*nu)); %������÷����
mu = E/(2*(1+nu)); %������÷����
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
Sx = diff(S,x); %��Ԫ�κ����ֱ��x��ƫ��������
Sy = diff(S,y); %��Ԫ�κ����ֱ��y��ƫ��������
Sz = diff(S,z); %��Ԫ�κ����ֱ��y��ƫ��������
S11 = Sx'*Sx; 
S22 = Sy'*Sy;
S33 = Sz'*Sz;
S12 = Sx'*Sy;
S21 = Sy'*Sx;
S13 = Sx'*Sz;
S31 = Sz'*Sx;
S23 = Sy'*Sz;
S32 = Sz'*Sy;
CK1int = zeros(24,24,24,24); %��ʼ��CK1int
CK2int = zeros(24,24,24,24); %��ʼ��CK2int
CK3int = zeros(24,24,24,24); %��ʼ��CK3int
for i = 1:24
    for j = 1:24
        %----------����CK1int----------------------------------------------
        CK1int_elm_xyz = S11(i,:)'*S11(j,:) + S22(i,:)'*S22(j,:) + S33(i,:)'*S33(j,:); %���ɱ�������
        CK1int_elm_yz = int(CK1int_elm_xyz, x, 0, elmlen); %�Ա���x���л���
        CK1int_elm_z = int(CK1int_elm_yz, y, -halfelmhgt, halfelmhgt); %�Ա���y���л���
        CK1int_elm = int(CK1int_elm_z, z, -halfelmwid, halfelmwid); %�Ա���z���л���
        CK1int(:,:,i,j) = eval(CK1int_elm);
        %----------����CK2int----------------------------------------------
        CK2int_elm_xyz = S11(i,:)'*S22(j,:) + S22(i,:)'*S11(j,:) + S11(i,:)'*S33(j,:) + S33(i,:)'*S11(j,:) + S33(i,:)'*S22(j,:) + S22(i,:)'*S33(j,:); %���ɱ�������
        CK2int_elm_yz = int(CK2int_elm_xyz, x, 0, elmlen); %�Ա���x���л���
        CK2int_elm_z = int(CK2int_elm_yz, y, -halfelmhgt, halfelmhgt); %�Ա���y���л���
        CK2int_elm = int(CK2int_elm_z, z, -halfelmwid, halfelmwid); %�Ա���z���л���
        CK2int(:,:,i,j) = eval(CK2int_elm);
        %----------����CK3int----------------------------------------------
        CK3int_elm_xyz = S12(i,:)'*S21(j,:) + S21(i,:)'*S12(j,:) + S13(i,:)'*S31(j,:) + S31(i,:)'*S13(j,:) + S32(i,:)'*S23(j,:) + S23(i,:)'*S32(j,:); %���ɱ�������
        CK3int_elm_yz = int(CK3int_elm_xyz, x, 0, elmlen); %�Ա���x���л���
        CK3int_elm_z = int(CK3int_elm_yz, y, -halfelmhgt, halfelmhgt); %�Ա���y���л���
        CK3int_elm = int(CK3int_elm_z, z, -halfelmwid, halfelmwid); %�Ա���z���л���
        CK3int(:,:,i,j) = eval(CK3int_elm);
    end
end
CK_elm = (lamda + 2*mu)/2*CK1int + lamda/2*CK2int + mu*CK3int; %���㵥Ԫά���µĵ�Ԫ�նȾ�������Բ��ֶ�Ӧ�Ĳ������

end