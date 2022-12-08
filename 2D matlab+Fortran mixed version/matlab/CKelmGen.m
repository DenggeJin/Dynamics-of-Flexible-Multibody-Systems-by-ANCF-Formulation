function CK_elm = CKelmGen(E,nu,elmlen,elmwid,elmhgt)


lamda = nu*E/((1+nu)*(1-2*nu)); %计算拉梅常数
mu = E/(2*(1+nu)); %计算拉梅常数
halfelmhgt = elmhgt/2; %计算单元的半高度
syms x y
s1 = 1 - 3*(x/elmlen)^2 + 2*(x/elmlen)^3;
s2 = elmlen*(x/elmlen - 2*(x/elmlen)^2 + (x/elmlen)^3);
s3 = elmlen*(y/elmlen - (x*y)/(elmlen^2));
s4 = 3*(x/elmlen)^2 - 2*(x/elmlen)^3;
s5 = elmlen*(-(x/elmlen)^2 + (x/elmlen)^3);
s6 = elmlen*((x*y)/(elmlen^2));
I = eye(2);
S = [s1*I, s2*I, s3*I, s4*I, s5*I, s6*I]; %生成单元的形函数矩阵
Sx = diff(S,x); %单元形函数分别对x的偏导数矩阵
Sy = diff(S,y); %单元形函数分别对y的偏导数矩阵
S11 = Sx'*Sx; 
S22 = Sy'*Sy;
S12 = Sx'*Sy;
S21 = Sy'*Sx;
CK1int = zeros(12,12,12,12); %初始化CK1int
CK2int = zeros(12,12,12,12); %初始化CK2int
CK3int = zeros(12,12,12,12); %初始化CK3int
for i = 1:12
    for j = 1:12
        %----------计算CK1int----------------------------------------------
        CK1int_elm_xyz = S11(i,:)'*S11(j,:) + S22(i,:)'*S22(j,:); %生成被积函数
        CK1int_elm_yz = int(CK1int_elm_xyz, x, 0, elmlen); %对变量x进行积分
        CK1int_elm_z = int(CK1int_elm_yz, y, -halfelmhgt, halfelmhgt); %对变量y进行积分
        CK1int_elm = elmwid * CK1int_elm_z; %对变量z进行积分
        CK1int(:,:,i,j) = eval(CK1int_elm);
        %----------计算CK2int----------------------------------------------
        CK2int_elm_xyz = S11(i,:)'*S22(j,:) + S22(i,:)'*S11(j,:); %生成被积函数
        CK2int_elm_yz = int(CK2int_elm_xyz, x, 0, elmlen); %对变量x进行积分
        CK2int_elm_z = int(CK2int_elm_yz, y, -halfelmhgt, halfelmhgt); %对变量y进行积分
        CK2int_elm = elmwid * CK2int_elm_z; %对变量z进行积分
        CK2int(:,:,i,j) = eval(CK2int_elm);
        %----------计算CK3int----------------------------------------------
        CK3int_elm_xyz = S12(i,:)'*S21(j,:) + S21(i,:)'*S12(j,:); %生成被积函数
        CK3int_elm_yz = int(CK3int_elm_xyz, x, 0, elmlen); %对变量x进行积分
        CK3int_elm_z = int(CK3int_elm_yz, y, -halfelmhgt, halfelmhgt); %对变量y进行积分
        CK3int_elm = elmwid * CK3int_elm_z; %对变量z进行积分
        CK3int(:,:,i,j) = eval(CK3int_elm);
    end
end
CK_elm = (lamda + 2*mu)/2*CK1int + lamda/2*CK2int + mu*CK3int; %计算单元维度下的单元刚度矩阵非线性部分对应的不变矩阵

end