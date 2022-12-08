function Svalue = ShapeFunc(xi,eta,elmlen)


x = xi*elmlen;
y = eta*elmlen;
z = zeta*elmlen;
s1 = 1 - 3*(x/elmlen)^2 + 2*(x/elmlen)^3;
s2 = elmlen*(x/elmlen - 2*(x/elmlen)^2 + (x/elmlen)^3);
s3 = elmlen*(y/elmlen - (x*y)/(elmlen^2));
s4 = elmlen*(z/elmlen - (x*z)/(elmlen^2));
s5 = 3*(x/elmlen)^2 - 2*(x/elmlen)^3;
s6 = elmlen*(-(x/elmlen)^2 + (x/elmlen)^3);
s7 = elmlen*((x*y)/(elmlen^2));
s8 = elmlen*((x*z)/(elmlen^2));
I = eye(3);
Svalue  = [s1*I, s2*I, s3*I, s4*I, s5*I, s6*I, s7*I, s8*I];%生成单元的形函数矩阵

end