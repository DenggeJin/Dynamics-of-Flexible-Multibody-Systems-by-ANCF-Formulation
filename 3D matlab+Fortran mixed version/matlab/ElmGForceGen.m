function Qg_elm = ElmGForceGen(rho,grav,elmlen,elmwid,elmhgt)
%--------------------------------------------------------------------------
% Copyright (c) 2014, Zhao Rui
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
%--------------------------------------------------------------------------
% This version was developed and run on Matlab 8.0.0.783 (R2012b)
% Date: March 20,2014
% Author: Zhao Rui   ( Harbin Institute of Technology, China )
% For problems with the code, please contact the author:
% Email: stevezhao1987@outlook.com
%        stevezhao1987@gmail.com
%        ruizhao@hit.edu.cn
%--------------------------------------------------------------------------

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


Qgint_elm_xyz = S(2,:)'; 
Qgint_elm_yz = int(Qgint_elm_xyz, x, 0, elmlen); 
Qgint_elm_z = int(Qgint_elm_yz, y, -halfelmhgt, halfelmhgt); 
Qgint_elm =  int(Qgint_elm_z, z, -halfelmwid, halfelmwid);
Qg_elm = -rho * grav * eval(Qgint_elm); %计算单位维度下的单元广义外力向量

end