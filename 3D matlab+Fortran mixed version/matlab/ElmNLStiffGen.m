function NLK_elm = ElmNLStiffGen(CK_elm,e_elm)


NLK_elm = zeros(24,24); %初始化单元维度下的单元刚度矩阵非线性部分
for i = 1:24
    for j = 1:24
        NLK_elm(i,j) = e_elm'*CK_elm(:,:,i,j)*e_elm; %计算单元维度下的单元刚度矩阵非线性部分
    end
end

end