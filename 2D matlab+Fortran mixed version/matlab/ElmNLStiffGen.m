function NLK_elm = ElmNLStiffGen(CK_elm,e_elm)


NLK_elm = zeros(12,12); %��ʼ����Ԫά���µĵ�Ԫ�նȾ�������Բ���
for i = 1:12
    for j = 1:12
        NLK_elm(i,j) = e_elm'*CK_elm(:,:,i,j)*e_elm; %���㵥Ԫά���µĵ�Ԫ�նȾ�������Բ���
    end
end

end