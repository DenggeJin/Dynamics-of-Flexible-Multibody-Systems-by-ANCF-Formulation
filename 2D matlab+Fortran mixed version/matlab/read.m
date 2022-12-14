function aa=read(elmlen_1,elmhgt_1,elmnum_1,elmnum_2)
%% open files
fprintf('open files\n');


out = fopen('output.txt', 'r');


%% read coordinates
fprintf('read coordinates\n');
tline = fgetl(out); 
A=str2num(tline);
DOF_SYS=A(1);
TIMENUM=A(2);
i=1;
while(i<TIMENUM+1)
    tline = fgetl(out);                   
    q_sys(:,i)=str2num(tline);
    i=i+1;
end
fclose(out);
% xynum=DOF_SYS/3;
% xy=zeros(xynum,TIMENUM);
% x=zeros(xynum/2,TIMENUM);
% y=zeros(xynum/2,TIMENUM);
% for i=1:TIMENUM
%     for j=1:2:xynum-1
%         xy(j,i)=q(j+4*(j-1)/2,i);
%         x((j+1)/2,i)=xy(j,i);
%         xy(j+1,i)=q(j+1+4*(j-1)/2,i);
%         y((j+1)/2,i)=xy(j+1,i);
%     end
% end


%%
disp([datestr(now,'HH:MM:SS.FFF  '),'开始后处理................']);
xi_1=0:0.003:elmlen_1;
ita_1=-elmhgt_1:0.003:elmhgt_1;
[x_1,y_1]=meshgrid(xi_1,ita_1);
elmlen=elmlen_1;
s1_1 = 1 - 3*(x_1/elmlen).^2 + 2*(x_1/elmlen).^3;
s2_1 = elmlen*(x_1/elmlen - 2*(x_1/elmlen).^2 + (x_1/elmlen).^3);
s3_1 = elmlen*(y_1/elmlen - (x_1.*y_1)/(elmlen^2));
s4_1 = 3*(x_1/elmlen).^2 - 2*(x_1/elmlen).^3;
s5_1 = elmlen*(-(x_1/elmlen).^2 + (x_1/elmlen).^3);
s6_1 = elmlen*((x_1.*y_1)/(elmlen^2));
lenth_2=size(xi_1,2)*size(ita_1,2);
X_1=0;
Y_1=0;

for nn=1:size(q_sys,2)
    mm=0;
for i=1:elmnum_1
    for j=1:size(s1_1,1)
        for k=1:size(s1_1,2) 
            mm=mm+1;
            X_temp=q_sys(6*(i-1)+1,nn)*s1_1(j,k)+q_sys(6*(i-1)+3,nn)*s2_1(j,k)+q_sys(6*(i-1)+5,nn)*s3_1(j,k)+q_sys(6*(i)+1,nn)*s4_1(j,k)+q_sys(6*(i)+3,nn)*s5_1(j,k)+q_sys(6*(i)+5,nn)*s6_1(j,k);
            Y_temp=q_sys(6*(i-1)+2,nn)*s1_1(j,k)+q_sys(6*(i-1)+4,nn)*s2_1(j,k)+q_sys(6*(i-1)+6,nn)*s3_1(j,k)+q_sys(6*(i)+2,nn)*s4_1(j,k)+q_sys(6*(i)+4,nn)*s5_1(j,k)+q_sys(6*(i)+6,nn)*s6_1(j,k);           
            X_1(mm,nn)=X_temp;
            Y_1(mm,nn)=Y_temp;
        end
    end
end
end
mmm=mm;



for nn=1:size(q_sys,2)
    mm=mmm;
for ii=elmnum_1+1:elmnum_1+elmnum_2
    i=ii+1;%第三个和第四个节点是铰点，中间没有单元
    for j=1:size(s1_1,1)
        for k=1:size(s1_1,2) 
             mm=mm+1;
            X_temp=q_sys(6*(i-1)+1,nn)*s1_1(j,k)+q_sys(6*(i-1)+3,nn)*s2_1(j,k)+q_sys(6*(i-1)+5,nn)*s3_1(j,k)+q_sys(6*(i)+1,nn)*s4_1(j,k)+q_sys(6*(i)+3,nn)*s5_1(j,k)+q_sys(6*(i)+5,nn)*s6_1(j,k);
            Y_temp=q_sys(6*(i-1)+2,nn)*s1_1(j,k)+q_sys(6*(i-1)+4,nn)*s2_1(j,k)+q_sys(6*(i-1)+6,nn)*s3_1(j,k)+q_sys(6*(i)+2,nn)*s4_1(j,k)+q_sys(6*(i)+4,nn)*s5_1(j,k)+q_sys(6*(i)+6,nn)*s6_1(j,k);           
            X_1(mm,nn)=X_temp;
            Y_1(mm,nn)=Y_temp;
        end
    end
end
end
disp([datestr(now,'HH:MM:SS.FFF  '),'后处理完成']);

%%
disp([datestr(now,'HH:MM:SS.FFF  '),'开始绘图..................']);
for i=1:size(q_sys,2)%共有8个单元，为了给不同的单元颜色区分，此处没有做通用化
plot(X_1(1:lenth_2,i),Y_1(1:lenth_2,i),'OM');axis([-2,2.5,-2.5,1]);hold on
plot(X_1(lenth_2+1:2*lenth_2,i),Y_1(lenth_2+1:2*lenth_2,i),'OC');axis([-2,2.5,-2.5,1]);hold on
plot(X_1(2*lenth_2+1:3*lenth_2,i),Y_1(2*lenth_2+1:3*lenth_2,i),'OY');axis([-2,2.5,-2.5,1]);hold on
plot(X_1(3*lenth_2+1:4*lenth_2,i),Y_1(3*lenth_2+1:4*lenth_2,i),'OM');axis([-2,2.5,-2.5,1]);hold on
plot(X_1(4*lenth_2+1:5*lenth_2,i),Y_1(4*lenth_2+1:5*lenth_2,i),'OC');axis([-2,2.5,-2.5,1]);hold on
plot(X_1(5*lenth_2+1:6*lenth_2,i),Y_1(5*lenth_2+1:6*lenth_2,i),'OY');axis([-2,2.5,-2.5,1]);hold on
plot(X_1(6*lenth_2+1:7*lenth_2,i),Y_1(6*lenth_2+1:7*lenth_2,i),'OM');axis([-2,2.5,-2.5,1]);hold on
plot(X_1(7*lenth_2+1:8*lenth_2,i),Y_1(7*lenth_2+1:8*lenth_2,i),'OC');axis([-2,2.5,-2.5,1]);hold off
pause(0.1)
end

aa=1;
end

% %还需要用插值形函数插值
% 
% for i=1:TIMENUM
%     plot(x(:,i),y(:,i))
%     axis ([0 10 0 10])
%     pause(0.2)
% end
