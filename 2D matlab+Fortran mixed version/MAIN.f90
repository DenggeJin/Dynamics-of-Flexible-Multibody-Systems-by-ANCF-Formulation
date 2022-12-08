
    
    
    
    
    
  
PROGRAM MAIN

USE BASIC_DATA

IMPLICIT NONE
 
INTEGER I

!-----------------------读取matlab前处理参数以及初始化-----------------------------
CALL  READ_INFO
PRINT*,'READ_INFO已完成'
PRINT*,''

!-------------------------时间推进计算------------------------------------------
!*********开始计时
PRINT*,'下面开始进行时间推进计算'
WRITE(*,"('共计算 ', I0, '步，每一步推进时间为:',F6.3,'S')")TIMENUM,TIMESTEP
PRINT*,''
CALL SYSTEM_CLOCK(START_CLK, COUNT_RATE)   !计算开始时间


!*********开始计算
DO I=2,TIMENUM
    
    T=I
    
    !非线性方程组的初值
    Q_PRIDICT=Q(:,I-1)+TIMESTEP*DQ(:,I-1)+TIMESTEP*TIMESTEP/2*DDQ(:,I-1)
    QLAMDA(1:DOF_SYS)=Q_PRIDICT
    QLAMDA(DOF_SYS+1:DOF)=LAMDA(:,I-1)
    
    !调用NEWTONSOLVER求解
    CALL NEWTONSOLVER(QLAMDA)
    
    !计算这一时间步的q,lamda,dq,ddq
    Q(:,I)=X(1:DOF_SYS)
    LAMDA(:,I)=X(DOF_SYS+1:DOF)
    DDQ(:,I)=1/(BETA*TIMESTEP*TIMESTEP)*(Q(:,I)-Q(:,I-1))-1/(BETA*TIMESTEP)*DQ(:,I-1)-(1/(2*BETA)-1)*DDQ(:,I-1)
    DQ(:,I)=DQ(:,I-1)+(1-GAMMA)*TIMESTEP*DDQ(:,I-1)+GAMMA*TIMESTEP*DDQ(:,I)
    
    PRINT*,''
    WRITE(*,"('第', I0, '步计算完毕')"),T
    
ENDDO


!*********结束计时
CALL SYSTEM_CLOCK(END_CLK, COUNT_RATE)   !结束时间
TIME_DIF = (END_CLK - START_CLK) / COUNT_RATE
PRINT*,''
WRITE(*,"('计算总时长为   : ', I0, ' hour ', I0, ' minutes ', I0, ' seconds.')")&
TIME_DIF / 3600, MOD(TIME_DIF, 3600) / 60, MOD(TIME_DIF, 60)
PRINT*,''


!----------------------释放内存，输出结果并保存为output.txt文件---------------------------------
CALL SUMMARY


END PROGRAM MAIN