
    
    
    
    
    
  
PROGRAM MAIN

USE BASIC_DATA

IMPLICIT NONE
 
INTEGER I

!-----------------------��ȡmatlabǰ��������Լ���ʼ��-----------------------------
CALL  READ_INFO
PRINT*,'READ_INFO�����'
PRINT*,''

!-------------------------ʱ���ƽ�����------------------------------------------
!*********��ʼ��ʱ
PRINT*,'���濪ʼ����ʱ���ƽ�����'
WRITE(*,"('������ ', I0, '����ÿһ���ƽ�ʱ��Ϊ:',F6.3,'S')")TIMENUM,TIMESTEP
PRINT*,''
CALL SYSTEM_CLOCK(START_CLK, COUNT_RATE)   !���㿪ʼʱ��


!*********��ʼ����
DO I=2,TIMENUM
    
    T=I
    
    !�����Է�����ĳ�ֵ
    Q_PRIDICT=Q(:,I-1)+TIMESTEP*DQ(:,I-1)+TIMESTEP*TIMESTEP/2*DDQ(:,I-1)
    QLAMDA(1:DOF_SYS)=Q_PRIDICT
    QLAMDA(DOF_SYS+1:DOF)=LAMDA(:,I-1)
    
    !����NEWTONSOLVER���
    CALL NEWTONSOLVER(QLAMDA)
    
    !������һʱ�䲽��q,lamda,dq,ddq
    Q(:,I)=X(1:DOF_SYS)
    LAMDA(:,I)=X(DOF_SYS+1:DOF)
    DDQ(:,I)=1/(BETA*TIMESTEP*TIMESTEP)*(Q(:,I)-Q(:,I-1))-1/(BETA*TIMESTEP)*DQ(:,I-1)-(1/(2*BETA)-1)*DDQ(:,I-1)
    DQ(:,I)=DQ(:,I-1)+(1-GAMMA)*TIMESTEP*DDQ(:,I-1)+GAMMA*TIMESTEP*DDQ(:,I)
    
    PRINT*,''
    WRITE(*,"('��', I0, '���������')"),T
    
ENDDO


!*********������ʱ
CALL SYSTEM_CLOCK(END_CLK, COUNT_RATE)   !����ʱ��
TIME_DIF = (END_CLK - START_CLK) / COUNT_RATE
PRINT*,''
WRITE(*,"('������ʱ��Ϊ   : ', I0, ' hour ', I0, ' minutes ', I0, ' seconds.')")&
TIME_DIF / 3600, MOD(TIME_DIF, 3600) / 60, MOD(TIME_DIF, 60)
PRINT*,''


!----------------------�ͷ��ڴ棬������������Ϊoutput.txt�ļ�---------------------------------
CALL SUMMARY


END PROGRAM MAIN