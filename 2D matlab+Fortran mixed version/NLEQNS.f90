!***************************************************************************************************
!- PURPORSE:
!  DEFINES BASIC DATA TYPES
!
!- CONSTANTS
! 
!
!- PROGAMMED BY:
!  
!***************************************************************************************************
    
SUBROUTINE NLEQNS(QLAMDA)

USE BASIC_DATA,ONLY:INT_KIND,REAL_KIND,Q,DQ,DDQ,T,BETA,TIMESTEP,DOF,DOF_SYS,LAMDANUM,M,CQ,QG,LK,NLK,RES
    
IMPLICIT NONE

REAL(REAL_KIND)    QLAMDA(DOF)

REAL(REAL_KIND)    Q_1(DOF_SYS)
REAL(REAL_KIND)    DDQ_1(DOF_SYS)
REAL(REAL_KIND)    L_1(LAMDANUM) 

REAL(REAL_KIND)    Q_0(DOF_SYS)    !��һʱ�䲽������
REAL(REAL_KIND)    DQ_0(DOF_SYS)   !
REAL(REAL_KIND)    DDQ_0(DOF_SYS)  !

REAL(REAL_KIND)    TEMP_1(DOF_SYS,1)    
REAL(REAL_KIND)    TEMP_2(DOF_SYS,1)   !
REAL(REAL_KIND)    TEMP_3(DOF_SYS,1)  !
REAL(REAL_KIND)    TEMP_4(DOF_SYS,1)  !

Q_1= QLAMDA(1:DOF_SYS)
L_1= QLAMDA(DOF_SYS+1:DOF)

Q_0=Q(1:DOF_SYS,T-1)!Q��T��BASIC_DATA���
DQ_0=DQ(1:DOF_SYS,T-1)
DDQ_0=DDQ(1:DOF_SYS,T-1)

DDQ_1=1/(BETA*TIMESTEP*TIMESTEP)*(Q_1-Q_0) - 1/(BETA*TIMESTEP)*DQ_0 - (1/(2*BETA)-1)*DDQ_0

CALL NLSTIFF(Q_1)    !����NLSTIFF�������ɷ����ԸնȾ���

TEMP_1(:,1)=DDQ_1
TEMP_2(:,1)=L_1
TEMP_3(:,1)=Q_1

TEMP_4=MATMUL(M,TEMP_1)+ MATMUL(TRANSPOSE(CQ),TEMP_2)+MATMUL((LK+NLK),TEMP_3)
RES(1:DOF_SYS)=TEMP_4(:,1)-QG(:,1)

RES(DOF_SYS+1:DOF) = MATMUL(CQ,DDQ_1)
    
END SUBROUTINE NLEQNS