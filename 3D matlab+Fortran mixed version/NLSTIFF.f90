!***************************************************************************************************
!- PURPORSE:
!  CALCUTE SYSTEM NONLINER STIFFNESS MATRIX
!  
!- INPUT ARGUMENTS:
!   QN 
!
!- OUTPUT ARGUMENTS:
!  NLK
!
!- CALL PROCEDURES:
!  ELENLSTIFF
!
!- CALLED BY:
!  SOLVER_1
!
!- PROGAMMED BY:
!  
!***************************************************************************************************
    
    
SUBROUTINE NLSTIFF(QN)

USE BASIC_DATA,ONLY:INT_KIND,REAL_KIND,ELENUM1,ELENUM2,DOF_SYS,CK_1,CK_2,BEQ,ELENLK,NLK

REAL(REAL_KIND) QN(DOF_SYS)
INTEGER(INT_KIND) A

NLK=0.0

!-----------------------���㲿��1������նȾ���---------------------------------------------
DO A = 1,ELENUM1
    CALL ELENLSTIFF(CK_1,MATMUL(BEQ(:,:,A) , QN))                           !���㵥Ԫ�����ԸնȾ���
    NLK = NLK + MATMUL(MATMUL(TRANSPOSE(BEQ(:,:,A)) , ELENLK),BEQ(:,:,A))   !�鼯Ϊ����նȾ���
ENDDO

!-----------------------���㲿��2������նȾ���---------------------------------------------
DO A = ELENUM1+1,ELENUM1+ELENUM2
    CALL ELENLSTIFF(CK_2,MATMUL(BEQ(:,:,A), QN))                            !���㵥Ԫ�����ԸնȾ���
    NLK = NLK +MATMUL( MATMUL(TRANSPOSE(BEQ(:,:,A)) , ELENLK) ,BEQ(:,:,A))  !�鼯Ϊ����նȾ���
ENDDO

END SUBROUTINE NLSTIFF