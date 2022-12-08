!***************************************************************************************************
!- PURPORSE:
!  CALCUTE ELEMENT NONLINER STIFFNESS MATRIX
!  
!- INPUT ARGUMENTS:
!   CK,E_ELM 
!
!- OUTPUT ARGUMENTS:
!  ELENLK
!
!- CALL PROCEDURES:
!  NONE
!
!- CALLED BY:
!  NLSTIFF
!
!- PROGAMMED BY:
!  
!***************************************************************************************************
    
SUBROUTINE ELENLSTIFF(CK,E_ELM)

USE BASIC_DATA

REAL(REAL_KIND) CK(2*POINTDOF,2*POINTDOF,2*POINTDOF,2*POINTDOF)                                 !单元刚度矩阵非线性部分的不变矩阵
REAL(REAL_KIND) E_ELM( 2*POINTDOF , 1)
REAL(REAL_KIND) TEMP(1,1)
INTEGER(INT_KIND) I
INTEGER(INT_KIND) J 

DO I=1,2*POINTDOF
    DO J=1,2*POINTDOF
        TEMP=MATMUL(MATMUL(TRANSPOSE(E_ELM) ,CK(:,:,I,J)),E_ELM)!此处要注意一个1*1的矩阵和一个数不一样
        ELENLK(I,J)=TEMP(1,1)                                   !单元刚度矩阵非线性部分
    ENDDO
ENDDO

END SUBROUTINE ELENLSTIFF