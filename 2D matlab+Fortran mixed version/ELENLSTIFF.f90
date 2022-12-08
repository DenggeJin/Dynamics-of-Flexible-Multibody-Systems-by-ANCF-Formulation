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

REAL(REAL_KIND) CK(2*POINTDOF,2*POINTDOF,2*POINTDOF,2*POINTDOF)                                 !��Ԫ�նȾ�������Բ��ֵĲ������
REAL(REAL_KIND) E_ELM( 2*POINTDOF , 1)
REAL(REAL_KIND) TEMP(1,1)
INTEGER(INT_KIND) I
INTEGER(INT_KIND) J 

DO I=1,2*POINTDOF
    DO J=1,2*POINTDOF
        TEMP=MATMUL(MATMUL(TRANSPOSE(E_ELM) ,CK(:,:,I,J)),E_ELM)!�˴�Ҫע��һ��1*1�ľ����һ������һ��
        ELENLK(I,J)=TEMP(1,1)                                   !��Ԫ�նȾ�������Բ���
    ENDDO
ENDDO

END SUBROUTINE ELENLSTIFF