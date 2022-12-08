!***************************************************************************************************
!- PURPORSE:
!     ����ţ�ٵ��������ÿһʱ���ƽ����еõ��ķ����Է�����
!
!- CALL PROCEDURES:
!     NLEQNS,���ӳ����������ɷ����Է����飬����õ��µķ��̲в�RES
!     SPARSE_SOLVER,���ӳ����������ɷ����Է��̵��ſɱȾ��󣬲����ţ�ٵ��������е����Է�����
!
!- CALLED BY:
!     MAIN.F90
! 
!- PROGAMMED BY:
!  
!***************************************************************************************************

SUBROUTINE NEWTONSOLVER(QLAMDA)

!----INPUT:QLAMDA����ǰʱ�̵ġ�ϵͳ����λ��+�������ճ���������
!----OUTPUT:X��X������BASIC_DATA�ʱ���ƽ�һ�����µġ�ϵͳ����λ��+�������ճ���������

USE BASIC_DATA,ONLY:INT_KIND,REAL_KIND,DOF,DOF_SYS,RES,X

IMPLICIT NONE

REAL(REAL_KIND) DET1
REAL(REAL_KIND) DET2
REAL(REAL_KIND) EPS1
REAL(REAL_KIND) EPS2
REAL(REAL_KIND) QQQ(DOF_SYS)
REAL(REAL_KIND) QLAMDA(DOF)
REAL(REAL_KIND) QLAMDATHEN(DOF)
REAL(REAL_KIND) JACOBI(DOF,DOF)

double precision    XX(DOF)
double precision    RESS(DOF)
double precision    INVJACOBIRES(DOF)

DET1=1.0
DET2=1.0
EPS1=1.D-5
EPS2=1.D-5


DO WHILE(ABS(DET2)>EPS2)
    
    XX=QLAMDA
    
    CALL NLEQNS(QLAMDA)!�õ�RES
    
    RESS=RES

    CALL SPARSE_SOLVER(XX,RESS,INVJACOBIRES)
    
    QLAMDATHEN=QLAMDA-INVJACOBIRES
    
    DET1=MAXVAL(QLAMDATHEN-QLAMDA)
    
    DET2=MAXVAL(RES)
    
    QLAMDA=QLAMDATHEN
    
END DO

X=QLAMDATHEN

END SUBROUTINE NEWTONSOLVER