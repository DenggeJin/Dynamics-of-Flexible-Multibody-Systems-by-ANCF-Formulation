!***************************************************************************************************
!- PURPORSE:
!     采用牛顿迭代法求解每一时间推进步中得到的非线性方程组
!
!- CALL PROCEDURES:
!     NLEQNS,该子程序用于生成非线性方程组，并求得当下的方程残差RES
!     SPARSE_SOLVER,该子程序用于生成非线性方程的雅可比矩阵，并求解牛顿迭代过程中的线性方程组
!
!- CALLED BY:
!     MAIN.F90
! 
!- PROGAMMED BY:
!  
!***************************************************************************************************

SUBROUTINE NEWTONSOLVER(QLAMDA)

!----INPUT:QLAMDA，当前时刻的“系统广义位移+拉格朗日乘子向量”
!----OUTPUT:X，X定义在BASIC_DATA里，时间推进一步后新的“系统广义位移+拉格朗日乘子向量”

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
    
    CALL NLEQNS(QLAMDA)!得到RES
    
    RESS=RES

    CALL SPARSE_SOLVER(XX,RESS,INVJACOBIRES)
    
    QLAMDATHEN=QLAMDA-INVJACOBIRES
    
    DET1=MAXVAL(QLAMDATHEN-QLAMDA)
    
    DET2=MAXVAL(RES)
    
    QLAMDA=QLAMDATHEN
    
END DO

X=QLAMDATHEN

END SUBROUTINE NEWTONSOLVER