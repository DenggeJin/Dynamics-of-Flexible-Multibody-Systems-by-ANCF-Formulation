!***************************************************************************************************
!- PURPORSE:
!  定义：（1）模型基本变量；（2）求解过程中间变量；（3）计算结果变量
!
!- CONSTANTS：
!  NONE
!
!- PROGAMMED BY:
!  
!***************************************************************************************************


MODULE BASIC_DATA

IMPLICIT NONE

!*****************************PART1:   基本求解控制参数***************************
!*****************************此部分变量用户可自定义********************************

!--------------------------定义变量类型------------------------------
INTEGER(4), PARAMETER::        REAL_KIND    = SELECTED_REAL_KIND(P=15)
INTEGER(4), PARAMETER::        INT_KIND     = SELECTED_INT_KIND(8)
INTEGER(4), PARAMETER::        WORD_KIND    = 32  


!---------------------------定义所求解问题的基本变量变量--------------------------------
INTEGER(INT_KIND):: ELENUM1 = 2     !部件1的单元数量
INTEGER(INT_KIND):: ELENUM2 = 6     !部件2的单元数量
INTEGER(INT_KIND):: POINTDOF = 6    !节点自由度
INTEGER(INT_KIND):: LAMDANUM = 4    !约束方程个数
INTEGER(INT_KIND):: POINTNUM1       !部件1的节点数
INTEGER(INT_KIND):: POINTNUM2       !部件2的节点数
INTEGER(INT_KIND):: DOF_SYS         !系统自由度数
INTEGER(INT_KIND):: DOF             !总自由度数

!---------------------------定义时间推进所需变量--------------------------------
REAL(REAL_KIND):: TIMESTEP=0.001    !时间推进步长
INTEGER(INT_KIND):: TIMENUM = 1000  !时间推进步数
REAL(REAL_KIND):: GAMMA=0.5         !newmark方法的参数
REAL(REAL_KIND):: BETA=0.25         !newmark方法的参数




!****************************PART2:  求解过程变量***********************************

!定义记录程序计算时间所需变量
INTEGER(INT_KIND) START_CLK     ! THE START CLOCK OF THE PROGRAM
INTEGER(INT_KIND) END_CLK       ! THE END CLOCK OF THE PROGRAM
INTEGER(INT_KIND) TIME_DIF        ! THE TIME SPAN
INTEGER(INT_KIND) COUNT_RATE    ! THE NUMBER OF CLICKS PER SECOND

INTEGER(INT_KIND)::T=0          !时间推进步数

DOUBLE PRECISION, ALLOCATABLE:: X(:)!迭代求解时的中间变量
REAL(REAL_KIND),ALLOCATABLE::RES(:) !迭代求解时的中间变量

REAL(REAL_KIND), ALLOCATABLE:: M(:,:)!系统定常质量矩阵
INTEGER(INT_KIND) DOF_M1             !系统质量矩阵第一维度长度
INTEGER(INT_KIND) DOF_M2             !系统质量矩阵第二维度长度
CHARACTER(WORD_KIND):: FILENAME_M    !存储系统质量矩阵的txt文件名

REAL(REAL_KIND), ALLOCATABLE:: LAMDA0(:,:)!拉格朗日乘子初始值
INTEGER(INT_KIND) DOF_LAMDA01             !拉格朗日乘子初始值第一维度长度
INTEGER(INT_KIND) DOF_LAMDA02             !拉格朗日乘子初始值第二维度长度
CHARACTER(WORD_KIND):: FILENAME_LAMDA0    !存储拉格朗日乘子初始值的txt文件名

REAL(REAL_KIND), ALLOCATABLE:: Q0(:,:)    !广义节点坐标初始值
INTEGER(INT_KIND) DOF_Q01
INTEGER(INT_KIND) DOF_Q02
CHARACTER(WORD_KIND):: FILENAME_Q0

REAL(REAL_KIND), ALLOCATABLE:: DQ0(:,:)    !广义节点速度初始值
INTEGER(INT_KIND) DOF_DQ01
INTEGER(INT_KIND) DOF_DQ02
CHARACTER(WORD_KIND):: FILENAME_DQ0

REAL(REAL_KIND), ALLOCATABLE:: DDQ0(:,:)   !广义节点加速度初始值
INTEGER(INT_KIND) DOF_DDQ01
INTEGER(INT_KIND) DOF_DDQ02
CHARACTER(WORD_KIND):: FILENAME_DDQ0

REAL(REAL_KIND), ALLOCATABLE:: CQ(:,:)     !约束方程矩阵
INTEGER(INT_KIND) DOF_CQ1
INTEGER(INT_KIND) DOF_CQ2
CHARACTER(WORD_KIND):: FILENAME_CQ

REAL(REAL_KIND), ALLOCATABLE:: LK(:,:)     !系统线性刚度矩阵
INTEGER(INT_KIND) DOF_LK1
INTEGER(INT_KIND) DOF_LK2
CHARACTER(WORD_KIND):: FILENAME_LK

REAL(REAL_KIND), ALLOCATABLE:: QG(:,:)     !外加载荷矩阵（向量）
INTEGER(INT_KIND) DOF_QG1
INTEGER(INT_KIND) DOF_QG2
CHARACTER(WORD_KIND):: FILENAME_QG

REAL(REAL_KIND), ALLOCATABLE:: LK_1(:,:)     !部件1单元线性刚度矩阵
INTEGER(INT_KIND) DOF_LK_11
INTEGER(INT_KIND) DOF_LK_12
CHARACTER(WORD_KIND):: FILENAME_LK_1

REAL(REAL_KIND), ALLOCATABLE:: LK_2(:,:)     !部件2单元线性刚度矩阵
INTEGER(INT_KIND) DOF_LK_21
INTEGER(INT_KIND) DOF_LK_22
CHARACTER(WORD_KIND):: FILENAME_LK_2

REAL(REAL_KIND),ALLOCATABLE:: ELENLK(:,:)    !单元刚度矩阵非线性部分
REAL(REAL_KIND), ALLOCATABLE:: NLK(:,:)      !非线性总体刚度矩阵
!NLK矩阵的DOF与LK矩阵一致，因此不再重复定义

REAL(REAL_KIND), ALLOCATABLE:: BEQ(:,:,:)    !系统坐标变换矩阵
INTEGER(INT_KIND) DOF_BEQ1
INTEGER(INT_KIND) DOF_BEQ2
INTEGER(INT_KIND) DOF_BEQ3
CHARACTER(WORD_KIND):: FILENAME_BEQ

REAL(REAL_KIND), ALLOCATABLE:: CK_1(:,:,:,:)!部件1的单元维度下的单元刚度矩阵非线性部分的不变矩阵
INTEGER(INT_KIND) DOF_CK_11
INTEGER(INT_KIND) DOF_CK_12
INTEGER(INT_KIND) DOF_CK_13
INTEGER(INT_KIND) DOF_CK_14
CHARACTER(WORD_KIND):: FILENAME_CK_1
REAL(REAL_KIND), ALLOCATABLE:: CK_2(:,:,:,:)!部件2的单元维度下的单元刚度矩阵非线性部分的不变矩阵
INTEGER(INT_KIND) DOF_CK_21
INTEGER(INT_KIND) DOF_CK_22
INTEGER(INT_KIND) DOF_CK_23
INTEGER(INT_KIND) DOF_CK_24
CHARACTER(WORD_KIND):: FILENAME_CK_2

REAL(REAL_KIND),ALLOCATABLE:: Q(:,:)        !广义节点坐标矩阵
REAL(REAL_KIND),ALLOCATABLE:: DQ(:,:)       !广义节点速度矩阵
REAL(REAL_KIND),ALLOCATABLE:: DDQ(:,:)      !广义节点加速度矩阵
REAL(REAL_KIND),ALLOCATABLE:: LAMDA(:,:)    !拉格朗日乘子矩阵
REAL(REAL_KIND),ALLOCATABLE:: Q_PRIDICT(:)  !由前一时刻的位移、速度、加速度预测当前时刻的位移
REAL(REAL_KIND),ALLOCATABLE:: QLAMDA(:)     !当前时刻的位移和拉格朗日乘子


END MODULE