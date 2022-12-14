!***************************************************************************************************
!- PURPORSE:
!     读取matlab输出的质量阵、刚度阵、初始参数等，对BSIC_DATA里定义的变量进行初始化操作
!
!- CALL PROCEDURES:
!     NONE
!
!- CALLED BY:
!     MAIN.F90
!
!- PROGAMMED BY:
!  
!***************************************************************************************************
SUBROUTINE READ_INFO

USE BASIC_DATA

IMPLICIT NONE

POINTNUM1 = ELENUM1 + 1                         !部件1的节点数
POINTNUM2 = ELENUM2 + 1                         !部件2的节点数
DOF_SYS = (POINTNUM1 + POINTNUM2) * POINTDOF    !系统自由度数
DOF=DOF_SYS+LAMDANUM                            !总自由度数


!--------------------------------------为矩阵分配内存---------------------------------------
!------随时间推进变化的矩阵
ALLOCATE(X(DOF))
ALLOCATE(RES(DOF))

ALLOCATE(Q(DOF_SYS,TIMENUM))
ALLOCATE(DQ(DOF_SYS,TIMENUM))
ALLOCATE(DDQ(DOF_SYS,TIMENUM))
ALLOCATE(LAMDA(LAMDANUM,TIMENUM))
ALLOCATE(Q_PRIDICT(DOF_SYS))
ALLOCATE(QLAMDA(DOF))

!------定常矩阵
!二维矩阵（把向量看作第二维的维度是1）
DOF_M1=DOF_SYS
DOF_M2=DOF_SYS
DOF_LAMDA01=LAMDANUM
DOF_LAMDA02=1
DOF_Q01=DOF_SYS
DOF_Q02=1
DOF_DQ01=DOF_SYS
DOF_DQ02=1
DOF_DDQ01=DOF_SYS
DOF_DDQ02=1
DOF_CQ2=DOF_SYS
DOF_CQ1=LAMDANUM
DOF_LK1=DOF_SYS
DOF_LK2=DOF_SYS
DOF_QG1=DOF_SYS
DOF_QG2=1
DOF_LK_11=2*POINTDOF
DOF_LK_12=2*POINTDOF
DOF_LK_21=2*POINTDOF
DOF_LK_22=2*POINTDOF


!三维矩阵
DOF_BEQ1=2*POINTDOF
DOF_BEQ2=DOF_SYS
DOF_BEQ3=ELENUM1+ELENUM2

!四维矩阵
DOF_CK_11=2*POINTDOF
DOF_CK_12=2*POINTDOF
DOF_CK_13=2*POINTDOF
DOF_CK_14=2*POINTDOF
DOF_CK_21=2*POINTDOF
DOF_CK_22=2*POINTDOF
DOF_CK_23=2*POINTDOF
DOF_CK_24=2*POINTDOF

ALLOCATE(M(DOF_M1,DOF_M2))
FILENAME_M='M'
ALLOCATE(LAMDA0(DOF_LAMDA01,DOF_LAMDA02))
FILENAME_LAMDA0='lamda0'
ALLOCATE(Q0(DOF_Q01,DOF_Q02))
FILENAME_Q0='q0'
ALLOCATE(DQ0(DOF_DQ01,DOF_DQ02))
FILENAME_DQ0='dq0'
ALLOCATE(DDQ0(DOF_DDQ01,DOF_DDQ02))
FILENAME_DDQ0='ddq0'
ALLOCATE(CQ(DOF_CQ1,DOF_CQ2))
FILENAME_CQ='Cq'
ALLOCATE(LK(DOF_LK1,DOF_LK2))
FILENAME_LK='LK'
ALLOCATE(QG(DOF_QG1,DOF_QG2))
FILENAME_QG='Qg'
ALLOCATE(LK_1(DOF_LK_11,DOF_LK_12))
FILENAME_LK_1='LK_1'
ALLOCATE(LK_2(DOF_LK_21,DOF_LK_22))
FILENAME_LK_2='LK_2'

ALLOCATE(NLK(DOF_LK1,DOF_LK2))
ALLOCATE(ELENLK(12,12))

ALLOCATE(BEQ(DOF_BEQ1,DOF_BEQ2,DOF_BEQ3))
FILENAME_BEQ='Beq'

ALLOCATE(CK_1(DOF_CK_11,DOF_CK_12,DOF_CK_13,DOF_CK_14))
FILENAME_CK_1='CK_1'
ALLOCATE(CK_2(DOF_CK_21,DOF_CK_22,DOF_CK_23,DOF_CK_24))
FILENAME_CK_2='CK_2'


!----------------------------------------从txt文档中读入参数------------------------------------
M=TXT2MAT(DOF_M1,DOF_M2,FILENAME_M)                     !调用函数TXT2MAX函数读取系统定常质量矩阵
LAMDA0=TXT2MAT(DOF_LAMDA01,DOF_LAMDA02,FILENAME_LAMDA0) !调用函数TXT2MAX函数读取初始拉格朗日乘子
Q0=TXT2MAT(DOF_Q01,DOF_Q02,FILENAME_Q0)                 !调用函数TXT2MAX函数读取初始广义节点坐标
DQ0=TXT2MAT(DOF_DQ01,DOF_DQ02,FILENAME_DQ0)             !调用函数TXT2MAX函数读取初始广义节点速度
DDQ0=TXT2MAT(DOF_DDQ01,DOF_DDQ02,FILENAME_DDQ0)         !调用函数TXT2MAX函数读取初始广义节点加速度
CQ=TXT2MAT(DOF_CQ1,DOF_CQ2,FILENAME_CQ)                 !调用函数TXT2MAX函数读取约束矩阵
LK=TXT2MAT(DOF_LK1,DOF_LK2,FILENAME_LK)                 !调用函数TXT2MAX函数读取系统线性刚度矩阵
QG=TXT2MAT(DOF_QG1,DOF_QG2,FILENAME_QG)                 !调用函数TXT2MAX函数读取外加载荷矩阵（向量）
LK_1=TXT2MAT(DOF_LK_11,DOF_LK_12,FILENAME_LK_1)         !调用函数TXT2MAX函数读取部件1单元线性刚度矩阵
LK_2=TXT2MAT(DOF_LK_21,DOF_LK_22,FILENAME_LK_2)         !调用函数TXT2MAX函数读取部件2单元线性刚度矩阵
 
BEQ=TXT3MAT(DOF_BEQ1,DOF_BEQ2,DOF_BEQ3,FILENAME_BEQ)    !调用函数TXT3MAX函数读取系统坐标转换矩阵

CK_1=TXT4MAT(DOF_CK_11,DOF_CK_12,DOF_CK_13,DOF_CK_14,FILENAME_CK_1)!调用函数TXT4MAX函数读取部件1的单元维度下的单元刚度矩阵非线性部分的不变矩阵
CK_2=TXT4MAT(DOF_CK_21,DOF_CK_22,DOF_CK_23,DOF_CK_24,FILENAME_CK_2)!调用函数TXT4MAX函数读取部件2的单元维度下的单元刚度矩阵非线性部分的不变矩阵


!为Q,DQ,DDQ和LAMDA赋初值
Q(:,1)=Q0(:,1)
DQ(:,1)=DQ0(:,1)
DDQ(:,1)=DDQ0(:,1)
LAMDA(:,1)=LAMDA0(:,1)

NLK(:,:)=0.0


CONTAINS
 
!-------------------内置函数TXT2MAT，用于读取二维和一维矩阵----------------------------
FUNCTION TXT2MAT(DOF1,DOF2,FILENAME)

IMPLICIT NONE

INTEGER(4), PARAMETER::        INT_KIND     = SELECTED_INT_KIND(8)
INTEGER(4), PARAMETER::        REAL_KIND    = SELECTED_REAL_KIND(P=15)
INTEGER(4), PARAMETER::        WORD_KIND    = 32  

INTEGER(INT_KIND) DOF1
INTEGER(INT_KIND) DOF2
REAL(REAL_KIND) TXT2MAT(DOF1,DOF2)
INTEGER(INT_KIND) I
INTEGER(INT_KIND) J   
CHARACTER(WORD_KIND):: FILENAME 
CHARACTER(WORD_KIND):: FILE

FILE=TRIM(FILENAME)//'.txt'                    !把字符串连起来

OPEN(10,file=FILE,action='read')
DO I = 1, DOF2
    READ(10,*) (TXT2MAT(J,I), J = 1, DOF1)  !此处需要注意读取的行列，因为矩阵在fortran中是按列储存，为读取快速，matlab输出的是矩阵的转置
ENDDO

END FUNCTION

!-------------------内置函数TXT3MAT，用于读取三维矩阵----------------------------------
FUNCTION TXT3MAT(DOF1,DOF2,DOF3,FILENAME)

IMPLICIT NONE

INTEGER(4), PARAMETER::        INT_KIND     = SELECTED_INT_KIND(8)
INTEGER(4), PARAMETER::        REAL_KIND    = SELECTED_REAL_KIND(P=15)
INTEGER(4), PARAMETER::        WORD_KIND    = 32  

INTEGER(INT_KIND) DOF1
INTEGER(INT_KIND) DOF2
INTEGER(INT_KIND) DOF3
REAL(REAL_KIND) TXT3MAT(DOF1,DOF2,DOF3)
INTEGER(INT_KIND) I
INTEGER(INT_KIND) J   
INTEGER(INT_KIND) K   
CHARACTER(WORD_KIND):: FILENAME 
CHARACTER(WORD_KIND):: FILE

FILE=TRIM(FILENAME)//'.txt'

OPEN(10,file=FILE,action='read')
DO K = 1, DOF3
DO I = 1, DOF2
    READ(10,*) (TXT3MAT(J,I,K), J = 1, DOF1)!此处需要注意读取的行列，因为矩阵在fortran中是按列储存，为读取快速，matlab输出的是矩阵的转置
    !此处读取时的维度顺序需要和matlab输出时的顺序对应
ENDDO
ENDDO

END FUNCTION

!-------------------内置函数TXT4MAT，用于读取四维矩阵----------------------------------
FUNCTION TXT4MAT(DOF1,DOF2,DOF3,DOF4,FILENAME)

IMPLICIT NONE

INTEGER(4), PARAMETER::        INT_KIND     = SELECTED_INT_KIND(8)
INTEGER(4), PARAMETER::        REAL_KIND    = SELECTED_REAL_KIND(P=15)
INTEGER(4), PARAMETER::        WORD_KIND    = 32  

INTEGER(INT_KIND) DOF1
INTEGER(INT_KIND) DOF2
INTEGER(INT_KIND) DOF3
INTEGER(INT_KIND) DOF4
REAL(REAL_KIND) TXT4MAT(DOF1,DOF2,DOF3,DOF4)
INTEGER(INT_KIND) I
INTEGER(INT_KIND) J   
INTEGER(INT_KIND) K   
INTEGER(INT_KIND) L   
CHARACTER(WORD_KIND):: FILENAME 
CHARACTER(WORD_KIND):: FILE

FILE=TRIM(FILENAME)//'.txt'

OPEN(10,file=FILE,action='read')
DO L = 1, DOF3
DO K = 1, DOF4
DO I = 1, DOF2
    READ(10,*) (TXT4MAT(J,I,L,K), J = 1, DOF1)!此处需要注意读取的行列，因为矩阵在fortran中是按列储存，为读取快速，matlab输出的是矩阵的转置
    !此处读取时的维度顺序需要和matlab输出时的顺序对应
ENDDO
ENDDO
ENDDO

END FUNCTION

!end program Console2
END SUBROUTINE READ_INFO