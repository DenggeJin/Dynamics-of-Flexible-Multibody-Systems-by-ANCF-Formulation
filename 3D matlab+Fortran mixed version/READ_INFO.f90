!***************************************************************************************************
!- PURPORSE:
!     ��ȡmatlab����������󡢸ն��󡢳�ʼ�����ȣ���BSIC_DATA�ﶨ��ı������г�ʼ������
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

POINTNUM1 = ELENUM1 + 1                         !����1�Ľڵ���
POINTNUM2 = ELENUM2 + 1                         !����2�Ľڵ���
DOF_SYS = (POINTNUM1 + POINTNUM2) * POINTDOF    !ϵͳ���ɶ���
DOF=DOF_SYS+LAMDANUM                            !�����ɶ���


!--------------------------------------Ϊ��������ڴ�---------------------------------------
!------��ʱ���ƽ��仯�ľ���
ALLOCATE(X(DOF))
ALLOCATE(RES(DOF))

ALLOCATE(Q(DOF_SYS,TIMENUM))
ALLOCATE(DQ(DOF_SYS,TIMENUM))
ALLOCATE(DDQ(DOF_SYS,TIMENUM))
ALLOCATE(LAMDA(LAMDANUM,TIMENUM))
ALLOCATE(Q_PRIDICT(DOF_SYS))
ALLOCATE(QLAMDA(DOF))

!------��������
!��ά���󣨰����������ڶ�ά��ά����1��
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


!��ά����
DOF_BEQ1=2*POINTDOF
DOF_BEQ2=DOF_SYS
DOF_BEQ3=ELENUM1+ELENUM2

!��ά����
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
ALLOCATE(ELENLK(2*POINTDOF,2*POINTDOF))

ALLOCATE(BEQ(DOF_BEQ1,DOF_BEQ2,DOF_BEQ3))
FILENAME_BEQ='Beq'

ALLOCATE(CK_1(DOF_CK_11,DOF_CK_12,DOF_CK_13,DOF_CK_14))
FILENAME_CK_1='CK_1'
ALLOCATE(CK_2(DOF_CK_21,DOF_CK_22,DOF_CK_23,DOF_CK_24))
FILENAME_CK_2='CK_2'


!----------------------------------------��txt�ĵ��ж������------------------------------------
M=TXT2MAT(DOF_M1,DOF_M2,FILENAME_M)                     !���ú���TXT2MAX������ȡϵͳ������������

LAMDA0=TXT2MAT(DOF_LAMDA01,DOF_LAMDA02,FILENAME_LAMDA0) !���ú���TXT2MAX������ȡ��ʼ�������ճ���
Q0=TXT2MAT(DOF_Q01,DOF_Q02,FILENAME_Q0)                 !���ú���TXT2MAX������ȡ��ʼ����ڵ�����
DQ0=TXT2MAT(DOF_DQ01,DOF_DQ02,FILENAME_DQ0)             !���ú���TXT2MAX������ȡ��ʼ����ڵ��ٶ�
DDQ0=TXT2MAT(DOF_DDQ01,DOF_DDQ02,FILENAME_DDQ0)         !���ú���TXT2MAX������ȡ��ʼ����ڵ���ٶ�
CQ=TXT2MAT(DOF_CQ1,DOF_CQ2,FILENAME_CQ)                 !���ú���TXT2MAX������ȡԼ������
LK=TXT2MAT(DOF_LK1,DOF_LK2,FILENAME_LK)                 !���ú���TXT2MAX������ȡϵͳ���ԸնȾ���
QG=TXT2MAT(DOF_QG1,DOF_QG2,FILENAME_QG)                 !���ú���TXT2MAX������ȡ����غɾ���������
LK_1=TXT2MAT(DOF_LK_11,DOF_LK_12,FILENAME_LK_1)         !���ú���TXT2MAX������ȡ����1��Ԫ���ԸնȾ���
LK_2=TXT2MAT(DOF_LK_21,DOF_LK_22,FILENAME_LK_2)         !���ú���TXT2MAX������ȡ����2��Ԫ���ԸնȾ���
BEQ=TXT3MAT(DOF_BEQ1,DOF_BEQ2,DOF_BEQ3,FILENAME_BEQ)    !���ú���TXT3MAX������ȡϵͳ����ת������

CK_1=TXT4MAT(DOF_CK_11,DOF_CK_12,DOF_CK_13,DOF_CK_14,FILENAME_CK_1)!���ú���TXT4MAX������ȡ����1�ĵ�Ԫά���µĵ�Ԫ�նȾ�������Բ��ֵĲ������
CK_2=TXT4MAT(DOF_CK_21,DOF_CK_22,DOF_CK_23,DOF_CK_24,FILENAME_CK_2)!���ú���TXT4MAX������ȡ����2�ĵ�Ԫά���µĵ�Ԫ�նȾ�������Բ��ֵĲ������

!ΪQ,DQ,DDQ��LAMDA����ֵ
Q(:,1)=Q0(:,1)
DQ(:,1)=DQ0(:,1)
DDQ(:,1)=DDQ0(:,1)
LAMDA(:,1)=LAMDA0(:,1)

NLK(:,:)=0.0


CONTAINS
 
!-------------------���ú���TXT2MAT�����ڶ�ȡ��ά��һά����----------------------------
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

FILE=TRIM(FILENAME)//'.txt'                    !���ַ���������

OPEN(10,file=FILE,action='read')
DO I = 1, DOF2
    READ(10,*) (TXT2MAT(J,I), J = 1, DOF1)  !�˴���Ҫע���ȡ�����У���Ϊ������fortran���ǰ��д��棬Ϊ��ȡ���٣�matlab������Ǿ����ת��
ENDDO

END FUNCTION

!-------------------���ú���TXT3MAT�����ڶ�ȡ��ά����----------------------------------
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
    READ(10,*) (TXT3MAT(J,I,K), J = 1, DOF1)!�˴���Ҫע���ȡ�����У���Ϊ������fortran���ǰ��д��棬Ϊ��ȡ���٣�matlab������Ǿ����ת��
    !�˴���ȡʱ��ά��˳����Ҫ��matlab���ʱ��˳���Ӧ
ENDDO
ENDDO

END FUNCTION

!-------------------���ú���TXT4MAT�����ڶ�ȡ��ά����----------------------------------
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
    READ(10,*) (TXT4MAT(J,I,L,K), J = 1, DOF1)!�˴���Ҫע���ȡ�����У���Ϊ������fortran���ǰ��д��棬Ϊ��ȡ���٣�matlab������Ǿ����ת��
    !�˴���ȡʱ��ά��˳����Ҫ��matlab���ʱ��˳���Ӧ
ENDDO
ENDDO
ENDDO

END FUNCTION

END SUBROUTINE READ_INFO