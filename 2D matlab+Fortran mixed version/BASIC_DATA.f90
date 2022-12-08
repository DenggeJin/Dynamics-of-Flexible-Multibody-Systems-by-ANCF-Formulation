!***************************************************************************************************
!- PURPORSE:
!  ���壺��1��ģ�ͻ�����������2���������м��������3������������
!
!- CONSTANTS��
!  NONE
!
!- PROGAMMED BY:
!  
!***************************************************************************************************


MODULE BASIC_DATA

IMPLICIT NONE

!*****************************PART1:   ���������Ʋ���***************************
!*****************************�˲��ֱ����û����Զ���********************************

!--------------------------�����������------------------------------
INTEGER(4), PARAMETER::        REAL_KIND    = SELECTED_REAL_KIND(P=15)
INTEGER(4), PARAMETER::        INT_KIND     = SELECTED_INT_KIND(8)
INTEGER(4), PARAMETER::        WORD_KIND    = 32  


!---------------------------�������������Ļ�����������--------------------------------
INTEGER(INT_KIND):: ELENUM1 = 2     !����1�ĵ�Ԫ����
INTEGER(INT_KIND):: ELENUM2 = 6     !����2�ĵ�Ԫ����
INTEGER(INT_KIND):: POINTDOF = 6    !�ڵ����ɶ�
INTEGER(INT_KIND):: LAMDANUM = 4    !Լ�����̸���
INTEGER(INT_KIND):: POINTNUM1       !����1�Ľڵ���
INTEGER(INT_KIND):: POINTNUM2       !����2�Ľڵ���
INTEGER(INT_KIND):: DOF_SYS         !ϵͳ���ɶ���
INTEGER(INT_KIND):: DOF             !�����ɶ���

!---------------------------����ʱ���ƽ��������--------------------------------
REAL(REAL_KIND):: TIMESTEP=0.001    !ʱ���ƽ�����
INTEGER(INT_KIND):: TIMENUM = 1000  !ʱ���ƽ�����
REAL(REAL_KIND):: GAMMA=0.5         !newmark�����Ĳ���
REAL(REAL_KIND):: BETA=0.25         !newmark�����Ĳ���




!****************************PART2:  �����̱���***********************************

!�����¼�������ʱ���������
INTEGER(INT_KIND) START_CLK     ! THE START CLOCK OF THE PROGRAM
INTEGER(INT_KIND) END_CLK       ! THE END CLOCK OF THE PROGRAM
INTEGER(INT_KIND) TIME_DIF        ! THE TIME SPAN
INTEGER(INT_KIND) COUNT_RATE    ! THE NUMBER OF CLICKS PER SECOND

INTEGER(INT_KIND)::T=0          !ʱ���ƽ�����

DOUBLE PRECISION, ALLOCATABLE:: X(:)!�������ʱ���м����
REAL(REAL_KIND),ALLOCATABLE::RES(:) !�������ʱ���м����

REAL(REAL_KIND), ALLOCATABLE:: M(:,:)!ϵͳ������������
INTEGER(INT_KIND) DOF_M1             !ϵͳ���������һά�ȳ���
INTEGER(INT_KIND) DOF_M2             !ϵͳ��������ڶ�ά�ȳ���
CHARACTER(WORD_KIND):: FILENAME_M    !�洢ϵͳ���������txt�ļ���

REAL(REAL_KIND), ALLOCATABLE:: LAMDA0(:,:)!�������ճ��ӳ�ʼֵ
INTEGER(INT_KIND) DOF_LAMDA01             !�������ճ��ӳ�ʼֵ��һά�ȳ���
INTEGER(INT_KIND) DOF_LAMDA02             !�������ճ��ӳ�ʼֵ�ڶ�ά�ȳ���
CHARACTER(WORD_KIND):: FILENAME_LAMDA0    !�洢�������ճ��ӳ�ʼֵ��txt�ļ���

REAL(REAL_KIND), ALLOCATABLE:: Q0(:,:)    !����ڵ������ʼֵ
INTEGER(INT_KIND) DOF_Q01
INTEGER(INT_KIND) DOF_Q02
CHARACTER(WORD_KIND):: FILENAME_Q0

REAL(REAL_KIND), ALLOCATABLE:: DQ0(:,:)    !����ڵ��ٶȳ�ʼֵ
INTEGER(INT_KIND) DOF_DQ01
INTEGER(INT_KIND) DOF_DQ02
CHARACTER(WORD_KIND):: FILENAME_DQ0

REAL(REAL_KIND), ALLOCATABLE:: DDQ0(:,:)   !����ڵ���ٶȳ�ʼֵ
INTEGER(INT_KIND) DOF_DDQ01
INTEGER(INT_KIND) DOF_DDQ02
CHARACTER(WORD_KIND):: FILENAME_DDQ0

REAL(REAL_KIND), ALLOCATABLE:: CQ(:,:)     !Լ�����̾���
INTEGER(INT_KIND) DOF_CQ1
INTEGER(INT_KIND) DOF_CQ2
CHARACTER(WORD_KIND):: FILENAME_CQ

REAL(REAL_KIND), ALLOCATABLE:: LK(:,:)     !ϵͳ���ԸնȾ���
INTEGER(INT_KIND) DOF_LK1
INTEGER(INT_KIND) DOF_LK2
CHARACTER(WORD_KIND):: FILENAME_LK

REAL(REAL_KIND), ALLOCATABLE:: QG(:,:)     !����غɾ���������
INTEGER(INT_KIND) DOF_QG1
INTEGER(INT_KIND) DOF_QG2
CHARACTER(WORD_KIND):: FILENAME_QG

REAL(REAL_KIND), ALLOCATABLE:: LK_1(:,:)     !����1��Ԫ���ԸնȾ���
INTEGER(INT_KIND) DOF_LK_11
INTEGER(INT_KIND) DOF_LK_12
CHARACTER(WORD_KIND):: FILENAME_LK_1

REAL(REAL_KIND), ALLOCATABLE:: LK_2(:,:)     !����2��Ԫ���ԸնȾ���
INTEGER(INT_KIND) DOF_LK_21
INTEGER(INT_KIND) DOF_LK_22
CHARACTER(WORD_KIND):: FILENAME_LK_2

REAL(REAL_KIND),ALLOCATABLE:: ELENLK(:,:)    !��Ԫ�նȾ�������Բ���
REAL(REAL_KIND), ALLOCATABLE:: NLK(:,:)      !����������նȾ���
!NLK�����DOF��LK����һ�£���˲����ظ�����

REAL(REAL_KIND), ALLOCATABLE:: BEQ(:,:,:)    !ϵͳ����任����
INTEGER(INT_KIND) DOF_BEQ1
INTEGER(INT_KIND) DOF_BEQ2
INTEGER(INT_KIND) DOF_BEQ3
CHARACTER(WORD_KIND):: FILENAME_BEQ

REAL(REAL_KIND), ALLOCATABLE:: CK_1(:,:,:,:)!����1�ĵ�Ԫά���µĵ�Ԫ�նȾ�������Բ��ֵĲ������
INTEGER(INT_KIND) DOF_CK_11
INTEGER(INT_KIND) DOF_CK_12
INTEGER(INT_KIND) DOF_CK_13
INTEGER(INT_KIND) DOF_CK_14
CHARACTER(WORD_KIND):: FILENAME_CK_1
REAL(REAL_KIND), ALLOCATABLE:: CK_2(:,:,:,:)!����2�ĵ�Ԫά���µĵ�Ԫ�նȾ�������Բ��ֵĲ������
INTEGER(INT_KIND) DOF_CK_21
INTEGER(INT_KIND) DOF_CK_22
INTEGER(INT_KIND) DOF_CK_23
INTEGER(INT_KIND) DOF_CK_24
CHARACTER(WORD_KIND):: FILENAME_CK_2

REAL(REAL_KIND),ALLOCATABLE:: Q(:,:)        !����ڵ��������
REAL(REAL_KIND),ALLOCATABLE:: DQ(:,:)       !����ڵ��ٶȾ���
REAL(REAL_KIND),ALLOCATABLE:: DDQ(:,:)      !����ڵ���ٶȾ���
REAL(REAL_KIND),ALLOCATABLE:: LAMDA(:,:)    !�������ճ��Ӿ���
REAL(REAL_KIND),ALLOCATABLE:: Q_PRIDICT(:)  !��ǰһʱ�̵�λ�ơ��ٶȡ����ٶ�Ԥ�⵱ǰʱ�̵�λ��
REAL(REAL_KIND),ALLOCATABLE:: QLAMDA(:)     !��ǰʱ�̵�λ�ƺ��������ճ���


END MODULE