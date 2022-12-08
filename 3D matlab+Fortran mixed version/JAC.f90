

SUBROUTINE JAC(XX,JACOBI)
    
!--------INPUT:XX
!--------OUTPUT:JACOBI

    USE BASIC_DATA,ONLY:REAL_KIND,INT_KIND ,DOF,DOF_SYS,LAMDANUM,M,CQ,BETA,TIMESTEP
    
    IMPLICIT NONE
    
    !----------------------------------此部分用来变量声明---------------------------
    REAL(REAL_KIND)::TEMP1(DOF,DOF)
    REAL(REAL_KIND)::TEMP2(DOF,DOF)
    REAL(REAL_KIND)  COE
    REAL(REAL_KIND)::EYE1(DOF_SYS,DOF_SYS)    !需要人为事先定义单位矩阵
    REAL(REAL_KIND)::EYE2(LAMDANUM,LAMDANUM)  !EYE1:DOF_SYS    EYE2:LAMDANUM=4
    
    DOUBLE PRECISION::JACOBI(DOF,DOF)                
    REAL(REAL_KIND)::JACOBI_1(DOF,DOF)
    REAL(REAL_KIND)::JACOBI_2(DOF,DOF)
    REAL(REAL_KIND)::FJAC(DOF_SYS,DOF_SYS)
    
    DOUBLE PRECISION::XX(DOF)                     
    DOUBLE PRECISION::Q(DOF_SYS,1)
    INTEGER(INT_KIND) I
    
    EYE1=0.0
    EYE2=0.0
    DO I=1,DOF_SYS
        EYE1(I,I)=1
    END DO
    DO I=1,LAMDANUM
        EYE2(I,I)=1
    END DO
    Q(:,1)=XX(1:DOF_SYS)
    
    !-------------------------------计算系统方程的雅可比矩阵----------------------------
    !计算惯性力及约束力的雅可比矩阵
    COE=1/(BETA*TIMESTEP*TIMESTEP)
    TEMP2(1:DOF_SYS,1:DOF_SYS)=COE*EYE1
    TEMP2(1:DOF_SYS,DOF_SYS+1:DOF)=0.0
    TEMP2(DOF_SYS+1:DOF,1:DOF_SYS)=0.0
    TEMP2(DOF_SYS+1:DOF,DOF_SYS+1:DOF)=EYE2
    
    TEMP1(1:DOF_SYS,1:DOF_SYS)=M
    TEMP1(DOF_SYS+1:DOF,DOF_SYS+1:DOF)=0.0
    TEMP1(DOF_SYS+1:DOF,1:DOF_SYS)=CQ
    TEMP1(1:DOF_SYS,DOF_SYS+1:DOF)=TRANSPOSE(CQ)
    
    JACOBI_1=MATMUL(TEMP1,TEMP2)
    
    !计算弹性力的雅可比矩阵
    CALL F_JAC(Q,FJAC)   !注意：此处调用F_JAC函数计算弹性力对广义坐标q（不包括lamda）的雅可比矩阵FJAC
    JACOBI_2(1:DOF_SYS,1:DOF_SYS)=FJAC
    JACOBI_2(1:DOF_SYS,DOF_SYS+1:DOF)=0.0
    JACOBI_2(DOF_SYS+1:DOF,1:DOF)=0.0
    
    !计算总的雅可比矩阵
    JACOBI=JACOBI_1+JACOBI_2
    
END SUBROUTINE JAC
    


SUBROUTINE F_JAC(Q,FJAC)

!------------本函数用于计算弹性力部分的雅可比矩阵------------
!-------INPUT:Q，当前的广义坐标向量，元素个数为DOF_SYS
!-------OUTPUT:FJAC,弹性力的雅可比矩阵部分


    USE BASIC_DATA,ONLY:REAL_KIND,INT_KIND,DOF_SYS,POINTDOF, &
                        & LK_1,LK_2,CK_1,CK_2,ELENUM1,ELENUM2,BEQ
    
    IMPLICIT NONE
    
    REAL(REAL_KIND)::FJAC(DOF_SYS,DOF_SYS)                             !系统弹性力雅可比矩阵
    REAL(REAL_KIND)::FELM_JAC(2*POINTDOF,2*POINTDOF,ELENUM1+ELENUM2)   !单元弹性力雅可比数组
    
    DOUBLE PRECISION::Q(DOF_SYS,1)     !当前的系统广义坐标
    REAL(REAL_KIND)::E(2*POINTDOF,1)   !两节点单元广义坐标
    INTEGER(INT_KIND) I                !循环指标
    INTEGER(INT_KIND) J                !循环指标
    INTEGER(INT_KIND) K                !循环指标
    INTEGER(INT_KIND) MM               !循环指标
    
    REAL(REAL_KIND)::JAC1=0.0
    REAL(REAL_KIND)::JAC2=0.0    
    REAL(REAL_KIND)::JAC3=0.0    

    REAL(REAL_KIND)::TEMP(1,2*POINTDOF)!下面5个变量都是计算过程中的中间变量，仅仅起到数据类型统一的作用
    REAL(REAL_KIND)::TEMPX(2*POINTDOF,1)
    REAL(REAL_KIND)::TEMP1(1,1)
    REAL(REAL_KIND)::TEMP2(1,1)
    REAL(REAL_KIND)::TEMP3(1,1)
    
    FJAC = 0.0
    FELM_JAC = 0.0

    !对第一个部件：
    DO MM=1,ELENUM1
        DO I=1,2*POINTDOF
            DO K=1,2*POINTDOF
                
                TEMP1=LK_1(I,K)
                JAC1=TEMP1(1,1)
                
                E=MATMUL(BEQ(:,:,MM),Q)
                TEMP2=MATMUL(MATMUL(TRANSPOSE(E),CK_1(:,:,I,K)),E)
                JAC2=TEMP2(1,1)
                
                JAC3=0.0
                TEMP3=0.0
                DO J=1,2*POINTDOF
                    TEMP(1,:)=CK_1(K,:,I,J)
                    TEMPX(:,1)=CK_1(:,K,I,J)
                    TEMP3=TEMP3+(MATMUL(TEMP,E)+MATMUL(TRANSPOSE(E),TEMPX))*E(J,1)
                END DO
                JAC3=TEMP3(1,1)
                
                FELM_JAC(I,K,MM)=JAC1+JAC2+JAC3
                
            END DO
        END DO
    END DO
    
    !对第二个部件：
    DO MM=ELENUM1+1,ELENUM1+ELENUM2
        DO I=1,2*POINTDOF
            DO K=1,2*POINTDOF
                
                TEMP1=LK_2(I,K)
                JAC1=TEMP1(1,1)
                
                E=MATMUL(BEQ(:,:,MM),Q)
                TEMP2=MATMUL(MATMUL(TRANSPOSE(E),CK_2(:,:,I,K)),E)
                JAC2=TEMP2(1,1)
                
                JAC3=0.0
                TEMP3=0.0
                DO J=1,2*POINTDOF
                    TEMP(1,:)=CK_2(K,:,I,J)
                    TEMPX(:,1)=CK_2(:,K,I,J)
                    TEMP3=TEMP3+(MATMUL(TEMP,E)+MATMUL(TRANSPOSE(E),TEMPX))*E(J,1)
                END DO
                JAC3=TEMP3(1,1)
                
                FELM_JAC(I,K,MM)=JAC1+JAC2+JAC3   
            END DO
        END DO
    END DO    
    
    !将两个部件的单元弹性力雅可比矩阵组集起来
    DO MM=1,ELENUM1+ELENUM2
        FJAC=FJAC+MATMUL(MATMUL(TRANSPOSE(BEQ(:,:,MM)),FELM_JAC(:,:,MM)),BEQ(:,:,MM))
    END DO
    
END SUBROUTINE F_JAC 