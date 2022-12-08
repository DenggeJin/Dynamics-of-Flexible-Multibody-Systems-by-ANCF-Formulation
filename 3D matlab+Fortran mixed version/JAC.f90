

SUBROUTINE JAC(XX,JACOBI)
    
!--------INPUT:XX
!--------OUTPUT:JACOBI

    USE BASIC_DATA,ONLY:REAL_KIND,INT_KIND ,DOF,DOF_SYS,LAMDANUM,M,CQ,BETA,TIMESTEP
    
    IMPLICIT NONE
    
    !----------------------------------�˲���������������---------------------------
    REAL(REAL_KIND)::TEMP1(DOF,DOF)
    REAL(REAL_KIND)::TEMP2(DOF,DOF)
    REAL(REAL_KIND)  COE
    REAL(REAL_KIND)::EYE1(DOF_SYS,DOF_SYS)    !��Ҫ��Ϊ���ȶ��嵥λ����
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
    
    !-------------------------------����ϵͳ���̵��ſɱȾ���----------------------------
    !�����������Լ�������ſɱȾ���
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
    
    !���㵯�������ſɱȾ���
    CALL F_JAC(Q,FJAC)   !ע�⣺�˴�����F_JAC�������㵯�����Թ�������q��������lamda�����ſɱȾ���FJAC
    JACOBI_2(1:DOF_SYS,1:DOF_SYS)=FJAC
    JACOBI_2(1:DOF_SYS,DOF_SYS+1:DOF)=0.0
    JACOBI_2(DOF_SYS+1:DOF,1:DOF)=0.0
    
    !�����ܵ��ſɱȾ���
    JACOBI=JACOBI_1+JACOBI_2
    
END SUBROUTINE JAC
    


SUBROUTINE F_JAC(Q,FJAC)

!------------���������ڼ��㵯�������ֵ��ſɱȾ���------------
!-------INPUT:Q����ǰ�Ĺ�������������Ԫ�ظ���ΪDOF_SYS
!-------OUTPUT:FJAC,���������ſɱȾ��󲿷�


    USE BASIC_DATA,ONLY:REAL_KIND,INT_KIND,DOF_SYS,POINTDOF, &
                        & LK_1,LK_2,CK_1,CK_2,ELENUM1,ELENUM2,BEQ
    
    IMPLICIT NONE
    
    REAL(REAL_KIND)::FJAC(DOF_SYS,DOF_SYS)                             !ϵͳ�������ſɱȾ���
    REAL(REAL_KIND)::FELM_JAC(2*POINTDOF,2*POINTDOF,ELENUM1+ELENUM2)   !��Ԫ�������ſɱ�����
    
    DOUBLE PRECISION::Q(DOF_SYS,1)     !��ǰ��ϵͳ��������
    REAL(REAL_KIND)::E(2*POINTDOF,1)   !���ڵ㵥Ԫ��������
    INTEGER(INT_KIND) I                !ѭ��ָ��
    INTEGER(INT_KIND) J                !ѭ��ָ��
    INTEGER(INT_KIND) K                !ѭ��ָ��
    INTEGER(INT_KIND) MM               !ѭ��ָ��
    
    REAL(REAL_KIND)::JAC1=0.0
    REAL(REAL_KIND)::JAC2=0.0    
    REAL(REAL_KIND)::JAC3=0.0    

    REAL(REAL_KIND)::TEMP(1,2*POINTDOF)!����5���������Ǽ�������е��м��������������������ͳһ������
    REAL(REAL_KIND)::TEMPX(2*POINTDOF,1)
    REAL(REAL_KIND)::TEMP1(1,1)
    REAL(REAL_KIND)::TEMP2(1,1)
    REAL(REAL_KIND)::TEMP3(1,1)
    
    FJAC = 0.0
    FELM_JAC = 0.0

    !�Ե�һ��������
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
    
    !�Եڶ���������
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
    
    !�����������ĵ�Ԫ�������ſɱȾ����鼯����
    DO MM=1,ELENUM1+ELENUM2
        FJAC=FJAC+MATMUL(MATMUL(TRANSPOSE(BEQ(:,:,MM)),FELM_JAC(:,:,MM)),BEQ(:,:,MM))
    END DO
    
END SUBROUTINE F_JAC 