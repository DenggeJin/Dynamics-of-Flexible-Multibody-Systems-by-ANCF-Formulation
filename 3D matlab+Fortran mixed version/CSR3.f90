
    
    SUBROUTINE CSR3(JACOBI,NZNUM,A,IA,JA)
    !---INPUT:
    !         JACOBI：雅可比矩阵
    !---OUTPUT:
    !          A:VALUE
    !          IA:ROWINDEX
    !          JA:COLUMS      此处参见MKL函数库官方说明文档：p3352
    

    USE BASIC_DATA,ONLY:REAL_KIND,INT_KIND,DOF          !NZNUM：JAC中非零元素个数
    IMPLICIT NONE
    
    INTEGER::NZNUM
    DOUBLE PRECISION::JACOBI(DOF,DOF)
    
    DOUBLE PRECISION::A(NZNUM)
    INTEGER(INT_KIND)::JA(NZNUM)
    INTEGER(INT_KIND)::IA(DOF+1)
    
    INTEGER(INT_KIND)::I
    INTEGER(INT_KIND)::J
    INTEGER(INT_KIND)::COUNT
    
    IA(1)=1
    COUNT=0
    DO I=1,DOF
        DO J=1,DOF
            IF(JACOBI(I,J).NE.0.0)THEN
                COUNT=COUNT+1
                A(COUNT)=JACOBI(I,J)
                JA(COUNT)=J
            END IF
        END DO
        IA(I+1)=COUNT+1
    END DO
    
    END SUBROUTINE CSR3
    