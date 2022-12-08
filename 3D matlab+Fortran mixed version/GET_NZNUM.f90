
    SUBROUTINE GET_NZNUM(JACOBI,NZNUM)
    !本程序用于得到雅可比矩阵中非零元素的个数NZNUM
    !--------INPUT:JACOBI,DOF
    !--------OUTPUT:NZNUM

    USE BASIC_DATA,ONLY:DOF
    IMPLICIT NONE
    
    INTEGER::I
    INTEGER::J
    INTEGER::NZNUM
    DOUBLE PRECISION::JACOBI(DOF,DOF)
    
    NZNUM=0
    
    DO I=1,DOF
        DO J=1,DOF
            IF(JACOBI(I,J).NE.0.0) THEN
                NZNUM=NZNUM+1
            END IF
        END DO
    END DO
    
    END SUBROUTINE GET_NZNUM           
    