
SUBROUTINE SUMMARY

USE BASIC_DATA

IMPLICIT NONE
INTEGER J

!---------------------------输出计算得到的广义节点坐标到output.txt文件------------------------
OPEN(10,file='output.txt',action='write')
 WRITE(10,'( 2(2X, I7))')  DOF_SYS, TIMENUM
!  WRITE(10,*) '*** DOF_SYS ***'
!  WRITE(10,*) 60
!  WRITE(10,*) '*** TIMENUM ***'
!  WRITE(10,*) 20
!  WRITE(10,*) '*** COORDINATE ***'
    DO J = 1, TIMENUM        
            WRITE(10,'( 60(2X, ES13.6))')  Q(:,J)
    ENDDO
CLOSE(10)


!----------------------------释放已经分配的内存--------------------------------
DEALLOCATE(M)!
DEALLOCATE(LAMDA0)!
DEALLOCATE(Q0)!
DEALLOCATE(DQ0)!
DEALLOCATE(DDQ0)!
DEALLOCATE(CQ)!
DEALLOCATE(LK)!
DEALLOCATE(QG)!

DEALLOCATE(BEQ)!
DEALLOCATE(CK_1)!
DEALLOCATE(CK_2)!
DEALLOCATE(LK_1)!
DEALLOCATE(LK_2)!

DEALLOCATE (Q)!
DEALLOCATE(DQ)!
DEALLOCATE(DDQ)!
DEALLOCATE(LAMDA)!
DEALLOCATE(Q_PRIDICT)!
DEALLOCATE(QLAMDA)!
DEALLOCATE(X)!

DEALLOCATE(NLK)
DEALLOCATE(ELENLK)

END SUBROUTINE SUMMARY