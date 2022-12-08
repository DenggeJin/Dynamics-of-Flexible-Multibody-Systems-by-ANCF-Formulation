

INCLUDE 'mkl_pardiso.f90'

SUBROUTINE SPARSE_SOLVER(XX,b,x)
!-----INPUT:XX
!-----OUTPUT:b,x

USE mkl_pardiso
USE BASIC_DATA,ONLY:DOF

IMPLICIT NONE

INTEGER, PARAMETER :: dp = KIND(1.0D0)
!Internal solver memory pointer 
TYPE(MKL_PARDISO_HANDLE), ALLOCATABLE  :: pt(:)
!All other variables
INTEGER maxfct, mnum, mtype, phase, n, nnz, nrhs, error, msglvl
INTEGER error1
INTEGER i, idum(1)
REAL(KIND=DP) ddum(1)
INTEGER, ALLOCATABLE :: iparm( : )

REAL(KIND=DP)::XX(DOF)
REAL(KIND=DP)::JACOBI(DOF,DOF)

REAL(KIND=DP), ALLOCATABLE :: a( : )
INTEGER, ALLOCATABLE :: ia( : )
INTEGER, ALLOCATABLE :: ja( : )
REAL(KIND=DP):: b( DOF )
REAL(KIND=DP):: x( DOF )


!--------------------Fill all arrays containing matrix data.---------------------------------
n=DOF

CALL JAC(XX,JACOBI)               !调用SUBROUTINE:JAC,计算在XX处的雅可比矩阵JACOBI
    
CALL GET_NZNUM(JACOBI,nnz)        !调用SUBROUTINE:GET_NZNUM,得到JACOBI矩阵中非零元素个数nnz
    
ALLOCATE(a(nnz))
ALLOCATE(ja(nnz))
ALLOCATE(ia(n+1)) 

CALL CSR3(JACOBI,nnz,a,ia,ja)     !调用SUBROUTINE:CSR3,利用压缩存储技术
                                  !将JACOBI矩阵进行一维变带宽存储,得到a,ia,ja三个一维数组，以便调用MKL库

!-----------------------Set up PARDISO control parameter--------------------
maxfct = 1 
mnum = 1
mtype  = 11!非对称
nrhs = 1 
error  = 0 ! initialize error flag
msglvl = 0 ! not print statistical information


ALLOCATE(iparm(64))
DO i = 1, 64
   iparm(i) = 0
END DO
iparm(1) = 1 ! no solver default
iparm(2) = 2 ! fill-in reordering from METIS
iparm(4) = 0 ! no iterative-direct algorithm
iparm(5) = 0 ! no user fill-in reducing permutation
iparm(6) = 0 ! =0 solution on the first n components of x
iparm(8) = 2 ! numbers of iterative refinement steps
iparm(10) = 13 ! perturb the pivot elements with 1E-13
iparm(11) = 1 ! use nonsymmetric permutation and scaling MPS
iparm(13) = 0 ! maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm(13) = 1 in case of inappropriate accuracy
iparm(14) = 0 ! Output: number of perturbed pivots
iparm(18) = -1 ! Output: number of nonzeros in the factor LU
iparm(19) = -1 ! Output: Mflops for LU factorization
iparm(20) = 0 ! Output: Numbers of CG Iterations

!--------------------------------------------------------------------------------------
!Initialize the internal solver memory pointer. This is only
! necessary for the FIRST call of the PARDISO solver.

ALLOCATE (pt(64))
DO i = 1, 64
   pt(i)%DUMMY =  0 
END DO


!-------------------------------------------------------------------------------------------
!Reordering and Symbolic Factorization, This step also allocates
! all memory that is necessary for the factorization
phase = 11 ! only reordering and symbolic factorization
CALL pardiso (pt, maxfct, mnum, mtype, phase, n, a, ia, ja, &
              idum, nrhs, iparm, msglvl, ddum, ddum, error)
    
WRITE(*,*) 'Reordering completed ... '
IF (error /= 0) THEN
   WRITE(*,*) 'The following ERROR was detected: ', error
   GOTO 1000
END IF
!WRITE(*,*) 'Number of nonzeros in factors = ',iparm(18)
!WRITE(*,*) 'Number of factorization MFLOPS = ',iparm(19)

!----------------------------------------------------------------------------
!Factorization.
phase = 22 ! only factorization
CALL pardiso (pt, maxfct, mnum, mtype, phase, n, a, ia, ja, &
              idum, nrhs, iparm, msglvl, ddum, ddum, error)
WRITE(*,*) 'Factorization completed ... '
IF (error /= 0) THEN
   WRITE(*,*) 'The following ERROR was detected: ', error
   GOTO 1000
ENDIF

!---------------------------------------------------------------------------
!Back substitution and iterative refinement
iparm(8) = 2 ! max numbers of iterative refinement steps
phase = 33 ! only solving
CALL pardiso (pt, maxfct, mnum, mtype, phase, n, a, ia, ja, &
              idum, nrhs, iparm, msglvl, b, x, error)
WRITE(*,*) 'Solve completed ... '
IF (error /= 0) THEN
   WRITE(*,*) 'The following ERROR was detected: ', error
   GOTO 1000
ENDIF
WRITE(*,*) 'The solution of the system is '
DO i = 1, n
   WRITE(*,*) ' x(',i,') = ', x(i)
END DO
      
1000 CONTINUE
!Termination and release of memory
phase = -1 ! release internal memory
CALL pardiso (pt, maxfct, mnum, mtype, phase, n, ddum, idum, idum, &
              idum, nrhs, iparm, msglvl, ddum, ddum, error1)

IF (ALLOCATED(ia))      DEALLOCATE(ia)
IF (ALLOCATED(ja))      DEALLOCATE(ja)
IF (ALLOCATED(a))       DEALLOCATE(a)
IF (ALLOCATED(b))       DEALLOCATE(b)
IF (ALLOCATED(x))       DEALLOCATE(x)
IF (ALLOCATED(iparm))   DEALLOCATE(iparm)

IF (error1 /= 0) THEN
   WRITE(*,*) 'The following ERROR on release stage was detected: ', error1
   STOP 1
ENDIF

IF (error /= 0) STOP 1
END SUBROUTINE SPARSE_SOLVER
    
    