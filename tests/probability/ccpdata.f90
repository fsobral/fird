module ccpdata

  implicit none

  ! PARAMETERS

  ! Number of observations in order to build the correlation matrix
  integer, parameter :: NOBS = 5
  ! Sparsity rate
  real(8), parameter :: SPRATE = 1.0D-01
  ! Penalization term for the chance constraint
  real(8), parameter :: PEN = 1.0D+04
  
  ! COMMON ARRAYS

  integer, allocatable :: cA(:), pA(:)
  real(8), allocatable :: A(:), b(:), c(:), corr(:), mu(:)

  ! COMMON SCALARS

  integer :: m, np
  real(8) :: plim

contains

  subroutine initialize(n, np_, mi, x, l, u, plim_)

    ! This subroutine randomly initializes the structure needed by the
    ! Chance Constrained Problem.

    implicit none

    ! SCALAR ARGUMENTS
    integer :: n, np_, mi
    real(8) :: plim_

    ! ARRAY ARGUMENTS
    real(8) :: l(n), u(n), x(n)

    intent(in ) :: n, np_, mi, plim_
    intent(out) :: l, u, x

    ! INTERFACES

    interface

       SUBROUTINE UNBIASED_CORRELATION( N, A, B, msg, lag, r, t, m)

         INTEGER, INTENT(IN)                :: N
         REAL(8), DIMENSION(N), INTENT(IN)  :: A 
         REAL(8), DIMENSION(N), INTENT(IN)  :: B           
         REAL(8), INTENT(IN)                :: msg 
         INTEGER, INTENT(IN)                :: lag 
         REAL(8), INTENT(OUT)               :: r
         REAL(8), INTENT(OUT)               :: t
         INTEGER, INTENT(OUT)               :: m

       END SUBROUTINE UNBIASED_CORRELATION

    end interface


    ! LOCAL SCALARS

    integer :: i, j, k, nover, nseed, pos
    real(8) :: pi, randnum, t, u1, u2

    ! LOCAL ARRAYS

    integer, allocatable :: seed(:)
    real(8), target, allocatable :: mObs(:,:)

    ! POINTERS

    real(8), pointer :: tmpA(:), tmpB(:)

    m  = mi

    np = np_

    plim = plim_

    ! Initialize the random structure

    call random_seed(SIZE=nseed)

    allocate(seed(nseed))

    do i = 1,nseed
       seed(i) = 123456789
    end do

    call random_seed(PUT=seed)

    ! Initialize module objects

    allocate(b(m), c(n), corr((np - 1) * np / 2), mu(np), &
             mObs(np, np), A(n * m), pA(m + 1), cA(n * m))

    ! Randomly generates vector 'c'
    
    do i = 1, np

       call random_number(randnum)

       c(i) = 1.0D0 + 1.0D+01 * randnum

    end do

    do i = np + 1, n

       call random_number(randnum)

       c(i) = - 1.0D+01 + 2.0D+01 * randnum

    end do

    ! Randomly generates sparse matrix 'A' in vector form

    pA(1) = 1

    do i = 1, m

       pos = pA(i)

       do j = 1, n

          call random_number(randnum)

          if ( randnum .le. SPRATE ) then

             call random_number(A(pos))

             if ( j .le. np ) A(pos) = - abs(A(pos))

             cA(pos) = j

             pos = pos + 1

          end if

       end do

       pA(i + 1) = pos

    end do

    ! Randomly generates mean 'mu'

    do i = 1, np

       call random_number(randnum)

       mu(i) = 5.0D0 * randnum

    end do

    ! Randomly generates the initial guess, lower and upper bounds

    do i = 1, np

       x(i) = 0.0D0

       l(i) = mu(i)

       u(i) = 1.0D+20

    end do

    do i = np + 1, n

       call random_number(randnum)

       x(i) = - 1.0D+01 + 2.0D+01 * randnum

       l(i) = - 1.0D+01

       u(i) =   1.0D+01

    end do

    ! Generate feasible right hand side 'b'

    do i = 1, m

       b(i) = 0.0D0

       do j = pA(i), pA(i + 1) - 1

          b(i) = b(i) + x(cA(j)) * A(j)

       end do

    end do

    ! Randomly generates correlation matrix 'corr', which is simmetric
    ! and has 1's in the diagonal (in vector form)

    call genRandHouseDP(np,mObs,0.0D0,1.0D+01)

    do j = 1, np

       do i = 1, j - 1

          corr(i + (j - 1) * (j - 2) / 2) = mObs(i,j)

          write(*,*) i, j, corr(i + (j - 1) * (j - 2) / 2)

       end do

    end do

!!$    pi  = acos(-1.0D0)
!!$
!!$    do j = 1, np
!!$
!!$       do i = 1, NOBS, 2
!!$
!!$          call random_number(u1)
!!$          
!!$          call random_number(u2)
!!$
!!$          mObs(i,j) = sqrt(- 2.0D0 * log(u1)) * cos(2.0D0 * pi * u2)
!!$          
!!$!          if ( j .gt. 1 ) mObs(i,j) = mObs(i,1) * mObs(i,j)
!!$
!!$          if ( i + 1 .le. NOBS ) then
!!$                
!!$             mObs(i + 1,j) = sqrt(- 2.0D0 * log(u1)) * sin(2.0D0 * pi * u2)
!!$
!!$          end if
!!$
!!$       end do
!!$
!!$    end do
!!$
!!$    do j = 1, np
!!$
!!$       tmpA => mObs(:,j)
!!$
!!$       do i = 1, j - 1
!!$
!!$          tmpB => mObs(:,i)
!!$
!!$          call unbiased_correlation(NOBS, tmpA, tmpB, 1.0D0, 0, &
!!$               corr(i + (j - 1) * (j - 2) / 2), t, nover)
!!$
!!$          !write(*,*) (tmpA(k), k = 1,NOBS)
!!$          write(*,*) i, j, corr(i + (j - 1) * (j - 2) / 2)
!!$
!!$       end do
!!$
!!$    end do

    deallocate(mObs, seed)
    
  end subroutine initialize

  ! ---------------------------------------------------------- !
  ! ---------------------------------------------------------- !
  
  subroutine genRandHouseDP(n,A,prob,magnitude)
    ! This subroutine generates a random symmetric matrix A with
    ! dimension n.
    !
    ! prob: is the probability of a eigenvalue to be very near to zero
    !
    ! magnitude: if is positive, defines the range of random
    ! values. If is negative, also defines this range, but the matrix
    ! has chance to be indefinite


    ! PARAMETERS
    real(8),parameter :: epsdiag = 1.0D-5

    ! ARRAY ARGUMENTS
    real(8) :: A(n,n)

    ! SCALAR ARGUMENTS
    integer :: n
    real(8) :: prob,magnitude

    intent( in) :: prob,n,magnitude
    intent(out) :: A

    ! LOCAL SCALARS
    integer :: i,j,k
    real(8) :: rnumber

    ! LOCAL ARRAYS
    real(8) :: Q(n,n),d(n)

    call genHouseholderOrt(n,Q,magnitude)

    ! Generates a random diagonal matrix D and calculates Q * D * Q^t

    do i = 1,n
       CALL RANDOM_NUMBER(rnumber)
       if (rnumber .lt. prob) then
          ! Semi-definite positiveness probability
          CALL RANDOM_NUMBER(rnumber)
          d(i) = (1.0D0 + rnumber) * epsdiag
       else if (magnitude .gt. 0.0D0) then
          ! Matrix is positive definite
          CALL RANDOM_NUMBER(rnumber)
          d(i) = epsdiag + magnitude * rnumber
       else
          ! Matrix is indefinite
          CALL RANDOM_NUMBER(rnumber)
          d(i) = - abs(magnitude) + 2.0D0 * abs(magnitude) * rnumber          
       end if
    end do

    do j = 1,n
       do i = 1,n
          A(i,j) = 0.0D0
       end do
    end do

    do j = 1,n
       do i = 1,n
          do k = 1,n
             A(k,j) = A(k,j) + Q(k,i) * Q(j,i) * d(i)
          end do
       end do
    end do

  end subroutine genRandHouseDP

  ! ---------------------------------------------------------- !
  ! ---------------------------------------------------------- !
  
  subroutine genHouseholderOrt(n,Q,magnitude)
    ! This subroutine generates an randomly generated orthogonal
    ! Householder matrix:
    !
    ! I - 2 * (u * u^T) / (u^T * u)
    !
    ! and returns such matrix in Q

    ! SCALAR ARGUMENTS
    integer :: n
    real(8) :: magnitude

    ! ARRAY ARGUMENTS
    real(8) :: Q(n,n)

    intent( in) :: n,magnitude
    intent(out) :: Q

    ! LOCAL SCALARS
    integer :: i,j
    real(8) :: rnumber,dp

    ! LOCAL ARRAYS
    real(8) :: u(n)

    dp = 0.0D0
    do i = 1,n
       CALL RANDOM_NUMBER(rnumber)
       u(i) = abs(magnitude) + 2.0D0 * abs(magnitude) * rnumber
       dp = dp + u(i) ** 2.0D0
    end do

    do j = 1,n
       do i = 1,n
          Q(i,j) = - 2.0D0 * u(i) * u(j) / dp
       end do
       Q(j,j) = Q(j,j) + 1.0D0
    end do

  end subroutine genHouseholderOrt

  ! ---------------------------------------------------------- !
  ! ---------------------------------------------------------- !
  
  subroutine destroy()

    deallocate(A, b, c, corr, mu, pA, cA)

  end subroutine destroy

  ! ---------------------------------------------------------- !
  ! ---------------------------------------------------------- !

  subroutine tofile(n, np, mi, x, l, u, filename)

    implicit none

    ! SCALAR ARGUMENTS
    integer :: n, np, mi

    ! ARRAY ARGUMENTS
    real(8)       :: l(n), u(n), x(n)
    character(80) :: filename

    intent(in ) :: filename, n, np, mi, l, u, x

    ! LOCAL SCALARS
    integer :: i, j

    ! LOCAL ARRAYS
    real(8) :: tmpA(n)

    open(99, FILE=filename)

    write(99, FMT=010) n, np, mi

    do i = 1, mi

       do j = 1, n

          tmpA(j) = 0.0D0

       end do

       do j = pA(i), pA(i + 1) - 1

          tmpA(cA(j)) = A(j)

       end do

       do j = 1, n
       
          write(99,FMT=020) tmpA(j)

       end do

       write(99,*)

    end do

    write(99,*)

    close(99)

    ! NON-EXECUTABLE STATEMENTS

010 FORMAT(I5,1X,I5,1X,I5,/)
020 FORMAT(1P,E15.8)

  end subroutine tofile

end module ccpdata
