module hsccpdata

  implicit none

  ! PARAMETERS

  ! Sparsity rate
  real(8), parameter :: SPRATE = 1.0D-01
  ! Penalization term for the chance constraint
  real(8), parameter :: PEN = 1.0D+04

  ! Relative error for subroutine MVNDST
  real(8), parameter :: RELERR = 1.0D-05
  ! Absolute error for subroutine MVNDST
  real(8), parameter :: ABSERR = 1.0D-08
  
  ! COMMON ARRAYS

  integer, allocatable :: cA(:), pA(:)
  real(8), allocatable :: A(:), c(:), corr(:), mu(:)

  ! COMMON SCALARS

  integer :: m, np
  real(8) :: plim

  ! PRIVATE SCALARS

  integer :: initial_seed

  private :: initial_seed

contains

  ! ---------------------------------------------------------- !
  ! ---------------------------------------------------------- !
  
  subroutine set_seed(is)

    ! This subroutine sets the initial seed. It is usefull for massive
    ! random tests.

    ! SCALAR ARGUMENTS
    integer :: is

    initial_seed = is

  end subroutine set_seed

  ! ---------------------------------------------------------- !
  ! ---------------------------------------------------------- !
  
  subroutine initialize(np_, mi, x, l, u, plim_)

    ! This subroutine randomly initializes the structure needed by the
    ! Chance Constrained Problem.

    implicit none

    ! SCALAR ARGUMENTS
    integer :: np_, mi
    real(8) :: plim_

    ! ARRAY ARGUMENTS
    real(8) :: l(np_), u(np_), x(np_)

    intent(in ) :: np_, mi, plim_
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
       seed(i) = initial_seed
    end do

    call random_seed(PUT=seed)

    ! Initialize module objects

    allocate(c(np), corr((np - 1) * np / 2), mu(np), &
             mObs(np, np), A(np * m), pA(m + 1), cA(np * m))

    ! Randomly generates vector 'c'
    
    do i = 1, np

       call random_number(randnum)

       c(i) = 1.0D0 + 1.0D+01 * randnum

    end do

    ! Randomly generates sparse matrix 'A' in vector form

    pA(1) = 1

    do i = 1, m

       pos = pA(i)

       do j = 1, np

          call random_number(randnum)

          if ( randnum .le. SPRATE ) then

             call random_number(A(pos))

             A(pos) = - abs(A(pos))

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

    ! Randomly generates correlation matrix 'corr', which is
    ! simmetric, definite positive and has 1's in the diagonal (in
    ! vector form)

    call genRandHouseDP(np,mObs,0.0D0,1.0D+02)

    do j = 1, np

       do i = 1, j - 1

          pos = i + (j - 1) * (j - 2) / 2

          corr(pos) = mObs(i,j) / (sqrt(mObs(i,i)) * sqrt(mObs(j,j)))

          ! write(*,*) i, j, corr(pos)

       end do

    end do

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

    implicit none

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

    implicit none

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
       u(i) = - abs(magnitude) + 2.0D0 * abs(magnitude) * rnumber
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

    deallocate(A, c, corr, mu, pA, cA)

  end subroutine destroy

end module hsccpdata
