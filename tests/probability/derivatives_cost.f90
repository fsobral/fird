program derivatives_cost

  integer :: np, flag
  real(8) :: start, finish, p, ABSERR, RELERR

  real(8), allocatable :: MU(:), CORR(:), x(:)

  write(*,*) "This program tests the cost of derivatives"
  write(*,*) "of CCP functions"

  np = 10

  allocate(x(np))

  call CPU_TIME(start)

  call evalprob(np, x, MU, CORR, ABSERR, RELERR, p, flag)

  call CPU_TIME(finish)

contains

  function fgauss(x, m, s)

    implicit none

    ! PARAMETERS
    real(8), parameter :: PI = asin(1.0)

    ! SCALAR ARGUMENTS
    real(8) :: m, s, x

    intent(in) :: m, s, x

    ! RETURN VALUE
    real(8) :: fgauss

    fgauss = 1 / (sqrt(2 * PI) * s) * &
             exp(- 5.0D-1 * (x - m) ** 2.0D0)
    
  end function fgauss

  !------------------------------------------------------------!
  ! SUBROUTINE EVALDPROB                                       !
  !                                                            !
  !                                                            !
  !------------------------------------------------------------!

  subroutine evaldprob(n, x, MU, CORR, ABSERR, RELERR, dp, flag)

    implicit none

    ! SCALAR ARGUMENTS
    integer :: flag, n
    real(8) :: ABSERR, RELERR

    ! ARRAY ARGUMENTS
    real(8) :: CORR((n - 1) * (n - 2) / 2), dp(n), MU(n), x(n)

    intent(in ) :: ABSERR, CORR, MU, n, RELERR, x
    intent(out) :: dp, flag

    ! Interface

    interface
       SUBROUTINE MVNDST( N, LOWER, UPPER, INFIN, CORREL, MAXPTS, &
            ABSEPS, RELEPS, ERROR, VAL, INFORM )

         integer :: N, INFORM, MAXPTS
         real(8) :: ABSEPS, RELEPS, VAL, ERROR

         integer :: INFIN(N)
         real(8) :: LOWER(N), UPPER(N), CORREL((N - 1) * N / 2)

         intent(in ) :: N, MAXPTS, INFIN, LOWER, UPPER, CORREL, &
                        ABSEPS, RELEPS
         intent(out) :: INFORM, VAL, ERROR

       end SUBROUTINE MVNDST
    end interface

    ! LOCAL SCALARS
    integer :: i, j, k, pos
    real(8) :: error, p, ff

    ! LOCAL ARRAYS
    real(8) :: l(n - 1), u(n - 1), ncorr((n - 2) * (n - 3) / 2)
    integer :: infty(n - 1)

    do k = 1,n

       ff = fgauss(x(k), MU(k), 1.0D0)

       if ( ff .lt. 1.0D-60 ) then

          pos = 1

          do j = 1, n

             if ( j .eq. k ) continue

             do i = 1, j - 1

                if ( i .eq. k ) continue

                ncorr(pos) = corr(i + (j - 1) * (j - 2) / 2)

                pos = pos + 1

             end do

          end do

          ! Affine transformation.
          ! Only need to shift 'mu'

          pos = 1
          
          do i = 1, n

             if ( j .eq. k ) continue

             u(pos) = x(i) - MU(i)

             pos = pos + 1

          end do

          do i = 1, n - 1

             infty(i) = 0

          end do

          flag = 0

          call mvndst(n - 1, l, u, infty, ncorr, &
               5000 * (n - 1) * (n - 1) * (n - 1), &
               ABSERR, RELERR, error, p, flag)

       end if

    end do

  end subroutine evaldprob

  !------------------------------------------------------------!
  ! SUBROUTINE EVALPROB                                        !
  !                                                            !
  ! This subroutine evaluates the probability P(\x_i \le x_i)  !
  ! for i \in \{1, \dots, np\}                                 !
  !                                                            !
  !------------------------------------------------------------!

  subroutine evalprob(n, x, MU, CORR, ABSERR, RELERR, p, flag)

    implicit none

    ! SCALAR ARGUMENTS
    integer :: flag, n
    real(8) :: ABSERR, p, RELERR

    ! ARRAY ARGUMENTS
    real(8) :: CORR((n - 1) * (n - 2) / 2), MU(n), x(n)

    intent(in ) :: ABSERR, CORR, MU, n, RELERR, x
    intent(out) :: p, flag

    ! Interface

    interface
       SUBROUTINE MVNDST( N, LOWER, UPPER, INFIN, CORREL, MAXPTS, &
            ABSEPS, RELEPS, ERROR, VAL, INFORM )

         integer :: N, INFORM, MAXPTS
         real(8) :: ABSEPS, RELEPS, VAL, ERROR

         integer :: INFIN(N)
         real(8) :: LOWER(N), UPPER(N), CORREL((N - 1) * N / 2)

         intent(in ) :: N, MAXPTS, INFIN, LOWER, UPPER, CORREL, &
                        ABSEPS, RELEPS
         intent(out) :: INFORM, VAL, ERROR

       end SUBROUTINE MVNDST
    end interface

    ! LOCAL SCALARS
    integer :: i
    real(8) :: error

    ! LOCAL ARRAYS
    real(8) :: l(n), u(n)
    integer :: infty(n)

    ! Affine transformation.
    ! Only need to shift 'mu'

    do i = 1, n

       u(i) = x(i) - MU(i)

    end do

    do i = 1, n

       infty(i) = 0

    end do

    flag = 0

    call mvndst(n, l, u, infty, corr, 5000 * n * n * n, ABSERR, &
                RELERR, error, p, flag)

  end subroutine evalprob


end program derivatives_cost
