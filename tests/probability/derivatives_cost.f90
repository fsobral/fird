program derivatives_cost

  use ccpdata, only: MU, CORR, ABSERR, RELERR, set_seed, initialize, &
                     destroy

  integer :: np, flag, nruns
  real(8) :: start, finish, p, dfttime, fttime

  real(8), allocatable :: x(:), l(:), u(:), dp(:)

  write(*,*) "This program tests the cost of derivatives"
  write(*,*) "of CCP functions"

  nruns = 10

  write(*,FMT=0002) 'DIM', 'F TIME', 'G TIME', 'RELATION'

  do i = 1, 5

     fttime = 0.0D0

     dfttime = 0.0D0

     do j = 1, nruns

        np = 10 * i

        allocate(x(np),l(np),u(np),dp(np))

        call set_seed(12345678 + i * j)

        ! Function

        call initialize(np, np, 0, x, l, u, 0.0D0)

        call CPU_TIME(start)

        call evalprob(np, x, MU, CORR, ABSERR, RELERR, p, flag)
        
        call CPU_TIME(finish)

        fttime = fttime + (finish - start)

        ! Derivatives

        call CPU_TIME(start)

        call evaldprob(np, x, MU, CORR, ABSERR, RELERR, dp, flag)
        
        call CPU_TIME(finish)

        dfttime = dfttime + (finish - start)

        deallocate(x, l, u, dp)

        call destroy()

     end do

     write(*,FMT=0001) np, (fttime / nruns), (dfttime / nruns), &
                       (dfttime / fttime)

  end do

  0002 FORMAT(/,A5,1X,A11,1X,A11,1X,A11)
  0001 FORMAT(I5,1X,F11.8,1X,F11.8,1X,F11.8)

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
    real(8) :: CORR(n * (n - 1) / 2), dp(n), MU(n), x(n)

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
    integer :: i, j, k, pos, pi, pj
    real(8) :: error, ff

    ! LOCAL ARRAYS
    real(8) :: l(n - 1), u(n - 1), ncorr((n - 1) * (n - 2) / 2)
    integer :: infty(n - 1)

    flag = 0

    do k = 1,n

       ff = fgauss(x(k), MU(k), 1.0D0)

       dp(k) = 0.0D0

       if ( ff .gt. 1.0D-60 ) then

          pos = 1

          do j = 1, n

             if ( j .eq. k ) cycle

             do i = 1, j - 1

                if ( i .eq. k ) cycle

                if ( i .lt. k ) then

                   pi = i + (k - 1) * (k - 2) / 2

                else

                   pi = k + (i - 1) * (i - 2) / 2

                end if

                if ( j .lt. k ) then

                   pj = j + (k - 1) * (k - 2) / 2

                else

                   pj = k + (j - 1) * (j - 2) / 2

                end if

                ncorr(pos) = corr(i + (j - 1) * (j - 2) / 2)

                ncorr(pos) = ncorr(pos) - corr(pi) * corr(pj)

                pos = pos + 1

             end do

          end do

          ! Affine transformation.
          ! Only need to shift 'mu'

          do i = 1, n

             if ( i .lt. k ) then

                u(i) = x(i) - (MU(i) + (x(k) - MU(k)) * &
                                 corr(i + (k - 1) * (k - 2) / 2))

             else if ( i .gt. k ) then

                u(i - 1) = x(i) - (MU(i) + (x(k) - MU(k)) * &
                                   corr(k + (i - 1) * (i - 2) / 2))

             end if

          end do

          do i = 1, n - 1

             infty(i) = 0

          end do

          flag = 0

          call mvndst(n - 1, l, u, infty, ncorr, &
               5000 * (n - 1) * (n - 1) * (n - 1), &
               ABSERR, RELERR, error, dp(k), flag)

          dp(k) = dp(k) * ff

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
    real(8) :: CORR(n * (n - 1) / 2), MU(n), x(n)

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
