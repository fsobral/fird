program CCP

  use ccpdata, only: ABSERR, initialize, destroy, MU, CORR, RELERR, &
                     set_seed, set_epsfeas, fromfile, PEN

  implicit none

  ! LOCAL SCALARS
  integer :: flag, ftype, i, me, mi, n, np, fcnt, prob, stat, ftcnt
  logical :: verbose
  real(8) :: epsfeas, epsopt, f, p, feas, plim, npfrac
  character(200) :: filename

  ! LOCAL ARRAYS
  real(8), allocatable :: l(:), u(:), x(:)

  call GETARG(1, filename)

  ! Verify if the argument is a filename or a problem number

  open(99, FILE = filename, STATUS = "old", IOSTAT = stat)

  if ( stat .eq. 0 ) then

     close(99)

     prob = -1

     me = 0

     call set_seed(12345678)

     call fromfile(n, np, mi, x, l, u, plim, filename)

  else

     ! Read the seed

     write(*, *) 'Write the number of a problem or just a seed'

     read(*, *) prob

     npfrac = 1.0D0 / 2.0D0

     n  = 4

     np = INT(n * npfrac)

     me = 0

     mi = 10

     plim = 8.0D-01
     
     allocate(x(n),l(n),u(n))
  
     call set_seed(12345678 + prob)

     call initialize(n, np, mi, x, l, u, plim)

  end if

  ! Create basic output files, in case of external failure

  open(99, FILE='ccp.out')

  if ( prob .eq. -1 ) then

     write(99, FMT=021) filename, n, np, mi, 1.0D+20, 1.0D+20, &
          1.0D+20, 1.0D+20, -1, -1

  else

     write(99, FMT=020) prob, n, np, mi, 1.0D+20, 1.0D+20, &
          1.0D+20, 1.0D+20, -1, -1

  end if

  close(99)

  open(99, FILE = 'ccp.sol')

  write(99, FMT=022) (x(i), i = 1, n)

  close(99)

  ! Call the solver - First run

  verbose = .false.

  epsfeas = 1.0D-04

  epsopt  = 1.0D-01

  call set_epsfeas(epsfeas)

  PEN = 1.0D+03

  ftype = 2

  call fird(n,x,l,u,me,mi,evalf,evalc,evaljac,verbose,ftype, &
       epsfeas,epsopt,f,feas,fcnt,flag)

  ftcnt = fcnt

  ! Call the solver - Second run

  verbose = .false.

  epsfeas = 1.0D-04

  epsopt  = 1.0D-04

  call set_epsfeas(epsfeas)

  PEN = 1.0D+10

  ftype = 2

  call fird(n,x,l,u,me,mi,evalf,evalc,evaljac,verbose,ftype, &
       epsfeas,epsopt,f,feas,fcnt,flag)

  ftcnt = ftcnt + fcnt

  ! Print status

  call evalprob(np, x, MU, CORR, ABSERR, RELERR, p, flag)

  if ( verbose ) write(*, FMT=020) prob, n, np, mi, f, feas, &
       f - PEN * max(0.0D0, plim - p - epsfeas) ** 2.0D0, p, ftcnt, flag

  open(99, FILE='ccp.out')

  if ( prob .eq. -1 ) then

     write(99, FMT=021) filename, n, np, mi, f, feas, &
          f - PEN * max(0.0D0, plim - p - epsfeas) ** 2.0D0, &
          p, ftcnt, flag

  else

     write(99, FMT=020) prob, n, np, mi, f, feas, &
          f - PEN * max(0.0D0, plim - p - epsfeas) ** 2.0D0, &
          p, ftcnt, flag

  end if

  close(99)

  ! Save solution

  open(99, FILE = 'ccp.sol')

  write(99, FMT=022) (x(i), i = 1, n)

  close(99)

  deallocate(x,l,u)

  call destroy()

  ! NON-EXECUTABLE STATEMENTS

020 FORMAT(I10,1X,I5,1X,I5,1X,I5,1X,E15.8,1X,E15.8,1X,1PE15.8,1X, &
         0PF10.8,1X,I10,I3)
021 FORMAT(A100,1X,I5,1X,I5,1X,I5,1X,E15.8,1X,E15.8,1X,1PE15.8,1X, &
         0PF10.8,1X,I10,I3)
022 FORMAT((1PE23.16))

contains

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
    real(8) :: CORR((n - 1) * n / 2), MU(n), x(n)

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

  !------------------------------------------------------------!
  ! SUBROUTINE EVALF                                           !
  !                                                            !
  ! Defines the objetive function.                             !
  !                                                            !
  !------------------------------------------------------------!

  subroutine evalf(n,x,f,flag)

    use ccpdata

    implicit none

    ! SCALAR ARGUMENTS
    integer :: flag,n
    real(8) :: f

    ! ARRAY ARGUMENTS
    real(8) :: x(n)

    intent(in ) :: n,x
    intent(out) :: f,flag

    ! LOCAL SCALARS
    integer :: i

    flag = - 1

    ! Penalize the probability

    call evalprob(np, x, mu, corr, ABSERR, RELERR, f, flag)

    if ( flag .ne. 0 ) then

       write(*,*) 'Error: MVNDST returned', flag

       flag = - 1

       return

    end if

    f = PEN * max(0.0D0, plim - f - epsfeas) ** 2.0D0

    ! Quadratic term

    do i = 1, n

       f = f + c(i) * x(i) + 5.0D-01 * x(i) ** 2.0D0

    end do

  end subroutine evalf

  !------------------------------------------------------------!
  ! SUBROUTINE EVALC                                           !
  !                                                            !
  ! Defines the contraints.                                    !
  !                                                            !
  !------------------------------------------------------------!

  subroutine evalc(n,x,ind,c,flag)

    use ccpdata, only: m, cA, pA, A, b

    implicit none

    ! SCALAR ARGUMENTS
    integer :: flag,ind,n
    real(8) :: c
    
    ! ARRAY ARGUMENTS
    real(8) :: x(n)
    
    intent(in ) :: ind,n,x
    intent(out) :: c,flag

    ! LOCAL SCALARS
    integer :: j

    flag = 0

    if ( ind .ge. 1 .and. ind .le. m ) then

       c = - b(ind)

       do j = pA(ind), pA(ind + 1) - 1

          c = c + x(cA(j)) * A(j)

       end do

    else

       flag = 1

    end if

  end subroutine evalc

  !------------------------------------------------------------!
  ! SUBROUTINE EVALJAC                                         !
  !                                                            !
  ! Defines the Jacobian of the constraints.                   !
  !                                                            !
  !------------------------------------------------------------!
  
  subroutine evaljac(n,x,ind,jcvar,jcval,jcnnz,flag)

    use ccpdata, only: m, cA, pA, A

    implicit none

    ! SCALAR ARGUMENTS
    integer :: flag,ind,jcnnz,n

    ! ARRAY ARGUMENTS
    integer :: jcvar(n)
    real(8) :: jcval(n),x(n)

    intent(in ) :: ind,n,x
    intent(out) :: flag,jcnnz,jcval,jcvar

    ! LOCAL SCALARS
    integer :: j

    flag = 0

    if ( ind .ge. 1 .and. ind .le. m ) then

       jcnnz = 0

       do j = pA(ind), pA(ind + 1) - 1

          jcnnz = jcnnz + 1

          jcvar(jcnnz) = cA(j)

          jcval(jcnnz) =  A(j)

       end do

    else

       flag = 1

    end if

  end subroutine evaljac

end program CCP
