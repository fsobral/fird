program CCP

  use ccpdata, only: initialize, destroy, MU, CORR

  implicit none

  ! LOCAL SCALARS
  integer :: flag,ftype,i,me,mi,n,np,fcnt
  logical :: verbose
  real(8) :: epsfeas,epsopt,f,p,feas,plim,npfrac

  ! LOCAL ARRAYS
  real(8), allocatable :: l(:), u(:), x(:)

  npfrac = 3.0D0 / 4.0D0

  n  = 10

  np = INT(n * npfrac)

  me = 0

  mi = 10

  plim = 8.0D-01

  allocate(x(n),l(n),u(n))

  call initialize(n, np, mi, x, l, u, plim)

  ! Call the solver

  verbose = .true.

  epsfeas = 1.0D-8
  epsopt  = 1.0D-4

  ftype = 2

  call fird(n,x,l,u,me,mi,evalf,evalc,evaljac,verbose,ftype, &
       epsfeas,epsopt,f,feas,fcnt,flag)

  call evalprob(np, x, MU, CORR, p, flag)

  write(*, FMT=010) p

  deallocate(x,l,u)

  call destroy()

  ! NON-EXECUTABLE STATEMENTS

010 FORMAT(/,'Probability:',1X,1P,E15.8,/)

contains

  !------------------------------------------------------------!
  ! SUBROUTINE EVALPROB                                        !
  !                                                            !
  ! This subroutine evaluates the probability P(\x_i \le x_i)  !
  ! for i \in \{1, \dots, np\}                                 !
  !                                                            !
  !------------------------------------------------------------!

  subroutine evalprob(n, x, MU, CORR, p, flag)

    implicit none

    ! SCALAR ARGUMENTS
    integer :: flag, n
    real(8) :: p

    ! ARRAY ARGUMENTS
    real(8) :: CORR((n - 1) * (n - 2) / 2), MU(n), x(n)

    intent(in ) :: CORR, MU, n, x
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
    real(8) :: error,releps,abseps

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

    releps = 5.0D-05
    
    abseps = 0.0D0

    flag = 0

    call mvndst(n, l, u, infty, corr, 5000 * n * n * n, abseps, &
                releps, error, p, flag)

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

    flag = 1

    ! Penalize the probability

    call evalprob(np, x, mu, corr, f, flag)

    if ( flag .ne. 0 ) then

       return

    end if

    f = PEN * (plim - f)

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
