program CCP

  use ccpstruct, only: initialize, destroy

  implicit none

  ! LOCAL SCALARS
  integer :: flag,ftype,i,me,mi,n,np,fcnt
  logical :: verbose
  real(8) :: epsfeas,epsopt,f,feas

  ! LOCAL ARRAYS
  real(8), allocatable :: l(:), u(:), x(:)

  n  = 2

  np = n / 2

  me = 0

  mi = 1

  allocate(x(n),l(n),u(n))

  call initialize(n, np, mi, x, l, u)

  ! Call the solver

  verbose = .true.

  epsfeas = 1.0D-8
  epsopt  = 1.0D-4

  ftype = 2

  call fird(n,x,l,u,me,mi,evalf,evalc,evaljac,verbose,ftype, &
       epsfeas,epsopt,f,feas,fcnt,flag)

  deallocate(x,l,u)

  call destroy()

contains

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
                releps, error, f, flag)

    f = - PEN * f

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

    use ccpdata, only: m, cA, pA

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

       c = 0.0

       do j = pA(ind), pA(ind + 1)

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

    use ccpdata, only: m, cA, pA

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

       do j = pA(ind), pA(ind + 1)

          jcnnz = jcnnz + 1

          jcvar(jcnnz) = cA(j)

          jcval(jcnnz) =  A(j)

       end do

    else

       flag = 1

    end if

  end subroutine evaljac

end program CCP