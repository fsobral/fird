program FIRDMAIN

  implicit none

  ! LOCAL SCALARS
  integer :: fcnt,flag,ftype,i,me,mi,n
  logical :: verbose
  real(8) :: epsfeas,epsopt,f,feas

  ! LOCAL ARRAYS
  real(8), allocatable :: l(:), u(:), x(:)

  n  = 2
  me = 1
  mi = 1
  
  allocate(x(n),l(n),u(n))

  x(1) = 2.0D0
  x(2) = 2.0D0

  l(1) = - 1.0D+20
  l(2) = - 1.0D+20
  u(1) =   1.0D+20
  u(2) =   1.0D+20

  verbose = .true.

  epsfeas = 1.0D-8
  epsopt  = 1.0D-4

  ftype = 1 ! Flat Filter

  call fird(n,x,l,u,me,mi,evalf,evalc,evaljac,verbose,ftype, &
       epsfeas,epsopt,f,feas,fcnt,flag)

  deallocate(x,l,u)

contains

  !------------------------------------------------------------!
  ! SUBROUTINE EVALF                                           !
  !                                                            !
  ! Defines the objetive function.                             !
  !                                                            !
  !------------------------------------------------------------!

  subroutine evalf(n,x,f,flag)

    ! SCALAR ARGUMENTS
    integer :: flag,n
    real(8) :: f

    ! ARRAY ARGUMENTS
    real(8) :: x(n)

    intent(in ) :: n,x
    intent(out) :: f,flag

    flag = 0

    f = (x(1) - 2.0D0) ** 2.0D0 + (x(2) - 1.0D0) ** 2.0D0

  end subroutine evalf

  !------------------------------------------------------------!
  ! SUBROUTINE EVALC                                           !
  !                                                            !
  ! Defines the contraints.                                    !
  !                                                            !
  !------------------------------------------------------------!

  subroutine evalc(n,x,ind,c,flag)

  ! SCALAR ARGUMENTS
  integer :: flag,ind,n
  real(8) :: c

  ! ARRAY ARGUMENTS
  real(8) :: x(n)

  intent(in ) :: ind,n,x
  intent(out) :: c,flag

  flag = 0

  if ( ind .eq. 1 ) then
     c = x(1) - 2.0D0 * x(2) + 1.0D0
  else if ( ind .eq. 2 ) then
     c = 0.25D0 * x(1) ** 2.0D0 + x(2) ** 2.0D0 - 1.0D0
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

    ! SCALAR ARGUMENTS
    integer :: flag,ind,jcnnz,n

    ! ARRAY ARGUMENTS
    integer :: jcvar(n)
    real(8) :: jcval(n),x(n)

    intent(in ) :: ind,n,x
    intent(out) :: flag,jcnnz,jcval,jcvar

    flag = 0

    if ( ind .eq. 1 ) then
       jcnnz = 2
       jcvar(1) = 1
       jcval(1) = 1.0D0
       jcvar(2) = 2
       jcval(2) = - 2.0D0
    else if ( ind .eq. 2 ) then
       jcnnz = 2
       jcvar(1) = 1
       jcval(1) = 2.0 * 0.25D0 * x(1)
       jcvar(2) = 2
       jcval(2) = 2.0D0 * x(2)
    else
       flag = 1
    end if

  end subroutine evaljac

end program FIRDMAIN
