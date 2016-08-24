program FIRDMAIN

  implicit none

  ! LOCAL SCALARS
  integer :: flag,ftype,i,me,mi,n,fcnt
  logical :: verbose
  real(8) :: epsfeas,epsopt,f,feas

  ! LOCAL ARRAYS
  real(8), allocatable :: l(:), u(:), x(:)

  n  = 2
  me = 0
  mi = 1
  
  allocate(x(n),l(n),u(n))

  x(1) = 30.0D0
  x(2) = 30.0D0

  l(1) = 0.0D+00
  l(2) = 0.0D+00
  u(1) = 1.0D+20
  u(2) = 1.0D+20

  verbose = .true.

  epsfeas = 1.0D-8
  epsopt  = 1.0D-4

  ftype = 2

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
    integer :: i,j
    real(8) :: error,val,releps,abseps

    ! LOCAL ARRAYS
    real(8) :: correl((n - 1) * n / 2),l(n),u(n)
    integer :: infty(n)
    real(8) :: SIG(n,n), A(n,n), B(n), MU(n)

    ! Mu

    MU(:) = (/ 180.0D0, 162.0D0 /)

    ! Matrix A

    A(:,1) = (/ 2.0D0, 3.0D0 /)
    A(:,2) = (/ 6.0D0, 3.0D0 /)

    ! Covariance Matrix
    
    SIG(:,1) = (/ 144.0D0, 80.0D0 /)
    SIG(:,2) = (/  80.0D0, 81.0D0 /)

    ! Affine transformation

    B(:) = (/ (1.0D0 / 12.0D0), (1.0D0 / 9.0D0) /)

    do i = 1,n
       u(i) = - MU(i)
       do j = 1,n
          u(i) = u(i) + A(i,j) * x(j)
       end do
       u(i) = B(i) * u(i)
    end do

    do i = 1, n
       infty(i) = 0
    end do

    ! This multiplication is for inserting ones in the diagonal

    do j = 1, n
       do i = j + 1, n
          correl(j + (i - 1) * (i - 2) / 2) = SIG(i,j) * B(i) * B(j)
       end do
    end do

    releps = 5.0D-05
    
    abseps = 0.0D0

    flag = 0

    call mvndst(n,l,u,infty,correl,5000 * n * n * n,abseps,releps, &
         error,f,flag)

    f = - f

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
     c = 2.0D0 * x(1) + 3.0D0 * x(2) - 140.0D0
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
       jcval(1) = 2.0D0
       jcvar(2) = 2
       jcval(2) = 3.0D0
    else
       flag = 1
    end if

  end subroutine evaljac

end program FIRDMAIN
