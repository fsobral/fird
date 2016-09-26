program HSCCP

  use hsccpdata, only: ABSERR, initialize, destroy, MU, CORR, RELERR, &
       set_seed

  implicit none

  ! INTERFACES

  interface

     subroutine hsgetdim(nprob, n_, me_, mi_)
       ! SCALAR ARGUMENTS
       integer nprob, n_, me_, mi_

       intent(in) :: nprob
       intent(out) :: n_, me_, mi_
     end subroutine hsgetdim

     subroutine hsinitp(n_, x_, l_, u_)
       ! SCALAR ARGUMENTS
       integer :: n_
       ! ARRAY ARGUMENTS
       real(8) :: x_(n_), l_(n_), u_(n_)

       intent(in) :: n_
       intent(out) :: x_, l_, u_
     end subroutine hsinitp

  end interface

  ! LOCAL SCALARS
  integer :: flag,ftype,i,me,mi,n,nt,np,fcnt,prob
  logical :: verbose
  real(8) :: epsfeas,epsopt,f,p,feas,plim,npfrac

  ! LOCAL ARRAYS
  real(8), allocatable :: l(:), u(:), x(:)

  ! Read the seed

  read(*,*) prob

  npfrac = 1.0D0 / 2.0D0

  call hsgetdim(nprob, n, me, mi)

  np = INT(n * npfrac)

  n = n + np

  plim = 8.0D-01

  allocate(x(n),l(n),u(n))

  call hsinitp(n - np, x(np + 1), l(np + 1), u(np + 1))

  call set_seed(12345678 + prob)

  call initialize(np, mi, x, l, u, plim)

  ! Call the solver

  verbose = .true.

  epsfeas = 1.0D-8
  epsopt  = 1.0D-4

  ftype = 2

  call fird(n,x,l,u,me,mi,evalf,evalc,evaljac,verbose,ftype, &
       epsfeas,epsopt,f,feas,fcnt,flag)

  call evalprob(np, x, MU, CORR, ABSERR, RELERR, p, flag)

  if ( verbose ) write(*, FMT=020) prob, n, np, mi, f, feas, p, &
       fcnt, flag

  open(99, FILE='ccp.out')

  write(99, FMT=020) prob, n, np, mi, f, feas, p, fcnt, flag

  close(99)

  deallocate(x,l,u)

  call destroy()

  ! NON-EXECUTABLE STATEMENTS

020 FORMAT(I10,1X,I5,1X,I5,1X,I5,1X,E15.8,1X,E15.8,1X,F10.8,1X,I10, &
       I3)

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

  !------------------------------------------------------------!
  ! SUBROUTINE EVALF                                           !
  !                                                            !
  ! Defines the objetive function.                             !
  !                                                            !
  !------------------------------------------------------------!

  subroutine evalf(n,x,f,flag)

    use hsccpdata

    implicit none

    ! SCALAR ARGUMENTS
    integer :: flag,n
    real(8) :: f

    ! ARRAY ARGUMENTS
    real(8) :: x(n)

    intent(in ) :: n,x
    intent(out) :: f,flag

    ! INTERFACES

    interface

       subroutine hsobjf(n,x_,f,flag)
         ! SCALAR ARGUMENTS
         integer :: flag,n
         real(8) :: f
         ! ARRAY ARGUMENTS
         real(8) :: x_(n)

         intent(in) :: n, x_
         intent(out) :: f, flag
       end subroutine hsobjf

    end interface

    ! LOCAL SCALARS
    integer :: i
    real(8) :: hsf, pf

    flag = 1

    ! Penalize the probability

    call evalprob(np, x, mu, corr, ABSERR, RELERR, hsf, flag)

    if ( flag .ne. 0 ) then

       return

    end if

    ! Evaluate the HS objective function

    call hsobjf(n - np, x(np + 1), pf, flag)

    if ( flag .ne. 0 ) then

       return

    end if

    f = hsf + PEN * (plim - f)

    ! Linear term in the probability variables

    do i = 1, np

       f = f + c(i) * x(i)

    end do

  end subroutine evalf

  !------------------------------------------------------------!
  ! SUBROUTINE EVALC                                           !
  !                                                            !
  ! Defines the contraints.                                    !
  !                                                            !
  !------------------------------------------------------------!

  subroutine evalc(n,x,ind,c,flag)

    use hsccpdata, only: m, cA, pA, A

    implicit none

    ! SCALAR ARGUMENTS
    integer :: flag,ind,n
    real(8) :: c

    ! ARRAY ARGUMENTS
    real(8) :: x(n)

    intent(in ) :: ind,n,x
    intent(out) :: c,flag

    ! INTERFACES

    interface

       subroutine hscon(nn,x_,ind,c,flag)
         ! SCALAR ARGUMENTS
         integer :: flag,ind,nn
         real(8) :: c
         ! ARRAY ARGUMENTS
         real(8) :: x_(nn)

         intent(in) :: nn, x_, ind
         intent(out) :: c, flag
       end subroutine hscon

    end interface

    ! LOCAL SCALARS
    integer :: j

    flag = 0

    if ( ind .ge. 1 .and. ind .le. m ) then

       ! Evaluate HS constraint

       call hscon(n - np, x(np + 1), ind, c, flag)

       if ( flag .ne. 0 ) then

          flag = 1

          return

       end if

       ! Append the linear terms, related to probability variables

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

    use hsccpdata, only: np, m, cA, pA, A

    implicit none

    ! SCALAR ARGUMENTS
    integer :: flag,ind,jcnnz,n

    ! ARRAY ARGUMENTS
    integer :: jcvar(n)
    real(8) :: jcval(n),x(n)

    intent(in ) :: ind,n,x
    intent(out) :: flag,jcnnz,jcval,jcvar

    ! INTERFACES

    interface

       subroutine hsjac(nn,x_,ind,jcvar,jcval,jcnnz,flag)
         ! SCALAR ARGUMENTS
         integer :: flag,ind,jcnnz,nn
         ! ARRAY ARGUMENTS
         integer :: jcvar(:)
         real(8) :: x_(nn), jcval(:)

         intent(in) :: nn, x_, ind
         intent(out) :: jcvar, jcval, jcnnz, flag
       end subroutine hsjac

    end interface
    
    ! LOCAL SCALARS
    integer :: j

    flag = 0
    
    if ( ind .ge. 1 .and. ind .le. m ) then

       call hsjac(n - np, x_(np + 1), ind, jcvar, jcval, jcnnz, flag)

       if ( flag .ne. 0 ) then

          flag = 1

          return

       end if

       do j = 1, jcnnz

          jcvar(j) = jcvar(j) + np

       end do

       ! Append the linear terms, related to probability variables

       do j = pA(ind), pA(ind + 1) - 1

          jcnnz = jcnnz + 1

          jcvar(jcnnz) = cA(j)

          jcval(jcnnz) =  A(j)

       end do

    else

       flag = 1

    end if

  end subroutine evaljac

end program HSCCP
