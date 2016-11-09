program HSCCP

  use hsccpdata, only: ABSERR, initialize, destroy, MU, CORR, RELERR, &
       set_seed

  implicit none

  ! INTERFACES

  interface

     subroutine hsgetdim(nprob, n_, me_, mi_)
       ! SCALAR ARGUMENTS
       integer :: nprob, n_, me_, mi_

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
  real(8) :: epsfeas,epsopt,f,p,feas,plim

  ! LOCAL ARRAYS
  real(8), allocatable :: l(:), u(:), x(:)

  ! Read the seed

  read(*,*) prob

  if ( .not. verifyProb(prob) ) stop

  call hsgetdim(prob, n, me, mi)

  np = n

  n = n + np

  plim = 8.0D-01

  allocate(x(n),l(n),u(n))

  call hsinitp(n - np, x(np + 1), l(np + 1), u(np + 1))

  call set_seed(12345678 + prob)

  call initialize(np, mi, x, l, u, plim)

  ! Call the solver

  verbose = .false.

  epsfeas = 1.0D-8
  epsopt  = 1.0D-4

  ftype = 2

  call fird(n,x,l,u,me,mi,evalf,evalc,evaljac,verbose,ftype, &
       epsfeas,epsopt,f,feas,fcnt,flag)

  call evalprob(np, x, MU, CORR, ABSERR, RELERR, p, flag)

  write(*, FMT=021) 'PRB', 'NO', 'NP', 'INEQ', 'EQ', 'F', 'FEAS', &
       'PROB', 'FEVAL', 'FLG'
     
  write(*, FMT=022) prob, n - np, np, mi, me, f, feas, &
       p, fcnt, flag
 
  open(99, FILE='ccp.out')

  write(99, FMT=020) prob, n - np, np, mi, me, f, feas, p, fcnt, flag

  close(99)

  deallocate(x,l,u)

  call destroy()

  ! NON-EXECUTABLE STATEMENTS

020 FORMAT(I10,1X,I5,1X,I5,1X,I5,1X,I5,1X,E15.8,1X,E15.8,1X,F10.8,1X, &
         I10,1X,I3)
021 FORMAT(/,A3,1X,A4,1X,A4,1X,A4,1X,A4,1X,A12,1X,A12,1X,A6,1X, &
         A9,1X,A3)
022 FORMAT(I3,1X,I4,1X,I4,1X,I4,1X,I4,1X,1P,E12.4,1X,1P,E12.4,1X, &
         0P,F6.4,1X,I9,1X,I3)

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

    ! Evaluate the probability

    call evalprob(np, x, mu, corr, ABSERR, RELERR, pf, flag)

    if ( flag .ne. 0 ) then

       return

    end if

    ! Evaluate the HS objective function

    call hsobjf(n - np, x(np + 1), hsf, flag)

    if ( flag .ne. 0 ) then

       return

    end if

    f = hsf + PEN * (plim - pf)

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
         ! PARAMETERS
         integer, parameter :: MMAX = 1000
         integer, parameter :: NMAX = 1000
         ! SCALAR ARGUMENTS
         integer :: flag,ind,jcnnz,nn
         ! ARRAY ARGUMENTS
         integer :: jcvar(MMAX * NMAX)
         real(8) :: x_(nn), jcval(MMAX * NMAX)

         intent(in) :: nn, x_, ind
         intent(out) :: jcvar, jcval, jcnnz, flag
       end subroutine hsjac

    end interface
    
    ! LOCAL SCALARS
    integer :: j

    flag = 0
    
    if ( ind .ge. 1 .and. ind .le. m ) then

       call hsjac(n - np, x(np + 1), ind, jcvar, jcval, jcnnz, flag)

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

  !------------------------------------------------------------!
  ! SUBROUTINE VERIFYPROB                                      !
  !                                                            !
  ! Tests whether the HS problem can be selected or not        !
  !                                                            !
  !------------------------------------------------------------!

  function verifyProb(prob)

    ! SCALAR ARGUMENTS

    integer :: prob

    ! RETURN VALUE

    logical :: verifyProb

    verifyProb = .false.

    if ( prob .eq. 10 .or. & 
         prob .eq. 11 .or. & 
         prob .eq. 12 .or. & 
         prob .eq. 13 .or. & 
         prob .eq. 15 .or. & 
         prob .eq. 16 .or. & 
         prob .eq. 17 .or. & 
         prob .eq. 18 .or. & 
         prob .eq. 19 .or. & 
         prob .eq. 20 .or. & 
         prob .eq. 21 .or. & 
         prob .eq. 22 .or. & 
         prob .eq. 23 .or. & 
         prob .eq. 24 .or. & 
         prob .eq. 29 .or. & 
         prob .eq. 30 .or. & 
         prob .eq. 31 .or. & 
         prob .eq. 33 .or. & 
         prob .eq. 34 .or. & 
         prob .eq. 35 .or. & 
         prob .eq. 36 .or. & 
         prob .eq. 37 .or. & 
         prob .eq. 43 .or. & 
         prob .eq. 44 .or. & 
         prob .eq. 57 .or. & 
         prob .eq. 58 .or. & 
         prob .eq. 59 .or. & 
         prob .eq. 64 .or. & 
         prob .eq. 65 .or. & 
         prob .eq. 66 .or. & 
         prob .eq. 67 .or. & 
         prob .eq. 70 .or. & 
         prob .eq. 72 .or. & 
         prob .eq. 76 .or. & 
         prob .eq. 83 .or. & 
         prob .eq. 84 .or. & 
         prob .eq. 85 .or. & 
         prob .eq. 86 .or. & 
         prob .eq. 88 .or. & 
         prob .eq. 89 .or. & 
         prob .eq. 90 .or. & 
         prob .eq. 91 .or. & 
         prob .eq. 92 .or. & 
         prob .eq. 93 .or. & 
         prob .eq. 95 .or. & 
         prob .eq. 96 .or. & 
         prob .eq. 97 .or. & 
         prob .eq. 98 .or. & 
         prob .eq. 100 .or. & 
         prob .eq. 101 .or. & 
         prob .eq. 102 .or. & 
         prob .eq. 103 .or. & 
         prob .eq. 104 .or. & 
         prob .eq. 105 .or. & 
         prob .eq. 106 .or. & 
         prob .eq. 108 .or. & 
         prob .eq. 113 .or. & 
         prob .eq. 116 .or. & 
         prob .eq. 117 .or. & 
         prob .eq. 118 .or. & 
         prob .eq. 215 .or. & 
         prob .eq. 218 .or. & 
         prob .eq. 220 .or. & 
         prob .eq. 221 .or. & 
         prob .eq. 222 .or. & 
         prob .eq. 223 .or. & 
         prob .eq. 224 .or. & 
         prob .eq. 225 .or. & 
         prob .eq. 226 .or. & 
         prob .eq. 227 .or. & 
         prob .eq. 228 .or. & 
         prob .eq. 230 .or. & 
         prob .eq. 231 .or. & 
         prob .eq. 232 .or. & 
         prob .eq. 233 .or. & 
         prob .eq. 234 .or. & 
         prob .eq. 236 .or. & 
         prob .eq. 237 .or. & 
         prob .eq. 238 .or. & 
         prob .eq. 239 .or. & 
         prob .eq. 249 .or. & 
         prob .eq. 250 .or. & 
         prob .eq. 251 .or. & 
         prob .eq. 253 .or. & 
         prob .eq. 264 .or. & 
         prob .eq. 268 .or. & 
         prob .eq. 270 .or. & 
         prob .eq. 277 .or. & 
         prob .eq. 278 .or. & 
         prob .eq. 279 .or. & 
         prob .eq. 280 .or. & 
         prob .eq. 284 .or. & 
         prob .eq. 285 .or. & 
         prob .eq. 315 .or. & 
         prob .eq. 323 .or. & 
         prob .eq. 324 .or. & 
         prob .eq. 326 .or. & 
         prob .eq. 327 .or. & 
         prob .eq. 329 .or. & 
         prob .eq. 330 .or. & 
         prob .eq. 331 .or. & 
         prob .eq. 332 .or. & 
         prob .eq. 337 .or. & 
         prob .eq. 339 .or. & 
         prob .eq. 340 .or. & 
         prob .eq. 341 .or. & 
         prob .eq. 342 .or. & 
         prob .eq. 343 .or. & 
         prob .eq. 346 .or. & 
         prob .eq. 347 .or. & 
         prob .eq. 349 .or. & 
         prob .eq. 354 .or. & 
         prob .eq. 356 .or. & 
         prob .eq. 359 .or. & 
         prob .eq. 360 .or. & 
         prob .eq. 361 .or. & 
         prob .eq. 362 .or. & 
         prob .eq. 363 .or. & 
         prob .eq. 364 .or. & 
         prob .eq. 365 .or. & 
         prob .eq. 366 .or. & 
         prob .eq. 369 .or. & 
         prob .eq. 372 .or. & 
         prob .eq. 374 .or. & 
         prob .eq. 380 .or. & 
         prob .eq. 384 .or. & 
         prob .eq. 385 .or. & 
         prob .eq. 386 .or. & 
         prob .eq. 387 .or. & 
         prob .eq. 388 .or. & 
         prob .eq. 389 .or. & 
         prob .eq. 392 ) then

       verifyProb = .true.

    end if

    return

  end function verifyProb

end program HSCCP
