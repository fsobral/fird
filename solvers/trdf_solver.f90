module trdf_solver

  use userinterface

  implicit none

  ! GLOBAL USER-DEFINED SUBROUTINES
  procedure(evalf  ), pointer :: uevalf
  procedure(evalc  ), pointer :: uevallc,uevalc
  procedure(evaljac), pointer :: uevalljac

  private

  public :: solver

contains

  ! Uses the adapted TRDF algorithm for solving the optimality phase

  subroutine solver(n,y,l,u,me,mi,uevalf_,uevalc_,uevallc_,uevalljac_, &
       nf,alpha,ffilter,hfilter,outiter,epsfeas,epsopt,verbose, &
       delta,fy,hynorm,rho,flag)

    use trdf

    implicit none

    ! SCALAR ARGUMENTS
    logical :: verbose
    integer :: flag,me,mi,n,nf,outiter
    real(8) :: alpha,delta,epsfeas,epsopt,fy,hynorm,rho

    ! ARRAY ARGUMENTS
    real(8) :: ffilter(nf),hfilter(nf),l(n),u(n),y(n)

    ! EXTERNAL SUBROUTINES
    external :: uevalf_,uevallc_,uevalljac_,uevalc_

    ! LOCAL SCALARS
    integer :: i,m,maxfcnt,npt,fcnt
    real(8) :: rbeg,rend,xeps

    ! LOCAL ARRAYS
    logical :: ccoded(2),equatn(me + mi),linear(me + mi)

    uevalf   => uevalf_
    uevallc   => uevallc_
    uevalljac => uevalljac_
    uevalc    => uevalc_

    m = me + mi

    do i = 1,me
       equatn(i) = .true.
       linear(i) = .true.
    end do
    do i = me + 1,m
       equatn(i) = .false.
       linear(i) = .true.
    end do

    if ( N .le. 2 ) then
       NPT = 2 * N + 1
    else
       NPT = 2 * N + 3
    end if

    ccoded(1) = .true.
    ccoded(2) = .true.

    maxfcnt = 1000 * n

    rbeg = rho !max(10.0D0 * epsopt, rho)

    rend = epsopt

    xeps = 1.0D-08

    !delta = max(DELMIN, rbeg, delta)

    call TRDFSUB(N,NPT,Y,L,U,M,EQUATN,LINEAR,CCODED,UEVALF,UEVALLC, &
         TRDF_EVALJAC,TRDF_EVALHC,UEVALC,MAXFCNT,RBEG,REND,XEPS,VERBOSE, &
         NF,ALPHA,FFILTER,HFILTER,OUTITER,DELTA,EPSFEAS,FY, &
         HYNORM,FCNT,RHO,FLAG)

  end subroutine solver

  ! ******************************************************************
  ! ******************************************************************

  subroutine trdf_evaljac(n,x,ind,jcvar,jcval,jcnnz,lim,lmem,flag)

    implicit none

    ! SCALAR ARGUMENTS
    logical :: lmem
    integer :: flag,ind,jcnnz,lim,n

    ! ARRAY ARGUMENTS
    integer :: jcvar(lim)
    real(8) :: x(n),jcval(lim)

    flag = -1

    lmem = .false.

    call uevalljac(n,x,ind,jcvar,jcval,jcnnz,flag)

  end subroutine trdf_evaljac

  ! ******************************************************************
  ! ******************************************************************

  subroutine trdf_evalhc(n,x,ind,hcrow,hccol,hcval,hcnnz,lim,lmem,flag)

    implicit none

    ! SCALAR ARGUMENTS
    logical :: lmem
    integer :: flag,hcnnz,ind,lim,n

    ! ARRAY ARGUMENTS
    integer :: hccol(lim),hcrow(lim)
    real(8) :: hcval(lim),x(n)

    flag = 0

    lmem = .false.

    hcnnz = 0

  end subroutine trdf_evalhc
  
end module trdf_solver
