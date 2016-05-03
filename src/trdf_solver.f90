module trdf_solver

  use userinterface

  implicit none

  ! GLOBAL USER-DEFINED SUBROUTINES
  procedure(evalf  ), pointer :: uevalf
  procedure(evalc  ), pointer :: uevalc
  procedure(evaljac), pointer :: uevaljac

  private

  public :: solver

contains

  ! Uses the adapted TRDF algorithm for solving the optimality phase

  subroutine solver(n,y,l,u,me,mi,uevalf_,uevalc_,uevaljac_, &
       nf,alpha,ffilter,hfilter,epsopt,fy,flag)

    use trdf

    implicit none

    ! SCALAR ARGUMENTS
    integer :: flag,me,mi,n,nf
    real(8) :: alpha,epsopt,fy

    ! ARRAY ARGUMENTS
    real(8) :: ffilter(nf),hfilter(nf),l(n),u(n),y(n)

    ! EXTERNAL SUBROUTINES
    external :: uevalf_,uevalc_,uevaljac_

    ! LOCAL SCALARS
    integer :: i,m,maxfcnt,npt,fcnt
    real(8) :: feas,rbeg,rend,xeps

    ! LOCAL ARRAYS
    logical :: ccoded(2),equatn(me + mi),linear(me + mi)

    uevalf   => uevalf_
    uevalc   => uevalc_
    uevaljac => uevaljac_

    m = me + mi

    do i = 1,me
       equatn(i) = .true.
       linear(i) = .true.
    end do
    do i = me + 1,m
       equatn(i) = .false.
       linear(i) = .true.
    end do

    NPT = 2 * N + 3

    ccoded(1) = .true.
    ccoded(2) = .true.

    maxfcnt = 1000 * n

    rbeg = 1.0D-1

    rend = epsopt

    xeps = 1.0D-8

    call TRDFSUB(N,NPT,Y,L,U,M,EQUATN,LINEAR,CCODED,UEVALF,UEVALC, &
         TRDF_EVALJAC,TRDF_EVALHC,MAXFCNT,RBEG,REND,XEPS, &
         NF,ALPHA,FFILTER,HFILTER,FY,FEAS,FCNT)     

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

    call uevaljac(n,x,ind,jcvar,jcval,jcnnz,flag)

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
