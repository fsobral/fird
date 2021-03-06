module algencan_restoration

  use userinterface, only: evalc,evaljac

  implicit none

  ! EXTERNAL SUBROUTINES
  procedure(evalc)  , pointer :: evalc_
  procedure(evaljac), pointer :: evaljac_

  ! WORK ARRAYS
  real(8), dimension(:), pointer :: xprev

  private

  public :: algrestore

contains

  subroutine algrestore(n,x,l,u,me,mi,uevalc,uevaljac,epsfeas, &
       verbose,infeas,flag)

    implicit none

    ! SCALAR ARGUMENTS
    integer :: flag,me,mi,n
    real(8) :: epsfeas,infeas
    logical :: verbose

    ! ARRAY ARGUMENTS
    real(8) :: l(n),u(n),x(n)

    intent(in   ) :: epsfeas,l,me,mi,n,u,verbose
    intent(out  ) :: flag,infeas
    intent(inout) :: x

    ! EXTERNAL SUBROUTINES
    procedure(evalc)   :: uevalc
    procedure(evaljac) :: uevaljac

    ! INTERFACES
    interface
       subroutine algencan(fsub,gsub,hsub,csub,jacsub,hcsub,fcsub,gjacsub, &
            gjacpsub,hlsub,hlpsub,jcnnzmax,hnnzmax,epsfeas,epsopt,efstain, &
            eostain,efacc,eoacc,outputfnm,specfnm,nvparam,vparam,n,x,l,u,m, &
            lambda,equatn,linear,coded,checkder,fu,cnormu,snorm,nlpsupn,inform)

         logical,      intent(in)    :: checkder
         integer,      intent(in)    :: hnnzmax,jcnnzmax,m,n,nvparam
         integer,      intent(out)   :: inform
         real(kind=8), intent(inout) :: cnormu,efacc,efstain,eoacc,eostain, &
              epsfeas,epsopt,fu,nlpsupn,snorm

         character(len=80), intent(in)    :: specfnm,outputfnm,vparam(nvparam)
         logical,           intent(in)    :: coded(11),linear(m)
         logical,           intent(inout) :: equatn(m)
         real(kind=8),      intent(in)    :: l(n),u(n)
         real(kind=8),      intent(inout) :: lambda(m),x(n)

         external :: fsub,gsub,hsub,csub,jacsub,hcsub,fcsub,gjacsub,gjacpsub, &
              hlsub,hlpsub
       end subroutine algencan
    end interface

    ! LOCAL SCALARS
    logical :: checkder
    integer :: hnnzmax,jcnnzmax,m,nvparam
    real(8) :: aepsfeas,efacc,efstain,eoacc,eostain,epsopt,f,nlpsupn,snorm

    ! LOCAL ARRAYS
    integer           :: i
    logical           :: coded(11),equatn(me + mi),linear(me + mi)
    real(8)           :: lambda(me + mi),xp(n)
    character(len=15) :: strtmp
    character(len=80) :: specfnm,outputfnm,vparam(10)

    target :: xp

    evalc_   => uevalc
    evaljac_ => uevaljac

    xprev => xp
    do i = 1,n
       xprev(i) = x(i)
    end do

    aepsfeas =  epsfeas
    epsopt   =  1.0D-01

    checkder = .false.

    m = me + mi

    hnnzmax  = n
    jcnnzmax =  n * m

    coded(1:11) = .false.
    coded(1: 3) =  .true.
    coded(4: 5) =  .true.

    lambda(     1: m) =  0.0D0
    equatn(     1:me) =  .true.
    equatn(me + 1: m) = .false.
    linear(     1: m) = .false.

    efstain   = sqrt(aepsfeas)
    ! Disable early stopping criterium
    eostain   = - 1.0D0

    efacc     = sqrt(aepsfeas)
    eoacc     = sqrt(epsopt  )

    outputfnm = ''
    specfnm   = ''

    nvparam   = 0
!    vparam(1) = 'IGNORE-OBJECTIVE-FUNCTION'

    ! Optimize

    call algencan(r_evalf,r_evalg,r_evalh,r_evalc,r_evaljac,r_evalhc, &
         r_evalfc,r_evalgjac,r_evalgjacp,r_evalhl,r_evalhlp,jcnnzmax, &
         hnnzmax,aepsfeas,epsopt,efstain,eostain,efacc,eoacc,outputfnm, &
         specfnm,nvparam,vparam,n,x,l,u,m,lambda,equatn,linear,coded, &
         checkder,f,infeas,snorm,nlpsupn,flag)

  end subroutine algrestore

  ! ******************************************************************
  ! ******************************************************************

  subroutine r_evalf(n,x,f,flag)

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

    f = 0.0D0

    do i = 1,n
       f = f + (x(i) - xprev(i)) ** 2
    end do

    f = 5.0D-01 * f

    flag = 0

  end subroutine r_evalf

  ! ******************************************************************
  ! ******************************************************************

  subroutine r_evalg(n,x,g,flag)

    implicit none

    ! SCALAR ARGUMENTS
    integer :: flag,n

    ! ARRAY ARGUMENTS
    real(8) :: g(n),x(n)

    intent(in ) :: n,x
    intent(out) :: flag,g

    ! LOCAL SCALARS
    integer :: i

    do i = 1,n
       g(i) = x(i) - xprev(i)
    end do

    flag = 0

  end subroutine r_evalg

  subroutine r_evalh(n,x,hrow,hcol,hval,hnnz,lim,lmem,flag)

    implicit none

    ! SCALAR ARGUMENTS
    integer :: flag,hnnz,lim,n

    ! ARRAY ARGUMENTS
    logical :: lmem
    integer :: hcol(*),hrow(*)
    real(8) :: hval(*),x(n)

    intent(in ) :: lim,n,x
    intent(out) :: flag,lmem,hcol,hrow,hnnz,hval

    ! LOCAL SCALARS
    integer :: i

    lmem = .false.

    hnnz = n

    do i = 1,n
       hcol(i) = i
       hrow(i) = i
       hval(i) = 1.0D0
    end do

    flag = 0

  end subroutine r_evalh

  ! ******************************************************************
  ! ******************************************************************

  subroutine r_evalc(n,x,ind,c,flag)

    implicit none

    ! SCALAR ARGUMENTS
    integer :: flag,ind,n
    real(8) :: c

    ! ARRAY ARGUMENTS
    real(8) :: x(n)

    intent(in ) :: ind,n,x
    intent(out) :: c,flag

    ! LOCAL SCALARS
    integer         :: i
    real(8),pointer :: gamma(:)

    flag  =  0

    call evalc_(n,x,ind,c,flag)

    if ( flag .ne. 0 ) then
       return
    end if

  end subroutine r_evalc

  ! ******************************************************************
  ! ******************************************************************

  subroutine r_evaljac(n,x,ind,jcvar,jcval,jcnnz,lim,lmem,flag)

    implicit none

    ! SCALAR ARGUMENTS
    logical :: lmem
    integer :: flag,ind,jcnnz,lim,n

    ! ARRAY ARGUMENTS
    integer :: jcvar(n)
    real(8) :: jcval(n),x(n)

    intent(in ) :: ind,lim,n,x
    intent(out) :: flag,jcnnz,jcval,jcvar,lmem

    ! LOCAL SCALARS
    integer :: i,m

    lmem = .false.
    flag =       0

    call evaljac_(n,x,ind,jcvar,jcval,jcnnz,flag)

  end subroutine r_evaljac

  ! ******************************************************************
  ! ******************************************************************

  subroutine r_evalhc(n,x,ind,hcrow,hccol,hcval,hcnnz,lim,lmem,flag)

    implicit none

    ! SCALAR ARGUMENTS
    logical :: lmem
    integer :: flag,hcnnz,ind,lim,n

    ! ARRAY ARGUMENTS
    integer :: hccol(:),hcrow(:)
    real(8) :: hcval(:),x(n)

    lmem = .false.
    flag = - 1

  end subroutine r_evalhc

  ! ******************************************************************
  ! ******************************************************************

  subroutine r_evalfc(n,x,f,m,c,flag)

    implicit none

    ! SCALAR ARGUMENTS
    integer :: flag,m,n
    real(8) :: f

    ! ARRAY ARGUMENTS
    real(8) :: c(m),x(n)

    flag = - 1

  end subroutine r_evalfc

  ! ******************************************************************
  ! ******************************************************************

  subroutine r_evalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,lim,lmem,flag)

    implicit none

    ! SCALAR ARGUMENTS
    logical :: lmem
    integer :: flag,jcnnz,lim,m,n

    ! ARRAY ARGUMENTS
    integer :: jcfun(:),jcvar(:)
    real(8) :: g(n),jcval(:),x(n)

    lmem = .false.
    flag = - 1

  end subroutine r_evalgjac

  ! ******************************************************************
  ! ******************************************************************

  subroutine r_evalgjacp(n,x,g,m,p,q,work,gotj,flag)

    implicit none

    ! SCALAR ARGUMENTS
    logical,   intent(inout) :: gotj
    integer,   intent(in)    :: m,n
    integer,   intent(out)   :: flag
    character, intent(in)    :: work

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in)    :: x(n)
    real(kind=8), intent(inout) :: p(m),q(n)
    real(kind=8), intent(out)   :: g(n)

    flag = - 1

  end subroutine r_evalgjacp

  ! ******************************************************************
  ! ******************************************************************

  subroutine r_evalhl(n,x,m,lambda,sf,sc,hlrow,hlcol,hlval,hlnnz,lim,&
       lmem,flag)

    implicit none

    ! SCALAR ARGUMENTS
    logical :: lmem
    integer :: flag,hlnnz,lim,m,n
    real(8) :: sf

    ! ARRAY ARGUMENTS
    integer :: hlcol(:),hlrow(:)
    real(8) :: hlval(:),lambda(m),sc(m),x(n)

    lmem = .false.

    flag = - 1

  end subroutine r_evalhl

  ! ******************************************************************
  ! ******************************************************************

  subroutine r_evalhlp(n,x,m,lambda,sf,sc,p,hp,goth,flag)

    implicit none

    ! SCALAR ARGUMENTS
    logical ::goth
    integer :: flag,m,n
    real(8) :: sf

    ! ARRAY ARGUMENTS
    real(8) :: hp(n),lambda(m),p(n),sc(m),x(n)

    flag = - 1

  end subroutine r_evalhlp

end module algencan_restoration
