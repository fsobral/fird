module restoration

  ! EXTERNAL SUBROUTINES
  pointer :: evalc,evaljac

  ! WORK ARRAYS
  real(8), allocatable, dimension(:) :: xprev

  ! INTERFACES

  interface
     subroutine evalc(n,x,ind,c,flag)
       ! SCALAR ARGUMENTS
       integer :: flag,ind,n
       real(8) :: c
       ! ARRAY ARGUMENTS
       real(8) :: x(n)

       intent(in ) :: ind,n,x
       intent(out) :: c,flag
     end subroutine evalc

     subroutine evaljac(n,x,ind,jcvar,jcval,jcnnz,flag)
       ! SCALAR ARGUMENTS
       integer :: flag,ind,jcnnz,n
       ! ARRAY ARGUMENTS
       integer :: jcvar(n)
       real(8) :: jcval(n),x(n)

       intent(in ) :: ind,n,x
       intent(out) :: flag,jcnnz,jcval,jcvar
     end subroutine evaljac
  end interface

  private

contains

  !------------------------------------------------------------!
  ! SUBROUTINE RESTORATION                                     !
  !                                                            !
  ! This subroutine returns a feasible point 'x'. Up to now it !
  ! uses ALGENCAN to solve the 'projection'                    !
  !                                                            !
  !                min  ||x - x_k||_2^2                        !
  !                s.t. g(x) <= gamma_k                        !
  !                     l <= x <= u                            !
  !                                                            !
  ! ARGUMENTS                                                  !
  !                                                            !
  ! NOR: integer,scalar,INPUT                                  !
  !      Number of variables                                   !
  !                                                            !
  ! X(NOR) : real(8),array,INPUT/OUTPUT                        !
  !          On INPUT is the current point, on OUTPUT is the   !
  !          (probably) feasible point                         !
  !                                                            !
  ! L(NOR) : real(8),array,INPUT                               !
  !          Lower bounds for the variables                    !
  !                                                            !
  ! U(NOR) : real(8),array,INPUT                               !
  !          Upper bounds for the variables                    !
  !                                                            !
  ! MOR: integer,scalar,INPUT                                  !
  !      Number of constraints                                 !
  !                                                            !
  ! INFEAS: real(8),scalar,OUTPUT                              !
  !         2 norm of infeasibility                            !
  !                                                            !
  ! FLAG: integer,scalar,OUTPUT                                !
  !       Status returned by the method:                       !
  !       = 0 - solution found                                 !
  !       > 0 - something happend                              !
  !                                                            !
  !                                                            !
  ! This subroutine uses the following modules                 !
  !                                                            !
  ! - skinny (SK_PRINTE,SK_RESTTYPE)                           !
  ! - engdata (engXPrev,ENG_FORRES)                            !
  !                                                            !
  !------------------------------------------------------------!
  
  subroutine restore(n,x,l,u,me,mi,uevalc,uevaljac,epsfeas,verbose, &
                     infeas,flag)

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
    external uevalc,uevaljac

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
    real(8) :: aepsfeas,cnorm,efacc,efstain,eoacc,eostain,epsopt,nlpsupn,snorm

    ! LOCAL ARRAYS
    integer           :: i
    logical           :: coded(11),equatn(me + mi),linear(me + mi)
    real(8)           :: lambda(me + mi)
    character(len=15) :: strtmp
    character(len=80) :: specfnm,outputfnm,vparam(10)

    ! EXTERNAL SUBROUTINES
!    external :: r_evalf,r_evalg,r_evalh,r_evalc,r_evaljac,r_evalhc, &
!         r_evalfc,r_evalgjac,r_evalgjacp,r_evalhl,r_evalhlp

    evalc   => uevalc
    evaljac => uevaljac

    aepsfeas =  epsfeas
    epsopt   =  1.0D-08

    checkder = .false.

    m = me + mi

    hnnzmax  = n ** 2
    jcnnzmax =  n * m

    coded(1:11) = .false.
    coded(4: 5) =  .true.

    lambda(     1: m) =  0.0D0
    equatn(     1:me) = .false.
    equatn(me + 1: m) = .false.
    linear(     1: m) = .false.

    efstain   = sqrt(aepsfeas)
    eostain   = epsopt ** 1.5d0

    efacc     = sqrt(aepsfeas)
    eoacc     = sqrt(epsopt  )

    outputfnm = ''
    specfnm   = ''

    nvparam   = 0

    ! Optimize

    call algencan(r_evalf,r_evalg,r_evalh,r_evalc,r_evaljac,r_evalhc, &
         r_evalfc,r_evalgjac,r_evalgjacp,r_evalhl,r_evalhlp,jcnnzmax, &
         hnnzmax,aepsfeas,epsopt,efstain,eostain,efacc,eoacc,outputfnm, &
         specfnm,nvparam,vparam,n,x,l,u,m,lambda,equatn,linear,coded, &
         checkder,infeas,cnorm,snorm,nlpsupn,flag)


  end subroutine restore

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

    !  write(*,*) 'Entrou r_evalf',x

!!$    f = 0.0D0
!!$
!!$    do i = 1,n
!!$       f = f + (x(i) - xprev(i)) ** 2
!!$    end do
!!$
!!$    f = 5.0D-01 * f
!!$
!!$    flag = 0

    flag = -1

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

!!$    do i = 1,n
!!$       g(i) = x(i) - xprev(i)
!!$    end do
!!$
!!$    flag = 0

    flag = -1

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

    lmem = .false.
    flag = - 1

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

    ! INTERFACES
    interface
       subroutine evalc(n,x,ind,c,flag)
         integer :: flag,ind,n
         real(8) :: c
         real(8) :: x(n)

         intent(in ) :: ind,n,x
         intent(out) :: c,flag
       end subroutine evalc
    end interface

    ! LOCAL SCALARS
    integer         :: i
    real(8),pointer :: gamma(:)

    !  write(*,*) 'Entrou r_evalc'

    flag  =  0

    call evalc(n,x,ind,c,flag)

    if ( flag .ne. 0 ) then
       return
    end if

    !  write(*,*) 'Saiu r_evalc'

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

    ! INTERFACES
    interface
       subroutine evaljac(n,x,ind,jcvar,jcval,jcnnz,flag)
         integer :: flag,ind,jcnnz,jcvar(n),n
         real(8) :: jcval(n),x(n)

         intent(in ) :: ind,n,x
         intent(out) :: flag,jcnnz,jcval,jcvar
       end subroutine evaljac
    end interface

    ! LOCAL SCALARS
    integer :: i,m

    lmem = .false.
    flag =       0

    !  write(*,*) 'Entrou r_evaljac'

    call evaljac(n,x,ind,jcvar,jcval,jcnnz,flag)

    !  write(*,*) 'Saiu r_evaljac'

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

end module restoration
