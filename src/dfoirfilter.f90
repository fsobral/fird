module dfoirfilter

  implicit none

  ! PARAMETERS
  real(8), parameter :: BETA = 1.0D-4
  real(8), parameter :: ALPHA = 1.0D-1
  real(8), parameter :: DELMIN = 1.0D-12
  real(8), parameter :: MU = 1.0D-1
  real(8), parameter :: ETA = 2.5D-1
  ! Maximum number of filter elements
  integer, parameter :: MAXNF = 1000

  ! SCALARS
  integer :: ncev,nfev,njev

  ! EXTERNAL SUBROUTINES
  pointer :: evalf,evalc,evaljac

  ! INTERFACES

  interface
     subroutine evalf(n,x,f,flag)
       ! SCALAR ARGUMENTS
       integer :: flag,n
       real(8) :: f
       ! ARRAY ARGUMENTS
       real(8) :: x(n)

       intent(in ) :: n,x
       intent(out) :: f,flag
     end subroutine evalf

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

  public :: dfoirfalg

contains

  subroutine dfoirfalg(n,x,l,u,me,mi,evalf_,evalc_,evaljac_,verbose)

    use restoration

    implicit none

    ! SCALAR ARGUMENTS
    integer :: me,mi,n
    logical :: verbose

    ! ARRAY ARGUMENTS
    real(8) :: l(n),u(n),x(n)

    ! EXTERNAL SUBROUTINES
    external :: evalf_,evalc_,evaljac_

    ! LOCAL ARRAYS
    real(8) :: ffilter(MAXNF),hfilter(MACNF),rl(n),ru(n),y(n),z(n)

    ! LOCAL SCALARS
    integer :: flag,i,m,nf
    real(8) :: c,cfeas,hxnorm,hznorm,rinfeas

    nf = 1

    nfev = 0
    ncev = 0
    njev = 0

    evalf   => evalf_
    evalc   => evalc_
    evaljac => evaljac_

    m = me + mi

    ! Initialization
    ! TODO: check for errors when allocating

    hxnorm = 0.0D0
    do i = 1,me
       call uevalc(n,x,i,c,flag)
       hxnorm = max(hxnorm, abs(c))
    end do
    do i = me + 1,m
       call uevalc(n,x,i,c,flag)
       hxnorm = max(hxnorm, max(0.0D0, c))
    end do

    write(*,*) 'HXNORM=',hxnorm
    write(*,*) 'BETA=',BETA

    ! ----------------- !
    ! Feasibility phase !
    ! ----------------- !

    if ( hxnorm .gt. epsfeas ) then

       do i = 1,n
          rl(i) = max(l(i),x(i) - BETA * hxnorm)
          ru(i) = min(u(i),x(i) + BETA * hxnorm)
       end do

       cfeas = 9.5D-1 * (1.0D0 - ALPHA) * hxnorm

       call restore(n,x,rl,ru,me,mi,uevalc,uevaljac,cfeas,verbose,hznorm,flag)

       write(*,*) 'HZNORM=',hznorm

       ! TODO: Test alpha and filter conditions. In case of failure,
       ! decrease feasibility tolerance.

    end if

!!$    hznorm = 0.0D0
!!$    do i = 1,m
!!$       call uevalc(n,x,i,c,flag)
!!$       hznorm = max(hznorm, c)
!!$    end do

    ! Verify convergence conditions

    ! NON-EXECUTABLE STATEMENTS
    
600 FORMAT('Iteration',1X,I10,/)

  end subroutine dfoirfalg

  !----------------------------------------------------------!
  ! SUBROUTINE UEVALF                                        !
  !----------------------------------------------------------!

  subroutine uevalf(n,x,f,flag)
    
    ! SCALAR ARGUMENTS
    integer :: flag,n
    real(8) :: f
    
    ! ARRAY ARGUMENTS
    real(8) :: x(n)
    
    intent(in ) :: n,x
    intent(out) :: f,flag

    call evalf(n,x,f,flag)

    nfev = nfev + 1

  end subroutine uevalf

  !----------------------------------------------------------!
  ! SUBROUTINE UEVALC                                        !
  !----------------------------------------------------------!

  subroutine uevalc(n,x,ind,c,flag)

    ! SCALAR ARGUMENTS
    integer :: flag,ind,n
    real(8) :: c
    
    ! ARRAY ARGUMENTS
    real(8) :: x(n)
    
    intent(in ) :: ind,n,x
    intent(out) :: c,flag
    
    call evalc(n,x,ind,c,flag)
    
    ncev = ncev + 1
    
  end subroutine uevalc

  !----------------------------------------------------------!
  ! SUBROUTINE UEVALJAC                                      !
  !----------------------------------------------------------!
  
  subroutine uevaljac(n,x,ind,jcvar,jcval,jcnnz,flag)

    ! SCALAR ARGUMENTS
    integer :: flag,ind,jcnnz,n
    
    ! ARRAY ARGUMENTS
    integer :: jcvar(n)
    real(8) :: jcval(n),x(n)

    intent(in ) :: ind,n,x
    intent(out) :: flag,jcnnz,jcval,jcvar

    call evaljac(n,x,ind,jcvar,jcval,jcnnz,flag)

    njev = njev + 1

  end subroutine uevaljac


end module dfoirfilter
