module dfoirfilter

  implicit none

  ! PARAMETERS
  real(8), parameter :: BETA = 1.0D-4

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
    real(8) :: rl(n),ru(n),y(n),z(n)

    ! LOCAL SCALARS
    integer :: flag,i,m
    real(8) :: beta,c,cfeas,hxnorm,hznorm,rinfeas

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
    do i = 1,m
       call uevalc(n,x,i,c,flag)
       hxnorm = max(hxnorm, c)
    end do

    ! ----------------- !
    ! Feasibility phase !
    ! ----------------- !

    do i = 1,n
       rl(i) = max(l(i),x(i) - BETA * hxnorm)
       ru(i) = min(u(i),x(i) + BETA * hxnorm)
    end do

    call restore(n,x,rl,ru,me,mi,evalc,evaljac,cfeas,verbose,BETA * hxnorm,flag)

    hznorm = 0.0D0
    do i = 1,m
       call uevalc(n,x,i,c,flag)
       hznorm = max(hznorm, c)
    end do

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
