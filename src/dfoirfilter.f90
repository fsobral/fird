module dfoirfilter

  use userinterface

  implicit none

  ! PARAMETERS
  real(8), parameter :: BETA = 1.0D-4
  real(8), parameter :: ALPHA = 1.0D-1
  real(8), parameter :: DELMIN = 1.0D-12
  real(8), parameter :: MU = 1.0D-1
  real(8), parameter :: ETA = 2.5D-1
  ! Maximum number of iterations
  integer, parameter :: MAXITER = 3

  ! ARRAYS
  integer, allocatable :: linpos(:),linvar(:)
  real(8), allocatable :: linrhs(:),linval(:)

  ! SCALARS
  integer :: ncev,nfev,njev

  ! EXTERNAL SUBROUTINES
  procedure(evalf  ), pointer :: uevalf
  procedure(evalc  ), pointer :: uevalc
  procedure(evaljac), pointer :: uevaljac

  private

  public :: dfoirfalg

contains

  subroutine dfoirfalg(n,x,l,u,me,mi,evalf_,evalc_,evaljac_,verbose,epsfeas,epsopt,flag)

    use restoration

    implicit none

    ! SCALAR ARGUMENTS
    integer :: flag,me,mi,n
    logical :: verbose
    real(8) :: epsfeas,epsopt

    ! ARRAY ARGUMENTS
    real(8) :: l(n),u(n),x(n)

    ! EXTERNAL SUBROUTINES
    external :: evalf_,evalc_,evaljac_

    ! LOCAL ARRAYS
    real(8) :: ffilter(MAXITER),hfilter(MAXITER),rl(n),ru(n),y(n),z(n),zd(n)

    ! LOCAL SCALARS
    integer :: i,j,k,iter,jcnnz,m,nf
    real(8) :: c,cfeas,dnorm,fx,fy,fz,hxnorm,hznorm,hynorm

    iter = 1

    nfev = 0
    ncev = 0
    njev = 0

    uevalf   => evalf_
    uevalc   => evalc_
    uevaljac => evaljac_

    m = me + mi

    ! Initialization
    ! TODO: check for errors when allocating

    allocate(linrhs(m), linpos(m  + 1), linvar(m * n), linval(m * n))

    hxnorm = evalinfeas(n,x,me,mi,flag)

    call aevalf(n,x,fx,flag)

    ! Filter initialization

    nf = 1

    do while ( .true. )

       ! Creates the temporary filter \hat F
       ! TODO: Change this!
       ffilter(nf) = fx
       hfilter(nf) = hxnorm

       if ( verbose ) write(*,900) iter

       write(*,*) 'FX=',fx
       write(*,*) 'HXNORM=',hxnorm
       write(*,*) 'BETA=',BETA
       write(*,904) (x(i), i = 1,n)

       ! ----------------- !
       ! Feasibility phase !
       ! ----------------- !

       cfeas = (1.0D0 - ALPHA) * hxnorm

       if ( hxnorm .gt. epsfeas ) then

!!$          do i = 1,n
!!$             rl(i) = max(l(i),x(i) - BETA * hxnorm)
!!$             ru(i) = min(u(i),x(i) + BETA * hxnorm)
!!$          end do

010       cfeas = 9.5D-1 * cfeas

          call restore(n,x,l,u,me,mi,uevalc,uevaljac,cfeas,verbose,hznorm,flag)

          call aevalf(n,x,fz,flag)

          write(*,*) 'HZNORM=',hznorm

          ! TODO: Test alpha and filter conditions. In case of failure,
          ! decrease feasibility tolerance.

          do i = 1,nf
             if ( hznorm .ge. (1.0D0 - ALPHA) * hfilter(i) .and. &
                  fz .ge. ffilter(i) - ALPHA * hfilter(i) ) goto 010
          end do

       end if

       write(*,903) (x(i), i = 1,n)

       ! ---------------- !
       ! Optimality phase !
       ! ---------------- !

       ! Construct the linear system

       k = 1

       do i = 1,m
          call aevalc(n,x,i,linrhs(i),flag)
          linrhs(i) = - linrhs(i)

          call aevaljac(n,x,i,linvar(k),linval(k),jcnnz,flag)
          linpos(i) = k

          k = k + jcnnz
       end do

       linpos(m + 1) = k

       ! Solve the subproblem

       ! Here we can use 'x', since its old value is not necessary anymore.
       do i = 1,n
          x(i) = z(i)
       end do

       call qpsolver(n,x,l,u,me,mi,aevalf,levalc,levaljac, &
            nf,ALPHA,ffilter,hfilter,cfeas,fy,flag)

       ! Verify convergence conditions

       dnorm = 0.0D0
       do i = 1,n
          dnorm = max(dnorm, abs(z(i) - y(i)))
       end do
       
       hynorm = evalinfeas(n,x,me,mi,flag)

!!$       if ( hynorm .le. epsfeas .and. &
!!$            dnorm .le. epsopt ) then
!!$          flag = 0
!!$          exit
!!$       end if
!!$
!!$       if ( dnorm .le. epsopt ) then
!!$          flag = 1
!!$          exit
!!$       end if

       ! ------------- !
       ! Filter Update !
       ! ------------- !

       if ( fy .ge. fx ) then

          ! This is an h-iteration
          ! The temporary filter turns effective
          nf = nf + 1
          if ( verbose ) write(*,901) fx, hxnorm

       else

          if ( verbose ) write(*,902)

       end if

       ! Prepare for next iteration !

       fx = fy

       hxnorm = evalinfeas(n,x,me,mi,flag)
       
       iter = iter + 1

       if ( iter .gt. MAXITER ) then
          flag = 2
          exit
       end if

    end do

    deallocate(linrhs,linpos,linvar,linval)


    ! NON-EXECUTABLE STATEMENTS
    
900 FORMAT('Iteration',1X,I10,/)
901 FORMAT(1X,'H-iteration: the pair (',E9.1,',',E9.1,') was added.',/)
902 FORMAT(1X,'F-iteration.',/)
903 FORMAT(1X,'Restored point:',/6X,3(1X,D21.8))
904 FORMAT(1X,'Current point:',/6X,3(1X,D21.8))

  end subroutine dfoirfalg

  !----------------------------------------------------------!
  ! FUNCTION EVALINFEAS                                      !
  !                                                          !
  ! This function evaluates the sup-norm of the              !
  ! infeasibilities.                                         !
  !                                                          !
  !----------------------------------------------------------!

  function evalinfeas(n,x,me,mi,flag)

    integer :: flag,me,mi,n
    real(8) :: x(n)

    real(8) :: evalinfeas

    real(8) :: c
    integer :: i

    evalinfeas = 0.0D0
    do i = 1,me
       call uevalc(n,x,i,c,flag)
       evalinfeas = max(evalinfeas, abs(c))
    end do
    do i = me + 1,me + mi
       call uevalc(n,x,i,c,flag)
       evalinfeas = max(evalinfeas, max(0.0D0, c))
    end do

  end function evalinfeas

  !----------------------------------------------------------!
  ! SUBROUTINE AEVALF                                        !
  !----------------------------------------------------------!

  subroutine aevalf(n,x,f,flag)
    
    ! SCALAR ARGUMENTS
    integer :: flag,n
    real(8) :: f
    
    ! ARRAY ARGUMENTS
    real(8) :: x(n)
    
    intent(in ) :: n,x
    intent(out) :: f,flag

    call uevalf(n,x,f,flag)

    nfev = nfev + 1

  end subroutine aevalf

  !----------------------------------------------------------!
  ! SUBROUTINE AEVALC                                        !
  !----------------------------------------------------------!

  subroutine aevalc(n,x,ind,c,flag)

    ! SCALAR ARGUMENTS
    integer :: flag,ind,n
    real(8) :: c
    
    ! ARRAY ARGUMENTS
    real(8) :: x(n)
    
    intent(in ) :: ind,n,x
    intent(out) :: c,flag
    
    call uevalc(n,x,ind,c,flag)
    
    ncev = ncev + 1
    
  end subroutine aevalc

  !----------------------------------------------------------!
  ! SUBROUTINE AEVALJAC                                      !
  !----------------------------------------------------------!
  
  subroutine aevaljac(n,x,ind,jcvar,jcval,jcnnz,flag)

    ! SCALAR ARGUMENTS
    integer :: flag,ind,jcnnz,n
    
    ! ARRAY ARGUMENTS
    integer :: jcvar(n)
    real(8) :: jcval(n),x(n)

    intent(in ) :: ind,n,x
    intent(out) :: flag,jcnnz,jcval,jcvar

    call uevaljac(n,x,ind,jcvar,jcval,jcnnz,flag)

    njev = njev + 1

  end subroutine aevaljac

  !----------------------------------------------------------!
  ! SUBROUTINE LEVALC                                        !
  !                                                          !
  ! This subroutine evaluates the linearized constraints of  !
  ! the optimality subproblem.                               !
  !                                                          !
  !----------------------------------------------------------!

  subroutine levalc(n,x,ind,c,flag)

    ! SCALAR ARGUMENTS
    integer :: flag,ind,n
    real(8) :: c
    
    ! ARRAY ARGUMENTS
    real(8) :: x(n)
    
    intent(in ) :: ind,n,x
    intent(out) :: c,flag
    
    ! LOCAL SCALARS
    integer :: end,i,start

    start = linpos(ind)
    end   = linpos(ind + 1) - 1

    c = linrhs(ind)

    do i = start,end
       c = c + linval(i) * x(linvar(i))
    end do

    flag = 0
    
  end subroutine levalc
  
  !----------------------------------------------------------!
  ! SUBROUTINE LEVALJAC                                      !
  !                                                          !
  ! This subroutine evaluates the gradients of the           !
  ! linearized constraints of the optimality subproblem.     !
  !----------------------------------------------------------!
  
  subroutine levaljac(n,x,ind,jcvar,jcval,jcnnz,flag)

    ! SCALAR ARGUMENTS
    integer :: flag,ind,jcnnz,n
    
    ! ARRAY ARGUMENTS
    integer :: jcvar(n)
    real(8) :: jcval(n),x(n)

    intent(in ) :: ind,n,x
    intent(out) :: flag,jcnnz,jcval,jcvar

    ! LOCAL SCALARS
    integer :: end,i,start

    start = linpos(ind)
    end   = linpos(ind + 1) - 1

    jcnnz = end - start + 1

    jcvar(1:jcnnz) = linvar(start:end)
    jcval(1:jcnnz) = linval(start:end)

    flag = 0

  end subroutine levaljac

end module dfoirfilter
