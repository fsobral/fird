module dfoirfilter

  use userinterface

  implicit none

  ! PARAMETERS
  real(8), parameter :: BETA = 1.0D+03
  real(8), parameter :: ALPHA = 1.0D-1
  real(8), parameter :: MU = 1.0D-1
  real(8), parameter :: ETA = 2.5D-1
  ! Restoration reduction factor
  real(8), parameter :: RESRFAC = 9.5D-01
  ! Maximum number of iterations
  integer, parameter :: MAXITER = 10000
  ! Maximum number of printing elements
  integer, parameter :: MAXNEL  = 20
  ! Minimum trust region radius
  real(8), parameter :: DELMIN = 1.0D-30
  ! Initial value of RHO
  real(8), parameter :: RHOINI = 1.0D0

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

  subroutine dfoirfalg(n,x,l,u,me,mi,evalf_,evalc_,evaljac_, &
       verbose,epsfeas,epsopt,f,feas,fcnt,flag)

    use restoration

    implicit none

    ! SCALAR ARGUMENTS
    integer :: fcnt,flag,me,mi,n
    logical :: verbose
    real(8) :: epsfeas,epsopt,f,feas

    ! ARRAY ARGUMENTS
    real(8) :: l(n),u(n),x(n)

    ! EXTERNAL SUBROUTINES
    external :: evalf_,evalc_,evaljac_

    ! LOCAL ARRAYS
    real(8) :: ffilter(MAXITER),hfilter(MAXITER),rl(n),ru(n),xp(n)

    ! LOCAL SCALARS
    logical :: isforb,isalph,isbeta
    integer :: i,j,k,iter,jcnnz,m,nf
    real(8) :: c,currfeas,delta,dzynorm,dxznorm,fx,fy,fz,hxnorm,&
         hznorm,hynorm,rho

    iter = 1

    nfev = 0
    ncev = 0
    njev = 0

    uevalf   => evalf_
    uevalc   => evalc_
    uevaljac => evaljac_

    m = me + mi

    rho = RHOINI

    delta = rho

    ! Initialization
    ! TODO: check for errors when allocating

    allocate(linrhs(m), linpos(m  + 1), linvar(m * n), linval(m * n))

    hxnorm = evalinfeas(n,x,l,u,me,mi,flag)

    currfeas = hxnorm

    dzynorm = 1.0D+20

    call aevalf(n,x,fx,flag)

    ! Filter initialization

    nf = 1

    do while ( .true. )

       ! Creates the temporary filter \hat F
       ! TODO: Change this!
       ffilter(nf) = fx
       hfilter(nf) = hxnorm

       if ( verbose ) write(*,900) iter,fx,hxnorm
       if ( verbose ) write(*,904) min(MAXNEL,n),(x(i), i = 1,min(MAXNEL,n))

       ! ----------------- !
       ! Feasibility phase !
       ! ----------------- !

       xp(1:n) = x(1:n)

       ! Update current feasibility
       
!       currfeas = max(epsfeas, min((1.0D0 - ALPHA) * hxnorm, hxnorm * dzynorm))
       currfeas = max(epsfeas, (1.0D0 - ALPHA) * hxnorm)

       if ( verbose ) write(*,905)

       isforb = .true.
       isalph = .false.
       isbeta = .false.

       fz     = fx
       hznorm = hxnorm
       dxznorm = 0.0D0

       do while ( hxnorm .gt. epsfeas .and.           &
                  ( isforb .or. ( .not. isalph ) ) )

          isforb = .false.
          isalph = .true.
          isbeta = .true.

!!$          do i = 1,n
!!$             rl(i) = max(l(i),x(i) - BETA * hxnorm)
!!$             ru(i) = min(u(i),x(i) + BETA * hxnorm)
!!$          end do

          call restore(n,x,l,u,me,mi,aevalc,aevaljac,currfeas,verbose, &
               hznorm,flag)

          if ( flag .ne. 0 ) then
             Write(*,*) 'Problems in the restoration solver.'
             exit
          end if
          
          if ( hznorm .gt. currfeas ) then
             write(*,*) 'Unable to reach desired feasibility.'
             flag = 4
             exit
          end if

          call aevalf(n,x,fz,flag)

          ! Test alpha and filter conditions. In case of failure,
          ! decrease feasibility tolerance.

          do i = 1,nf
             if ( hznorm .ge. (1.0D0 - ALPHA) * hfilter(i) .and. &
                  fz .ge. ffilter(i) - ALPHA * hfilter(i) ) then
                isforb = .true.
                exit
             end if
          end do

          dxznorm = evalDist(n,xp,x)

          if ( hznorm .ge. (1 - ALPHA) * hxnorm ) isalph = .false.
          if ( dxznorm .gt. BETA * hxnorm ) isbeta = .false.

          if ( verbose ) write(*,906) isforb,isalph,isbeta
          
          if ( isforb .or. .not. isalph ) then 

             currfeas = max(epsfeas, RESRFAC * currfeas)
             if ( verbose ) write(*,907) currfeas

          end if

       end do

       if ( flag .ne. 0 ) exit

       if ( verbose ) WRITE(*,908) fz,hznorm,min(n,MAXNEL), &
                      (x(i), i = 1,min(MAXNEL,n))

       ! ---------------- !
       ! Optimality phase !
       ! ---------------- !

       ! Construct the linear system

       if ( verbose ) write(*,909)

       k = 1

       do i = 1,m

          linrhs(i) = 0.0D0

          ! Just for inequalities
          if ( i .gt. me ) then
             call aevalc(n,x,i,c,flag)
             linrhs(i) = max(0.0D0, c) - c
          end if

          call aevaljac(n,x,i,linvar(k),linval(k),jcnnz,flag)
          linpos(i) = k

          do j = k,k + jcnnz - 1
             linrhs(i) = linrhs(i) + linval(j) * x(linvar(j))
          end do

          k = k + jcnnz
       end do

       linpos(m + 1) = k

       ! Solve the subproblem

       do i = 1,n
          xp(i) = x(i)
       end do

       delta = max(DELMIN, rho, delta)

       rho = max(10.0D0 * epsopt, delta)
!       rho = max(10.0D0 * epsopt, min(delta, dzynorm))

       call qpsolver(n,x,l,u,me,mi,aevalf,aevalc,levalc,levaljac, &
            nf,ALPHA,ffilter,hfilter,currfeas,epsopt,verbose,delta, &
            fy,hynorm,rho,flag)

       dzynorm = evalDist(n,xp,x)

       if ( flag .ne. 0 ) write(*,912) flag

       if ( verbose ) write(*,910) fy,hynorm,delta,rho,dzynorm, &
            nfev,min(n,MAXNEL),(x(i), i=1,min(n,MAXNEL))

       ! ------------- !
       ! Filter Update !
       ! ------------- !

       if ( verbose ) write(*,903)
       
       if ( fy .ge. fx ) then

          ! This is an h-iteration
          ! The temporary filter turns effective
          nf = nf + 1
          if ( verbose ) write(*,901) fx, hxnorm

       else

          if ( verbose ) write(*,902)

       end if

       ! -------------------------- !
       ! Prepare for next iteration !
       ! -------------------------- !

       fx = fy

       hxnorm = hynorm
       
       iter = iter + 1

       ! Verify convergence conditions

       if ( hynorm .le. epsfeas .and. &
            dzynorm .le. epsopt ) then
          flag = 0
          exit
       end if

       if ( iter .gt. MAXITER ) then
          flag = 2
          exit
       end if

    end do

    deallocate(linrhs,linpos,linvar,linval)

    if ( verbose ) write(*,911) fx,hxnorm,nfev,min(n,MAXNEL), &
         (x(i),i=1,min(n,MAXNEL))

    f = fy

    feas = hynorm

    fcnt = nfev

    ! NON-EXECUTABLE STATEMENTS
    
900 FORMAT(/,70('-'),/,'Iteration',I61,/,70('-'),/,/,'F(X) = ', &
           40X,1PD23.8,/,'H(X) = ',40X,1PD23.8)
901 FORMAT('H-iteration: the pair (',1PD17.8,',',1PD17.8,') was added.',/)
902 FORMAT('F-iteration.',/)
903 FORMAT(/,'Filter update',/,13('-'))
904 FORMAT('Current point (first',1X,I5,' elements):',/, &
           (2X,4(1X,1PD16.8)))
905 FORMAT(/,'Restoration Phase',/,17('-'),/)
906 FORMAT(3X,'Forbidden?',56X,L,/,3X,'Alpha?',60X,L,/,3X,'Beta?',61X,L,/)
907 FORMAT(3X,'Updated feasibility requirement to',27X,E9.1,/)
908 FORMAT(3X,'F(Z) =',45X,1PD16.8,/,3X,'H(Z) =',45X,D16.8,/,3X, &
           'Restored point (first',1X,I5,' elements):',/, &
           (3X,16X,3(1X,1PD16.8)))
909 FORMAT(/,'Optimization Phase',/,18('-'),/)
910 FORMAT(3X,'F(Z+D) =',43X,1PD16.8,/,3X,'H(Z+D) =',43X,D16.8,/,   &
           3X,'DELTA =',44X,1PD16.8,/,3X,'RHO =',46X,1PD16.8,/,     &
           3X,'||D|| =',44X,D16.8,/,3X,'Number of f evaluat. =',1X, &
           I44,/,3X,'Optimized point (first',1X, &
           I5,' elements):',/,(3X,16X,3(1X,1PD16.8)))
911 FORMAT(/,70('-'),/,'Final Iteration',/,70('-'),/,'F(X) =',48X, &
           1PD16.8,/,'H(X) =',48X,D16.8,/,   &
           'Function Evaluations =',28X,I20,/,'Solution (first',1X, &
           I5,' elements):',/,(2X,4(1X,1PD16.8)))
912 FORMAT(3X,'WARNING! Solver return FLAG ',I3,'!',/)

  end subroutine dfoirfalg

  function evaldist(n,x,y)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: n

    ! ARRAY ARGUMENTS
    real(8), intent(in) :: x(n),y(n)

    real(8) :: evaldist

    ! LOCAL SCALARS
    integer :: i

    evaldist = 0.0D0
    do i = 1,n
       evaldist = max(evaldist, abs(x(i) - y(i)))
    end do

  end function evaldist

  !----------------------------------------------------------!
  ! FUNCTION EVALINFEAS                                      !
  !                                                          !
  ! This function evaluates the sup-norm of the              !
  ! infeasibilities.                                         !
  !                                                          !
  !----------------------------------------------------------!

  function evalinfeas(n,x,l,u,me,mi,flag)

    implicit none

    integer :: flag,me,mi,n
    real(8) :: l(n),u(n),x(n)

    real(8) :: evalinfeas

    real(8) :: c
    integer :: i

    evalinfeas = 0.0D0
    do i = 1,n
       evalinfeas = max(evalinfeas, x(i) - u(i), l(i) - x(i))
    end do
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

    c = - linrhs(ind)

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
