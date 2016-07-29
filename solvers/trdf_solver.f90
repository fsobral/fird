module trdf_solver

  use userinterface

  implicit none


  ! PARAMETERS
  ! Maximum number of elements to print
  integer, parameter :: MAXXEL = 30
  ! Relation between AREd and PRED
  real(8), parameter :: ETA  = 1.0D-1
  real(8), parameter :: ETA1 = 7.0D-1
  ! Stationarity parameter
  real(8), parameter :: MU = 5.0D-1

  ! COMMON SCALARS

  ! Number of function evaluations
  integer :: IC
  ! Maximum number of function evaluations
  integer :: MAXIC
  ! Index in Y of the best point of model Q
  integer :: TBAR_

  real(8) :: VQUAD_A

  ! COMMON ARRAYS
  real(8), allocatable :: XBASE_A(:),GOPT_A(:),HQ_A(:)
  real(8), allocatable :: FF_(:),Q_(:),H_(:,:),Y_(:,:)

  ! GLOBAL USER-DEFINED SUBROUTINES
  procedure(evalf  ), pointer :: uevalf
  procedure(evalc  ), pointer :: uevallc,uevalc
  procedure(evaljac), pointer :: uevalljac

  private

  public :: trdfsolver

contains

  ! Uses the adapted TRDF algorithm for solving the optimality phase
  
  subroutine trdfsolver(n,y,l,u,me,mi,uevalf_,uevalc_,uevallc_,uevalljac_, &
       nf,alpha,ffilter,hfilter,filterTest,outiter,epsfeas,epsopt,verbose, &
       delta,fy,hynorm,rho,flag)

    use filters, only: absfilter

    implicit none

    ! SCALAR ARGUMENTS
    logical :: verbose
    integer :: flag,me,mi,n,nf,outiter
    real(8) :: alpha,delta,epsfeas,epsopt,fy,hynorm,rho

    ! ARRAY ARGUMENTS
    real(8) :: ffilter(nf),hfilter(nf),l(n),u(n),y(n)

    ! EXTERNAL SUBROUTINES
    procedure(evalf  )   :: uevalf_
    procedure(evalc  )   :: uevallc_,uevalc_
    procedure(evaljac)   :: uevalljac_
    procedure(absfilter) :: filterTest

    ! LOCAL SCALARS
    integer :: i,m,maxfcnt,npt,fcnt
    real(8) :: rbeg,rend,xeps

    ! LOCAL ARRAYS
    logical :: ccoded(2),equatn(me + mi),linear(me + mi)

    uevalf    => uevalf_
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

    rbeg = rho

    rend = epsopt

    xeps = 1.0D-08

    call TRDFSUB(N,NPT,Y,L,U,M,EQUATN,LINEAR,CCODED,MAXFCNT,RBEG, &
         REND,XEPS,VERBOSE,NF,ALPHA,FFILTER,HFILTER,FILTERTEST,   &
         OUTITER,DELTA,EPSFEAS,FY,HYNORM,FCNT,RHO,FLAG)

  end subroutine trdfsolver

  ! ******************************************************************
  ! ******************************************************************

  SUBROUTINE TRDFSUB(N,NPT,X,XL,XU,M,EQUATN,LINEAR,CCODED,MAXFCNT, &
       RBEG,REND,XEPS,OUTPUT,NF,ALPHA,FFILTER,HFILTER,FILTERTEST,  &
       OUTITER,DELTA,EPSFEAS,F,FEAS,FCNT,RHO,FLAG)

    ! This subroutine is strongly based on the implementation of the
    ! Derivative-free Trust-region algorithm for constrained
    ! optimization described in
    !
    !     P. D. Conejo, E. W. Karas, L. G. Pedroso, "A trust-region
    !     derivative-free algorithm for constrained optimization".
    !     Optimization Methods & Software, v. 30 (6), p. 1126-1145,
    !     2015.
    !
    ! For more information on configurations, recommended values of
    ! the parameters and examples for the original method, please read
    ! the documentation provided with the method.
    !
    ! The modifications performed in the algorithm were based in the
    ! "Optimality Phase" described in
    ! 
    ! "Global convergence of a derivative-free inexact restoration
    ! filter algorithm for nonlinear programming", P.S. Ferreira,
    ! E.W. Karas, M. Sachine and F.N.C. Sobral, Submitted, 2016.

    use filters, only: absfilter

    IMPLICIT NONE

    ! SCALAR ARGUMENTS
    logical :: OUTPUT
    integer :: flag,m,maxfcnt,N,NF,NPT,FCNT,outiter
    real(8) :: ALPHA,F,DELTA,EPSFEAS,FEAS,RBEG,REND,RHO,XEPS

    ! ARRAY ARGUMENTS
    REAL(8) :: FFILTER(NF),HFILTER(NF),X(N),XL(N),XU(N)
    logical :: ccoded(2),equatn(m),linear(m)

    ! EXTERNAL SUBROUTINES
    procedure(absfilter) :: filterTest

    intent(in   ) :: m,maxfcnt,n,npt,rbeg,rend,xeps,xl,xu,ccoded, &
         equatn,linear,alpha,nf,ffilter,hfilter,epsfeas, &
         outiter
    intent(out  ) :: feas,fcnt,flag,rho
    intent(inout) :: delta,f,x

    ! LOCAL ARRAYS
    REAL(8) :: D(n),XNOVO(n),VETOR1(NPT + N + 1),Z(n)

    ! LOCAL SCALARS
    logical :: forbidden
    integer :: i,j,k,kn,previt,t
    real(8) :: alfa,beta,c,cnorm,distsq,dsq,gamma, &
         mindelta,rhobeg,rhoend,sigm,sum,tau,tempofinal, &
         tempoinicial,fz,distz,qx,qz

    ! -------------- !
    ! Initialization !
    ! -------------- !

    flag = 0

    ! Set the user-defined functions

    FZ      = F
    F       = 1.0D+300
    IC      = 0
    MAXIC   = maxfcnt ! MAXIMUM NUMBER OF FUNCTION EVALUATIONS
    RHOBEG  = RBEG
    RHOEND  = REND
    RHO     = RHOBEG
    DELTA   = MAX(DELTA, RHO) ! Correcting delta if necessary
    GAMMA   = 0.1D0

    DO I=1,N
       Z(I)     = X(I)
       XNOVO(I) = X(I)
    END DO

    if ( outiter .gt. 1 ) then

       ! Repurposing the interpolation points from the previous
       ! iterations

       GOTO 11

    else

       ! First outer iteration. Allocates the whole structure.
       ! TODO: Maybe we have to deallocate it?
       ! TODO: Test allocation errors
       allocate(Y_(NPT,N),FF_(NPT),Q_(1+N+N*(N+1)/2),H_(NPT+N+1,NPT+N+1))
       allocate(XBASE_A(N),GOPT_A(N),HQ_A(N * (N + 1) / 2))

    end if

    DO I=1,N
       XBASE_A(I) = X(I) 
    END DO

    GO TO 5

4   CONTINUE

    DO I=1,N
       X(I)       = XNOVO(I)
       XBASE_A(I) = XNOVO(I)
    END DO

5   continue

    ! Since we are rebuilding, z_k is the center of the model
    tbar_ = 1

    CALL  PRIMEIROMODELO1 (N,Z,FZ,Q_,H_,NPT,RHO,Y_,FF_,FLAG) 

    IF ( OUTPUT ) WRITE(*,1002) RHO,DELTA,FF_(1),IC,MIN(N,MAXXEL), &
         (X(I), I=1,MIN(N,MAXXEL))
    IF ( FLAG .NE. 0 ) GOTO 31

    FZ = FF_(1)         

11  CALL SUBPROBLEMA(N,NPT,Q_,DELTA,D,X,XL,XU,DSQ,M,EQUATN,LINEAR, &
         CCODED,XEPS,FLAG)

    ! Actually, we should sum Q(1) to have the correct value
    ! of the model at the points. But, in order to calculate
    ! the difference, we can omit Q(1), since it will be
    ! canceled.

    call mevalf(N,d,QX,flag)

    do i = 1,n
       d(i) = z(i) - xbase_a(i)
    end do

    call mevalf(N,d,QZ,flag)

    IF ( OUTPUT ) WRITE(*,1003) RHO,DELTA,QX + Q_(1),FZ,IC

    IF ( FLAG .NE. 0 ) THEN

       IF ( RHO .LE. RHOEND ) THEN
          GOTO 31
       ELSE
          RHO = GAMMA * RHO
          GOTO 4
       END IF

    END IF

    ! Evaluate the distance of the solution to z^k

    DISTSQ = (10.0D0 * RHO) ** 2.0D0                 

    DISTZ = 0.0D0
    DO I = 1,N
       DISTZ = MAX(DISTZ, ABS(X(I) - Z(I)))
    END DO

    ! If the distance to z^k is small, then verify if it is not time
    ! to stop the whole algorithm. Otherwise, try to build the model
    ! in a smaller radius.

    IF ( DISTZ .LT. MU * RHO ) THEN 

       FEAS = 0.0D0
       do I = 1,M
          CALL UEVALC(N,Z,I,C,FLAG)
          IF ( FLAG .NE. 0 ) GOTO 31          
          IF ( EQUATN(I) ) THEN
             FEAS = MAX(FEAS,ABS(C))
          ELSE
             FEAS = MAX(FEAS,MAX(0.0D0,C))
          END IF
       end do

       ! If z^k is feasible, then stop then exit (and stop the
       ! algorithm)
       
       IF (RHO .LE. RHOEND .AND. FEAS .LE. EPSFEAS) GO TO 31

       ! Test the distance of the points in Y to the feasibility point
       ! z^k

       KN = 0
       DO K = 1,NPT
          SUM = 0.0D0
          DO J = 1,N
             SUM = SUM + (Y_(K,J) - Z(J)) ** 2.0D0
          END DO
          IF ( SUM .GT. DISTSQ ) THEN
             KN = K
             exit
          END IF
       END DO

       IF ( KN .EQ. 0 ) RHO = GAMMA * RHO      

       GO TO 4
    END IF

    CALL CALFUN(N,X,F,FLAG)
    IF ( FLAG .NE. 0 ) GOTO 31

    IF ( OUTPUT ) WRITE(*,1006) F

    ! CHOOSE WHO LEAVE Y CALCULATING THE VALUE OF SIGMA. THE VARIABLE
    ! 'T' IS CHOOSEN FOR DEFINING WHO LEAVES. WHEN THE SET Y IS BEING
    ! REPURPOSED, THEN TBAR_ CONTAINS THE BEST POINT OF THE PREVIOUS
    ! ITERATION.    

    t = tbar_
    CALL SIGMA(H_,N,NPT,Y_,X,VETOR1,SIGM,ALFA,BETA,TAU,t,DELTA)

    ! IF ANY REDUCTION IN F, PUT X IN INTERPOLATION SET.

    IF ( F .LE. FZ ) THEN  
       IF ( OUTPUT ) WRITE(*,1005) t
       DO I = 1,N
          Y_(t,I) = X(I) 
       END DO
    ELSE
       GO TO 23
    END IF

    ! UPDATE H              
    CALL INVERSAH(H_, N, NPT,VETOR1, SIGM, t, ALFA, BETA,TAU)

    CALL ATUALIZAQ(H_, N, NPT, Q_, DELTA, Y_, X, F, t) 

    ! Test sufficient reduction of the objective function and filter

23  IF ( F .LE. FZ + ETA * (QX - QZ) ) THEN

       FEAS = 0.0D0

       do I = 1,M

          CALL UEVALC(N,X,I,C,FLAG)
          IF ( FLAG .NE. 0 ) GOTO 31          

          IF ( EQUATN(I) ) THEN
             FEAS = MAX(FEAS,ABS(C))
          ELSE
             FEAS = MAX(FEAS,MAX(0.0D0,C))
          END IF

       end do

       ! Filter test

       forbidden = filterTest(F,FEAS,ALPHA,nf,ffilter,hfilter)

       if ( .not. forbidden ) then

          ! Increase TR radius in case of high decrease

          IF ( F - FZ .GE. ETA1 * (QX - QZ) .OR. &
               DISTZ .LT. DELTA ) THEN
             DELTA = DELTA  
          ELSE
             DELTA = DELTA + DELTA  
          END IF

          FZ    = F
          tbar_ = t

          DO I=1, N
             XNOVO(I) = X(I)
          END DO

          FLAG = 0

          GO TO 31 
       end if

    END IF

    IF (sigm .le. 0d0) THEN
       DELTA = 5.0D-1 * DELTA
       GOTO 4
    END IF

    IF (IC == MAXIC) THEN
       FLAG = 3
       GO TO 31
    END IF

    ! If RHO is too low and no point is found not belonging to the
    ! filter, then declare failure in the optimization phase.

    IF ( RHO .LE. 1.0D-15 ) THEN
       FLAG = 4
       GOTO 31
    END IF

    ! Poisedness test

    KN = 0
    DO K = 1,NPT
       SUM = 0.0D0
       DO J = 1,N
          SUM = SUM + (Y_(K,J) - Z(J)) ** 2.0D0
       END DO
       IF ( SUM .GT. DISTSQ ) THEN
          KN = K
          exit
       END IF
    END DO

    IF ( KN .GT. 0 ) THEN
       DELTA = 5.0D-1 * DELTA
       GOTO 4
    END IF

    DELTA = RHO  
    RHO   = GAMMA * RHO      

    GO TO 11    

    ! Finish iterations and return

31  continue           

    if ( OUTPUT .and. flag .eq.  0 ) write(*,1020)
    if ( OUTPUT .and. flag .eq. -1 ) write(*,1021)
    if ( OUTPUT .and. flag .eq.  2 ) write(*,1022)
    if ( OUTPUT .and. flag .eq.  3 ) write(*,1023) MAXIC
    if ( OUTPUT .and. flag .eq.  4 ) write(*,1024)

    F = FZ

    do i = 1,n
       x(i) = XNOVO(i)
    end do

    FEAS = 0.0D0
    do I = 1,M

       CALL UEVALC(N,X,I,C,FLAG)
       IF ( FLAG .NE. 0 ) GOTO 31

       IF ( EQUATN(I) ) THEN
          FEAS = MAX( FEAS,ABS(C) )
       ELSE
          FEAS = MAX( FEAS,MAX(0.0D0,C) )
       END IF
    end do

    FCNT = IC

    ! NON-EXECUTABLE STATEMENTS

1000 FORMAT(/,'PHASE 0',/,7('-'),/,/,'FEASIBILITY =',36X,D23.8,/, &
         'NEW POINT',/,3(1X,D23.8))
1001 FORMAT(/,'PHASE 1',/,7('-'),/)
1002 FORMAT(3X,'(RE)BUILDING MODEL from scratch.',/, &
         3X,5X,'RHO =',45X,D12.5,/,               &
         3X,5X,'Delta =',43X,D12.5,/,             &
         3X,5X,'Objective function =',19X,D23.8,/,&
         3X,5X,'Function evaluations =',30X,I10,/,&
         3X,5X,'Current model center (first ',I3, &
         ' elements)',/,3X,6X,3(1X,D21.8))
1003 FORMAT(/,3X,'SOLVED TR SUBPROBLEM.',/,                 &
         3X,5X,'RHO =',45X,D12.5,/,                      &
         3X,5X,'Delta =',43X,D12.5,/,                    &
         3X,5X,'Model value =',26X,D23.8,/,              &
         3X,5X,'Objective function (at Z) =',12X,D23.8,/,&
         3X,5X,'Function evaluations =',30X,I10)
1004 FORMAT(5X,'Objective function =',24X,D23.8)
1005 FORMAT(/,'REMOVING sampling point',1X,I4,'.')
1006 FORMAT(5X,'Objective function (at trial point) =',&
         7X,D23.8)


1020 FORMAT(/,3X,'Solution was found!',/)
1021 FORMAT(/,3X,'Flag -1: Error while evaluating functions.',/)
1022 FORMAT(/,3X,'Flag 2: Error in the internal solver.',/)
1023 FORMAT(/,3X,'Flag 3: Reached the maximum of',1X,I10,1X, &
         'function evaluations.',/)
1024 FORMAT(/,3X,'Flag 4: Rho smaller than tolerance.',/)

3000 FORMAT(/,'Welcome to TRDF Algorithm!',/,                   &
         'This algorithm was based on paper',/,                &
         'P.D. Conejo, E.W. Karas, and L.G. Pedroso',/,        &
         '"A trust-region derivative-free algorithm for',/,    &
         'constrained problems", Optimization Methods ',/,     &
         'and Software, v. 30 (6), p. 1126-1145, 2015]',/)
  END SUBROUTINE TRDFSUB

  ! ******************************************************************
  ! ******************************************************************

  ! ********************************  FIRST MODEL  *******************************
  SUBROUTINE  PRIMEIROMODELO1 (N,X,FX,Q,H,NPT,DELTA,Y,FF,FLAG)

    implicit none

    ! SCALAR ARGUMENTS
    integer :: n,npt,flag
    real(8) :: delta,fx

    ! ARRAY ARGUMENTS
    real(8) :: Q(*), FF(*), x(n), H(NPT+N+1,NPT+N+1),YY(N)

    ! NPT IS THE NUMBER INTERPOLATION POINTS.
    ! Y IS THE INTERPOLATION SET.
    ! FF KEEP IN VALUES OF F IN Y. 
    ! Q STORES THE HESSIAN AND GRADIENT OF MODEL.
    ! H  IS THE INVERSE ASSOCIATED WITH SYSTEM.
    ! YY STORES EACH VECTOR OF Y. 
    ! HQ IS THE HESSIAN IN MATRIX FORMAT.  

    ! LOCAL SCALARS
    integer :: i,ii,j,k
    real(8) :: ACUMULADOR

    ! LOCAL ARRAYS
    real(8) ::  E(N+1,NPT),OMEGA(NPT,NPT),Y(NPT,N),GAMMA(N+1,N+1), &
         Z(NPT,NPT-N-1),HQ(n,n),FFTEMP(npt)
    INTEGER :: IP(npt), IQ(npt)

    DO I=1, 1+N+N*(N+1)/2
       Q(I)=0.0D0
    END DO ! START Q
    DO I=1, NPT+N+1
       DO J=1, NPT+N+1
          H(I,J)=0.0D0
       END DO
    END DO

    ! START HQ
    DO I=1, N
       DO J=1, N
          HQ(I,J)=0.0D0
       END DO
    END DO

    ! NPT2N = 2*N+1
    DO I=1, N
       Y(1,I)=X(I) 
    END DO
    DO I=1, N
       DO J=1, N 
          Y(I+1,J)= X(J)
          Y(I+1,I )= X(I )+DELTA 
          Y(I+N+1,J)= X(J)
          Y(I+N+1,I)= X(I )-DELTA                 
       END DO
    END DO

    ! It is possible to use the value of the objective function at the
    ! center
    FF(1) = FX
    DO I = 2,2 * N + 1 
       DO J = 1,N
          YY(J) = Y(I,J) 
       END DO
       CALL CALFUN(N,YY,FF(I),FLAG)

       IF ( FLAG .NE. 0 ) RETURN

    END DO

    ! ******************* MODEL ***************************

    Q(1)=FF(1)          
    ! DEFINE THE GRADIENT GOPT OF THE FIRST MODEL
    DO I=1, N 
       Q(I+1)=(1D0/(2*DELTA)) * (FF(I+1)-FF(I+1+N))      
    END DO
    ! DEFINE THE DIAGONAL OF THE HESSIAN MODEL   
    DO I=1, N 
       HQ(I,I)=(1D0/(DELTA**2))*(FF(I+1)+FF(I+1+N)-2*FF(1)) 
    END DO

    ! NPT >= 2N+1       

    IF (NPT .GT. 2*N+1) THEN             
       ! SETTING THE POITS M-2N+1
       IF (NPT .GT. 2*N+1) THEN
          DO J= 2*N+2, NPT
             IF (J .LE. 3*N+1)  IP(J) = J-2*N-1
             IF (J .GE. 3*N+2)  IP(J) = IP(J-N) 
          END DO
       END IF

       II =1 
       DO J= 2*N+2, NPT 
          IF (IP(J) + II .LE. N) THEN
             IQ(J) = IP(J) + II
          ELSE 
             IQ(J) = IP(J) + II - N 
          END IF
          IF (MOD(IP(J)  ,N) .EQ. 0) II = II+1 
       END DO

       ! OBTAIN THE POINTS Y OF 2N+1 TO NPT.
       DO I=2*N+2, NPT
          DO J= 1, N 
             Y(I,J) = Y(IP(I)+1, J) + Y(IQ(I)+1, J) - Y(1,J)                
          END DO
       END DO
       DO I=2*N+2, NPT
          DO J=1, N
             YY(J) = Y(I,J) 
          END DO
          CALL CALFUN(N,YY,FF(I),FLAG)   

          IF ( FLAG .NE. 0 ) RETURN

       END DO

       ! DEFINE OTHERS INPUTS OF HESSIAN FOR OVER 2N+1.     
       DO J=2*N+2, NPT
          HQ(IP(J),IQ(J))=(1.0D0/(DELTA**2))*(FF(J)-FF(IP(J)+1) &
               -FF(IQ(J)+1)+ FF(1))  
          HQ(IQ(J),IP(J)) =  HQ(IP(J),IQ(J)) 
       END DO
    END IF
    K=1
    DO I=1 ,N
       DO J=1, I
          Q(K+N+1) = HQ(I,J)
          K = K+1
       END DO
    END DO

    ! UPDATE THE GRADIENT AND HESSIAN FOR THE FIRST MODEL. 

    DO I=1, N
       GOPT_A(I) = Q(I+1)       
    END DO
    DO J=1, (N+1)*N/2                  
       HQ_A(J) =  Q(1+N+J)  

    END DO

    ! ******************* FIRST INVERSE***************************

    DO I=1, NPT
       DO J=1, NPT
          OMEGA(I,J)=0D0
       END DO
    END DO
    DO I=1, N+1
       DO J=1, N+1
          GAMMA(I,J)=0D0
       END DO
    END DO
    DO I=1, NPT  
       DO J=1, NPT-N -1
          Z(I,J)=0D0
       END DO
    END DO
    ! MATRIX E( N+1 X NPT)
    DO I=1, N+1
       DO J=1,NPT  
          E(I,J)=0D0
       END DO
    END DO
    E(1,1)=1D0     
    DO I=2, N+1
       E(I,I)= 1.0D0/(2.0D0*DELTA)
       E(I,I+N)= - 1.0D0/(2.0D0*DELTA)
    END DO
    ! MATRIX Z(NPT X NPT-N-1)            
    DO I=1, N 
       Z(1,I)= -SQRT(2.0D0)/(DELTA**2)
       Z(I+1,I)=  SQRT(2.0D0)/(2.0D0*DELTA**2)
       Z(N+I+1,I)= SQRT(2.0D0)/(2.0D0*DELTA**2)
    END DO

    ! THE NEW INVERSE FOR MORE OF 2N+1 POINTS IN Y
    IF (NPT .GT.  2*N+1) THEN             
       DO I=N+1, NPT-N-1 
          Z(1,I)= 1.0D0/(DELTA**2)
          Z(N+I+1, I) = Z(1,I) 
          Z(IP(N+I+1)+1,I) = -1.0D0/(DELTA**2)
          Z(IQ(N+I+1)+1,I) = -1.0D0/(DELTA**2)        
       END DO
    END IF

    ! MULTIPLYING ZZ^T FOR DETERMINE OMEGA          
    ACUMULADOR=0D0
    DO I=1, NPT   
       DO K=1, NPT              
          DO J=1, NPT-N-1
             ACUMULADOR =  ACUMULADOR + Z(I,J)*Z(K,J)
          END DO
          OMEGA(I,K) = ACUMULADOR               
          ACUMULADOR = 0D0              
       END DO
    END DO

    ! THE MATRIX INVERSE H     
    DO I=1, NPT
       DO J=1, NPT
          H(I,J)=OMEGA(I,J)         
       END DO
    END DO

    ! THE N+1 LINES OF H                 
    DO I=NPT+1, NPT+N+1
       DO J= 1, NPT 
          H(I,J) = E(I-NPT,J)         
       END DO
    END DO
    ! THE N+1 COLUMNS OF H 

    DO I=1, NPT
       DO J= NPT+1, NPT+N+1
          H(I,J) = H(J,I)        
       END DO
    END DO

    RETURN
  END SUBROUTINE PRIMEIROMODELO1

  ! ******************************************************************
  ! ******************************************************************

  SUBROUTINE SIGMA(H,N,NPT,Y,X,VETOR1,SIGM,ALFA,BETA,TAU,IT,DELTA)

    implicit none

    ! SCALAR ARGUMENTS
    integer :: IT,N,NPT
    real(8) :: alfa,beta,delta,sigm,tau

    ! ARRAY ARGUMENTS
    real(8) :: X(*), VETOR1(*), H(NPT+N+1,NPT+N+1), Y(NPT,N)

    ! WW STORAGE THE VETOR1 IN W TO PRODUCE ALFA BETA TAU HOW IN
    ! DEFINITION.  SIGMA = ALFA BETA + TAU**2. ALFA = ET^T H ET,
    ! NAMELY, HTT (T = IT) IT INDICATE THAT THE VETOR1 IN POSITION TWO
    ! (SECOND LINE OF Y) LEAVE OF Y CHOOSE THE LARGEST SIGMA (AND
    ! THEREFORE IT) UNLESS DE CURRENT POINT (IT CURRENT)

    ! LOCAL SCALARS
    integer :: i,IAUXILIAR,ITT,j,k,kkk
    real(8) :: AGUARD,CONT,SIGMI

    ! LOCAL ARRAYS
    real(8) :: WW(NPT+N+1,1),AUXILIAR(4)           

    ITT=IT    ! STORAGE THE POSITION OF BEST ITERATING YET.          
    IT =1
    SIGMI = -1.0D100
    WW(NPT+1,1) = 1.0D0
    DO I=1, N
       WW(I+NPT+1, 1) =  X(I)- XBASE_A(I)
    END DO
    DO I=1, NPT
       CONT = 0.0D0
       DO K=1, N         
          CONT = CONT +  (Y(I,K) - XBASE_A(K)) * (X(K) - XBASE_A(K))          
       END DO
       WW(I, 1) = 0.5D0*CONT**2
    END DO

    ! CALCULATING ALL SIGMA AND CHOOSE IT FOR THE GREATER SIGMA. MULTIPLY T-th LINE OF H FOR W FOR OBTAIN TAU. 
    DO KKK = 1, NPT-1            
       CONT=0.0D0
       DO I=1, NPT+ N +1
          CONT = CONT +  H(IT, I) * WW(I,1)
       END DO
       TAU = CONT
       ! CALCULUS OF BETA = 0.5||X^+-XBASE||-WHW                  
       DO I=1, NPT+N+1
          CONT = 0.0D0 
          DO K=1, NPT+N+1
             CONT = CONT + H(I,K) * WW(K,1)
             VETOR1(I) = CONT                  
          END DO
       END DO

       CONT = 0.0D0
       DO I=1,  NPT+N+1
          CONT = CONT+WW(I,1)*VETOR1(I)
          AGUARD = CONT
       END DO

       ! CALCULUS  OF X-XB^4
       CONT = 0.0D0
       DO I=1, N
          CONT = CONT + (X(I)-XBASE_A(I))**2
       END DO
       BETA =   0.5D0 * CONT**2 - AGUARD
       ALFA = H(IT,IT) 
       SIGM = ALFA * BETA + TAU**2  
       IF (SIGM .GE. SIGMI .AND. IT .NE. ITT) THEN                      
          SIGMI = SIGM
          IAUXILIAR  =  IT                
          AUXILIAR(1) = ALFA
          AUXILIAR(2) = BETA
          AUXILIAR(3) = TAU
          AUXILIAR(4) = SIGM                           
       END IF
       IT = IT + 1
    END DO

    ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

    IT   = IAUXILIAR  
    ALFA = AUXILIAR(1) 
    BETA = AUXILIAR(2)
    TAU  = AUXILIAR(3) 
    SIGM = AUXILIAR(4)

    CONT=0.0D0
    DO I=1, NPT+ N +1
       CONT = CONT +  H(IT, I) * WW(I,1)
    END DO
    TAU = CONT

    RETURN
  END SUBROUTINE SIGMA

  ! ******************************************************************
  ! ******************************************************************

  ! **************** UPDAT THE INVERSE H ******************************
  SUBROUTINE INVERSAH(H, N, NPT,VETOR1,SIGM,IT,ALFA,BETA,TAU)
!!$    IMPLICIT REAL*8 (A-H,O-Z)

    ! SCALAR ARGUMENTS
    integer :: IT,N,NPT
    real(8) :: alfa,beta,sigm,tau

    ! ARRAY ARGUMENTS
    real(8) :: VETOR1(*),H(NPT+N+1,NPT+N+1)

    !#include "tr_params.par"

    ! ALFA*(E-MM) * (E-MM)'- BETA* H * E * E'*H+TAU*H*E*(E-MM)'+ (E-MM)* E'* H
    ! MM = H*WW THAT IS STORED IN VETOR1 

    ! LOCAL SCALARS
    integer :: i,j

    ! LOCAL ARRAYS
    real(8) :: P1(NPT+N+1,NPT+N+1),P2(NPT+N+1,NPT+N+1),P3(NPT+N+1,NPT+N+1)

    VETOR1(IT) = VETOR1(IT)-1.0D0              
    DO I=1, N+NPT+1
       DO J=1, N+NPT+1
          P1(I,J)= VETOR1(I) * VETOR1(J)
          P2(I,J)= H(I, IT) * H(J, IT)
          P3(I,J)= (H(IT, I)*(-VETOR1(J)))+(H(IT,J)*(-VETOR1(I)))
       END DO
    END DO
    DO I=1, N+NPT+1
       DO J=1, N+NPT+1
          if (sigm .eq. 0d0) return
          H(I,J)=H(I,J)+(1/SIGM)*(ALFA*P1(I,J)-BETA*P2(I,J)+TAU*P3(I,J)) 
       END DO
    END DO
    RETURN
  END SUBROUTINE INVERSAH

  ! ******************************************************************
  ! ******************************************************************

  SUBROUTINE  SUBPROBLEMA(N, NPT, Q, DELTA, D, X, XL, XU, DSQ, M, &
       EQUATN, LINEAR, CCODED, XEPS, FLAG)

    implicit none

    ! SCALAR ARGUMENTS
    integer :: flag,m,N,NPT
    real(8) :: delta,dsq,xeps

    ! ARRAY ARGUMENTS
    real(8) :: Q(*), X(*), XL(*),XU(*),D(n)
    logical :: ccoded(2),equatn(m), linear(m)

    ! LOCAL SCALARS
    integer :: i
    real(8) :: cnorm,f,qxCur,qxNew,sum

    ! LOCAL ARRAYS
    real(8) ::  XANTIGO(n),L(N),H(NPT + N + 1,NPT + N + 1),U(N) 

    DO I = 1,N
       L(I) = DMAX1(XL(I) - XBASE_A(I),X(I) - XBASE_A(I) - DELTA)
       U(I) = DMIN1(XU(I) - XBASE_A(I),X(I) - XBASE_A(I) + DELTA)
       XANTIGO(I) = X(I)       
    END DO

    DO I = 1,N    
       D(I) = X(I) - XBASE_A(I)
    END DO

    CALL MEVALF(N,D,F,FLAG)

    IF ( FLAG .NE. 0 ) RETURN

    qxCur = F + Q(1)

    CALL SOLVER(N, L, U, D, M, EQUATN, LINEAR, CCODED,      &
         mevalf, mevalg, mevalh, mevalc, mevaljac, mevalhc, &
         .false., XEPS, CNORM, FLAG)

    IF ( FLAG .NE. 0 ) RETURN

    ! CALCULUS THE STEP LENGTH 
    SUM = 0.0D0
    DO I=1, N
       SUM = SUM + (X(I)-(D(I)+XBASE_A(I)))**2
    END DO
    DSQ = SUM
    DO I=1, N
       X(I)= D(I) +  XBASE_A(I)   
    END DO

    CALL MEVALF(N,D,F,FLAG)

    IF ( FLAG .NE. 0 ) RETURN

    VQUAD_A = F + Q(1) ! MODEL IN XNOVO 

    qxNew = VQUAD_A

    IF ( qxNew - qxCur .GE. 0.D0 .or. cnorm .gt. xeps )  THEN
       DO I = 1,N
          X(I) = XANTIGO(I) 
       END DO
       DSQ = 0D0
    END IF

    RETURN
  END SUBROUTINE SUBPROBLEMA

  ! ******************************************************************
  ! ****************************************************************** 

  SUBROUTINE  ATUALIZAQ(H, N, NPT, Q, DELTA, Y, X, F, IT)

    implicit none

    ! SCALAR ARGUMENTS
    integer :: IT,N,NPT
    real(8) :: delta,f

    ! ARRAY ARGUMENTS
    REAL(8) :: Q(*),X(*),Y(NPT,N),H(NPT+N+1,NPT+N+1)

    ! LOCAL SCALARS
    integer :: i,ii,j,jj,k

    ! LOCAL ARRAYS
    real(8) :: VETORAUX(n),DD(N,N),QQ((N+1)*(N+2)/2),TEMP(1+N+NPT)

    DO I=1, 1+N+NPT
       TEMP(I) =  (F - VQUAD_A)* H(I, IT) ! IS LAMBDA 
    END DO
    DO I=1, N 
       DO J=1, N
          DD(I,J)=0.0D0
       END DO
    END DO
    ! M=M+     LAMBCG(J)*((Y(:,J)-XB') * ( Y(: ,J)-XB')')  
    DO I=1, NPT
       DO JJ=1, N
          VETORAUX(JJ) =  Y(I,JJ)-XBASE_A(JJ)  
       END DO
       DO K=1, N                    
          DO J=1, N
             DD(K,J) = DD(K,J)+ TEMP(I)* VETORAUX(K) * VETORAUX(J)  
             ! DDEH IS THE GRADIENT OF QUADRATIC D                               
          END DO
       END DO
    END DO
    ! N+1 FIRST ELEMENTS OF THE QQ (PARAMETER OF THE NEW MODEL)
    QQ(1) = TEMP(NPT+1)
    DO I=2, N+1 
       QQ(I) = TEMP(NPT+I)
    END DO
    ! PUT IN THE QQ FOR SYMMETRY OF DD
    II=1
    J=1
    ! DO WHILE (II .LE.  (N+1)*N/2 ) 
    DO WHILE (II .LE. N) 
       DO I=1, II 
          QQ(J+1+N) = DD(II,I) 
          J=J+1
       END DO
       II = II + 1
    END  DO
    ! ADD Q TO QQ.  Q IS THE OLD MODEL, QQ IS THE MODEL D, AND Q + D = Q+    
    DO I=1, 1+N+ (N+1)*N/2
       Q(I) = Q(I) + QQ(I)
    END DO
    ! UPDAT THE MODEL IN THE FILE INIP FOR COMMON FUNCTION 
    DO I=1, N      
       GOPT_A(I) = Q(I+1)
    END DO
    DO J=1, (N+1)*N/2 
       HQ_A(J) =  Q(1+N+J) 
    END DO
    RETURN
  END SUBROUTINE  ATUALIZAQ

  ! ******************************************************************
  ! ******************************************************************

  subroutine mvv(v1,v2,n,gradd)

    implicit none

    ! multiplica vetor por vetor

    ! SCALAR ARGUMENTS
    real(8) :: gradd

    ! ARRAY ARGUMENTS
    real(8) :: v1(n),v2(n)

    ! LOCAL SCALARS
    real(8) :: soma
    integer :: j,n

    soma = 0.0D0
    do j = 1,n
       soma  = soma + v1(j) * v2(j)
       gradd = soma
    end do
    return
  end subroutine mvv

  ! ******************************************************************
  ! ******************************************************************

  subroutine mmv(HQ,S,n,v)

    implicit none

    ! multiplica matriz simetrica (dada como vetor) por vetor
    ! estava com hs, e troquei para hss para nao atualizar hs desnec

    ! ARRAY ARGUMENTS
    real(8) :: S(N), HQ(N * (N + 1) / 2),v(n)

    ! LOCAL ARRAYS
    real(8) :: HSS(n)

    ! LOCAL SCALARS
    integer :: n,i,IH,j

    IH = 0
    DO J = 1,N
       HSS(J) = 0.0D0
       DO I = 1,J
          IH = IH + 1
          IF (I .LT. J) HSS(J) = HSS(J) + HQ(IH) * S(I)
          HSS(I) = HSS(I) + HQ(IH) * S(J)
          v(I) = HSS(I)
       end DO
    end DO
    return
  end subroutine mmv

  ! ******************************************************************
  ! ******************************************************************

  subroutine calfun(n,x,f,flag)

    implicit none

    ! SCALAR ARGUMENTS
    integer :: flag,n
    real(8) :: f

    ! ARRAY ARGUMENTS
    real(8) :: x(n)

    flag = 0

    ! Returns flag = 3 if reached the maximum of function evaluations
    if ( IC .eq. MAXIC ) then 
       flag = 3
       return
    end if

    ! TODO: add 'flag' to calobjf
    call uevalf(n,x,f,flag)

    IC = IC + 1

  end subroutine calfun

  subroutine mevalf(n,x,f,flag)

    implicit none

    ! SCALAR ARGUMENTS
    integer :: flag,n
    real(8) :: f

    ! ARRAY ARGUMENTS
    real(8) :: x(n)

    ! LOCAL ARRAYS
    real(8) :: hqd(n)

    ! LOCAL SCALARS
    real(8) :: gradd,dhqd

    flag = 0

    ! avalia o modelo quadratico f, no ponto d

    ! definindo o modelo quadratico f
    call mvv(GOPT_a,x, n, gradd)
    call mmv(HQ_a, x, n, hqd)
    call mvv(x, hqd, n, dhqd)

    f = gradd + dhqd / 2.0D0

  end subroutine mevalf

  !     ******************************************************************
  !     ******************************************************************

  subroutine mevalg(n,x,g,flag)

    implicit none

    ! SCALAR ARGUMENTS
    integer :: flag,n

    ! ARRAY ARGUMENTS
    real(8) :: g(n),x(n)

    ! LOCAL ARRAYS
    real(8) :: hqd(n)

    ! LOCAL SCALARS
    integer :: i

    flag = 0

    ! avalia o gradiente do modelo quadratico g, no ponto d
    call mmv(HQ_a, x, n, hqd)
    do i=1, n
       g(i)=  GOPT_a(i) + hqd(i)
    end do

  end subroutine mevalg

  ! ******************************************************************
  ! ******************************************************************

  subroutine mevalh(n,x,hrow,hcol,hval,hnnz,lim,lmem,flag)

    implicit none

    ! SCALAR ARGUMENTS
    logical :: lmem
    integer :: flag,n,hnnz,lim

    ! ARRAY ARGUMENTS
    integer :: hcol(lim),hrow(lim)
    real(8) :: hval(lim),x(n)

    ! LOCAL SCALARS
    integer :: i,iii,j

    flag = 0
    lmem = .false.

    hnnz = (N + 1) * N / 2

    ! cria hlin e hcol com as coordenadas da triangular inferior

    III=1
    DO WHILE (III .LE.  (N+1)*N/2 )
       do I=1, N
          DO J= 1, I
             hrow(III) = I
             hcol(III) = J
             III= III + 1
          end do
       end do
    END DO

    do i=1,(N+1)*N/2
       hval(i) = HQ_a(i)
    end do

  end subroutine mevalh

  !     ******************************************************************
  !     ******************************************************************

  subroutine mevalc(n,x,ind,c,flag)

    implicit none

    ! SCALAR ARGUMENTS
    integer :: ind,flag,n
    real(8) :: c

    ! ARRAY ARGUMENTS
    real(8) :: x(n)

    ! LOCAL ARRAYS
    real(8) :: XA(n)

    ! LOCAL SCALARS
    integer :: i

    flag = -1

    DO I = 1,N
       XA(I) = X(I) + XBASE_A(I)
    END DO

    call uevallc(n,XA,ind,c,flag)

  end subroutine mevalc

  ! ******************************************************************
  ! ******************************************************************

  subroutine mevaljac(n,x,ind,jcvar,jcval,jcnnz,lim,lmem,flag)

    implicit none

    ! SCALAR ARGUMENTS
    logical :: lmem
    integer :: flag,ind,jcnnz,lim,n

    ! ARRAY ARGUMENTS
    integer :: jcvar(lim)
    real(8) :: x(n),jcval(lim)

    ! LOCAL ARRAYS
    real(8) :: XA(n)

    ! LOCAL SCALARS
    integer :: i

    flag = -1

    lmem = .false.

    DO I = 1,N
       XA(I) = X(I) + XBASE_A(I)
    END DO

    call uevalljac(n,XA,ind,jcvar,jcval,jcnnz,flag)

  end subroutine mevaljac

  ! ******************************************************************
  ! ******************************************************************

  subroutine mevalhc(n,x,ind,hcrow,hccol,hcval,hcnnz,lim,lmem,flag)

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

  end subroutine mevalhc
  
end module trdf_solver
