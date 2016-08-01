subroutine fird(n,x,l,u,me,mi,uevalf,uevalc,uevaljac,verbose,ftype, &
     epsfeas,epsopt,f,feas,fcnt,flag)

  ! This subroutine uses the Derivative-free Inexact Restoration
  ! Filter algorithm described in
  !
  ! "Global convergence of a derivative-free inexact restoration
  ! filter algorithm for nonlinear programming" P.S. Ferreira,
  ! E.W. Karas, M. Sachine and F.N.C. Sobral, Submitted, 2016.
  !
  ! This subroutine is simply a wrapper in order to make the code
  ! portable to be easily used by F77 subroutines.
  !
  ! N, input: number of variables
  !
  ! X: on input is the initial point provided by the user. On output
  ! is the solution found by FIRD
  !
  ! L(N), U(N), input: lower and upper bound on the variables,
  ! respectively
  !
  ! ME, MI, input: the number of equality and inequality
  ! constrains. Note that the constraints have to be coded in that
  ! order
  !
  ! UEVALF, UEVALC, UEVALJAC, input: subroutines provided by the user
  ! for evaluating the objective function, constraints and derivatives
  ! of the constraints, respectively. They have to follow the same
  ! interface defined by module 'userinterface.f90'
  !
  ! VERBOSE, input: if '.true.', FIRD outputs information at each
  ! step, otherwise no output information is given
  !
  ! FTYPE, input: type of filter. 1 means the 'flat' filter and 2
  ! means the 'slanting' filter
  !
  ! EPSFEAS, input: feasibility tolerance
  !
  ! EPSOPT, input: optimality tolerance
  !
  ! F, output: objective function value at the solution
  !
  ! FEAS, output: sup norm of feasibility at the solution
  !
  ! FCNT, output: total number of objective function evaluations
  !
  ! FLAG, output: reason for stopping FIRD. The FLAG output means
  !
  ! 0  - Solution was found
  ! 1  - Maximum number of OUTER iterations was reached
  ! 2  - Maximum number of obj. function evaluations was reached
  ! 3x - Failure in the restoration phase
  ! 4x - Failure in the optimization phase
  ! 5  - Memory error
  ! 6  - User-defined function error

  use userinterface
  use dfoirfilter

  implicit none

  ! SCALAR ARGUMENTS
  integer :: fcnt,flag,ftype,me,mi,n
  logical :: verbose
  real(8) :: epsfeas,epsopt,f,feas

  ! ARRAY ARGUMENTS
  real(8) :: l(n),u(n),x(n)

  ! EXTERNAL SUBROUTINES
  procedure(evalf)   :: uevalf
  procedure(evalc)   :: uevalc
  procedure(evaljac) :: uevaljac

  call dfoirfalg(n,x,l,u,me,mi,uevalf,uevalc,uevaljac,verbose,ftype, &
       epsfeas,epsopt,f,feas,fcnt,flag)

end subroutine fird
