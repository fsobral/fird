module rinterface

  use userinterface

  abstract interface
     subroutine restoration(n,x,l,u,me,mi,uevalc,uevaljac, &
          epsfeas,verbose,infeas,flag)
       !
       ! Given an initial point, this subroutine returns an
       ! 'epsfeas'-feasible point.
       !
       ! N : integer,scalar,INPUT. Number of variables
       !                                                           
       ! X(N) : real(8),array,INPUT/OUTPUT. On INPUT is the current
       !        point, on OUTPUT is the (probably) feasible point
       !                                                           
       ! L(N) : real(8),array,INPUT. Lower bounds for the variables
       !                                                           
       ! U(N) : real(8),array,INPUT. Upper bounds for the variables
       !                                                           
       ! ME: integer,scalar,INPUT. Number of equality constraints
       !
       ! MI: integer,scalar,INPUT. Number of inequality constraints
       !
       ! UEVALC: procedure. User-defined subroutine for constraints
       !
       ! UEVALJAC: procedure. User-defined subroutine for the sparse
       !           Jacobian of the constraints
       !
       ! EPSFEAS: real,scalar,INPUT. Feasibility stopping criterium
       !
       ! VERBOSE: logical,scalar,INPUT. Toggle verbose mode ON (true)
       !          or OFF (false)
       !
       ! INFEAS: real(8),scalar,OUTPUT. Sup-norm of infeasibility
       !                                                           
       ! FLAG: integer,scalar,OUTPUT Status returned by the method:
       !       = 0 - solution found 
       !       > 0 - something happened
       
       ! SCALAR ARGUMENTS
       integer :: flag,me,mi,n
       real(8) :: epsfeas,infeas
       logical :: verbose

       ! ARRAY ARGUMENTS
       real(8) :: l(n),u(n),x(n)

       ! EXTERNAL SUBROUTINES
       procedure(evalc)   :: uevalc
       procedure(evaljac) :: uevaljac

       intent(in   ) :: epsfeas,l,me,mi,n,u,verbose
       intent(out  ) :: flag,infeas
       intent(inout) :: x
     end subroutine restoration
  end interface

contains

end module rinterface
