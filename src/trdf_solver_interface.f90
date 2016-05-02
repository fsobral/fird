! Uses the adapted TRDF algorithm for solving the optimality phase

subroutine qpsolver(n,y,l,u,me,mi,uevalf,uevalc,uevaljac, &
                    nf,ffilter,hfilter,epsopt,fy,hynorm,flag)

  use trdf_solver
  use userinterface

  implicit none

  ! SCALAR ARGUMENTS
  integer :: flag,me,mi,n,nf
  real(8) :: epsopt,fy,hynorm

  ! ARRAY ARGUMENTS
  real(8) :: ffilter(nf),hfilter(nf),l(n),u(n),y(n)

  ! EXTERNAL SUBROUTINES
  procedure(evalf) :: uevalf
  procedure(evalc) :: uevalc
  procedure(evaljac) :: uevaljac

  

end subroutine qpsolver
