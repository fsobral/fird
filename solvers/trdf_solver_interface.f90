! Uses the adapted TRDF algorithm for solving the optimality phase

subroutine qpsolver(n,y,l,u,me,mi,evalf_,evalc_,evallc_,evalljac_, &
                    nf,alpha,ffilter,hfilter,epsfeas,epsopt,       &
                    verbose,fy,hynorm,flag)

  use trdf_solver
  use userinterface

  implicit none

  ! SCALAR ARGUMENTS
  logical :: verbose
  integer :: flag,me,mi,n,nf
  real(8) :: alpha,epsfeas,epsopt,fy,hynorm

  ! ARRAY ARGUMENTS
  real(8) :: ffilter(nf),hfilter(nf),l(n),u(n),y(n)

  ! EXTERNAL SUBROUTINES
  procedure(evalf) :: evalf_
  procedure(evalc) :: evallc_,evalc_
  procedure(evaljac) :: evalljac_

  call solver(n,y,l,u,me,mi,evalf_,evalc_,evallc_,evalljac_, &
       nf,alpha,ffilter,hfilter,epsfeas,epsopt,verbose,fy,   &
       hynorm,flag)

end subroutine qpsolver
