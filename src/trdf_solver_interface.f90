! Uses the adapted TRDF algorithm for solving the optimality phase

subroutine qpsolver(n,y,l,u,me,mi,evalf_,evallc_,evalljac_,evalc_, &
                    nf,alpha,ffilter,hfilter,epsfeas,epsopt,fy,flag)

  use trdf_solver
  use userinterface

  implicit none

  ! SCALAR ARGUMENTS
  integer :: flag,me,mi,n,nf
  real(8) :: alpha,epsfeas,epsopt,fy

  ! ARRAY ARGUMENTS
  real(8) :: ffilter(nf),hfilter(nf),l(n),u(n),y(n)

  ! EXTERNAL SUBROUTINES
  procedure(evalf) :: evalf_
  procedure(evalc) :: evallc_,evalc_
  procedure(evaljac) :: evalljac_

  call solver(n,y,l,u,me,mi,evalf_,evallc_,evalljac_,evalc_, &
       nf,alpha,ffilter,hfilter,epsfeas,epsopt,fy,flag)

end subroutine qpsolver
