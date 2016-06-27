!-------------------------------------------------------------------!
! SUBROUTINE QPSOLVER                                               !
!-------------------------------------------------------------------!
!
! Users should modify this subroutine in order to have his/her own
! optimization method.
!
! This subroutine uses the adapted trust region derivative-free TRDF
! algorithm for finding a point not in the filter that reduces the
! functional value at restored point.
!

subroutine qpsolver(n,y,l,u,me,mi,evalf_,evalc_,evallc_,evalljac_,  &
                    nf,alpha,ffilter,hfilter,outiter,epsfeas,epsopt,&
                    verbose,delta,fy,hynorm,rho,flag)

  use trdf_solver
  use userinterface

  implicit none

  ! SCALAR ARGUMENTS
  logical :: verbose
  integer :: flag,me,mi,n,nf,outiter
  real(8) :: alpha,delta,epsfeas,epsopt,fy,hynorm,rho

  ! ARRAY ARGUMENTS
  real(8) :: ffilter(nf),hfilter(nf),l(n),u(n),y(n)

  ! EXTERNAL SUBROUTINES
  procedure(evalf) :: evalf_
  procedure(evalc) :: evallc_,evalc_
  procedure(evaljac) :: evalljac_

  intent(in   ) :: alpha,epsfeas,epsopt,ffilter,hfilter,l,me,mi,n, &
       nf,outiter,u,verbose
  intent(out  ) :: hynorm,flag
  intent(inout) :: delta,fy,rho,y

  call solver(n,y,l,u,me,mi,evalf_,evalc_,evallc_,evalljac_, &
       nf,alpha,ffilter,hfilter,outiter,epsfeas, &
       epsopt,verbose,delta,fy,hynorm,rho,flag)

end subroutine qpsolver
