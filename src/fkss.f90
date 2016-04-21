subroutine fkss(n,x,l,u,me,mi,evalf,evalc,evaljac,verbose,epsfeas,epsopt,flag)

  use dfoirfilter

  implicit none

  ! SCALAR ARGUMENTS
  integer :: flag,me,mi,n
  logical :: verbose
  real(8) :: epsfeas,epsopt

  ! ARRAY ARGUMENTS
  real(8) :: l(n),u(n),x(n)

  ! EXTERNAL SUBROUTINES
  external :: evalf,evalc,evaljac

  call dfoirfalg(n,x,l,u,me,mi,evalf,evalc,evaljac,verbose,epsfeas,epsopt,flag)

end subroutine fkss
