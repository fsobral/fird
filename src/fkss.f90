subroutine fkss(n,x,l,u,me,mi,evalf,evalc,evaljac,verbose,flag)

  use dfoirfilter

  implicit none

  ! SCALAR ARGUMENTS
  integer :: flag,me,mi,n
  logical :: verbose

  ! ARRAY ARGUMENTS
  real(8) :: l(n),u(n),x(n)

  ! EXTERNAL SUBROUTINES
  external :: evalf,evalc,evaljac

  call dfoirfalg(n,x,l,u,me,mi,evalf,evalc,evaljac,verbose)

end subroutine fkss
