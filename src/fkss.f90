subroutine fkss(n,x,l,u,me,mi,evalf,evalc,uevaljac,verbose,flag)

  use dfoirfilter

  implicit none

  ! SCALAR ARGUMENTS
  integer :: flag,me,mi,n
  logical :: verbose

  ! ARRAY ARGUMENTS
  real(8) :: x(n)

  ! EXTERNAL SUBROUTINES
  external :: evalf,evalc,evaljac

  ! LOCAL SCALARS
  integer :: flag
  real(1) :: beta,c,cfeas,rinfeas

  call dfoirfalg(n,x,l,u,me,mi,evalf,evalc,evaljac,verbose)

end subroutine fkss
