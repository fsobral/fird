subroutine fkss(n,x,l,u,me,mi,uevalf,uevalc,uevaljac,verbose,ftype, &
     epsfeas,epsopt,f,feas,fcnt,flag)

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

end subroutine fkss
