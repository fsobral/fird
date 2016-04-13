subroutine fkss(n,x,l,u,me,mi,uevalf,uevalc,uevaljac,verbose)

  use filterirdfo
  use restoration

  implicit none

  ! SCALAR ARGUMENTS
  integer :: me,mi,n
  logical :: verbose

  ! ARRAY ARGUMENTS
  real(8) :: x(n)

  ! EXTERNAL SUBROUTINES
  external :: uevalf,uevalc,uevaljac

  ! LOCAL ARRAYS
  real(8) :: y(n),z(n)

  ! LOCAL SCALARS
  integer :: flag
  real(1) :: cfeas,rinfeas

  ! ----------------- !
  ! Feasibility phase !
  ! ----------------- !

  restore(n,x,l,u,me,mi,uevalc,uevaljac,cfeas,verbose,rinfeas,flag)

  xinfeas = evalInfeas(n,x,me,mi,uevalc)

  ! Verify convergence conditions

end subroutine fkss
