!-------------------------------------------------------------------!
! SUBROUTINE RESTORATION                                            !
!-------------------------------------------------------------------!
!
! Users should modify this subroutine in order to have his/her own
! restoration method.
!
! In its current version, the subroutine uses the general nonlinear
! Augmented Lagrangian solver ALGENCAN to solve the orthogonal
! projection problem
!
!                                                            
!                min  ||x - x_k||_2^2                        
!                s.t. g(x) <= gamma_k                        
!                     l <= x <= u                            
!
! All the parameters of the subroutine are carefully explained in the
! interface file 'rinterface.f90' and in the document.

subroutine restore(n,x,l,u,me,mi,uevalc,uevaljac,epsfeas,verbose, &
     infeas,flag)

  use algencan_restoration

  implicit none

  ! SCALAR ARGUMENTS
  integer :: flag,me,mi,n
  real(8) :: epsfeas,infeas
  logical :: verbose

  ! ARRAY ARGUMENTS
  real(8) :: l(n),u(n),x(n)

  intent(in   ) :: epsfeas,l,me,mi,n,u,verbose
  intent(out  ) :: flag,infeas
  intent(inout) :: x

  ! EXTERNAL SUBROUTINES
  external :: uevalc,uevaljac

  call algrestore(n,x,l,u,me,mi,uevalc,uevaljac,epsfeas,verbose, &
       infeas,flag)

end subroutine restore
