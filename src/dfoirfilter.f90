module dfoirfilter

  implicit none

contains

  subroutine dfoirfalg(n,x,l,u,me,mi,evalf,evalc,evaljac,verbose)

    implicit none

    ! SCALAR ARGUMENTS
    integer :: me,mi,n
    logical :: verbose

    ! ARRAY ARGUMENTS
    real(8) :: x(n)

    ! EXTERNAL SUBROUTINES
    external :: evalf,evalc,evaljac

    ! LOCAL ARRAYS
    real(8) :: y(n),z(n)

    ! LOCAL SCALARS
    integer :: i,flag
    real(1) :: cfeas,rinfeas

    ! Initialization
    ! TODO: check for errors when allocating

    ! ----------------- !
    ! Feasibility phase !
    ! ----------------- !

    do i = 1,n
       xprev(i) = x(i)
    end do

    restore(n,x,l,u,me,mi,evalc,uevaljac,cfeas,verbose,rinfeas,flag)

    xinfeas = evalInfeas(n,x,me,mi,uevalc)

    ! Verify convergence conditions


    ! Finalization
    ! TODO: check for errors
    deallocate(xprev)

  end subroutine dfoirfalg

end module dfoirfilter
