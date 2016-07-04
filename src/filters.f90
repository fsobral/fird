module filters

  ! This module defines the prototype of a filter test used by the
  ! Filter DFO IR Algorithm.
  !
  ! It also defines two filters: flat filter and slanting filter.

  abstract interface
     function absfilter(fx,hxnorm,alpha_,nf,ffilter,hfilter)
       ! SCALAR ARGUMENTS
       integer :: nf
       real(8) :: alpha_,fx,hxnorm
       ! ARRAY ARGUMENTS
       real(8) :: ffilter(nf),hfilter(nf)
       intent(in) :: alpha_,ffilter,fx,hfilter,hxnorm,nf

       ! RETURN VALUE
       logical :: absfilter
     end function absfilter
  end interface

  public :: absfilter,flatFilter,slantingFilter

contains

  !----------------------------------------------------------!
  ! FUNCTION SLANTINGFILTER                                  !
  !                                                          !
  ! This function tests if the point for which the values fx !
  ! (objective function value) and hxnorm (norm of the       !
  ! infeasibilities) have been calculated belongs to the     !
  ! flat filter.                                             !
  !                                                          !
  ! Returns .true. if the point belongs to the filter or     !
  ! .false. otherwise. The considered filter is the          !
  ! 'slanting filter' defined by Chin and Fletcher (2003)    !
  !                                                          !
  !----------------------------------------------------------!
  
  function slantingFilter(fx,hxnorm,alpha_,nf,ffilter,hfilter)

    ! SCALAR ARGUMENTS
    integer :: nf
    real(8) :: alpha_,fx,hxnorm

    ! ARRAY ARGUMENTS
    real(8) :: ffilter(nf),hfilter(nf)

    intent(in) :: alpha_,ffilter,fx,hfilter,hxnorm,nf

    ! RETURN VALUE
    logical :: slantingFilter

    ! LOCAL SCALARS
    integer :: i

    slantingFilter = .false.

    do i = 1,nf
       if ( hxnorm .ge. (1.0D0 - alpha_) * hfilter(i) .and. &
            fx + alpha_ * hxnorm .ge. ffilter(i) ) then
          slantingFilter = .true.
          return
       end if
    end do

  end function slantingFilter

  !----------------------------------------------------------!
  ! FUNCTION FLATFILTER                                      !
  !                                                          !
  ! This function tests if the point for which the values fx !
  ! (objective function value) and hxnorm (norm of the       !
  ! infeasibilities) have been calculated belongs to the     !
  ! flat filter.                                             !
  !                                                          !
  ! Returns .true. if the point belongs to the filter or     !
  ! .false. otherwise. The considered filter is the 'flat    !
  ! filter' described in Fletcher and Leyffer (2002).        !
  !                                                          !
  !----------------------------------------------------------!

  function flatFilter(fx,hxnorm,alpha_,nf,ffilter,hfilter)

    ! SCALAR ARGUMENTS
    integer :: nf
    real(8) :: alpha_,fx,hxnorm

    ! ARRAY ARGUMENTS
    real(8) :: ffilter(nf),hfilter(nf)

    intent(in) :: alpha_,ffilter,fx,hfilter,hxnorm,nf

    ! RETURN VALUE
    logical :: flatFilter

    ! LOCAL SCALARS
    integer :: i

    flatFilter = .false.

    do i = 1,nf
       if ( hxnorm .ge. (1.0D0 - alpha_) * hfilter(i) .and. &
            fx .ge. ffilter(i) - alpha_ * hfilter(i) ) then
          flatFilter = .true.
          return
       end if
    end do

  end function flatFilter

end module filters
