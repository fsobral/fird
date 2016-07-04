module filters

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
  ! FUNCTION FLATFILTER                                      !
  !                                                          !
  ! This function tests if the point for which the values fx !
  ! (objective function value) and hxnorm (norm of the       !
  ! infeasibilities) have been calculated belongs to the     !
  ! flat filter.                                             !
  !                                                          !
  ! Returns .true. if the point belongs to the filter or     !
  ! .false. otherwise.                                       !
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

    Write(*,*)
    Write(*,*)
    Write(*,*)
    write(*,*) '----------> Slanting'
    Write(*,*)
    Write(*,*)
    Write(*,*)

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
  ! .false. otherwise.                                       !
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

    Write(*,*)
    Write(*,*)
    Write(*,*)
    write(*,*) '----------> Flat'
    Write(*,*)
    Write(*,*)
    Write(*,*)

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
