module ccpdata

  implicit none

  ! PARAMETERS

  ! Number of observations in order to build the correlation matrix
  integer, parameter :: NOBS = 5
  ! Sparsity rate
  real(8), parameter :: SPRATE = 1.0D-01
  ! Penalization term for the chance constraint
  real(8), parameter :: PEN = 1.0D+04
  
  ! COMMON ARRAYS

  integer, allocatable :: cA(:), pA(:)
  real(8), allocatable :: A(:), b(:), c(:), corr(:), mu(:)

  ! COMMON SCALARS

  integer :: m, np
  real(8) :: plim

contains

  subroutine initialize(n, np_, mi, x, l, u, plim_)

    ! This subroutine randomly initializes the structure needed by the
    ! Chance Constrained Problem.

    implicit none

    ! SCALAR ARGUMENTS
    integer :: n, np_, mi
    real(8) :: plim_

    ! ARRAY ARGUMENTS
    real(8) :: l(n), u(n), x(n)

    intent(in ) :: n, np_, mi, plim_
    intent(out) :: l, u, x

    ! LOCAL SCALARS

    integer :: i, j, k, nover, nseed, pos
    real(8) :: randnum, t

    ! LOCAL ARRAYS

    integer, allocatable :: seed(:)
    real(8), allocatable :: tmpA(:), tmpB(:)

    m  = mi

    np = np_

    plim = plim_

    ! Initialize the random structure

    call random_seed(SIZE=nseed)

    allocate(seed(nseed))

    do i = 1,nseed
       seed(i) = 123456789
    end do

    call random_seed(PUT=seed)

    ! Initialize module objects

    allocate(b(m), c(n), corr((np - 1) * np / 2), mu(np), &
             tmpA(NOBS), tmpB(NOBS), A(n * m), pA(m + 1), &
             cA(n * m))

    ! Randomly generates vector 'c'
    
    do i = 1, np

       call random_number(randnum)

       c(i) = 1.0D0 + 1.0D+01 * randnum

    end do

    do i = np + 1, n

       call random_number(randnum)

       c(i) = - 1.0D+01 + 2.0D+01 * randnum

    end do

    ! Randomly generates sparse matrix 'A' in vector form

    pA(1) = 1

    do i = 1, m

       pos = pA(i)

       do j = 1, n

          call random_number(randnum)

          if ( randnum .le. SPRATE ) then

             call random_number(A(pos))

             if ( j .le. np ) A(pos) = - abs(A(pos))

             cA(pos) = j

             pos = pos + 1

          end if

       end do

       pA(i + 1) = pos

    end do

    ! Randomly generates mean 'mu'

    do i = 1, np

       call random_number(randnum)

       mu(i) = 5.0D0 * randnum

    end do

    ! Randomly generates the initial guess, lower and upper bounds

    do i = 1, np

       x(i) = 0.0D0

       l(i) = mu(i)

       u(i) = 1.0D+20

    end do

    do i = np + 1, n

       call random_number(randnum)

       x(i) = - 1.0D+01 + 2.0D+01 * randnum

       l(i) = - 1.0D+01

       u(i) =   1.0D+01

    end do

    ! Randomly generates feasible right hand side 'b'

    do i = 1, m

       b(i) = 0.0D0

       do j = pA(i), pA(i + 1) - 1

          b(i) = b(i) + x(cA(j)) * A(j)

       end do

    end do

    ! Randomly generates correlation matrix 'corr', which is simmetric
    ! and has 1's in the diagonal (in vector form)

    do j = 1, np

       do i = 1, j - 1

          do k = 1, NOBS

             call random_number(tmpA)

             call random_number(tmpB)

             call unbiased_correlation(NOBS, tmpA, tmpB, 1.0D0, 0, &
                  corr(i + (j - 1) * (j - 2) / 2), t, nover)

          end do

       end do

    end do

    deallocate(tmpA, tmpB, seed)
    
  end subroutine initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine destroy()

    deallocate(A, b, c, corr, mu, pA, cA)

  end subroutine destroy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine tofile(n, np, mi, x, l, u, filename)

    implicit none

    ! SCALAR ARGUMENTS
    integer :: n, np, mi

    ! ARRAY ARGUMENTS
    real(8)       :: l(n), u(n), x(n)
    character(80) :: filename

    intent(in ) :: filename, n, np, mi, l, u, x

    ! LOCAL SCALARS
    integer :: i, j

    ! LOCAL ARRAYS
    real(8) :: tmpA(n)

    open(99, FILE=filename)

    write(99, FMT=010) n, np, mi

    do i = 1, mi

       do j = 1, n

          tmpA(j) = 0.0D0

       end do

       do j = pA(i), pA(i + 1) - 1

          tmpA(cA(j)) = A(j)

       end do

       do j = 1, n
       
          write(99,FMT=020) tmpA(j)

       end do

       write(99,*)

    end do

    write(99,*)

    close(99)

    ! NON-EXECUTABLE STATEMENTS

010 FORMAT(I5,1X,I5,1X,I5,/)
020 FORMAT(1P,E15.8)

  end subroutine tofile

end module ccpdata
