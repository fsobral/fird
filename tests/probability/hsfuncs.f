
C     ******************************************************************
C     ******************************************************************

      subroutine hsgetdim(nprob, n_, me_, mi_)

      IMPLICIT NONE

C     SCALAR ARGUMENTS

      integer nprob, n_, me_, mi_

C     COMMON SCALARS

      integer N,NILI,NINL,NELI,NENL,NTP

C     COMMON BLOCKS

      common/L1/N,NILI,NINL,NELI,NENL
      common/L8/NTP

      NTP = nprob

      call CONV(1)

      n_ = N

      me_ = NELI + NENL

      mi_ = NILI + NINL

      end

C     ******************************************************************
C     ******************************************************************

      subroutine hsinitp(n_, x_, l_, u_)

      IMPLICIT NONE

C     PARAMETERS
      integer MMAX,NMAX

      parameter(NMAX = 1000)
      parameter(MMAX = 1000)

C     SCALAR ARGUMENTS
      integer n_

C     ARRAY ARGUMENTS
      double precision x_(n_), l_(n_), u_(n_)

C     COMMON SCALARS

      integer N,NILI,NINL,NELI,NENL,NEX,NTP

C     COMMON ARRAYS

      logical INDEX1(MMAX),INDEX2(MMAX),LXL(NMAX),LXU(NMAX)
      double precision X(NMAX),XL(NMAX),XU(NMAX)

C     COMMON BLOCKS

      common/L1/N,NILI,NINL,NELI,NENL
      common/L2/X
      common/L8/NTP
      common/L9/INDEX1
      common/L10/INDEX2
      common/L11/LXL
      common/L12/LXU
      common/L13/XL
      common/L14/XU

C     LOCAL SCALARS

      integer i

C     INITIAL POINT AND BOX CONSTRAINTS.

      do i = 1,n_
         x_(i) = X(i)
      end do

      do i = 1,n_
         if ( LXL(i) ) then
            l_(i) = XL(i)
         else
            l_(i) = - 1.0D+20
         end if
      end do

      do i = 1,n_
         if ( LXU(i) ) then
            u_(i) = XU(i)
         else
            u_(i) = 1.0D+20
         end if
      end do

C     CONSTRAINTS

      do i = 1,NILI + NINL + NELI + NENL
         INDEX1(i) = .false.
         INDEX2(i) = .false.
      end do

      end

C     ******************************************************************
C     ******************************************************************

      subroutine hsobjf(n,x_,f,flag)

      implicit none

!     PARAMETERS
      integer MMAX,NMAX

      parameter(NMAX = 1000)
      parameter(MMAX = 1000)

!     SCALAR ARGUMENTS
      integer flag,n
      double precision f

!     ARRAY ARGUMENTS
      double precision x_(n)

!     COMMON SCALARS
      double precision FX

!     COMMON ARRAYS
      double precision X(NMAX)

!     COMMON BLOCKS
      common/L2/X
      common/L6/FX

!     LOCAL SCALARS
      integer i

      flag = 0

      do i = 1,n
         X(i) = x_(i)
      end do

      CALL CONV(2)

      f = FX

      end

C     ******************************************************************
C     ******************************************************************

      subroutine hscon(nn,x_,ind,c,flag)

      implicit none

!     PARAMETERS
      integer MMAX,NMAX

      parameter(NMAX = 1000)
      parameter(MMAX = 1000)

!     SCALAR ARGUMENTS
      integer flag,ind,nn
      double precision c

!     ARRAY ARGUMENTS
      double precision x_(nn)

!     COMMON SCALARS
      integer N,NILI,NINL,NELI,NENL

C     COMMON ARRAYS
      double precision G(MMAX),X(NMAX)
      logical INDEX1(MMAX)

C     COMMON BLOCKS
      common/L1/N,NILI,NINL,NELI,NENL
      common/L2/X
      common/L3/G
      common/L9/INDEX1

C     LOCAL SCALARS
      integer i,rind

      flag = 0

      do i = 1,nn
         X(i) = x_(i)
      end do
      
      if ( ind .le. NELI + NENL ) then
         rind = ind + NILI + NINL
      else
         rind = ind - NELI - NENL
      end if

      INDEX1(rind) = .true.
      
      call conv(4)

      INDEX1(rind) = .false.
      
      c = - G(rind)

      end

C     ******************************************************************
C     ******************************************************************

      subroutine hsjac(nn,x_,ind,jcvar,jcval,jcnnz,flag)

      implicit none

!     PARAMETERS
      integer MMAX,NMAX

      parameter(NMAX = 1000)
      parameter(MMAX = 1000)

C     SCALAR ARGUMENTS
      integer flag,ind,jcnnz,nn

C     ARRAY ARGUMENTS
      integer jcvar(NMAX * MMAX)
      double precision x_(nn),jcval(NMAX * MMAX)

C     COMMON SCALARS
      integer N,NILI,NINL,NELI,NENL

C     COMMON ARRAYS
      double precision GG(NMAX * MMAX),X(NMAX)
      logical INDEX2(MMAX)

C     COMMON BLOCKS
      common/L1/N,NILI,NINL,NELI,NENL
      common/L2/X
      common/L5/GG
      common/L10/INDEX2

C     LOCAL SCALARS
      integer i,m,rind
      
      flag = 0

      do i = 1,nn
         X(i) = x_(i)
      end do

      if ( ind .le. NELI + NENL ) then
         rind = ind + NILI + NINL
      else
         rind = ind - NELI - NENL
      end if

      INDEX2(rind) = .true.
      
      CALL CONV(5)

      INDEX2(rind) = .false.
      
      jcnnz = nn

      m = NILI + NINL + NELI + NENL

      do i = 1,nn
         jcvar(i) = i
         jcval(i) = - GG((i - 1) * m + rind)
      end do 

      end
