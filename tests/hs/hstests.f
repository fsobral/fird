      PROGRAM PRINCIPAL
      
      IMPLICIT NONE

C     PARAMETERS
      integer MMAX,NMAX

      parameter(NMAX = 1000)
      parameter(MMAX = 1000)

C     LOCAL SCALARS

      integer I,FLAG,ME,MI,M,NN
      double precision RELDIFF

C     LOCAL ARRAYS

      double precision X_(NMAX),L(NMAX),U(NMAX)

C     USER-DEFINED SUBROUTINES

      external calobjf,calcon,caljac

C     COMMON SCALARS

      integer N,NILI,NINL,NELI,NENL,NEX,NTP
      logical INDEX1(MMAX),INDEX2(MMAX),LXL(NMAX),LXU(NMAX)
      double precision G(MMAX),GG(NMAX * MMAX),X(NMAX),XL(NMAX),XU(NMAX)
      logical LEX
      double precision FEX

C     COMMON ARRAYS

      double precision XEX(NMAX * NMAX)

C     COMMON BLOCKS

      common/L1/N,NILI,NINL,NELI,NENL
      common/L2/X
      common/L3/G
      common/L5/GG
      common/L8/NTP
      common/L9/INDEX1
      common/L10/INDEX2
      common/L11/LXL
      common/L12/LXU
      common/L13/XL
      common/L14/XU
      COMMON/L20/LEX,NEX,FEX,XEX

C     LOCAL SCALARS

      character OPTM
      logical VERBOSE
      integer FCNT,FTYPE
      double precision C,F,FEAS,EPSFEAS,EPSOPT

C     WRITE(*,*) 'Number of the problem: '
      READ(*,*)NTP

      CALL CONV(1)

C     NUMBER OF VARIABLES
      nn = N

C     INITIAL POINT AND BOX CONSTRAINTS.

      do i = 1,nn
         x_(i) = X(i)
      end do

      do i = 1,nn
         if ( LXL(i) ) then
            l(i) = XL(i)
         else
            l(i) = - 1.0D+20
         end if
      end do

      do i = 1,nn
         if ( LXU(i) ) then
            u(i) = XU(i)
         else
            u(i) = 1.0D+20
         end if
      end do

C     NUMBER OF CONSTRAINTS
      ME = NELI + NENL
      MI = NILI + NINL
      M  = ME + MI

C     CONSTRAINTS

      do i = 1,M
         INDEX1(i) = .false.
         INDEX2(i) = .false.
      end do

C     CALLS THE ALGORITHM

      open(75,FILE='runhs.out')
      write(75,0020) NTP,N,NILI + NINL,NELI + NENL,1.0D20,1.0D20,
     +     1.0D20,-1,-1
      close(75)

C     Some HS problems do not have derivatives of the constraints
      if ( NTP .eq. 332 .or. NTP .eq. 348 .or. NTP .eq. 349 .or.
     +     NTP .eq. 356 .or. NTP .eq. 362 .or. NTP .eq. 363 .or.
     +     NTP .eq. 364 .or. NTP .eq. 365 .or. NTP .eq. 366 .or. 
     +     NTP .eq. 369 .or. NTP .eq. 390 .or. NTP .eq. 392 .or.
     +     NTP .eq. 393 ) then
         stop
      end if

      VERBOSE = .true.

      EPSFEAS = 1.0D-08
      
      EPSOPT = 1.0D-04

      FTYPE = 2

      call fird(n,x_,l,u,me,mi,calobjf,calcon,caljac,verbose,ftype,
     +     epsfeas,epsopt,f,feas,fcnt,flag)

      reldiff = (f - FEX) / max(1.0D0,abs(f),abs(FEX))

      optm = ' '
      if ( reldiff .le. 1.0D-01 .and. FEAS .le. EPSFEAS ) then
         optm = '*'
      end if

      open(75,FILE='runhs.out')
      write(75,0020) NTP,N,NILI + NINL,NELI + NENL,FEX,F,FEAS,FCNT,
     +     FLAG
      write(*,0021) NTP,N,NILI + NINL,NELI + NENL,FEX,F,FEAS,FCNT,
     +     FLAG,optm
      close(75)

!     NON-EXECUTABLE STATEMENTS

 0020 FORMAT(I4,1X,I4,1X,I4,1X,I4,5X,E15.8,1X,E15.8,1X,E15.8,1X,I15,
     +     1X,I4)
 0021 FORMAT(I4,1X,I4,1X,I4,1X,I4,5X,E15.8,1X,E15.8,1X,E15.8,1X,I15,
     +     1X,I4,1X,A1)

      END PROGRAM PRINCIPAL

C     ******************************************************************
C     ******************************************************************

      subroutine calobjf(n,x_,f,flag)

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

      end subroutine calobjf

C     ******************************************************************
C     ******************************************************************

      subroutine calcon(nn,x_,ind,c,flag)

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

      end subroutine calcon

C     ******************************************************************
C     ******************************************************************

      subroutine caljac(nn,x_,ind,jcvar,jcval,jcnnz,flag)

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

      end subroutine caljac
