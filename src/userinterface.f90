module userinterface

  abstract interface
     subroutine evalf(n,x,f,flag)
       ! SCALAR ARGUMENTS
       integer :: flag,n
       real(8) :: f
       ! ARRAY ARGUMENTS
       real(8) :: x(n)

       intent(in ) :: n,x
       intent(out) :: f,flag
     end subroutine evalf

     subroutine evalc(n,x,ind,c,flag)
       ! SCALAR ARGUMENTS
       integer :: flag,ind,n
       real(8) :: c
       ! ARRAY ARGUMENTS
       real(8) :: x(n)

       intent(in ) :: ind,n,x
       intent(out) :: c,flag
     end subroutine evalc

     subroutine evaljac(n,x,ind,jcvar,jcval,jcnnz,flag)
       ! SCALAR ARGUMENTS
       integer :: flag,ind,jcnnz,n
       ! ARRAY ARGUMENTS
       integer :: jcvar(n)
       real(8) :: jcval(n),x(n)

       intent(in ) :: ind,n,x
       intent(out) :: flag,jcnnz,jcval,jcvar
     end subroutine evaljac
  end interface

contains

end module userinterface
