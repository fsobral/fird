module ointerface

  use userinterface
  use filters      , only: absfilter

  abstract interface
     subroutine optimization(n,y,l,u,me,mi,uevalf,uevalc,uevallc,uevalljac, &
          nf,alpha,ffilter,hfilter,filterTest,outiter,epsfeas,epsopt,       &
          verbose,delta,fy,hynorm,rho,flag)
  
       ! TODO: Comment this interface!

       ! SCALAR ARGUMENTS
       logical :: verbose
       integer :: flag,me,mi,n,nf,outiter
       real(8) :: alpha,delta,epsfeas,epsopt,fy,hynorm,rho
       
       ! ARRAY ARGUMENTS
       real(8) :: ffilter(nf),hfilter(nf),l(n),u(n),y(n)
       
       ! EXTERNAL SUBROUTINES
       procedure(evalf)     :: uevalf
       procedure(evalc)     :: uevallc,uevalc
       procedure(evaljac)   :: uevalljac
       procedure(absfilter) :: filterTest

       intent(in   ) :: alpha,epsfeas,epsopt,ffilter,hfilter,l,me,mi,n, &
            nf,outiter,u,verbose
       intent(out  ) :: hynorm,flag
       intent(inout) :: delta,fy,rho,y
     end subroutine optimization
  end interface

contains

end module ointerface
