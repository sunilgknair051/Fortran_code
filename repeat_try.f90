program repeat
implicit none

module constants
real w0,ws,wd
parameter(w0=4./9., ws=1./9., wd= 1./36.)
integer j,x,y,nx,ny,ndir
parameter(nx=2,ny=2,ndir=8)
real f(1:nx,1:ny,0:ndir),fn(1:nx,1:ny,0:ndir),f_out(1:nx,1:ny,0:ndir)
real r(1:nx,1:ny),ux(1:nx,1:ny),uy(1:nx,1:ny)
integer :: dirx(0:8)=(/0,1,0,-1,0,1,-1,-1,1/)
integer ::diry(0:8)=(/0,0,1,0,-1,1,1,-1,-1/)
integer ::no_slip(0:8)=(/0,3,4,1,2,7,8,5,6/)
real ::wi(0:8)=(/w0,ws,ws,ws,ws,wd,wd,wd,wd/)
end module constants


use constants

 DO x= 1,nx
   DO y = 1,ny
     r(x,y) = 1
     ux(x,y) = 0
     uy(x,y) = 0
   ENDDO
 ENDDO

 call init_equilibrium(dirx,diry,wi,r,ux,uy,f)
 print*, f

end program repeat

subroutine init_equilibrium(dirx,diry,wi,r,ux,uy,f)
use constants
  !integer :: dirx(0:8)
  integer ::diry(0:8)
  !integer i,x,y,nx,ny,ndir
  !real cidotu, feq, omtauniv, tauinv
  !parameter(nx=2,ny=2,ndir=8)
  !real f(1:nx,1:ny,0:ndir),fn(1:nx,1:ny,0:ndir),f_out(1:nx,1:ny,0:ndir)
  !real r(1:nx,1:ny),ux(1:nx,1:ny),uy(1:nx,1:ny)
  !real ::wi(0:8)
  DO x= 1,nx
    DO y =1,ny
      DO i = 0,ndir
        cidotu = dirx(i)*ux(x,y)+diry(i)*uy(x,y)
        f(x,y,i) = wi(i)*r(x,y)*(1+3*cidotu+4.5*cidotu*cidotu-1.5*(ux(x,y)*ux(x,y)+uy(x,y)*uy(x,y)))
        !print *, f(x,y,i)
      END DO
    END DO
  END DO
END subroutine init_equilibrium


subroutine compute_rho_u(f,r,ux,uy)
use constants
!integer :: dirx(0:8)
!integer ::diry(0:8)
!integer i,x,y,nx,ny,ndir
!real cidotu, feq, omtauniv, tauinv
!parameter(nx=2,ny=2,ndir=8)
!real f(1:nx,1:ny,0:ndir),fn(1:nx,1:ny,0:ndir),f_out(1:nx,1:ny,0:ndir)
!real r(1:nx,1:ny),ux(1:nx,1:ny),uy(1:nx,1:ny)
!real ::wi(0:8)

DO x= 1,nx
  DO y =1,ny
            !real r(x,y) = 0
            !real ux(x,y) = 0
            !real uy(x,y) = 0
    !r(x,y)= 0, ux(x,y) = 0, uy(x,y)= 0
     DO i = 0,ndir
       r(x,y) = r(x,y)+f(x,y,i)
       ux(x,y) = ux(x,y)+f(x,y,i)*dirx(i)
       uy(x,y) = uy(x,y)+f(x,y,i)*diry(i)
     END DO
    ux(x,y) = ux(x,y)/r(x,y)
    uy(x,y) = uy(x,y)/r(x,y)
  END DO
END DO
end subroutine compute_rho_u


subroutine collide(f,r,ux,uy)
use constants
!integer :: dirx(0:8)
!integer ::diry(0:8)
!integer i,x,y,nx,ny,ndir
!real cidotu, feq, omtauniv, tauinv
!parameter(nx=2,ny=2,ndir=8)
!real f(1:nx,1:ny,0:ndir),fn(1:nx,1:ny,0:ndir),f_out(1:nx,1:ny,0:ndir)
!real r(1:nx,1:ny),ux(1:nx,1:ny),uy(1:nx,1:ny)
!real ::wi(0:8)

DO x= 1,nx
  DO y =1,ny
     DO i = 0,ndir
        cidotu = dirx(i)*ux(x,y)+diry(i)*uy(x,y)
        f(x,y,i) = wi(i)*r(x,y)*(1+3*cidotu+4.5*cidotu*cidotu-1.5*(ux(x,y)*ux(x,y)+uy(x,y)*uy(x,y)))
        f_out(x,y,i) = (omtauniv*f(x,y,i))+(tauinv*feq)
     END DO
  END DO
END DO
end subroutine collide
