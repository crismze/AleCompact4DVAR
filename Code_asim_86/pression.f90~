
!
! ****************************************************************
!
subroutine derivep(n1b,pr,n1) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
  USE paramx_m 
  USE derivex_m 
!...Translated by PSUITE Trans90                  4.3ZH 16:05:14   1/27/ 4  
!...Switches: -yv INDDO=0 -nbejkno
  implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  integer , intent(in) :: n1b 
  integer  :: n1 
  real(8) , intent(out) :: pr(n1b,11) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer :: n1bm1, n1bm2, n1bm3, n1bm4, n1bm5, n1bm6, n1bm7, i, j, ip 
!-----------------------------------------------
!
! ****************************************************************
!
   n1bm1=n1b - 1 
   n1bm2=n1b - 2 
   n1bm3=n1b - 3 
   n1bm4=n1b - 4 
   n1bm5=n1b - 5 
   n1bm6=n1b - 6 
   n1bm7=n1b - 7 
!
   do i=1,n1b 
   do j=1,11 
      pr(i,j)=0. 
   enddo
   enddo 
!
   pr(1,6)=af2x 
   pr(1,7)=alfa2x 
   pr(1,9)=1. 
   pr(1,10)=-af2x 
   pr(1,11)=alfa2x 
!
!     ------------------------------------------
!     pas de conditions aux limites sur p en i=1
!
!      pr(2,5 ) =-af1x
!      pr(2,6 ) = 1.
!      pr(2,7 ) =-bf1x
!      pr(2,8 ) = alfa1x
!      pr(2,9 ) =-cf1x
!     ------------------------------------------
!
!     --------------------------
!     condition dp/dy = ? en i=1
!
   pr(2,6)=1. 
!     --------------------------
!
   pr(3,9)=af2x 
   pr(4,6)=bf1x 
   pr(4,8)=cf1x 
   pr(5,5)=-afix 
   pr(5,9)=afix 
   pr(5,11)=bfix 
   pr(6,1)=bfix 
   pr(6,3)=afix 
   pr(6,4)=alfaix 
   pr(6,6)=1. 
   pr(6,7)=-afix 
   pr(6,8)=alfaix 
   pr(6,9)=-bfix 
!
!     ------------------------------------------------
!     coefficients conditions aux limites sur u en i=1
!     ------------------------------------------------
!
!      pr(3,5 ) =-af2x
!      pr(4,4 ) = af1x
!      pr(5,3 ) =-bfix
!
   do i=7,n1bm7,2 
      ip=i+1 
      pr(i,3)=-bfix 
      pr(i,5)=-afix 
      pr(i,9)=afix 
      pr(i,11)=bfix 
!
      pr(ip,1)=bfix 
      pr(ip,3)=afix 
      pr(ip,4)=alfaix 
      pr(ip,6)=1. 
      pr(ip,7)=-afix 
      pr(ip,8)=alfaix 
      pr(ip,9)=-bfix 
   enddo
!
   pr(n1bm5,3)=-bfix 
   pr(n1bm5,5)=-afix 
   pr(n1bm5,9)=afix 
   pr(n1bm4,1)=bfix 
   pr(n1bm4,3)=afix 
   pr(n1bm4,4)=alfaix 
   pr(n1bm4,6)=1. 
   pr(n1bm4,7)=-afix 
   pr(n1bm4,8)=alfaix 
   pr(n1bm4,9)=-bfix 
   pr(n1bm3,5)=-afmx 
   pr(n1bm2,4)=-cfnx 
   pr(n1bm2,6)=-bfnx 
   pr(n1bm1,2)=afmx 
   pr(n1bm1,3)=alfamx 
   pr(n1bm1,5)=1. 
   pr(n1bm1,6)=-afmx 
   pr(n1bm1,7)=alfamx 
!
!     -------------------------------------------
!     pas de conditions aux limites sur p en i=n1
!
!      pr(n1b  ,1 ) = cfnx
!      pr(n1b  ,3 ) = bfnx
!      pr(n1b  ,4 ) = alfanx
!      pr(n1b  ,5 ) = afnx
!      pr(n1b  ,6 ) = 1.
!     -------------------------------------------
!
!     ---------------------------
!     condition dp/dy = ? en i=n1
!
   pr(n1b,6)=1. 
!     ---------------------------
!
!
!     -------------------------------------------------
!     coefficients conditions aux limites sur u en i=n1
!     -------------------------------------------------
!
!      pr(n1bm5,11) = bfix
!      pr(n1bm3,9 ) = afmx
!      pr(n1bm2,8 ) =-afnx
!
   return  
end subroutine derivep
!
!********************************************************************
!
subroutine divpres(pp,ux,uy,uz,tx,ty,tz,tr,di,sx,sy,sz,ffx,fsx,fwx,ffy,&
     fsy,fwy,ffz,fsz,fwz,ffxp,fsxp,fwxp,ffyp,fsyp,fwyp,ffzp,fszp,fwzp,bxn,&
     byn,bzn,bx1,by1,bz1,nx,ny,nz) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
  USE paramx_m 
  USE derivex_m 
!...Translated by PSUITE Trans90                  4.3ZH 16:05:14   1/27/ 4  
!...Switches: -yv INDDO=0 -nbejkno
  implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  integer  :: nx 
  integer  :: ny 
  integer  :: nz 
  real(8) , intent(out) :: pp(nx,ny,nz) 
  real(8) , intent(inout) :: ux(nx,ny,nz) 
  real(8)  :: uy(nx,ny,nz) 
  real(8)  :: uz(nx,ny,nz) 
  real(8) , intent(inout) :: tx(nx,ny,nz) 
  real(8)  :: ty(nx,ny,nz) 
  real(8)  :: tz(nx,ny,nz) 
  real(8)  :: tr(nx,ny,nz) 
  real(8)  :: di(nx,ny,nz) 
  real(8)  :: sx(ny,nz) 
  real(8)  :: sy(nx,nz) 
  real(8)  :: sz(nx,ny) 
  real(8)  :: ffx(nx) 
  real(8)  :: fsx(nx) 
  real(8)  :: fwx(nx) 
  real(8)  :: ffy(ny) 
  real(8)  :: fsy(ny) 
  real(8)  :: fwy(ny) 
  real(8)  :: ffz(nz) 
  real(8)  :: fsz(nz) 
  real(8)  :: fwz(nz) 
  real(8)  :: ffxp(nx) 
  real(8)  :: fsxp(nx) 
  real(8)  :: fwxp(nx) 
  real(8)  :: ffyp(ny) 
  real(8)  :: fsyp(ny) 
  real(8)  :: fwyp(ny) 
  real(8)  :: ffzp(nz) 
  real(8)  :: fszp(nz) 
  real(8)  :: fwzp(nz) 
  real(8) , intent(in) :: bxn(ny,nz) 
  real(8) , intent(in) :: byn(ny,nz) 
  real(8) , intent(in) :: bzn(ny,nz) 
  real(8) , intent(in) :: bx1(ny,nz) 
  real(8) , intent(in) :: by1(ny,nz) 
  real(8) , intent(in) :: bz1(ny,nz) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer :: nyz, nxyz, k, j, ijk, jk, i 
!-----------------------------------------------
!
!********************************************************************
!
   nyz=ny*nz 
   nxyz=nx*ny*nz 
!
   do k=1,nz 
   do j=1,ny 
      ux(1,j,k)=bx1(j,k) 
      uy(1,j,k)=by1(j,k) 
      uz(1,j,k)=bz1(j,k) 
      ux(nx,j,k)=bxn(j,k) 
      uy(nx,j,k)=byn(j,k) 
      uz(nx,j,k)=bzn(j,k) 
   enddo 
   enddo 
!  
   stop 'pression.f90 ligne 226'
   call dery (ty,uy,tr,di,sy,ffy,fsy,fwy,nx,ny,nz,0) 
   if (nz>1) call derz (tz,uz,tr,sz,ffz,fsz,fwz,nx,ny,nz,0) 
!
   if (nz>1) then 
      do ijk=1,nxyz 
         ty(ijk,1,1)=ty(ijk,1,1)+tz(ijk,1,1) 
      enddo
   endif
!
   do jk=1,nyz 
      tx(1,jk,1)=af1x*ux(1,jk,1)+bf1x*ux(2,jk,1)+cf1x*ux(3,jk,1) 
      tx(2,jk,1)=af2x*(ux(3,jk,1)-ux(1,jk,1)) 
   enddo
!
   do i=3,nx-2 
   do jk=1,nyz 
      tx(i,jk,1)=afix*(ux(i+1,jk,1)-ux(i-1,jk,1))&
                +bfix*(ux(i+2,jk,1)-ux(i-2,jk,1)) 
   enddo
   enddo 
!
   do jk=1,nyz 
      tx(nx-1,jk,1)=afmx*(ux(nx,jk,1)-ux(nx-2,jk,1)) 
      tx(nx,jk,1)=(-afnx*ux(nx,jk,1))&
                 -bfnx*ux(nx-1,jk,1)-cfnx*ux(nx-2,jk,1) 
   enddo
!
   do jk=1,nyz 
      pp(1,jk,1)=tx(1,jk,1)+ty(1,jk,1)+alfa1x*ty(2,jk,1) 
      pp(2,jk,1)=tx(2,jk,1)+alfa2x*(ty(1,jk,1)+ty(3,jk,1))+ty(2,jk,1) 
   enddo
!
   do i=3,nx-2 
   do jk=1,nyz 
      pp(i,jk,1)=tx(i,jk,1)+alfaix*(ty(i-1,jk,1)+ty(i+1,jk,1))+ty(i,jk,1) 
   enddo
   enddo
!
   do jk=1,nyz 
      pp(nx-1,jk,1)=tx(nx-1,jk,1)&
                   +alfamx*(ty(nx-2,jk,1)+ty(nx,jk,1))+ty(nx-1,jk,1) 
      pp(nx,jk,1)=tx(nx,jk,1)+ty(nx,jk,1)+alfanx*ty(nx-1,jk,1) 
   enddo
!
   return  
end subroutine divpres
!
!********************************************************************
!
subroutine divbric2(pp,ux,uy,uz,tx,ty,tz,tr,di,sx,sy,sz,ffx,fsx,fwx,ffy,&
     fsy,fwy,ffz,fsz,fwz,ffxp,fsxp,fwxp,ffyp,fsyp,fwyp,ffzp,fszp,fwzp,bxn,&
     byn,bzn,bx1,by1,bz1,nx,ny,nz,esp,ax,ay,az) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
  USE paramx_m 
  USE paramy_m 
  USE paramz_m 
  USE derivex_m 
  USE param1_m 
  USE barreau_m 
!...Translated by PSUITE Trans90                  4.3ZH 16:05:14   1/27/ 4  
!...Switches: -yv INDDO=0 -nbejkno
  implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  integer  :: nx 
  integer  :: ny 
  integer  :: nz 
  real(8) , intent(out) :: pp(nx,ny,nz) 
  real(8) , intent(inout) :: ux(nx,ny,nz) 
  real(8) , intent(inout) :: uy(nx,ny,nz) 
  real(8) , intent(inout) :: uz(nx,ny,nz) 
  real(8)  :: tx(nx,ny,nz) 
  real(8)  :: ty(nx,ny,nz) 
  real(8)  :: tz(nx,ny,nz) 
  real(8)  :: tr(nx,ny,nz) 
  real(8)  :: di(nx,ny,nz) 
  real(8)  :: sx(ny,nz) 
  real(8)  :: sy(nx,nz) 
  real(8)  :: sz(nx,ny) 
  real(8)  :: ffx(nx) 
  real(8)  :: fsx(nx) 
  real(8)  :: fwx(nx) 
  real(8)  :: ffy(ny) 
  real(8)  :: fsy(ny) 
  real(8)  :: fwy(ny) 
  real(8)  :: ffz(nz) 
  real(8)  :: fsz(nz) 
  real(8)  :: fwz(nz) 
  real(8)  :: ffxp(nx) 
  real(8)  :: fsxp(nx) 
  real(8)  :: fwxp(nx) 
  real(8)  :: ffyp(ny) 
  real(8)  :: fsyp(ny) 
  real(8)  :: fwyp(ny) 
  real(8)  :: ffzp(nz) 
  real(8)  :: fszp(nz) 
  real(8)  :: fwzp(nz) 
  real(8) , intent(in) :: bxn(ny,nz) 
  real(8) , intent(in) :: byn(ny,nz) 
  real(8) , intent(in) :: bzn(ny,nz) 
  real(8) , intent(in) :: bx1(ny,nz) 
  real(8) , intent(in) :: by1(ny,nz) 
  real(8) , intent(in) :: bz1(ny,nz) 
  real(8) , intent(in) :: esp(nx,ny,nz) 
  real(8)  :: ax(nx,ny,nz) 
  real(8)  :: ay(nx,ny,nz) 
  real(8)  :: az(nx,ny,nz) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer :: nyz, nxyz, k, j, i, ijk, jk 
!-----------------------------------------------
!
!********************************************************************
!
   nyz=ny*nz 
   nxyz=nx*ny*nz 
!
   do k=1,nz 
   do j=1,ny 
      ux(1,j,k)=bx1(j,k) 
      uy(1,j,k)=by1(j,k) 
      uz(1,j,k)=bz1(j,k) 
      ux(nx,j,k)=bxn(j,k) 
      uy(nx,j,k)=byn(j,k) 
      uz(nx,j,k)=bzn(j,k) 
   enddo 
   enddo 
!
   do k=1,nz 
   do j=1,ny 
   do i=1,nx 
      ax(i,j,k)=(1.-esp(i,j,k))*ux(i,j,k) 
      ay(i,j,k)=(1.-esp(i,j,k))*uy(i,j,k) 
      az(i,j,k)=(1.-esp(i,j,k))*uz(i,j,k) 
   enddo 
   enddo 
   enddo 
!
   call derx (tx,ax,tr,sx,ffx,fsx,fwx,nx,ny,nz,0) 
   call dery (ty,ay,tr,di,sy,ffy,fsy,fwy,nx,ny,nz,0) 
   if (nz>1) call derz (tz,az,tr,sz,ffz,fsz,fwz,nx,ny,nz,0) 
!
   if (nz>1) then 
      do ijk=1,nxyz 
         tx(ijk,1,1)=tx(ijk,1,1)+ty(ijk,1,1)+tz(ijk,1,1) 
      enddo
   else 
      do ijk=1,nxyz 
         tx(ijk,1,1)=tx(ijk,1,1)+ty(ijk,1,1) 
      enddo
   endif
!
   do jk=1,nyz 
      pp(1,jk,1)=tx(1,jk,1)+alfa1x*tx(2,jk,1) 
      pp(2,jk,1)=alfa2x*(tx(1,jk,1)+tx(3,jk,1))+tx(2,jk,1) 
   enddo
!
   do i=3,nx-2 
   do jk=1,nyz 
      pp(i,jk,1)=alfaix*(tx(i-1,jk,1)+tx(i+1,jk,1))+tx(i,jk,1) 
   enddo
   enddo
!
   do jk=1,nyz 
      pp(nx-1,jk,1)=alfamx*(tx(nx-2,jk,1)+tx(nx,jk,1))+tx(nx-1,jk,1) 
      pp(nx,jk,1)=tx(nx,jk,1)+alfanx*tx(nx-1,jk,1) 
   enddo
!
   return  
end subroutine divbric2


!
!********************************************************************
!
subroutine divbric1(pp,ux,uy,uz,tx,ty,tz,tr,di,sx,sy,sz,ffx,fsx,fwx,ffy,&
     fsy,fwy,ffz,fsz,fwz,ffxp,fsxp,fwxp,ffyp,fsyp,fwyp,ffzp,fszp,fwzp,bxn,&
     byn,bzn,bx1,by1,bz1,nx,ny,nz,esp) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
  USE paramx_m 
  USE paramy_m 
  USE paramz_m 
  USE derivex_m 
  USE param1_m 
  USE barreau_m 
!...Translated by PSUITE Trans90                  4.3ZH 16:05:14   1/27/ 4  
!...Switches: -yv INDDO=0 -nbejkno
  implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  integer  :: nx 
  integer  :: ny 
  integer  :: nz 
  real(8) , intent(out) :: pp(nx,ny,nz) 
  real(8) , intent(in) :: ux(nx,ny,nz) 
  real(8) , intent(in) :: uy(nx,ny,nz) 
  real(8) , intent(in) :: uz(nx,ny,nz) 
  real(8)  :: tx(nx,ny,nz) 
  real(8)  :: ty(nx,ny,nz) 
  real(8)  :: tz(nx,ny,nz) 
  real(8)  :: tr(nx,ny,nz) 
  real(8)  :: di(nx,ny,nz) 
  real(8)  :: sx(ny,nz) 
  real(8)  :: sy(nx,nz) 
  real(8)  :: sz(nx,ny) 
  real(8)  :: ffx(nx) 
  real(8)  :: fsx(nx) 
  real(8)  :: fwx(nx) 
  real(8)  :: ffy(ny) 
  real(8)  :: fsy(ny) 
  real(8)  :: fwy(ny) 
  real(8)  :: ffz(nz) 
  real(8)  :: fsz(nz) 
  real(8)  :: fwz(nz) 
  real(8)  :: ffxp(nx) 
  real(8)  :: fsxp(nx) 
  real(8)  :: fwxp(nx) 
  real(8)  :: ffyp(ny) 
  real(8)  :: fsyp(ny) 
  real(8)  :: fwyp(ny) 
  real(8)  :: ffzp(nz) 
  real(8)  :: fszp(nz) 
  real(8)  :: fwzp(nz) 
  real(8)  :: bxn(ny,nz) 
  real(8)  :: byn(ny,nz) 
  real(8)  :: bzn(ny,nz) 
  real(8)  :: bx1(ny,nz) 
  real(8)  :: by1(ny,nz) 
  real(8)  :: bz1(ny,nz) 
  real(8) , intent(in) :: esp(nx,ny,nz) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer :: nyz, nxyz, k, j, i, ijk 
  real(8), dimension(nx,ny,nz) :: ax, ay, az 
!-----------------------------------------------
!
!********************************************************************
!
   nyz=ny*nz 
   nxyz=nx*ny*nz 
!
   do k=1,nz 
   do j=1,ny 
   do i=1,nx 
      ax(i,j,k)=(1.-esp(i,j,k))*ux(i,j,k) 
      ay(i,j,k)=(1.-esp(i,j,k))*uy(i,j,k) 
      az(i,j,k)=(1.-esp(i,j,k))*uz(i,j,k) 
   enddo 
   enddo 
   enddo 
!
   call derx (tx,ax,tr,sx,ffx,fsx,fwx,nx,ny,nz,0) 
   call dery (ty,ay,tr,di,sy,ffy,fsy,fwy,nx,ny,nz,0) 
   if (nz>1) call derz (tz,az,tr,sz,ffz,fsz,fwz,nx,ny,nz,0) 
!
   if (nz>1) then 
      do ijk=1,nxyz 
         tx(ijk,1,1)=tx(ijk,1,1)+ty(ijk,1,1)+tz(ijk,1,1) 
      enddo
   else 
      do ijk=1,nxyz 
         tx(ijk,1,1)=tx(ijk,1,1)+ty(ijk,1,1) 
      enddo
   endif
!
   do j=1,ny 
   do i=1,nx 
      pp(i,j,1)=tx(i,j,1) 
   enddo 
   enddo 
   return  
end subroutine divbric1
!
!****************************************************************
!
subroutine temporellebis(pp,ux,uy,uz,tr,us,ps,&
     xk2,yk2,zk2,xkx,yky,zkz,nx,ny,nz,&
     nxm,nym,nzm,mx,my,mz,nwork,ntrigsX,ntrigsY)
!
!****************************************************************
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
  USE paramx_m 
  USE paramy_m 
  USE paramz_m 
  USE derivex_m
  USE derivey_m
!...Translated by PSUITE Trans90                  4.3ZH 16:05:14   1/27/ 4  
!...Switches: -yv INDDO=0 -nbejkno
  implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  integer :: nx,nxm,nwork,ntrigsX,ntrigsY 
  integer :: ny,nym 
  integer :: nz,nzm 
  integer :: mx 
  integer :: my 
  integer :: mz 
  real(8) , intent(inout) :: pp(nx,ny,nz) 
  real(8) :: ux(nx,ny,nz) 
  real(8) :: uy(nx,ny,nz) 
  real(8) :: uz(nx,ny,nz) 
  real(8) :: tr(mx,my,mz) 
  real(8) :: us(mx,my,mz) 
  real(8) :: ps(mx,my,mz) 
  real(8) :: xk2(mx) 
  real(8) :: yk2(my) 
  real(8) :: zk2(mz) 
  real(8) :: xkx(mx) 
  real(8) :: yky(my) 
  real(8) :: zkz(mz) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer :: mxyz, mxy, j, i ,k
  real(8) , dimension(mx,my) :: ps2d, us2d 
  real(8) :: xyzk,pi, twopi,xxk1,w,wp,x,y 
!-----------------------------------------------!
!
  mxyz=mx*my*mz
  mxy=mx*my
  pi=acos(-1.)
  twopi=2.*pi
  xxk1=pi/xlx

  call annule (ps2d,mxy)
  call annule (ps,mxyz)
  call annule (us,mxyz)
  call annule (us2d,mxy)

  do k=1,nz
  do j=1,ny
     y=(j-1)*dy
     do i=1,nx
        x=(i-1)*dx
!      pp1(i,j,k)=cos(4.*xxk1*x)*cos(6.*xxk1*y)
!      pp(i,j,k)=-16.*xxk1*xxk1*cos(4.*xxk1*x)*cos(6.*xxk1*y)&
!              -36.*xxk1*xxk1*cos(6.*xxk1*y)*cos(4.*xxk1*x)
     enddo
   enddo
   enddo

   do j=1,my
   do i=1,mx
      ps2d(i,j)=0.
   enddo
   enddo
   do j=1,ny
   do i=1,nx
      ps2d(i,j)=pp(i,j,1)
   enddo
   enddo
!
   call SLFFT2D(ps2d,nx,ny,mx,my,nxm,nym,nwork,&
        ntrigsX,ntrigsY,nclx,ncly,-1)
   if ((nclx==1).or.(nclx==2)) then
      do i=1,mx-2
      w=twopi*0.5*(i-1)/(mx-2)
      wp=afix*2.*dx*sin(w)+(bfix*2.*dx)*sin(2.*w)
      wp=wp/(1.+2.*alfaix*cos(w))
!!         wp=w
      xkx(i)=(mx-2)*wp/xlx
      xk2(i)=xkx(i)*xkx(i)
   enddo
   xkx(mx-1)=0.
   xk2(mx-1)=0.
   xkx(mx  )=0.
   xk2(mx  )=0.
   endif
   if (nclx==0) then
      do i=1,mx-2,2
         w=twopi*0.5*(i-1)/(mx-2)
         wp=afix*2.*dx*sin(w)+bfix*2.*dx*sin(2.*w)
         wp=wp/(1.+2.*alfaix*cos(w))
         xkx(i)=(mx-2)*wp/xlx
         xkx(i+1)=xkx(i)
         xk2(i)=xkx(i)*xkx(i)
         xk2(i+1)=xk2(i) 
!         write(*,10) i,w,xkx(i),(mx-2)*w/xlx
      enddo
      xkx(mx-1)=0.
      xk2(mx-1)=0.
      xkx(mx  )=0.
      xk2(mx  )=0.
   endif
!          
   if (ncly==0) then
      do j=1,my-2,2
         w=twopi*0.5*(j-1)/(my-2)
         wp=afjy*2.*dy*sin(w)+bfjy*2.*dy*sin(2.*w)
         wp=wp/(1.+2.*alfajy*cos(w))
         yky(j)=(my-2)*wp/yly
         yky(j+1)=yky(j)
         yk2(j)=yky(j)*yky(j)
         yk2(j+1)=yk2(j) 
!         write(*,10) i,w,xkx(i),(mx-2)*w/xlx
      enddo
      yky(my-1)=0.
      yk2(my-1)=0.
      yky(my  )=0.
      yk2(my  )=0.
   endif 
   if ((ncly==1).or.(ncly==2)) then
      do j=1,my-2
         w=twopi*0.5*(j-1)/(my-2)
         wp=afjy*2.*dy*sin(w)+bfjy*2.*dy*sin(2.*w)
         wp=wp/(1.+2.*alfajy*cos(w))
!!         wp=w
         yky(j)=(my-2)*wp/yly
         yk2(j)=yky(j)*yky(j)
      enddo
      yky(my-1)=0.
      yk2(my-1)=0.
      yky(my)=0.
      yk2(my)=0.
   endif
!
   do j=1,ny
   do i=1,nx
      xyzk=xk2(i)+yk2(j)
      if (xyzk.eq.0.) then
!              print *,'mode nul : ',i,j,ps2d(i,j)
         us2d(i,j)=0.
      else
         us2d(i,j)=ps2d(i,j)/(-xyzk)
      endif
   enddo
   enddo
!
           
!      do j=1,ny
!      do i=1,nx
!         print *,i,j,us2d(i,j),ps2d(i,j)
!      enddo
!      print *,'***********************'
!      enddo
!      stop
!   

   call SLFFT2D(us2d,nx,ny,mx,my,nxm,nym,nwork,&
        ntrigsX,ntrigsY,nclx,ncly,1)
!
   do k=1,nz
   do j=1,ny
   do i=1,nx
      pp(i,j,k)=us2d(i,j)
   enddo
   enddo
   enddo

!===============================================================
   return
end subroutine temporellebis
!
!*******************************************************************
!
   subroutine bandec(a,n,m,m1,m2,al,indx,xd)
!
!*******************************************************************
!
     implicit none
!
     real(8),parameter :: tiny=1.e-20
     real(8),dimension(n/2,m1) :: al
     real(8),dimension(n/2,m1+m2+1) :: a
     real(8),dimension(n/2) :: indx     
     integer :: n,m1,m2,l,mm,i,k,j,m
     real(8) :: xd,dum
   
     
      mm=m1+m2+1
      l=m1
      do i=1,m1
         do j=m1+2-i,mm
            a(i,j-l)=a(i,j)
         enddo
         l=l-1
         do j=mm-l,mm
            a(i,j)=0.
         enddo
      enddo
      xd=1.
      l=m1
      do k=1,n/2
         dum=a(k,1)
!         write(*,*)'dum',dum
         i=k
         if (l.lt.n/2) l=l+1
         do j=k+1,l
            if (abs(a(j,1)).gt.abs(dum)) then
               dum=a(j,1)     
               i=j
            endif
         enddo
         indx(k)=i
         if (dum.eq.0.) a(k,1)=tiny
         if (i.ne.k) then
           xd=-xd
            do j=1,mm
               dum=a(k,j)
               
               a(k,j)=a(i,j)
               a(i,j)=dum
              
            enddo
         endif
         do i=k+1,l
            dum=a(i,1)/a(k,1)
            al(k,i-k)=dum
            do j=2,mm
               a(i,j-1)=a(i,j)-dum*a(k,j)
            enddo
            a(i,mm)=0.
         enddo
      enddo
      return
    end subroutine bandec
!*******************************************************
!*******************************************************
!
    subroutine banbks(a,n,m,m1,m2,al,indx,b)
!
!*******************************************************
!
      implicit none
!
     real(8),parameter :: tiny=1.e-20
     real(8),dimension(n/2,m1) :: al
     real(8),dimension(n/2,m1+m2+1) :: a
     real(8),dimension(n/2) :: indx
     real(8),dimension(n/2) :: b
     integer :: n,m1,m2,l,mm,m,i,k,j
     real(8) :: xd,dum
!      
      mm=m1+m2+1
      l=m1
      do k=1,n/2
         i=indx(k)
         if (i.ne.k) then
            dum=b(k)
            b(k)=b(i)
            b(i)=dum
         endif
         if (l.lt.n/2) l=l+1
         do i=k+1,l
            b(i)=b(i)-al(k,i-k)*b(k)
         enddo
      enddo
      l=1
      do i=n/2,1,-1
         dum=b(i)
         do k=2,l
            dum=dum-a(i,k)*b(k+i-1)
         enddo
         b(i)=dum/a(i,1)
         if (l.lt.mm) l=l+1
      enddo
      return
    end subroutine banbks
!*******************************************************
!*******************************************************
!
    subroutine banbks1(a,n,m,m1,m2,al,indx,b,x)
!
!*******************************************************
!
      implicit none
!
     real(8),dimension(n/2,m1) :: al
     real(8),dimension(n/2,m1+m2+1) :: a
     real(8),dimension(n/2) :: indx
     real(8),dimension(n/2) :: b,x
     integer,dimension(n/2) :: lv,lu
     integer :: n,m1,m2,l,mm,m,i,k,j
     real(8) :: piv1,piv2
!      
      mm=m1+m2+1
!
      do i=1,n/2
         x(i)=b(i)
      enddo
!
      l=m1
      do k=1,n/2
         if (l.lt.n/2) l=l+1
         lv(k)=l
      enddo
!
      do k=1,n/2
         j=indx(k) ! si (j.ne.k) se cumple (j-k)=1 
         piv1=(j-k)*x(j)+(k-j+1)*x(k)
         piv2=(j-k)*x(k)+(k-j+1)*x(j)
         x(k)=piv1
         x(j)=piv2
         do i=k+1,lv(k)
            x(i)=x(i)-al(k,i-k)*x(k)
         enddo
      enddo
!
      l=1
      do i=n/2,1,-1
         lu(i)=l
         if (l.lt.mm) l=l+1
      enddo
!
      do i=n/2,1,-1
         do k=2,lu(i)
            x(i)=x(i)-a(i,k)*x(k+i-1)
         enddo
         x(i)=x(i)/a(i,1)
      enddo
!
      return
    end subroutine banbks1
!****************************************************************
!****************************************************************
!
subroutine temporellebisSTR(pp,ux,uy,uz,tr,us,ps,&
     xk2,yk2,zk2,xkx,yky,zkz,&
     yp,ja,jb,d,d1,d5,a,al,indx,&
     a1,al1,indx1,e,c,a2,alpha,beta,m1,m2,nx,ny,nz,&
     nxm,nym,nzm,mx,my,mz,nwork,ntrigsX,ntrigsY)
!
!****************************************************************
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
  USE paramx_m 
  USE paramy_m 
  USE paramz_m 
  USE derivex_m
  USE derivey_m
!...Translated by PSUITE Trans90                  4.3ZH 16:05:14   1/27/ 4  
!...Switches: -yv INDDO=0 -nbejkno
  implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  integer :: nx,nxm,nwork,ntrigsX,ntrigsY,ny,nym,nz,nzm,mx,my,mz,m1,m2
  real(8),dimension(nx,ny,nz) :: pp,ux,uy,uz,tr,us,ps 
  real(8),dimension(mx) :: xkx,xk2 
  real(8),dimension(my) :: yky,yk2 
  real(8),dimension(mz) :: zkz,zk2 
  real(8),dimension(ny/2,m1+m2+1) :: a,a1,a2
  real(8),dimension(ny/2,m1) :: al,al1
  real(8),dimension(ny/2) :: indx,indx1,e,c
  real(8),dimension(nx) :: d5
  real(8),dimension(ny) :: d1,d
  real(8),dimension(2) :: ja,jb
  real(8) :: beta,den,xnum,alpha,xcx,yp
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer :: mxyz, mxy, j, i ,k,ip,npaire
  real(8) , dimension(mx,my) :: ps2d, us2d ,ps2d1
  real(8) :: xmo1,xd,xd1,xtoto
  real(8) :: xyzk,pi,twopi,xxk1,w,wp,x,y,xvar,xvar1,beta1,tw,xmo,xps 
!-----------------------------------------------!
!
!*****************ATTENTION*******************
  npaire=1
!*********************************************
  mxyz=mx*my*mz
  mxy=mx*my
  pi=acos(-1.)
  twopi=2.*pi
  xxk1=pi/xlx
!
  call annule (ps2d,mxy)
  call annule (ps,mxyz)
  call annule (us,mxyz)
  call annule (us2d,mxy)

   do j=1,my
   do i=1,mx
      ps2d(i,j)=0.
   enddo
   enddo
   do j=1,ny
   do i=1,nx
      ps2d(i,j)=pp(i,j,1)
   enddo
   enddo
!
   call slfft2d(ps2d,nx,ny,mx,my,nxm,nym,nwork,&
        ntrigsX,ntrigsY,nclx,ncly,-1)
!      do i=1,nx
!      do j=1,ny
!         print '(i4,i4,f15.8,f15.8)',i,j,ps2d(i,j)
!      enddo
!      pause
!      enddo
   if ((nclx==1).or.(nclx==2)) then
      do i=1,mx-2
      w=twopi*0.5*(i-1)/(mx-2)
      wp=afix*2.*dx*sin(w)+(bfix*2.*dx)*sin(2.*w)
      wp=wp/(1.+2.*alfaix*cos(w))
!!      wp=w
      xkx(i)=(mx-2)*wp/xlx
      xk2(i)=xkx(i)*xkx(i)
   enddo
   xkx(mx-1)=0.
   xk2(mx-1)=0.
   xkx(mx  )=0.
   xk2(mx  )=0.
   endif
!          
   if ((ncly==1).or.(ncly==2)) then
      do j=1,my-2
         w=twopi*0.5*(j-1)/(my-2)
         wp=afjy*2.*dy*sin(w)+bfjy*2.*dy*sin(2.*w)
         wp=wp/(1.+2.*alfajy*cos(w))
!!         wp=w
         yky(j)=(my-2)*wp
         yk2(j)=yky(j)*yky(j)
      enddo
      yky(my-1)=0.
      yk2(my-1)=0.
      yky(my)=0.
      yk2(my)=0.
   endif
!*******************STRET STRET STRET***********************************
!*****************petit tablo d onde dans b 
      do j=1,ny/2
         d(j)=yky(2*j-1)
         d1(j)=yky(2*j)
      enddo 
!      do j=1,ny/2
!         print *,'d',d(j),d1(j)
!      enddo
!      pause
!********************nombre d onde dans d       
   if (nclx==0) then
      do i=1,nx,2
         ip=i+1
         d5(i)=xkx(i)
         d5(ip)=d5(i)
      enddo
   endif
   if ((nclx==1).or.(nclx==2)) then
      do i=1,nx
         d5(i)=xkx(i)
      enddo
   endif
!      do i=1,nx
!         print *,'toto',d5(i)
!      enddo
!      pause
!*****************DEBUT DE LA BOUCLE EN XXXXXXXXXXXX******************
       do i=1,nx
!*********************************************************************
         xps=ps2d(i,ny)
         do j=1,ny/2
            e(j)=ps2d(i,2*j-1)
            c(j)=ps2d(i,2*j)
         enddo
!***************Construction de la matrice a**************************
      xvar=alpha+1./2./beta
      xmo=alpha+3./4./beta
      tw=twopi*twopi
      beta1=beta*beta
      xvar1=xvar*xvar
      a(1,m1+1)=-(d5(i)*d5(i))-d(1)*d(1)*xvar1/tw*4.&
           -(1./4./tw/beta1)*(d(1)*d(1+1))
      do j=2,ny/2-1
         a(j,m1+1)=-(d5(i)*d5(i))-d(j)*d(j)*xvar1/tw*4.&
              -(1./4./tw/beta1)*(d(j)*d(j-1)+d(j)*d(j+1))
      enddo
      a(ny/2,m1+1)=-(d5(i)*d5(i))-d(ny/2)*d(ny/2)*xvar1/tw*4.&
           -(1./4./tw/beta1)*(d(ny/2)*d(ny/2-1))
!****************element du cote des m1*******************************
      a(1,m1)=0.
      do j=2,ny/2
         a(j,m1)=xvar/tw/beta*(d(j)*d(j-1)+d(j-1)*d(j-1))
      enddo
!
      a(1,m1-1)=0.
      a(2,m1-1)=0.
      do j=3,ny/2
         a(j,m1-1)=-d(j-1)*d(j-2)/4./tw/beta1
      enddo
!****************element du cote des m2*******************************
      a(1,4)=xvar/tw/beta*(d(1)*d(2)+d(2)*d(2))
      do j=2,ny/2-1
         a(j,m1+m2)=xvar/tw/beta*(d(j)*d(j+1)+d(j+1)*d(j+1))
      enddo
      a(1,m1+m2)=a(1,m1+m2)*2.
!      a(ny/2-1,m1+m2)=xvar/tw/beta*(d(ny/2-1)*d(ny/2))
!      a(ny/2-1,m1+m2)=a(ny/2-1,m1+m2)+xmo/tw/beta*(d(ny/2)*d(ny/2))
      a(ny/2,m1+m2)=0.
!        
      do j=1,ny/2-2
         a(j,m1+m2+1)=-d(j+1)*d(j+2)/4./tw/beta1
      enddo
      a(1,m1+m2+1)=a(1,m1+m2+1)*2.
      a(ny/2-1,m1+m2+1)=0.
      a(ny/2,m1+m2+1)=0.
!********************************************************************
      if (d5(i).eq.0.) then 
         a(1,m1+1)=1.
         a(1,m1+m2)=0.
         a(1,m1+m2+1)=0.
!         a(ny/2,m1+1)=1.
!         a(ny/2,m1+m2)=0.
!         a(ny/2,m1+m2+1)=0.
      endif
!*********************************************************************
!*******construction de la matrice a2*********************************
      xvar=alpha+1./2./beta
      xmo=alpha+3./4./beta
      xmo1=alpha+1./4./beta
      tw=twopi*twopi
      beta1=beta*beta
      xvar1=xvar*xvar
      a2(1,m1+1)=-(d5(i)*d5(i))-d1(1)*d1(1)*xmo1*xmo/tw*4.&
           -(1./4./tw/beta1)*(d1(1)*d1(2))
      do j=2,ny/2-1
         a2(j,m1+1)=-(d5(i)*d5(i))-d1(j)*d1(j)*xvar1/tw*4.&
              -(1./4./tw/beta1)*(d1(j)*d1(j-1)+d1(j)*d1(j+1))
      enddo
      a2(ny/2,m1+1)=-(d5(i)*d5(i))-d1(ny/2)*d1(ny/2)*xmo1*xmo/tw*4.&
           -(1./4./tw/beta1)*(d1(ny/2)*d1(ny/2-1))
!***************element du cote des m1********************************
      a2(1,m1)=0.
      if (npaire.eq.0) then
         a2(2,m1)=xvar/tw/beta*(d1(2)*d1(1))&
              +xmo1/tw/beta*(d1(1)*d1(1))
      endif
      if (npaire.eq.1) then
         a2(2,m1)=xvar/tw/beta*(d1(2)*d1(1))&
              +xmo/tw/beta*(d1(1)*d1(1))
      endif
      do j=3,ny/2-1
         a2(j,m1)=xvar/tw/beta*(d1(j)*d1(j-1)+d1(j-1)*d1(j-1))
      enddo
      if (npaire.eq.0) then
         a2(ny/2,m1)=xmo/tw/beta*(d1(ny/2)*d1(ny/2-1))&
              +xvar/tw/beta*(d1(ny/2-1)*d1(ny/2-1))
      endif
      if (npaire.eq.1) then
         a2(ny/2,m1)=xmo1/tw/beta*(d1(ny/2)*d1(ny/2-1))&
              +xvar/tw/beta*(d1(ny/2-1)*d1(ny/2-1))
      endif
!
      a2(1,m1-1)=0.
      a2(2,m1-1)=0.
      do j=3,ny/2
         a2(j,m1-1)=-d1(j-1)*d1(j-2)/4./tw/beta1
      enddo
!****************element du cote des m2*********************************
      if (npaire.eq.0) then
         a2(1,m1+m2)=xmo/tw/beta*(d1(1)*d1(2))&
              +xvar/tw/beta*(d1(2)*d1(2))
      endif
      if (npaire.eq.1) then
         a2(1,m1+m2)=xmo1/tw/beta*(d1(1)*d1(2))&
              +xvar/tw/beta*(d1(2)*d1(2))
      endif
      do j=2,ny/2-2
         a2(j,m1+m2)=xvar/tw/beta*(d1(j)*d1(j+1)+d1(j+1)*d1(j+1))
      enddo
      if (npaire.eq.0) then
         a2(ny/2-1,m1+m2)=xvar/tw/beta*(d1(ny/2-1)*d1(ny/2))&
              +xmo1/tw/beta*(d1(ny/2)*d1(ny/2))
      endif
      if (npaire.eq.1) then
         a2(ny/2-1,m1+m2)=xvar/tw/beta*(d1(ny/2-1)*d1(ny/2))&
              +xmo/tw/beta*(d1(ny/2)*d1(ny/2))
      endif
      a2(ny/2,m1+m2)=0.
      do j=1,ny/2-2
         a2(j,m1+m2+1)=-d1(j+1)*d1(j+2)/4./tw/beta1
      enddo
      a2(ny/2-1,m1+m2+1)=0.
      a2(ny/2,m1+m2+1)=0.
!***********************************************************************
!***********************************************************************
       call bandec(a,ny,my,m1,m2,al,indx,xd)
       call bandec(a2,ny,my,m1,m2,al1,indx1,xd1)
       call banbks(a,ny,my,m1,m2,al,indx,e)
       call banbks(a2,ny,my,m1,m2,al1,indx1,c)
!********************on reconstruit les a+ib****************************
       do j=1,ny-1,2
          us2d(i,j)=e((j+1)/2)
       enddo
       do j=2,ny,2
          us2d(i,j)=c(j/2)
       enddo 
!***********************************************************************
       if (d5(i).ne.0.) then
          xtoto=xps+e(ny/2-1)*d(ny/2)*d(ny/2-1)/2./beta1/tw&
               -e(ny/2)*d(ny/2)*d(ny/2)*2.*xvar/tw/beta
          xtoto=-xtoto/d5(i)/d5(i)
       endif     
       if (d5(i).eq.0) then
          us2d(i,1)=0.
          us2d(i,ny)=0.
       endif
       us2d(i,ny)=xtoto
!       print *,xps
!       pause
!       if (i==1) us2d(1,ny)=-xps
!*******FIN BOUCLE EN XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXxxxxx*************
       enddo
!***********************************************************************
!***********************************************************************
!   do j=1,ny
!   do i=1,nx
!      xyzk=xk2(i)+yk2(j)
!      if (xyzk.eq.0.) then
!              print *,'mode nul : ',i,j,ps2d(i,j)
!         us2d(i,j)=0.
!      else
!         us2d(i,j)=ps2d(i,j)/(-xyzk)
!      endif
!   enddo
!   enddo
!
           !
!     do i=1,nx
!     do j=1,ny
!        print *,i,j,us2d(i,j)
!     enddo
 !    pause
 !    enddo
!      stop
!   
   call slfft2d(us2d,nx,ny,mx,my,nxm,nym,nwork,&
        ntrigsX,ntrigsY,nclx,ncly,1)
!
   do k=1,nz
   do j=1,ny
   do i=1,nx
      pp(i,j,k)=us2d(i,j)
   enddo
   enddo
   enddo

!===============================================================
   return
end subroutine temporellebisSTR
!
!*******************************************************************
