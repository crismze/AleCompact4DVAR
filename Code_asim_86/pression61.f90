!********************************************************************
!
subroutine divbis26 (ppm,ux,uy,uz,tx,ty,tz,tr,di,sx,sy,sz,t,r,&
     cfx6,csx6,cwx6,&
     cfxp6,csxp6,cwxp6,&
     nx,nxm,mx,ny,nym,my,nz,nvec,&
     cifx6,cisx6,ciwx6,cifxp6,cisxp6,ciwxp6,&
     cfy6,csy6,cwy6,cfyp6,csyp6,cwyp6,&
     cify6,cisy6,ciwy6,cifyp6,cisyp6,ciwyp6,yp,ypi)    
! 
!********************************************************************
!
   USE paramx_m
   USE paramy_m
   USE paramz_m
!
   implicit none
!
   integer :: i,j,k,nx,nxm,mx,ny,nym,my,nz,nvec
   real(8),dimension(nx,ny,nz) :: ux,uy,uz,tz,di,tr,df2
   real(8),dimension(nxm,nym,nz) :: tx,ppm,ty,tx2,ty2,ppm1
   real(8),dimension(nx,nym,nz) :: ty1,tx3
   real(8),dimension(nxm,ny,nz) :: tx1,df3
   real(8),dimension(ny,nz) :: sx
   real(8),dimension(nx,nz) :: sy
   real(8),dimension(nx,ny) :: sz
   real(8),dimension(nvec) :: t,r
   real(8),dimension(my) :: yp,ypi
!
   real(8),dimension(nxm) :: cfx6,csx6,cwx6,cifx6,cisx6,ciwx6
   real(8),dimension(nxm) :: cifxp6,cisxp6,ciwxp6,cfxp6,csxp6,cwxp6
   real(8),dimension(nym) :: cfy6,csy6,cwy6,cify6,cisy6,ciwy6
   real(8),dimension(nym) :: cifyp6,cisyp6,ciwyp6,cfyp6,csyp6,cwyp6   
!
   real(8) :: tmax,tmoy,pi,x,y,xl2,xxk1,xxk2
   integer :: imax,jmax,kmax,nxyz,ijk,l
!

!
   tmax=0.
   tmoy=0.
   imax=0
   jmax=0
   kmax=0
   do i=1,nvec
      t(i)=0.
      r(i)=0.
   enddo
!

   call intery6(ty1,ux,tr,di,sy,cifyp6,cisyp6,ciwyp6,nx,ny,nym,nz,1)
!   
   call decx6(tx,ty1,tr,sx,cfx6,csx6,cwx6,nx,nxm,ny,nym,nz,0)  
!     
   call inter6(tx1,uy,tr,sx,cifxp6,cisxp6,ciwxp6,nx,nxm,ny,nz,1)

!      
   call decy6(ty,tx1,tr,di,sy,cfy6,csy6,cwy6,nx,nxm,ny,nym,nz,0)
!


   do j=1,nym
   do i=1,nxm
      ppm(i,j,1)=tx(i,j,1)+ty(i,j,1)
   enddo
   enddo
!----------------------------------------------        
   do j=1,nym
   do i=1,nxm
      tmoy=tmoy+abs(ppm(i,j,1))
      tmax=(tmax+ppm(i,j,1))/2.+abs(tmax-ppm(i,j,1))/2.
   enddo
   enddo
   kmax=1
!
   tmoy=tmoy/(nxm*nym*nz)
   print *,'DIV U = ',tmax,tmoy
!
   return
end subroutine divbis26
!
!********************************************************************
!
      subroutine divbis36 (ppm,ux,uy,uz,tx,ty,tz,tr,di,sx,sy,sz,t,r,&
                     cfx6,csx6,cwx6,&
                     cfxp6,csxp6,cwxp6,&
                     nx,nxm,mx,ny,nym,my,nz,nvec,&
                     cifx6,cisx6,ciwx6,cifxp6,cisxp6,ciwxp6,&
                     cfy6,csy6,cwy6,cfyp6,csyp6,cwyp6,&
                     cify6,cisy6,ciwy6,cifyp6,cisyp6,ciwyp6,&
                     nwork,ntrigsX,ntrigsY,yp,ypi)    
! 
!********************************************************************
!
!
   USE paramx_m
   USE paramy_m
   USE paramz_m
!
   implicit none
!
   real(8),dimension(nx,ny,nz) :: ux,uy,uz,di,tr,tz,df2
   real(8),dimension(nxm,ny,nz) :: tx1,df3,tx3
   real(8),dimension(nx,nym,nz) :: ty1
   real(8),dimension(nxm,nym,nz) :: ppm,ppm1,ty,tx,tx2,ty2
   real(8),dimension(ny,nz) :: sx 
   real(8),dimension(nx,nz) :: sy
   real(8),dimension(nx,ny) :: sz
   real(8),dimension(nvec) :: t,r
   real(8),dimension(nxm) :: cfx6,ccx6,cwx6,cifx6,cisx6,ciwx6,csx6
   real(8),dimension(nxm) :: cfxp6,csxp6,cwxp6,cifxp6,cisxp6,ciwxp6 
   real(8),dimension(nym) :: cfy6,ccy6,cwy6,cify6,cisy6,ciwy6,csy6
   real(8),dimension(nym) :: cfyp6,csyp6,cwyp6,cifyp6,cisyp6,ciwyp6 
   real(8),dimension(mx,my) :: us2d,ps2d,ux2d,uy2d
   real(8),dimension(nx) :: xx1
   real(8),dimension(ny) :: yy1
   real(8),dimension(ny+1) :: yp,ypi
   real(8) ::tmax,tmoy,pi,x,y,xl2
   integer :: nx,ny,nz,nxm,nym,nvec,i,j,k,imax,jmax,kmax,ii,jj
   integer :: nyxz,ijk,mx,my,mxsz,nwork,ntrigsX,ntrigsY
!

   tmax=0.
   tmoy=0.
   imax=0
   jmax=0
   kmax=0
   do i=1,nvec
      t(i)=0.
      r(i)=0.
   enddo
!
   call intery6(ty1,ux,tr,di,sy,cifyp6,cisyp6,ciwyp6,nx,ny,nym,nz,1)
!  
   call decx6(tx,ty1,tr,sx,cfx6,csx6,cwx6,nx,nxm,ny,nym,nz,0) 
!     
   call inter6(tx1,uy,tr,sx,cifxp6,cisxp6,ciwxp6,nx,nxm,ny,nz,1)
!      
   call decy6(ty,tx1,tr,di,sy,cfy6,csy6,cwy6,nx,nxm,ny,nym,nz,0)
!

   do j=1,nym
   do i=1,nxm
      ppm(i,j,1)=tx(i,j,1)+ty(i,j,1)
   enddo
   enddo
!

!--------------------------------------------
!----------------------------------------------
   mxsz=mx*my
   call annule(us2d,mxsz)
   call annule(ps2d,mxsz)
   do j=1,nym
   do i=1,nxm
      us2d(i,j)=ppm(i,j,1)
      ps2d(i,j)=ppm1(i,j,1)
   enddo
   enddo
   call slfft2d(us2d,nxm,nym,mx,my,nxm,nym,nwork,&
        ntrigsX,ntrigsY,0,0,-1)
       x=0.
       do i=1,nx
          xx1(i)=x
          x=x+dx 
       enddo
       y=0.
       do j=1,ny
          yy1(j)=y
          y=y+dy
       enddo
!       open(30,file='vort',form='formatted',status='unknown')
!         do j=1,nym
!         do i=1,nxm
!         write(30,*)  xx1(i),yy1(j),ppm(i,j,1)
!         enddo
!         write(30,*)
!         enddo
!       close(30)
!      pause
!   do j=1,my
!   do i=1,mx
!      print *,i,j,us2d(i,j)
!   enddo
!   pause
!   enddo
 !  pause
!   call slfft2d(ps2d,nx,ny,mx,my,nxm,nym,nwork,&
!        ntrigsX,ntrigsY,nclx,ncly,-1)
!      open (10,file='courbe.dat',form='formatted',status='unknown')
!      do j=1,my
!         y=(j-1)*dy
!         write(10,*) y,us2d(1,j)
!      enddo
!      close(10)  
!      pause'11'
      us2d(1,1)=0.
!      us2d(nx,1)=0.
!      us2d(nx,ny)=0.
!      do i=1,mx
!         us2d(i,my-1)=0.
!        us2d(i,my)=0.
!      enddo


!   tmax=-1.e-14
!   ii=0.
!   jj=0.
!   do j=1,my
!   do i=1,mx
!      if (us2d(i,j)>tmax) then 
!         tmax=us2d(i,j) 
!         ii=i
!         jj=j
!      endif
!   enddo
!   enddo
!   print *,ii,jj,tmax
!   pause 'tot'
   
!      do i=1,mx
!         us2d(i,ny+1)=0.
!      enddo
!   do j=1,my
!   do i=1,mx
!      if (us2d(i,j)>1.e-13) print *,i,j,ppm(i,j,1)
!   enddo
!   enddo
!   pause

   call slfft2d(us2d,nxm,nym,mx,my,nxm,nym,nwork,&
        ntrigsX,ntrigsY,0,0,1)
   do k=1,nz
   do j=1,nym
   do i=1,nxm
!     ppm(i,j,k)=us2d(i,j)
   enddo
   enddo
   enddo
!
   do j=1,nym
   do i=1,nxm
      tmoy=tmoy+abs(ppm(i,j,1))
      tmax=(tmax+ppm(i,j,1))/2.+abs(tmax-ppm(i,j,1))/2.
   enddo
   enddo
   kmax=1
!
   tmoy=tmoy/(nxm*nym*nz)
   print *,'DIV U = ',tmax,tmoy
!
   return
end subroutine divbis36
!
!********************************************************************
!
subroutine gradpbis16(ppm,px,py,py1,pz,tr,di,sx,sy,sz,&
          cfx6,csx6,cwx6,&
          cfxp6,csxp6,cwxp6,&
          nx,nxm,mx,ny,nym,my,nz,gdt,cifx6,cisx6,ciwx6,&
          cifxp6,cisxp6,ciwxp6,&
          cfi6,cci6,cbi6,cfip6,ccip6,cbip6,&
          csip6,cwip6,csi6,cwi6,&
          cifi6,cici6,cibi6,cifip6,cicip6,&
          cibip6,cisip6,ciwip6,cisi6,ciwi6,cify6,cisy6,ciwy6,&
          cifyp6,cisyp6,ciwyp6,cfyp6,csyp6,cwyp6,&
          cfi6y,cci6y,cbi6y,cfip6y,ccip6y,cbip6y,&
          csip6y,cwip6y,csi6y,cwi6y,&
          cifi6y,cici6y,cibi6y,cifip6y,cicip6y,&
          cibip6y,cisip6y,ciwip6y,cisi6y,ciwi6y,cfy6,csy6,cwy6,yp,ypi)
!
!********************************************************************
!
   USE paramx_m
   USE paramy_m
!
   implicit none
!
   integer :: nx,nxm,mx,ny,nym,my,nz
   real(8),dimension(nx,ny,nz) :: px,py,pz,pp,tr,di,py1,uz,ux,uy
   real(8),dimension(nxm-1,ny,nz) :: df2
   real(8),dimension(nxm,nym,nz) :: ppm,ppm1
   real(8),dimension(ny,nz) :: sx
   real(8),dimension(nx,nz) :: sy
   real(8),dimension(nx,ny) :: sz
   real(8),dimension(nxm) :: cfx6,csx6,cwx6,cifx6,cisx6,ciwx6
   real(8),dimension(nxm) :: cifxp6,cisxp6,ciwxp6,cfxp6,csxp6,cwxp6
   real(8),dimension(nx) :: cfi6,cci6,cbi6,cfip6,ccip6,cbip6
   real(8),dimension(nx) :: csip6,cwip6,csi6,cwi6,cifi6,cici6,cibi6
   real(8),dimension(nx) :: cifip6,cicip6,cibip6,cisip6,ciwip6,cisi6,ciwi6
   real(8),dimension(ny) :: cfi6y,cci6y,cbi6y,cfip6y,ccip6y,cbip6y,cfy6,csy6,cwy6
   real(8),dimension(ny) :: csip6y,cwip6y,csi6y,cwi6y,cifi6y,cici6y,cibi6y
   real(8),dimension(ny) :: cifip6y,cicip6y,cibip6y,cisip6y,ciwip6y,cisi6y,ciwi6y
   real(8),dimension(nym) :: cify6,cisy6,ciwy6,cifyp6,cisyp6,ciwyp6
   real(8),dimension(nym) :: cfyp6,csyp6,cwyp6
   real(8),dimension(ny+1) :: yp,ypi
   real(8),dimension(mx,my,nz) :: wk1,wk2
   integer :: i,j,k
   real(8) :: gdt,pi,x,y,xl2,xxk1,xxk2
     
!   pi=acos(-1.)
!   xxk1=pi/xlx
!   xxk2=pi/yly
!   do k=1,nz
!   do j=1,ny
!      y=yp(j)
!      do i=1,nx
!         x=(i-1)*dx
!         ux(i,j,k)=-3.*xxk1*sin(3.*xxk1*x)*cos(3.*xxk2*y)
!         uy(i,j,k)=-3.*xxk2*sin(3.*xxk2*y)*cos(3.*xxk1*x)      
!      enddo
!   enddo
!   enddo
!   do k=1,nz
!   do j=1,nym
!      y=ypi(j)
!      do i=1,nxm
!         x=(i-0.5)*dx
!         ppm(i,j,k)=cos(3.*xxk2*y)*cos(3.*xxk1*x)
!      enddo
!   enddo
!   enddo

!   print *,nx,nxm,ny,nym,mx,my
!   pause

   if (ncly==1) then
      call interiy6(wk1,ppm,cifip6y,cisip6y,ciwip6y,nxm,nym,ny,nz,1)
   endif
   if ((nclx==1).or.(nclx==2)) then
      call deci6(px,wk1,tr,sx,cfip6,csip6,cwip6,nxm,nx,ny,nz,1)
   endif 
!************************************
  wk1(:,:,:)=0. ; wk2(:,:,:)=0.
!************************************
   if ((nclx==1).or.(nclx==2)) then
      call interi6(wk1,ppm,tr,sx,cifip6,cisip6,ciwip6,nxm,nx,nym,nz,1)
   endif
   if (ncly==1) then
      call deciy6(py,wk1,tr,di,sy,cfip6y,csip6y,cwip6y,nx,nym,ny,nz,1)
   endif
         do k=1,nz
         do j=1,ny
         do i=1,nx
            py1(i,j,k)=py(i,j,k)/gdt 
         enddo
         enddo
         enddo


!   xl2=0.
!   do j=1,ny
!   do i=1,nx
!      xl2=xl2+(py(i,j,1)-uy(i,j,1))*(py(i,j,1)-uy(i,j,1))
!   enddo
!   enddo
!   xl2=xl2/nx/ny
!   write(*,*)'erreur inter toto',sqrt(xl2) 
!      open (10,file='courbe.dat',form='formatted',status='unknown')
!      do j=1,ny
!         write(10,*) yp(j),uy(10,j,1),py(10,j,1)
!      enddo
!      close(10)
!   do i=1,nx
!   do j=1,ny
!     print '(i4,f15.8,f15.8,f15.8)',i,ux(i,j,1),px(i,j,1)
!   enddo
!   pause
!   enddo  
!!   stop

   return
end subroutine gradpbis16
!
!****************************************************************
!
subroutine temporellebis16STR(pp,ppm,ux,uy,uz,tr,us,ps,&
     xk2,yk2,zk2,xkx,yky,zkz,&
     yp,ypi,ja,jb,d,d1,d5,a,al,indx,&
     a1,al1,indx1,e,c,a2,alpha,beta,m1,m2,&
     nx,nxm,ny,nym,nz,mx,my,mz,&
     nwork,ntrigsX,ntrigsY)
!
!****************************************************************
!
   USE paramx_m 
   USE paramy_m
   USE paramz_m
   USE dericex6_m
   USE interpol6_m
   USE dericey6_m
   USE interpoly6_m
!
   implicit none
!
   real(8),dimension(nx,ny,nz) :: pp,ux,uy,uz
   real(8),dimension(nxm,nym,nz) :: ppm1,ppm,ppm2
   real(8),dimension(mx,my,mz) :: tr,us,ps,ps1
   real(8),dimension(mx,my) :: ps2d,ps2d1,ur1,ue1,w1,w2,w11,w22
   real(8),dimension(mx,my) :: us2d,ux2d,ur,ue,pp2d,ps2d2
   real(8),dimension(mx) :: exs,xkx,xk2
   real(8),dimension(my) :: eys,yky,yk2
   real(8),dimension(mz) :: zkz,zk2 
   real(8),dimension(ny/2,m1+m2+1) :: a,a1,a2
   real(8),dimension(ny/2,m1) :: al,al1
   real(8),dimension(ny/2) :: indx,indx1,e,c
   real(8),dimension(nx) :: d5
   real(8),dimension(ny) :: d1,d
   real(8),dimension(2) :: ja,jb
   real(8),dimension(ny+1) :: yp,ypi
   real(8),dimension(nx) :: transy
   real(8),dimension(ny) :: transx
   integer :: mxyz,mxy,i,j,k,nx,mx,nxm,ny,my,nym,nz,mz,mxxyyz,nwork,m1,m2
   integer :: ntrigsX,ntrigsY,ntoto,ip,npaire
   real(8) :: pi,twopi,xxk1,x,y,xx3,xx4,w,wp,ytt,xtt,yt,xt,yt1,xt1,xyzk
   real(8) :: xa,xb,xx1,xx2,xa1,xb1,xl2,xxa,alpha,beta
   real(8) :: xmo1,xd,xd1,xtoto
   real(8) :: xvar,xvar1,beta1,tw,xmo,xps 
!
!*****************ATTENTION*******************
   npaire=1
!*********************************************
   mxyz=mx*my*mz
   mxy=mx*my
   pi=acos(-1.)
   twopi=2.*pi
   xxk1=twopi/xlx
   xxa=200.
   mxxyyz=nxm*nym*nz
   ppm1(:,:,:)=0. ; ps(:,:,:)=0. ; ps1(:,:,:)=0. ; us(:,:,:)=0.
   us2d(:,:)=0. ; ps2d(:,:)=0.

!     
!   do k=1,nz
!   do j=1,nym
!      y=(j-0.5)*dy
!      do i=1,nxm
!         x=(i-0.5)*dx
!         ppm(i,j,k)=cos(6.*pi*x)!*cos(13.*pi*y)!
!
!   enddo
!   enddo
!   enddo
      
      do k=1,nz
      do j=1,nym
      do i=1,nxm/2
        ps(i,j,k)=ppm(2*(i-1)+1,j,k)
      enddo
      enddo
      enddo

      do k=1,nz
      do j=1,nym
      do i=nxm/2+1,nxm
         ps(i,j,k)=ppm(2*nxm-2*i+2,j,k)
      enddo
      enddo
      enddo
!
      do k=1,nz
      do i=1,nxm
      do j=1,nym/2
         ps1(i,j,k)=ps(i,2*(j-1)+1,k)
      enddo
      enddo
      enddo

      do k=1,nz
      do i=1,nxm
      do j=nym/2+1,nym
         ps1(i,j,k)=ps(i,2*nym-2*j+2,k)
      enddo
      enddo
      enddo
!***********TRAVAIL EN 2D********************
      do i=1,mx
      do j=1,my
         ps2d(i,j)=ps1(i,j,1)
      enddo
      enddo     
!
 print *,'toto'
!
   call SLFFT2D(ps2d,nxm,nym,mx,my,nxm,nym,nwork,&
        ntrigsX,ntrigsY,0,0,-1)

!      do j=1,ny
!      do i=1,nx
!         if (abs(ps2d(i,j)).gt.0.005) print '(i4,i4,f15.8)',i,j,ps2d(i,j)
!      enddo
!      enddo
!      stop
   ur(:,:)=0. ; ue(:,:)=0. ; ur1(:,:)=0. ; ue1(:,:)=0.

   do j=1,nym
   do i=1,nxm/2
      ur(i,j)=ps2d(2*i-1,j)
   enddo
   enddo
      
      do j=1,nym
      ur(nxm/2+1,j)=ps2d(nxm+1,j)
      enddo
      do j=1,nym
      do i=2,nxm/2
         ur(i+nxm/2,j)=ur(nxm/2-i+2,j)
      enddo
      enddo
      do j=1,nym
      do i=1,nxm/2
         ue(i,j)=ps2d(2*i,j)
      enddo
      enddo
      do j=1,nym
      ue(nxm/2+1,j)=ps2d(nx+1,j)
      enddo
      do j=1,nym
      do i=2,nxm/2
         ue(i+nxm/2,j)=-ue(nxm/2-i+2,j)
      enddo
      enddo

      do j=1,ny
      do i=1,nxm/2+1
         xa=cos((j-1)/2.*pi/(nym))
         xb=sin((j-1)/2.*pi/(nym))
         xx1=ur(i,j)*xa/2.+ue(i,j)*xb/2.
         xx2=ur(i,nym-j+2)*xa/2.-ue(i,nym-j+2)*xb/2.           
         ur1(i,j)=xx1+xx2      
      enddo
      enddo
      do j=1,ny
      do i=2,nxm/2
         ur1(i+nxm/2,j)=ur1(nxm/2-i+2,j)
      enddo
      enddo
!----------------------------------------------------
      do j=1,ny
      do i=1,nxm/2+1
         xa=cos((j-1)/2.*pi/(nym))
         xb=sin((j-1)/2.*pi/(nym))
         xa1=cos((nym-j+1)/2.*pi/(nym))
         xb1=sin((nym-j+1)/2.*pi/(nym))
         xx1=-ur(i,j)*xb/2.+ue(i,j)*xa/2.
         xx2=ur(i,nym-j+2)*xb/2.+ue(i,nym-j+2)*xa/2.           
         ue1(i,j)=xx1+xx2 
      enddo
      enddo
      do j=1,ny
      do i=2,nxm/2
         ue1(i+nxm/2,j)=-ue1(nxm/2-i+2,j)
      enddo
      enddo

      do i=1,nxm
         ur1(i,1)=ur1(i,1)*2.
         ue1(i,1)=ue1(i,1)*2.
      enddo
   
      do j=1,ny
      do i=1,nxm
         xb=sin((i-1)/2.*pi/(nxm))
         xa=cos((i-1)/2.*pi/(nxm))
         xx1=(ur1(i,j)+ur1(nxm-i+2,j))*xa/2
         xx2=(ue1(i,j)-ue1(nxm-i+2,j))*xb/2.
         ps2d(i,j)=xx1+xx2
!         print *,i,j,xx1,xx2
      enddo
!      pause
      enddo
      do j=1,ny
         ps2d(1,j)=2.*ps2d(1,j)
!!!!!!!!!!!!!!!!         ps2d(nxm+1,j)=0.
      enddo


      ps2d1=0. ; ps2d1=ps2d
    
!   do i=1,mx
!   do j=1,my
!      print '(i4,i4,f15.8)',i,j,ps2d(i,j)
!   enddo
!   pause
!   enddo      
!********************INTEGRATION******************************
           do i=1,mx
              exs(i)=0.
           enddo
!
           do j=1,my
              eys(j)=0.
           enddo
           do i=1,mx
              w=twopi*0.5*(i-1)/(mx-2)
              wp=acix6*2.*dx*sin(w/2.)+(bcix6*2.*dx)*sin(3./2.*w)
              wp=wp/(1.+2.*alcaix6*cos(w))
!!!               wp=w
              xkx(i)=(mx-2)*wp/xlx
              exs(i)=(mx-2)*w/xlx
              xk2(i)=xkx(i)*xkx(i)
           enddo 
           xkx(mx-1)=0.
           xk2(mx-1)=0.
           xkx(mx  )=0.
           xk2(mx  )=0.
!
           do j=1,my
              w=twopi*0.5*(j-1)/(my-2)
              wp=aciy6*2.*dy*sin(w/2.)+(bciy6*2.*dy)*sin(3./2.*w)
              wp=wp/(1.+2.*alcaiy6*cos(w))
!!!              wp=w
              yky(j)=(my-2)*wp
              eys(j)=(my-2)*w/yly
              yk2(j)=yky(j)*yky(j)
           enddo
           yky(my-1)=0.
           yk2(my-1)=0.
           xkx(my  )=0.
           yk2(my  )=0.
!******************la fonction de transfert*********************
!!           dy=dy*20.
           do j=1,my
              xtt=(biciy6*2.*cos(eys(j)*3.*dy/2.))
              xt=(aiciy6*2.*cos(eys(j)*dy/2.))
              xt1=(1.+2.*ailcaiy6*cos(eys(j)*dy))
              transx(j)=((xtt+xt)/xt1)
           enddo
           
           do i=1,mx
              ytt=(bicix6*2.*cos(exs(i)*3.*dx/2.))
              yt=(aicix6*2.*cos(exs(i)*dx/2.))
              yt1=(1.+2.*ailcaix6*cos(exs(i)*dx))
              transy(i)=((ytt+yt)/yt1)
           enddo
!           do i=1,nx
!              print *,transx(i),transy(i)
!           enddo
!           pause
!*******************STRET STRET STRET STRET********************

!********************nombre d onde dans d**********************       
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
!         print *,'toto',d5(i),transy(i)
!      enddo
!      pause
!*****************DEBUT DE LA BOUCLE EN XXXXXXXXXXXX******************
       do i=1,nxm
!*********************************************************************
!*****************petit tablo d onde dans b******************** 
      do j=1,ny/2
         d(j)=yky(2*j-1)*transy(i)
         d1(j)=yky(2*j)*transy(i)
      enddo 

      
      ps2d1=0. ; ps2d1=ps2d

         xps=ps2d(i,ny)
         e=0. ; c=0.
         do j=1,ny/2
            e(j)=ps2d(i,2*j-1)
            c(j)=ps2d(i,2*j)
!            print *,e(j),c(j)
         enddo
!         pause
!***************Construction de la matrice a**************************
      xvar=alpha+1./2./beta
      xmo=alpha+3./4./beta
      tw=twopi*twopi
      beta1=beta*beta
      xvar1=xvar*xvar
      a(1,m1+1)=-(d5(i)*d5(i)*transx(1)*transx(1))-d(1)*d(1)*xvar1/tw*4.&
           -(1./4./tw/beta1)*(d(1)*d(2))
      do j=2,ny/2-1
         a(j,m1+1)=-(d5(i)*d5(i)*transx(2*j-1)*transx(2*j-1))-d(j)*d(j)*xvar1/tw*4.&
              -(1./4./tw/beta1)*(d(j)*d(j-1)+d(j)*d(j+1))
      enddo
      a(ny/2,m1+1)=-(d5(i)*d5(i)*transx(ny-1)*transx(ny-1))-d(ny/2)*d(ny/2)*xvar1/tw*4.&
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
      a2(1,m1+1)=-(d5(i)*d5(i)*transx(2)*transx(2))-d1(1)*d1(1)*xmo*xmo1/tw*4.&
           -(1./4./tw/beta1)*(d1(1)*d1(2))
      do j=2,ny/2-1
         a2(j,m1+1)=-(d5(i)*d5(i)*transx(2*j)*transx(2*j))-d1(j)*d1(j)*xvar1/tw*4.&
              -(1./4./tw/beta1)*(d1(j)*d1(j-1)+d1(j)*d1(j+1))
      enddo
      a2(ny/2,m1+1)=-(d5(i)*d5(i)*transx(ny)*transx(ny))-d1(ny/2)*d1(ny/2)*xmo*xmo/tw*4.&
           -(1./4./tw/beta1)*(d1(ny/2)*d1(ny/2-1))
!***************element du cote des m1********************************
      a2(1,m1)=0.
      a2(2,m1)=xvar/tw/beta*(d1(2)*d1(1))&
           +xmo/tw/beta*(d1(1)*d1(1))
      do j=3,ny/2-1
         a2(j,m1)=xvar/tw/beta*(d1(j)*d1(j-1)+d1(j-1)*d1(j-1))
      enddo
      a2(ny/2,m1)=xmo/tw/beta*(d1(ny/2)*d1(ny/2-1))&
           +xvar/tw/beta*(d1(ny/2-1)*d1(ny/2-1))
      a2(1,m1-1)=0.
      a2(2,m1-1)=0.
      do j=3,ny/2
         a2(j,m1-1)=-d1(j-1)*d1(j-2)/4./tw/beta1
      enddo
!****************element du cote des m2*********************************
      a2(1,m1+m2)=xmo1/tw/beta*(d1(1)*d1(2))&
           +xvar/tw/beta*(d1(2)*d1(2))
      do j=2,ny/2-2
         a2(j,m1+m2)=xvar/tw/beta*(d1(j)*d1(j+1)+d1(j+1)*d1(j+1))
      enddo
      a2(ny/2-1,m1+m2)=xvar/tw/beta*(d1(ny/2-1)*d1(ny/2))&
           +xmo/tw/beta*(d1(ny/2)*d1(ny/2))
      a2(ny/2,m1+m2)=0.
      do j=1,ny/2-2
         a2(j,m1+m2+1)=-d1(j+1)*d1(j+2)/4./tw/beta1
      enddo
      a2(ny/2-1,m1+m2+1)=0.
      a2(ny/2,m1+m2+1)=0.
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
!       do j=1,my
!          print *,j,ps2d1(i,j)
!       enddo
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
!********ATTENTION********************
!*************************************       
     us2d(i,ny)=0.!xtoto
!       print *,xps
!       pause
!       if (i==1) us2d(1,ny)=-xps
!*******FIN BOUCLE EN XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXxxxxx*************
       enddo
!***********************************************************************bon
!**************************************************************
!   ur(:,:)=0. ; ue(:,:)=0. ; ur1(:,:)=0. ; ue1(:,:)=0.

!       do i=3,3
!       do j=1,my
!          print *,i,j,us2d(i,j)
!       enddo
!!       pause
!       enddo

!       us2d=0. ; us2d=ps2d1

      do j=1,ny
      do i=1,nxm/2+1
         w1(i,j)=us2d(i,j)
         w2(i,j)=0.
      enddo
      enddo

      do j=1,ny
      do i=1,nxm/2+1
         w11(i,j)=us2d(nxm-i+2,j)
         w22(i,j)=0.
      enddo
      enddo
!
      do j=1,my
      do i=1,mx
         ur1(i,j)=0.
         ue1(i,j)=0.
      enddo
      enddo


      do j=1,ny
      do i=1,nxm/2+1
         xb=sin((i-1)/2.*pi/(nxm))
         xa=cos((i-1)/2.*pi/(nxm))
         xb1=sin((nxm-i+1)/2.*pi/(nxm))
         xa1=cos((nxm-i+1)/2.*pi/(nxm))
         ur1(i,j)=w1(i,j)*xa+w11(i,j)*xa1
         ue1(i,j)=w1(i,j)*xb-w11(i,j)*xb1
      enddo
!      pause
      enddo
!
      do j=1,ny
      do i=2,nxm/2
         ur1(i+nxm/2,j)=ur1(nxm/2-i+2,j)
      enddo
      enddo
      do j=1,ny
      do i=2,nxm/2
         ue1(i+nxm/2,j)=-ue1(nxm/2-i+2,j)
      enddo
      enddo

!      do i=1,nx
!      do j=1,ny
!         print *,'W8',ur1(i,j),ue1(i,j)
!      enddo
!      pause
!      enddo
!
      do j=1,ny
      do i=1,nxm/2+1
         xa=cos((j-1)/2.*pi/(nym))
         xb=sin((j-1)/2.*pi/(nym))
         xa1=cos((nym-j+1)/2.*pi/(nym))
         xb1=sin((nym-j+1)/2.*pi/(nym))
         xx1=ur1(i,j)*xa-ue1(i,j)*xa1
         xx2=ur1(i,nym-j+2)*xb+ue1(i,nym-j+2)*xa 
         xx3=-ur1(i,nym-j+2)*xb1+ue1(i,nym-j+2)*xb
         xx4=ur1(i,j)*xb+ue1(i,j)*xa
         ps2d(2*i-1,j)=xx1+xx2
         ps2d(2*i,j)=xx3+xx4
!         print *,'XX',xx1,xx2
      enddo
!      pause
      enddo
      do i=1,mx
        ps2d(i,ny)=0.
      enddo
      do i=1,mx
         ps2d(i,ny)=ps2d(i,(ny+1)/2)
      enddo

!**********************************************************
!      print *,nxm,nym,mx,my,nxm,nym,nwork,ntrigsX,ntrigsY
!      pause
   call SLFFT2D(ps2d,nxm,nym,mx,my,nxm,nym,nwork,&
        ntrigsX,ntrigsY,0,0,1)

!
!      do j=1,ny
!         do i=1,nx
!            print *,i,j,ps2d(i,j)
!         enddo
!         pause
!      enddo 
      do j=1,ny
      do i=1,nxm/2
         pp2d(2*i-1,j)=ps2d(i,j)
      enddo
      enddo
      do j=1,ny
      do i=1,nxm/2
         pp2d(2*i,j)=ps2d(nxm-i+1,j) 
      enddo
      enddo
      
      do j=1,nym/2   
      do i=1,nxm
         ux2d(i,2*j-1)=pp2d(i,j)
      enddo
      enddo
      
      do j=1,nym/2
      do i=1,nxm
         ux2d(i,2*j)=pp2d(i,nym-j+1)
      enddo
      enddo 

      do k=1,nz
      do i=1,nxm
      do j=1,nym
         ppm(i,j,k)=ux2d(i,j)
      enddo
      enddo
      enddo
!      open (10,file='courbe.dat',form='formatted',status='unknown')
!      do j=1,nym
!         y=(j-0.5)*dy
!         write(10,*) y,ppm(10,j,1),ppm(100,j,1)
!      enddo
!      close(10)
!      pause
!      open (10,file='courbe.dat',form='formatted',status='unknown')
!      do i=1,nxm
!         x=(i-0.5)*dx
!         write(10,*) x,ppm(i,10,1),ppm(i,50,1)
!      enddo
!      close(10)
!      pause
!
!   xl2=0.
!   do j=1,nym
!   do i=1,nxm 
!     xl2=xl2+(ppm(i,j,1)-ppm1(i,j,1))*(ppm(i,j,1)-ppm1(i,j,1))
!   enddo
!   enddo
!   xl2=xl2/nxm/nym
!   write(*,*)'erreur finale ppm',sqrt(xl2)
!   stop
!      open (10,file='courbe.dat',form='formatted',status='unknown')
!      do j=1,nym
!         y=(j-0.5)*dy
!         write(10,*) y,ps1(15,j,1),ppm(5,j,1)
!      enddo
!      close(10)
!      pause
!   do j=1,nym
!   do i=1,nxm
!      print *,i,j,ppm(i,j,1)
!   enddo
!   pause
!   enddo  

  
   return
end subroutine temporellebis16STR

