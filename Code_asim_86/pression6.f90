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
!   pi=acos(-1.)
!   xxk1=pi/xlx
!   xxk2=pi/yly
!   do k=1,nz
!   do j=1,nym
!     y=ypi(j)!!INTERPOLE
!     y=(j-0.5)*dy
!      do i=1,nxm
!         x=(i-0.5)*dx
!         tx2(i,j,k)=5.*xxk1*cos(5.*xxk1*x)*cos(13.*xxk2*y)
!         ty2(i,j,k)=13.*xxk2*cos(13.*xxk2*y)*cos(15*xxk1*x)
!         ppm1(i,j,k)=tx2(i,j,k)+ty2(i,j,k)
!   enddo
!   enddo
!   enddo
!!
!   do k=1,nz
!   do j=1,ny
!      y=yp(j) !!!!!COLOC
!      y=(j-1)*dy
!      do i=1,nx
!         x=(i-1)*dx
!         ux(i,j,k)=sin(5.*xxk1*x)*cos(13.*xxk2*y)
!        uy(i,j,k)=sin(13.*xxk2*y)*cos(15.*xxk1*x)
!     enddo
!   enddo
!   enddo 
!
   tmax=0
   tmoy=0.
   imax=0
   jmax=0
   kmax=0
   do i=1,nvec
      t(i)=0.
      r(i)=0.
   enddo
!

   call intery6(ty1,ux,cifyp6,cisyp6,ciwyp6,nx,ny,nym,nz,1)
!   
   call decx6(tx,ty1,cfx6,csx6,cwx6,nx,nxm,ny,nym,nz,0)  
!     
   call inter6(tx1,uy,cifxp6,cisxp6,ciwxp6,nx,nxm,ny,nz,1)
!      
   call decy6(ty,tx1,cfy6,csy6,cwy6,nx,nxm,ny,nym,nz,0)
!

!  xl2=0.
!   do j=1,nym
!   do i=1,nxm
!      xl2=xl2+(ty(i,j,1)-ty2(i,j,1))*(ty(i,j,1)-ty2(i,j,1))
!   enddo
!   enddo
!   xl2=xl2/nxm/nym
!   write(*,*)'erreur inter ppm',sqrt(xl2)     
!!!   stop
!      open (10,file='courbe.dat',form='formatted',status='unknown')
!      do j=1,nym
!         write(10,*) ypi(j),ty(7,j,1),ty2(7,j,1)
!      enddo
!      close(10)
!!      pause
!   do i=1,nxm
!   do j=1,nym
!     print '(i4,f15.8,f15.8,f15.8)',i,ty(i,j,1),ty2(i,j,1)
!   enddo
!   pause
!   enddo

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
subroutine divU (ux,uy,&
             cfx6,csx6,cwx6,&
             cfxp6,csxp6,cwxp6,&
             cifx6,cisx6,ciwx6,cifxp6,cisxp6,ciwxp6,&
             cfy6,csy6,cwy6,cfyp6,csyp6,cwyp6,&
             cify6,cisy6,ciwy6,cifyp6,cisyp6,ciwyp6,&
             nx,ny,nz,nxm,nym)    
! 
!********************************************************************
!
   USE paramx_m
   USE paramy_m
   USE paramz_m
!
   implicit none
!
   real(8),dimension(nx,ny,nz) :: ux,uy
   real(8),dimension(nxm,ny,nz) :: tx1
   real(8),dimension(nx,nym,nz) :: ty1
   real(8),dimension(nxm,nym,nz) :: ppm,ty,tx
   real(8),dimension(nxm) :: cfx6,ccx6,cwx6,cifx6,cisx6,ciwx6,csx6
   real(8),dimension(nxm) :: cfxp6,csxp6,cwxp6,cifxp6,cisxp6,ciwxp6 
   real(8),dimension(nym) :: cfy6,ccy6,cwy6,cify6,cisy6,ciwy6,csy6
   real(8),dimension(nym) :: cfyp6,csyp6,cwyp6,cifyp6,cisyp6,ciwyp6 
   real(8) ::tmax,tmoy
   integer :: nx,ny,nz,nxm,nym,i,j,k,imax,jmax,kmax
!
   tmax=0.
   tmoy=0.
   imax=0
   jmax=0
   kmax=0
!
   call intery6(ty1,ux,cifyp6,cisyp6,ciwyp6,nx,ny,nym,nz,1)
   call decx6(tx,ty1,cfx6,csx6,cwx6,nx,nxm,ny,nym,nz,0) 
   call inter6(tx1,uy,cifxp6,cisxp6,ciwxp6,nx,nxm,ny,nz,1)
   call decy6(ty,tx1,cfy6,csy6,cwy6,nx,nxm,ny,nym,nz,0)
!
   do j=1,nym
   do i=1,nxm
      ppm(i,j,1)=tx(i,j,1)+ty(i,j,1)
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
end subroutine divU
!
!********************************************************************
!
subroutine gradpbis16(ps2d,u8x,u8y,uxn,uyn,&
          cfx6,csx6,cwx6,&
          cfxp6,csxp6,cwxp6,&
          cifx6,cisx6,ciwx6,&
          cifxp6,cisxp6,ciwxp6,&
          cfi6,cci6,cbi6,cfip6,ccip6,cbip6,&
          csip6,cwip6,csi6,cwi6,&
          cifi6,cici6,cibi6,cifip6,cicip6,&
          cibip6,cisip6,ciwip6,cisi6,ciwi6,cify6,cisy6,ciwy6,&
          cifyp6,cisyp6,ciwyp6,cfyp6,csyp6,cwyp6,&
          cfi6y,cci6y,cbi6y,cfip6y,ccip6y,cbip6y,&
          csip6y,cwip6y,csi6y,cwi6y,&
          cifi6y,cici6y,cibi6y,cifip6y,cicip6y,&
          cibip6y,cisip6y,ciwip6y,cisi6y,ciwi6y,cfy6,csy6,cwy6,&
          nx,nxm,mx,ny,nym,my,nz)
!
!********************************************************************
!
   USE paramx_m
   USE paramy_m
!
   implicit none
!
   real(8),dimension(nx,ny,nz) :: uxn,uyn,u8x,u8y,px,py
   real(8),dimension(nxm,nym,nz) :: ppm
   real(8),dimension(nxm,ny,nz) :: wk1
   real(8),dimension(nx,nym,nz) :: wk2
   real(8),dimension(mx,my) :: ps2d,ux2d,pp2d
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
   integer :: i,j,k,nx,mx,nxm,ny,my,nym,nz,mz
!
!  ******************* temporellebis16STR ********************
!  ************** p(mx,my) ---> p(nxm,nym,nz) ****************
   do j=1,ny
   do i=1,nxm/2
      pp2d(2*i-1,j)=ps2d(i,j)
   enddo
   enddo
!
   do j=1,ny
   do i=1,nxm/2
      pp2d(2*i,j)=ps2d(nxm-i+1,j) 
   enddo
   enddo
!      
   do j=1,nym/2   
   do i=1,nxm
      ux2d(i,2*j-1)=pp2d(i,j)
   enddo
   enddo
!     
   do j=1,nym/2
   do i=1,nxm
      ux2d(i,2*j)=pp2d(i,nym-j+1)
   enddo
   enddo 
!
   do k=1,nz
   do i=1,nxm
   do j=1,nym
      ppm(i,j,k)=ux2d(i,j)
   enddo
   enddo
   enddo
!
!  ********************* gradpbis16 **********************
   call interiy6(wk1,ppm,cifip6y,cisip6y,ciwip6y,nxm,nym,ny,nz,1)
   call deci6(px,wk1,cfip6,csip6,cwip6,nxm,nx,ny,nz,1)
   call interi6(wk2,ppm,cifip6,cisip6,ciwip6,nxm,nx,nym,nz,1)
   call deciy6(py,wk2,cfip6y,csip6y,cwip6y,nx,nym,ny,nz,1)
!
!  ********************** corgp ***********************
   do j=1,ny
   do i=1,nx
      uxn(i,j,1)=-px(i,j,1)+u8x(i,j,1)
      uyn(i,j,1)=-py(i,j,1)+u8y(i,j,1) 
   enddo
   enddo
!
   return
end subroutine gradpbis16
!
!****************************************************************
!
subroutine temporellebis16STR(ppm,ps2dre,nx,nxm,ny,nym,nz,mx,my,mz)
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
   real(8),dimension(nxm,nym,nz) :: ppm
   real(8),dimension(mx,my) :: ps2dre
   real(8),dimension(mx,my) :: ux2d,pp2d
   integer :: i,j,k,nx,mx,nxm,ny,my,nym,nz,mz
!
      do j=1,ny
      do i=1,nxm/2
         pp2d(2*i-1,j)=ps2dre(i,j)
      enddo
      enddo
      do j=1,ny
      do i=1,nxm/2
         pp2d(2*i,j)=ps2dre(nxm-i+1,j) 
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
!  
   return
end subroutine temporellebis16STR
!
!****************************************************************
!
subroutine construye_matriz_B(a_tot,al_tot,indx_tot,a2_tot,al1_tot,indx1_tot,d5,&
        alpha,beta,m1,m2,nx,nxm,ny,nym,nz,mx,my,mz)
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
   real(8),dimension(mx) :: exs,xkx,xk2
   real(8),dimension(my) :: eys,yky,yk2
   real(8),dimension(mz) :: zkz,zk2 
   real(8),dimension(ny/2,m1+m2+1) :: a,a1,a2
   real(8),dimension(ny/2,m1) :: al,al1
   real(8),dimension(ny/2) :: indx,indx1,e,c,e1,c1
   real(8),dimension(ny/2,m1+m2+1,nxm) :: a_tot,a2_tot
   real(8),dimension(ny/2,m1,nxm) :: al_tot,al1_tot
   real(8),dimension(ny/2,nxm) :: indx_tot,indx1_tot
   real(8),dimension(nx) :: d5
   real(8),dimension(ny) :: d1,d
   real(8),dimension(nx) :: transy
   real(8),dimension(ny) :: transx
   integer :: i,j,k,nx,mx,nxm,ny,my,nym,nz,mz,m1,m2
   real(8) :: pi,twopi,xxk1,x,y,w,wp,ytt,xtt,yt,xt,yt1,xt1
   real(8) :: alpha,beta
   real(8) :: xmo1,xd,xd1,xt2,yt2
   real(8) :: xvar,xvar1,beta1,tw,xmo 
!
   pi=acos(-1.)
   twopi=2.*pi
   xxk1=twopi/xlx
!
!  ****** Nombre d onde associe a la derivee premiere sur un maillage decale ******
   do i=1,mx
      exs(i)=0.
   enddo
!
   do j=1,my
      eys(j)=0.
   enddo
!
   do i=1,mx
      w=twopi*0.5*(i-1)/(mx-2)
      wp=acix6*2.*dx*sin(w/2.)+(bcix6*2.*dx)*sin(3./2.*w)
      wp=wp/(1.+2.*alcaix6*cos(w))
      xkx(i)=(mx-2)*wp/xlx
      exs(i)=(mx-2)*w/xlx
      xk2(i)=xkx(i)*xkx(i)
   enddo 
!         
   xkx(mx-1)=0.
   xk2(mx-1)=0.
   xkx(mx  )=0.
   xk2(mx  )=0.
!
   do j=1,my
      w=twopi*0.5*(j-1)/(my-2)
      wp=aciy6*2.*dy*sin(w/2.)+(bciy6*2.*dy)*sin(3./2.*w)
      wp=wp/(1.+2.*alcaiy6*cos(w))
      yky(j)=(my-2)*wp
      eys(j)=(my-2)*w/yly
      yk2(j)=yky(j)*yky(j)
   enddo
!
   yky(my-1)=0.
   yk2(my-1)=0.
   yky(my  )=0.
   yk2(my  )=0.

!  ****** Fonction de transfert associe au nombre d onde ******
   do j=1,ny
      xtt=(biciy6*2.*cos(eys(j)*3.*dy/2.))
      xt=(aiciy6*2.*cos(eys(j)*dy/2.))
      xt1=(1.+2.*ailcaiy6*cos(eys(j)*dy))
      transx(j)=(xtt+xt)/xt1  
   enddo
!           
   do i=1,nx
      ytt=(bicix6*2.*cos(exs(i)*3.*dx/2.))
      yt=(aicix6*2.*cos(exs(i)*dx/2.))
      yt1=(1.+2.*ailcaix6*cos(exs(i)*dx))
      transy(i)=(ytt+yt)/yt1
   enddo
!
!  ***************************************************************
!  ****************** Correccion por stretching ******************
!  ***************************************************************
!
!  nombre d onde dans D       
   if ((nclx==1).or.(nclx==2)) then
      do i=1,nx
         d5(i)=xkx(i)
      enddo
   endif
!
   do i=1,nxm ! Construye la matriz B=A^2
!
!     Petit tablo d onde dans b  
      do j=1,ny/2
         d(j)=yky(2*j-1)*transy(i)
         d1(j)=yky(2*j)*transy(i)
      enddo 
!
!     ****** Construye a: band diagonal matrix with m1 subdiagonal rows 
!     ****** and m2 superdiagonal rows, compactly stored
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
      a(ny/2,m1+1)=-(d5(i)*d5(i)*transx(ny-2)*transx(ny-2))-d(ny/2)*d(ny/2)*xvar1/tw*4.&
           -(1./4./tw/beta1)*(d(ny/2)*d(ny/2-1))
!
!     Element du cote des m1
      a(1,m1)=0.
      do j=2,ny/2
         a(j,m1)=xvar/tw/beta*(d(j)*d(j-1)+d(j-1)*d(j-1))
      enddo
      a(1,m1-1)=0.
      a(2,m1-1)=0.
      do j=3,ny/2
         a(j,m1-1)=-d(j-1)*d(j-2)/4./tw/beta1
      enddo
!
!     Element du cote des m2
      a(1,4)=xvar/tw/beta*(d(1)*d(2)+d(2)*d(2))
      do j=2,ny/2-1
         a(j,m1+m2)=xvar/tw/beta*(d(j)*d(j+1)+d(j+1)*d(j+1))
      enddo
      do j=1,ny/2-2
         a(j,m1+m2+1)=-d(j+1)*d(j+2)/4./tw/beta1
      enddo
      a(1,m1+m2+1)=a(1,m1+m2+1)*2.
      a(ny/2-1,m1+m2+1)=0.
      a(ny/2,m1+m2+1)=0.
!
      if (d5(i).eq.0.) then 
         a(1,m1+1)=1.
         a(1,m1+m2)=0.
         a(1,m1+m2+1)=0.
      endif
!
!     ****** Construye a2: band diagonal matrix with m1 subdiagonal rows 
!     ****** and m2 superdiagonal rows, compactly stored
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
      a2(ny/2,m1+1)=-(d5(i)*d5(i)*transx(ny-1)*transx(ny-1))-d1(ny/2)*d1(ny/2)*xmo*xmo/tw*4.&
           -(1./4./tw/beta1)*(d1(ny/2)*d1(ny/2-1))
!
!     Element du cote des m1
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
!
!     Element du cote des m2
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

!     Constructs an LU decomposition of a rowwise permutation of a. 
!     The upper triangular matrix replaces a, while the lower triangular matrix is returned in al.
!     The vector indx records the row permutation effected by the partial pivoting.
      call bandec(a,ny,my,m1,m2,al,indx,xd)
      call bandec(a2,ny,my,m1,m2,al1,indx1,xd1)
!
!     Guarda resultados
      a_tot(:,:,i) = a
      al_tot(:,:,i) = al
      indx_tot(:,i) = indx
      a2_tot(:,:,i) = a2
      al1_tot(:,:,i) = al1
      indx1_tot(:,i) = indx1
!
   enddo ! Construye la matriz B=A^2
!
   return
end subroutine construye_matriz_B
!
!********************************************************************
!
subroutine divbric16bis (ppm,ux,uy,&
     cfx6,csx6,cwx6,&
     cfxp6,csxp6,cwxp6,&
     nx,nxm,mx,ny,nym,my,nz,&
     cifx6,cisx6,ciwx6,cifxp6,cisxp6,ciwxp6,&
     cfy6,csy6,cwy6,cfyp6,csyp6,cwyp6,&
     cify6,cisy6,ciwy6,cifyp6,cisyp6,ciwyp6)    
! 
!********************************************************************
!
   USE paramx_m
   USE paramy_m
   USE paramz_m
   USE param1_m
   USE barreau_m
!
   implicit none
!
   integer :: nx,ny,nz,i,j,k,nxm,mx,nym,my
   real(8),dimension(nx,ny,nz) :: ux,uy
   real(8),dimension(nxm,nym,nz) :: tx,ty
   real(8),dimension(nx,nym,nz) :: ty1
   real(8),dimension(nxm,ny,nz) :: tx1
   real(8),dimension(nxm,nym,nz) :: ppm
   real(8),dimension(nxm) :: cfx6,csx6,cwx6,cifx6,cisx6,ciwx6 
   real(8),dimension(nxm) :: cifxp6,cisxp6,ciwxp6,cfxp6,csxp6,cwxp6
   real(8),dimension(nym) :: cfy6,csy6,cwy6,cify6,cisy6,ciwy6 
   real(8),dimension(nym) :: cifyp6,cisyp6,ciwyp6,cfyp6,csyp6,cwyp6
!
   call intery6(ty1,ux,cifyp6,cisyp6,ciwyp6,nx,ny,nym,nz,1)
!   
   call decx6(tx,ty1,cfx6,csx6,cwx6,nx,nxm,ny,nym,nz,0)  
!     
   call inter6(tx1,uy,cifxp6,cisxp6,ciwxp6,nx,nxm,ny,nz,1)
!      
   call decy6(ty,tx1,cfy6,csy6,cwy6,nx,nxm,ny,nym,nz,0)
!     
   do j=1,nym
   do i=1,nx
      ppm(i,j,1)=tx(i,j,1)+ty(i,j,1)
   enddo
   enddo
!
   return
end subroutine divbric16bis

