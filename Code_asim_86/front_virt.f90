!********************************************************************
!
subroutine virtuel(ux,uy,uz,nx,ny,nz,yp) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
  USE paramx_m 
  USE paramy_m 
  USE paramz_m 
  USE barreau_m 
!...Translated by PSUITE Trans90                  4.3ZH 16:58:09   1/27/ 4  
!...Switches: -yv INDDO=0 -nbejkno
  implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  integer , intent(in) :: nx 
  integer , intent(in) :: ny 
  integer , intent(in) :: nz 
  real(8) , intent(in) :: ux(nx,ny,nz) 
  real(8) , intent(in) :: uy(nx,ny,nz) 
  real(8) , intent(in) :: uz(nx,ny,nz) 
  real(8), dimension(ny+1) :: yp
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer :: np, k, j, i 
  real(8) :: dum, ux2, uy2, uz2, eps, uxmax, uymax, xm, ym, r 
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
  real(8) , external :: epmach 
!-----------------------------------------------
!
!********************************************************************
!
   dum=0. 
   ux2=0. 
   uy2=0. 
   uz2=0. 
   np=0 
   eps=epmach(dum) 
   uxmax=-1.E10 
   uymax=-1.E10 
!
   if (nz>1) then 
      do k=1,nz 
      do j=1,ny 
      do i=1,nx 
         xm=(i-1)*dx 
         ym=yp(j) 
         r=sqrt((xm-cex)*(xm-cex)+(ym-cey)*(ym-cey)) 
         if (r-ra >= eps) cycle  
         ux2=ux2+ux(i,j,k)*ux(i,j,k) 
         uy2=uy2+uy(i,j,k)*uy(i,j,k) 
         uz2=uz2+uz(i,j,k)*uz(i,j,k) 
         np=np+1 
!                     print *,'ijk,',i,j,k
      enddo
      enddo 
      enddo 
      ux2=ux2/np 
      uy2=uy2/np 
      uz2=uz2/np 
      ux2=sqrt(ux2) 
      uy2=sqrt(uy2) 
      uz2=sqrt(uz2) 
      print *,'Vitesses dans l''obstacle : ',ux2,uy2,uz2,np 
!         stop
   else 
      do j=1,ny 
      do i=1,nx 
         xm=(i-1)*dx 
         ym=yp(j) 
         r=sqrt((xm-cex)*(xm-cex)+(ym-cey)*(ym-cey)) 
         if (r-ra>=eps) cycle  
         ux2=ux2+ux(i,j,1)*ux(i,j,1) 
         uy2=uy2+uy(i,j,1)*uy(i,j,1) 
         if (abs(ux(i,j,1))>uxmax) uxmax=abs(ux(i,j,1))
         if (abs(uy(i,j,1))>uymax) uymax=abs(uy(i,j,1)) 
         np=np+1 
!                  print *,'ij',r
      enddo
!      pause
      enddo 
      ux2=ux2/np 
      uy2=uy2/np 
      ux2=sqrt(ux2) 
      uy2=sqrt(uy2) 
      print *,'vorms',ux2,uy2,np
      print *,'ux uy max ',uxmax,uymax 
   endif
! 
   return  
end subroutine virtuel
!
!********************************************************************
!
subroutine forcage(ux,uy,uz,nx,ny,nz,esp,qx,qy,qz,yp,&
           cexx,ceyy,ceyy2,ceyy3,ceyy21,ceyy31,&
           theta,h,iimin,i5,i55,jj,j1,j2,j3,j4,j5) 
!********************************************************************
!   M o d u l e s 
!-----------------------------------------------
  USE paramx_m 
  USE paramy_m 
  USE paramz_m 
  USE param1_m 
  USE barreau_m 
  USE controle_m !utile pour jets
  USE delta_m  !utile pour jets
  USE ecoulement_m !utile pour jets
  
  implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  integer :: nx,ny,nz 
  real(8),dimension(nx,ny) :: theta
  real(8),dimension(nx,ny,nz) :: ux,uy,uz,esp,qx,qy,qz
  real(8),dimension(ny+1) :: yp,y2a
  real(8),dimension(nx) :: cexx,ceyy,ceyy2,ceyy3,ceyy21,ceyy31
  integer :: ngxsym, ngysym, ndxsym, ndysym, j, i, k,np,jmax,jjmax,iimin,iimax
  integer :: ii,i21,i31,i5,i6,i7,i55
  integer :: j1,j2,j3,j4,j5,jj,jj1
  integer :: ixa1,ixb1,ncexx
  real(8) , dimension(1) :: uxs, uys, uzs, e 
  real(8) :: dum,pi,twopi,zero,xm,ym,zm,rsym,xsym,ysym,r,r2,t,u,dyb,r1,x11,xa1,xb1
  real(8) :: yl1,zl1,h,yy,angle3,distx,distx2,distx3,distx21,distx31 !variable ajoute pour jet de controle
  real(8) :: val                                    !variable poubelle a supprimer 
  real(8) , external :: epmach 
!********************************************************************
!
   dum=0. 
   pi=acos(-1.) 
   twopi=2*pi 
   zero=epmach(dum)*100. 
   r=0.
   r2=0.
   t=0. 
   u=0.
!
   do j=j5,ny
   do i=1,i5
      ym=yp(j)
      r=ym-ceyy(i)
      if (r.lt.zero) then
      ux(i,j,1)=0.
      uy(i,j,1)=0.
      esp(i,j,1)=1.
      endif 
   enddo
   enddo
!
   do j=j5,ny
   do i=1,i5
      ym=yp(j)
      r=ym-ceyy(i)
      r2=ym-ceyy21(i) 
      if ((r.gt.zero).and.(r2.lt.zero)) then
      ux(i,j,1)=0.
      uy(i,j,1)=0.
      esp(i,j,1)=1.
      endif 
   enddo
   enddo
!
   do j=j5,ny
   do i=1,i5
      ym=yp(j)
      r=ym-ceyy3(i)
      r2=ym-ceyy31(i) 
      if ((r.gt.zero).and.(r2.le.zero)) then
      ux(i,j,1)=0.
      uy(i,j,1)=0.
      esp(i,j,1)=1.
      endif 
   enddo
   enddo
!
   do j=j5,j4
   do i=1,i55
      ym=yp(j)
      r=ym-ceyy3(i)
      r2=(ym-ceyy31(iimin))
      if ((r.gt.zero).and.(r2.le.zero)) then
      ux(i,j,1)=0.
      uy(i,j,1)=0.
      esp(i,j,1)=1.
      endif
    enddo
    enddo
!
!1002 format(45F12.8)
!   open(10,file='forgnu.dat',form='formatted')
!   do j=1,ny
!   do i=1,nx
!    write (10,1002) yp(j),ux(1,j,1)-u1,((1.)*dx)+ux(1,j,1)-u1,((2.)*dx)+ux(2,j,1)-u1,(3.*dx)+ux(3,j,1)-u1,(4*dx)+ux(4,j,1)-u1,&
!                   ((5.)*dx)+ux(5,j,1)-u1,((6.)*dx)+ux(6,j,1)-u1,((7.)*dx)+ux(7,j,1)-u1,((8.)*dx)+ux(8,j,1)-u1,&
!                   ((9.)*dx)+ux(9,j,1)-u1,((10.)*dx)+ux(10,j,1)-u1,((11.)*dx)+ux(11,j,1)-u1,((12.)*dx)+ux(12,j,1)-u1, &
!                   ((12.)*dx)+ux(12,j,1)-u1,((13.)*dx)+ux(13,j,1)-u1,((14.)*dx)+ux(14,j,1)-u1,((15.)*dx)+ux(15,j,1)-u1, &
!                   ((16.)*dx)+ux(16,j,1)-u1,((17.)*dx)+ux(17,j,1)-u1,((18.)*dx)+ux(18,j,1)-u1,((19.)*dx)+ux(19,j,1)-u1, &
!                   ((20.)*dx)+ux(20,j,1)-u1,((21.)*dx)+ux(21,j,1)-u1,((22.)*dx)+ux(22,j,1)-u1,((23.)*dx)+ux(23,j,1)-u1, &
!                   ((24.)*dx)+ux(24,j,1)-u1,((25.)*dx)+ux(25,j,1)-u1,((26.)*dx)+ux(26,j,1)-u1,((27.)*dx)+ux(27,j,1)-u1, &
!                   ((28.)*dx)+ux(28,j,1)-u1,((29.)*dx)+ux(29,j,1)-u1,((30.)*dx)+ux(30,j,1)-u1,((31.)*dx)+ux(31,j,1)-u1, &
!                   ((32.)*dx)+ux(32,j,1)-u1,((33.)*dx)+ux(33,j,1)-u1,((34.)*dx)+ux(34,j,1)-u1,((35.)*dx)+ux(35,j,1)-u1, &
!                   ((36.)*dx)+ux(36,j,1)-u1,((37.)*dx)+ux(37,j,1)-u1,((38.)*dx)+ux(38,j,1)-u1,((39.)*dx)+ux(39,j,1)-u1, &
!                   ((40.)*dx)+ux(40,j,1)-u1,((41.)*dx)+ux(41,j,1)-u1,((42.)*dx)+ux(42,j,1)-u1,((43.)*dx)+ux(43,j,1)-u1
!    enddo
!    close(10)
!   write (10,*) 
!   enddo 
!   call system ('gnuplot doudou.gnu')
!   call system ('gnuplot plot.gp')
!
   return  
 end subroutine forcage
!*******************************************************************
!
 subroutine preparation_bruitvent (cexx,ceyy,ceyy2,ceyy21,ceyy3,ceyy31,&
            angle3,h,iimin,i5,i55,jj,j1,j2,j3,j4,j5,yp,nx,ny,nz)
!
!*******************************************************************
   USE paramx_m
   USE paramy_m
   USE paramz_m
   USE param1_m 
   USE barreau_m
   USE controle_m
   USE delta_m  
   USE ecoulement_m 
!
   implicit none
!
   integer :: nx,ny,nz 
   real(8),dimension(nx) :: cexx,ceyy,ceyy2,ceyy3,ceyy21,ceyy31
   real(8),dimension(ny+1) :: yp
   integer :: i,j,k,jmax,jjmax,iimin,iimax
   integer :: ii,i21,i31,i5,i6,i7,i55
   integer :: j1,j2,j3,j4,j5,jj,jj1
   real(8) :: dum,pi,twopi,zero,xm,ym,zm,rsym,xsym,ysym,r,r2,dyb,r1,x11,xa1,xb1
   real(8) :: yl1,zl1,h,yy,angle3,distx,distx2,distx3,distx21,distx31 
   real(8) :: val  
   real(8) , external :: epmach 
!  
   dum=0. 
   pi=acos(-1.) 
   twopi=2.*pi 
   zero=epmach(dum)*100. 
   r=0. 
   r2=0. 
! 
   h=rlat/5.
   angle3=-(0.5*pi-phi) 
   distx=sqrt((ra*tan(angle3))*(ra*tan(angle3))) !epaisseur inferieur externe du cylindre
   distx2=sqrt(((ra+h+rlat)*tan(angle3))*((ra+h+rlat)*tan(angle3))) !ligne definissant le centre du jet
   distx21=sqrt(((ra+h)*tan(angle3))*((ra+h)*tan(angle3))) !epaisseur inferieur interne du cylindre
   distx3=sqrt(((ra+h+2.*rlat)*tan(angle3))*((ra+h+2.*rlat)*tan(angle3))) !epaisseur superieur internet du cylindre 
   distx31=sqrt(((ra+2.*rlat+2.*h)*tan(angle3))*((ra+2.*rlat+2.*h)*tan(angle3))) !epaisseur superieur externe du cylindre
!   
   iimin=1
   i5=(distx/dx)+1
   i55=i5/4.
   i6=(distx2/dx)+1
   i7=(distx3/dx)+1
   i21=(distx21/dx)+1
   i31=(distx31/dx)+1
!
   cexx(iimin)=(iimin-1.)*dx   
   cexx(i5)=(i5-1.)*dx
   cexx(i6)=(i6-1.)*dx
   cexx(i7)=(i7-1.)*dx
   cexx(i21)=(i21-1.)*dx
   cexx(i31)=(i31-1.)*dx
!   
   do i=iimin,i5
   cexx(i)=(i-1.)*dx 
   if (tan(angle3).eq.0.) stop 'division par zero dans front_virt.f90'  
      ceyy(i)=(y1-ra)+(cexx(i5)-cexx(i))/abs(tan(angle3)) 
      ceyy2(i)=(y1-ra)+(cexx(i6)-cexx(i))/abs(tan(angle3)) 
   enddo
!
   do i=iimin,nx
      cexx(i)=(i-1.)*dx 
      ceyy21(i)=(y1-ra)+(cexx(i21)-cexx(i))/abs(tan(angle3)) 
      ceyy3(i)=(y1-ra)+(cexx(i7)-cexx(i))/abs(tan(angle3)) 
      ceyy31(i)=(y1-ra)+(cexx(i31)-cexx(i))/abs(tan(angle3)) 
   enddo
!
1000 format(6f12.8) 
   do i=iimin,nx
   write(10,1000) ((i-1)*dx),ceyy(i),ceyy2(i),ceyy21(i),ceyy3(i),ceyy31(i)
   enddo
!
   jj=int(cey/dy)+1
   j1=1
   do j=ny,jj,-1
      r=abs(yp(j)-cey)
      if (r.gt.(0.)) j1=j+1
   enddo
   j2=1
   do j=ny,jj,-1
      r=abs(yp(j)-cey)
      if (r.gt.(h)) j2=j+1
   enddo
   j3=1
   do j=ny,jj,-1
      r=abs(yp(j)-cey)
      if (r.gt.(h+2.*rlat)) j3=j+1
   enddo
   j4=1
   do j=ny,jj,-1
      r=abs(yp(j)-cey)
      if (r.gt.(2.*h+2.*rlat)) j4=j+1
   enddo
   j5=1
   do j=1,jj
      r=abs(yp(j)-cey)
      if (r.gt.(ra)) j5=j+1
   enddo
!
end subroutine preparation_bruitvent
 
