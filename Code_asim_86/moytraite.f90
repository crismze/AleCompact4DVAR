! ********************************************************
      program moytraite
! ********************************************************
  USE paramx_m 
  USE paramy_m 
  USE paramz_m 
  USE dericex_m 
  USE interpol4_m 
  USE dericex6_m 
  USE interpol6_m 
  USE dericey_m 
  USE interpoly4_m 
  USE dericey6_m 
  USE interpoly6_m 
  USE derivex_m 
  USE derivey_m 
  USE derivez_m 
  USE parfix_m 
  USE parfiy_m 
  USE parfiz_m 
  USE param1_m 
  USE param2_m 
  USE limit_m 
  USE ondelette_m 
  USE ecoulement_m 
  USE controle_m 
  USE sauvegarde_m 
  USE module_avs_m 
  USE delta_m 
  USE barreau_m 
!
   implicit none
!
  integer,parameter :: nx=601,ny=257,nz=1
  integer,parameter :: m1=2,m2=2;
  integer,parameter :: nxm=nx-1,nym=ny-1,nzm=nz
  integer,parameter :: mx=nx+1,my=ny+1,mz=nz
  integer,parameter :: nyp=257,nzp=1,mw=4*mx*my*mz,nvec=256
  integer,parameter :: mxf=2*nx,myf=2*ny,mzf=2*nz,mb=2*mx
  integer,parameter :: mxy=mx*my,myz=my*mz
  integer,parameter :: nwork=max((3*nym+3)*(3*nxm+3),2*ny*nx)
  integer,parameter :: ntrigsX=4*nx,ntrigsY=4*ny,nfax=19
! ****************************3D******************************
   real(8),dimension(nx,ny,nz) :: uxmt,uymt,uzmt,uxux,uyuy,uzuz,uxuy
   real(8),dimension(nx,ny,nz) :: ux,uy,uz,gx,gy,gz,tx,ty,tz,di,sc1
   real(8),dimension(ny,nz) :: sx,px,vx,wx,xx
   real(8),dimension(nx,nz) :: sy,py,vy,wy,xy
   real(8),dimension(nx,ny) :: sz,pz,vz,wz,xz
   real(8),dimension(nx) :: ffx,fcx,fbx,sfx,scx,sbx,fsx,fwx,ssx,swx,truc
   real(8),dimension(nx) :: ffxp,fsxp,fwxp,sfxp,ssxp,swxp
   real(8),dimension(ny) :: ffy,fcy,fby,sfy,scy,sby,fsy,fwy,ssy,swy
   real(8),dimension(ny) :: ffyp,fsyp,fwyp,sfyp,ssyp,swyp
   real(8),dimension(nz) :: ffz,fcz,fbz,sfz,scz,sbz,fsz,fwz,ssz,swz
   real(8),dimension(nz) :: ffzp,fszp,fwzp,sfzp,sszp,swzp
   real(8),dimension(nx) :: uaxe_y,uaxe_z
   real(8),dimension(nxm) :: cfx6,ccx6,cbx6,cfxp6,ccxp6,cbxp6,ciwxp6
   real(8),dimension(nxm) :: csxp6,cwxp6,csx6,cwx6,cifx6,cicx6,cisx6   
   real(8),dimension(nxm) :: cibx6,cifxp6,cicxp6,cibxp6,cisxp6,ciwx6
   real(8),dimension(nx) :: cfi6,cci6,cbi6,xx1 
   real(8),dimension(nx) :: cfip6,ccip6,cbip6,csip6,cwip6,csi6  
   real(8),dimension(nx) :: cwi6,cifi6,cici6,cibi6,cifip6,cicip6  
   real(8),dimension(nx) :: cibip6,cisip6,ciwip6,cisi6,ciwi6 
   real(8),dimension(nym) :: cfy6,ccy6,cby6 
   real(8),dimension(nym) :: cfyp6,ccyp6,cbyp6,csyp6,cwyp6,csy6 
   real(8),dimension(nym) :: cwy6,cify6,cicy6,ciby6,cifyp6  
   real(8),dimension(nym) :: cicyp6,cibyp6,cisyp6,ciwyp6,cisy6,ciwy6 
   real(8),dimension(ny) :: cfi6y,cci6y,cbi6y,cfip6y,ccip6y,cbip6y  
   real(8),dimension(ny) :: csip6y,cwip6y,csi6y,cwi6y,cifi6y,cici6y  
   real(8),dimension(ny) :: cibi6y,cifip6y,cicip6y,cibip6y  
   real(8),dimension(ny) :: cisip6y,ciwip6y,cisi6y,ciwi6y   
!**************************STRET*******************************************
   real(8),dimension(ny) :: ppy,pp2y,pp3y,pp4y,fux1,d1,d,hprime
   real(8),dimension(my) :: yp,yeta
   real(8),dimension(2) :: ja,jb
   real(8),dimension(nx) :: d5
   real(8),dimension(ny/2,m1+m2+1) :: a,a1,a2
   real(8),dimension(ny/2,m1) :: al,al1
   real(8),dimension(ny/2) :: indx,indx1,e,c
! ****************************2D******************************
! ****************************1D******************************
  real(8),dimension(nz) :: xlfux
! ************************************************************
  real(8) :: pi,twopi
  real(8) :: yinf,beta,den,xnum,alpha,xcx,den1,den2,den3,den4,xnum1,cst
  real(8) :: fourpi,r,x11,x,val,dec
! ************************************************************
  integer :: nxyz, ntime
  integer :: i,j,k,nxz,nxy
  integer :: longueur
  integer :: i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,&
             i15,i16,i17,i18,i19,i20
  integer :: j1,j2,j3,j4,jj
! ************************************************************
!  character(len=80) path_network,nom_network,path_files,&
!        nom_file,nom_script,nom_film,&
!        nom_fichier
!
!  character(len=80) lect_vort,avs_vort,lect_q,avs_q,lect_pdf,lect_moy,&
!                    lect_exp
  character(len=80) lect,fichier
  character(len=4)  suffixe

! 
  nxyz=nx*ny*nz
  pi=acos(-1.)
  twopi=2.*pi
!
   call parametre(ux,uy,uz,gx,gy,gz,sx,sy,sz,tx,ty,tz,di,&
        ffx,fcx,fbx,ffy,fcy,fby,ffz,fcz,fbz,&
        sfx,scx,sbx,sfy,scy,sby,sfz,scz,sbz,&
        fsx,fwx,fsy,fwy,fsz,fwz,&
        ssx,swx,ssy,swy,ssz,swz,&
        ffxp,fsxp,fwxp,&
        ffyp,fsyp,fwyp,&
        ffzp,fszp,fwzp,&
        sfxp,ssxp,swxp,&
        sfyp,ssyp,swyp,&
        sfzp,sszp,swzp,&
        nx,nxm,ny,nym,nz)
!
      yinf=-yly/2.
      beta=12.5
      den=2.*beta*yinf
      xnum=-yinf-sqrt(pi*pi*beta*beta+yinf*yinf)
      alpha=abs(xnum/den)
      xcx=1./beta/alpha
!     c'est ma fonction de stretching
      if (alpha.ne.0.) then 
      do j=1,ny
      if (ncly.eq.0) yeta(j)=(j-1.)*(1./ny)
      if (ncly.eq.1) yeta(j)=(j-1.)*(1./(ny-1.))
      den1=sqrt(alpha*beta+1.)
      xnum=den1/sqrt(alpha/pi)/sqrt(beta)/sqrt(pi)     
      den=2.*sqrt(alpha/pi)*sqrt(beta)*pi*sqrt(pi)
      den3=((sin(pi*yeta(j)))*(sin(pi*yeta(j)))/beta/pi)+alpha/pi
      den4=2.*alpha*beta-cos(2.*pi*yeta(j))+1.
      xnum1=(atan(xnum*tan(pi*yeta(j))))*den4/den1/den3/den
      cst=sqrt(beta)*pi/(2.*sqrt(alpha)*sqrt(alpha*beta+1.))
      if (yeta(j).lt.0.5) yp(j)=xnum1-cst-yinf
      if (yeta(j).eq.0.5) yp(j)=0.-yinf
      if (yeta(j).gt.0.5) yp(j)=xnum1+cst-yinf
      enddo
      endif
      if (alpha.eq.0.) then
         yp(1)=-1.e10
         do j=2,ny
            yeta(j)=(j-1.)*(1./ny)
            yp(j)=-beta*cos(pi*yeta(j))/sin(yeta(j)*pi)
         enddo
      endif
      do j=2,ny
         print *,'yy',yp(j),'dif',yp(j)-yp(j-1)
      enddo
!
      do j=1,ny
         ppy(j)=1./(alpha/pi+(1./pi/beta)*sin(pi*yeta(j))* &
            sin(pi*yeta(j)))
         pp2y(j)=ppy(j)*ppy(j)
         pp4y(j)=pp2y(j)*(-2./beta*cos(pi*yeta(j))*sin(pi*yeta(j)))
         hprime(j)=1./ppy(j)
      enddo
!
      open(10,file='moyt.dat',form='unformatted')
      read(10) uxmt
      read(10) uymt
      read(10) uzmt
      read(10) uxux
      read(10) uyuy
      read(10) uzuz
      read(10) uxuy
      close(10)                
!
      print *,'ntime attendu cms=',uxmt(1,1,1)
       print *,'ntime attendu cms=',uxmt(1,ny,1)
      print *,'ntime effectif =',uxmt(1,1,1)/100.
      print *,'ntime effectif =',uxmt(1,ny,1)/100.
      ntime=100
      print *,'ntime effectif =',ntime
!
      do i=1,nxyz
         uxmt(i,1,1)=uxmt(i,1,1)/ntime
         uymt(i,1,1)=uymt(i,1,1)/ntime
         uxux(i,1,1)=uxux(i,1,1)/ntime
         uyuy(i,1,1)=uyuy(i,1,1)/ntime
         uxuy(i,1,1)=uxuy(i,1,1)/ntime
         uxux(i,1,1)=uxux(i,1,1)-uxmt(i,1,1)*uxmt(i,1,1)
         uyuy(i,1,1)=uyuy(i,1,1)-uymt(i,1,1)*uymt(i,1,1)
         uxuy(i,1,1)=uxuy(i,1,1)-uxmt(i,1,1)*uymt(i,1,1)
         if (uxux(i,1,1).lt.0.) uxux(i,1,1)=0.
         uxux(i,1,1)=sqrt(uxux(i,1,1))
         if (uyuy(i,1,1).lt.0.) uyuy(i,1,1)=0.
         uyuy(i,1,1)=sqrt(uyuy(i,1,1))
      enddo
!      open(30,file='uxmt.dat',&
!              form='unformatted',status='unknown')
!      write(30) uxmt
!      close(30) 
!      stop           
!
      ra=1.
      cex=0.
      i1=nint((cex+30.*ra)/dx)+1  
      i2=nint((cex+45.*ra)/dx)+1  
      i3=nint((cex+60.*ra)/dx)+1  
      i4=nint((cex+75.*ra)/dx)+1 
!      i6=nint((cex+5.*ra)/dx)+1  
!      i7=nint((cex+6.*ra)/dx)+1  
!      i8=nint((cex+7.*ra)/dx)+1  
!      i9=nint((cex+8.*ra)/dx)+1 
!      i10=nint((cex+9.*ra)/dx)+1  
!      i11=nint((cex+10.*ra)/dx)+1  
!      i12=nint((cex+55.*ra)/dx)+1  
!      i13=nint((cex+60.*ra)/dx)+1  
!      i14=nint((cex+65.*ra)/dx)+1  
!      i15=nint((cex+70.*ra)/dx)+1  
!      i16=nint((cex+75.*ra)/dx)+1  
!      i17=nint((cex+80.*ra)/dx)+1  
!
      print *,'profils en i = ',i1,i2,i3,i4,i5,i6,i7,i8,i9
!
      open(20,file='uxmt.dat',form='formatted')
!      do i=1,nx
      do j=1,ny
         write(20,10) yp(j),uxmt(i1,j,1),uxmt(i2,j,1),&
         uxmt(i3,j,1),uxmt(i4,j,1)
!         ,uxmt(i5,j,1),&
!         uxmt(i6,j,1),uxmt(i7,j,1),uxmt(i8,j,1),&
!         uxmt(i9,j,1),uxmt(i10,j,1),uxmt(i11,j,1)
!         uxmt(i12,j,1),uxmt(i13,j,1),uxmt(i14,j,1),&
!         uxmt(i15,j,1),uxmt(i16,j,1),uxmt(i17,j,1)   
      enddo
         write(20,10) 
!      enddo
      close(20)
!
      open(20,file='uymt.dat',form='formatted')
!      do i=1,nx
      do j=1,ny
         write(20,10) yp(j),uymt(i1,j,1),uymt(i2,j,1),&
         uymt(i3,j,1),uymt(i4,j,1)
      enddo
         write(20,10) 
!      enddo
      close(20)
!
      open(20,file='uxux.dat',form='formatted')
!      do i=1,nx
      do j=1,ny
         write(20,10) yp(j),uxux(i1,j,1),uxux(i2,j,1),&
         uxux(i3,j,1),uxux(i4,j,1)
      enddo
         write(20,10) 
!      enddo
      close(20)
!
      open(20,file='uyuy.dat',form='formatted')
!      do i=1,nx
      do j=1,ny
         write(20,10) yp(j),uyuy(i1,j,1),uyuy(i2,j,1),&
         uyuy(i3,j,1),uyuy(i4,j,1)
      enddo
         write(20,10) 
!      enddo
      close(20)
 
     open(20,file='uxuy.dat',form='formatted')
!      do i=1,nx
      do j=1,ny
         write(20,10) yp(j),uxuy(i1,j,1),uxuy(i2,j,1),&
         uxuy(i3,j,1),uxuy(i4,j,1)
      enddo
         write(20,10) 
!      enddo
      close(20)
!
   10 format(6e30.15)
      end
!
!********************************************************************
!
subroutine parametre(ux,uy,uz,gx,gy,gz,&
     sx,sy,sz,tx,ty,tz,di,&
     ffx,fcx,fbx,ffy,fcy,fby,ffz,fcz,fbz,&
     sfx,scx,sbx,sfy,scy,sby,sfz,scz,sbz,&
     fsx,fwx,fsy,fwy,fsz,fwz,&
     ssx,swx,ssy,swy,ssz,swz,&
     ffxp,fsxp,fwxp,&
     ffyp,fsyp,fwyp,&
     ffzp,fszp,fwzp,&
     sfxp,ssxp,swxp,&
     sfyp,ssyp,swyp,&
     sfzp,sszp,swzp,&
     nx,nxm,ny,nym,nz)
!
!********************************************************************
!
   USE paramx_m 
   USE paramy_m 
   USE paramz_m 
   USE dericex_m 
   USE interpol4_m 
   USE dericex6_m 
   USE interpol6_m 
   USE dericey_m 
   USE interpoly4_m 
   USE dericey6_m 
   USE interpoly6_m 
   USE derivex_m 
   USE derivey_m 
   USE derivez_m 
   USE parfix_m 
   USE parfiy_m 
   USE parfiz_m 
   USE param1_m 
   USE param2_m 
   USE limit_m 
   USE ondelette_m 
   USE ecoulement_m 
   USE controle_m 
   USE sauvegarde_m 
   USE module_avs_m 
   USE delta_m 
   USE barreau_m 
!
   implicit none
!
   integer :: nx,ny,nz,nxm,nym
!
   real(8),dimension(nx,ny,nz) :: ux,uy,uz,gx,gy,gz,tx,ty,tz,di
   real(8),dimension(ny,nz) :: sx,px,vx,wx,xx,ppy,pp2y
   real(8),dimension(nx,nz) :: sy,py,vy,wy,xy
   real(8),dimension(nx,ny) :: sz,pz,vz,wz,xz
   real(8),dimension(nx) :: ffx,fcx,fbx,sfx,scx,sbx,fsx,fwx,ssx,swx
   real(8),dimension(nx) :: ffxp,fsxp,fwxp,sfxp,ssxp,swxp
   real(8),dimension(ny) :: ffy,fcy,fby,sfy,scy,sby,fsy,fwy,ssy,swy
   real(8),dimension(ny) :: ffyp,fsyp,fwyp,sfyp,ssyp,swyp
   real(8),dimension(nz) :: ffz,fcz,fbz,sfz,scz,sbz,fsz,fwz,ssz,swz
   real(8),dimension(nz) :: ffzp,fszp,fwzp,sfzp,sszp,swzp
   real(8),dimension(nx) :: uaxe_y,uaxe_z
   real(8),dimension(nxm) :: cfx6,ccx6,cbx6,cfxp6,ccxp6,cbxp6,ciwxp6
   real(8),dimension(nxm) :: csxp6,cwxp6,csx6,cwx6,cifx6,cicx6,cisx6   
   real(8),dimension(nxm) :: cibx6,cifxp6,cicxp6,cibxp6,cisxp6,ciwx6
   real(8),dimension(nx) :: cfi6,cci6,cbi6 
   real(8),dimension(nx) :: cfip6,ccip6,cbip6,csip6,cwip6,csi6  
   real(8),dimension(nx) :: cwi6,cifi6,cici6,cibi6,cifip6,cicip6  
   real(8),dimension(nx) :: cibip6,cisip6,ciwip6,cisi6,ciwi6 
   real(8),dimension(nym) :: cfy6,ccy6,cby6 
   real(8),dimension(nym) :: cfyp6,ccyp6,cbyp6,csyp6,cwyp6,csy6 
   real(8),dimension(nym) :: cwy6,cify6,cicy6,ciby6,cifyp6  
   real(8),dimension(nym) :: cicyp6,cibyp6,cisyp6,ciwyp6,cisy6,ciwy6 
   real(8),dimension(ny) :: cfi6y,cci6y,cbi6y,cfip6y,ccip6y,cbip6y  
   real(8),dimension(ny) :: csip6y,cwip6y,csi6y,cwi6y,cifi6y,cici6y  
   real(8),dimension(ny) :: cibi6y,cifip6y,cicip6y,cibip6y  
   real(8),dimension(ny) :: cisip6y,ciwip6y,cisi6y,ciwi6y   
   real(8) :: sc,theta,pi
   character*80 :: a
!
!   character(len=80) path_network,nom_network,path_files,&
!        nom_file,nom_script,nom_film,&
!        nom_fichier

   character(len=80) lect_vort,avs_vort,lect_q,avs_q,lect_pdf,&
                     lect_moy,lect_exp

   character(len=80) lect,fichier
   character(len=3)  suffixe
!
1000 format(A,80X)
1001 format(I6)
1002 format(F12.8)
!
   open (10,file='general.prm',status='unknown',form='formatted')
      read (10,1000) A
      read (10,1000) A
      read (10,1000) A
      read (10,1000) A
      read (10,1001) nclx
      read (10,1000) A
      read (10,1001) ncly
      read (10,1000) A
      read (10,1001) nclz
      read (10,1000) A
      read (10,1002) u1
      read (10,1000) A
      read (10,1002) u2
      read (10,1000) A
      read (10,1002) rayon
      read (10,1000) A
      read (10,1002) bruit
      read (10,1000) A
      read (10,1002) xlx
      read (10,1000) A
      read (10,1002) yly
      read (10,1000) A
      read (10,1002) zlz
      read (10,1000) A
      read (10,1002) re
      read (10,1000) A
      read (10,1002) sc
      read (10,1000) A
      read (10,1001) iecoule
      read (10,1000) A
      read (10,1001) iles
      read (10,1000) A
      read (10,1001) ifiltre
      read (10,1000) A
      read (10,1001) ivirtuel
      read (10,1000) A
      read (10,1001) iscalaire
      read (10,1000) A
      read (10,1001) ientree
      read (10,1000) A
      read (10,1001) idebut
      read (10,1000) A
      read (10,1001) ifin
      read (10,1000) A
      read (10,1001) nschema
      read (10,1000) A
      read (10,1001) nxboite
      read (10,1000) A
      read (10,1001) icontrole
      read (10,1000) A
      read (10,1001) iavs
   close(10)
!
   if (ivirtuel==1) then 
      open(10,file='barreau.prm',status='unknown',form='formatted') 
         read (10,1000) a 
         read (10,1000) a 
         read (10,1000) a 
         read (10,1000) a 
         read (10,1002) cex 
         read (10,1000) a 
         read (10,1002) cey 
         read (10,1000) a 
         read (10,1002) ra 
         read (10,1000) a 
         read (10,1002) alphav 
         read (10,1000) a 
         read (10,1002) betav 
         read (10,1000) a 
         read (10,1002) xkk 
      close(10) 
      print *,'  cex,cey=',cex,cey 
      print *,'  ra=',ra 
      print *,'  alphav=',alphav 
      print *,'  betav=',betav 
      print *,'  xkk=',xkk 
   endif
!
   if (icontrole==1) then 
      open(10,file='controle.prm',status='unknown',form='formatted') 
         read (10,1000) a 
         read (10,1000) a 
         read (10,1000) a 
         read (10,1000) a 
         read (10,1002) rlat 
         read (10,1000) a 
         read (10,1002) distance 
         read (10,1000) a 
         read (10,1002) phi 
         read (10,1000) a 
         read (10,1002) theta 
         read (10,1000) a 
         read (10,1002) ufreq 
         read (10,1000) a 
         read (10,1002) uampl 
         phi=(pi*phi)/180. 
         theta=(pi*theta)/180. 
      close(10) 
   endif
!

   print *,'nclx=',nclx
   print *,'ncly=',ncly
   print *,'nclz=',nclz
   print *,'u1,u2=',u1,u2
   print *,'rayon=',rayon
   print *,'bruit=',bruit
   print *,'Lx,Ly,Lz=',xlx,yly,zlz
   print *,'re,sc=',re,sc
   print *,'iecoule,iscalaire,ientree=',iecoule,iscalaire,ientree
   print *,'idebut,ifin=',idebut,ifin
   print *,'Nschema=',nschema
!
   if (nclx.eq.0) dx=xlx/nx
   if (nclx.eq.1) dx=xlx/(nx-1.)
   if (nclx.eq.2) dx=xlx/(nx-1.)
!
   if (ncly.eq.0) dy=yly/ny
   if (ncly.eq.1) dy=yly/(ny-1.)
   if (ncly.eq.2) dy=yly/(ny-1.)
!
   if (nz.gt.1) then
      if (nclz.eq.0) dz=zlz/nz
      if (nclz.eq.1) dz=zlz/(nz-1.)
      if (nclz.eq.2) dz=zlz/(nz-1.)
   else
      dz=zlz
   endif
!
   dx2=dx*dx
   dy2=dy*dy
   dz2=dz*dz
!
!
   return
end subroutine parametre

