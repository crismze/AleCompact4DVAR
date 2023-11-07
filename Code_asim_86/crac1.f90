!**********************************************************
!
program stats_cm
!
!***********************************************************
!
   USE paramx_m 
   USE paramy_m 
   USE paramz_m 
   USE derivex_m 
   USE derivey_m 
   USE derivez_m 
   USE param1_m 
   USE param2_m 
   USE limit_m 

!
   implicit none
!
   integer,parameter :: nx=181,ny=109,nz=1
!   integer,parameter :: nx=17,ny=33,nz=1
   integer,parameter :: m1=2,m2=2;
   integer,parameter :: nxm=nx-1,nym=ny-1,nzm=nz
   integer,parameter :: mx=nx+1,my=ny+1,mz=nz
   integer,parameter :: nyp=64,nzp=65,mw=4*mx*my*mz,nvec=256
   integer,parameter :: mxf=2*nx,myf=2*ny,mzf=2*nz,mb=2*mx
   integer,parameter :: mxy=mx*my,myz=my*mz
   integer,parameter :: nwork=max((3*nym+3)*(3*nxm+3),2*ny*nx)
   integer,parameter :: ntrigsX=4*nx,ntrigsY=4*ny,nfax=19
!
   real(8),dimension(nx,ny,nz) :: ux,uy,uz,gx,gy,gz,tx,ty,tz,di
   real(8),dimension(ny,nz) :: sx,px,vx,wx,xx
   real(8),dimension(nx,nz) :: sy,py,vy,wy,xy
   real(8),dimension(nx,ny) :: sz,pz,vz,wz,xz
   real(8),dimension(nx) :: ffx,fcx,fbx,sfx,scx,sbx,fsx,fwx,ssx,swx
   real(8),dimension(nx) :: ffxp,fsxp,fwxp,sfxp,ssxp,swxp
   real(8),dimension(ny) :: ffy,fcy,fby,sfy,scy,sby,fsy,fwy,ssy,swy
   real(8),dimension(ny) :: ffyp,fsyp,fwyp,sfyp,ssyp,swyp
   real(8),dimension(nz) :: ffz,fcz,fbz,sfz,scz,sbz,fsz,fwz,ssz,swz
   real(8),dimension(nz) :: ffzp,fszp,fwzp,sfzp,sszp,swzp
   real(8),dimension(nx) :: uaxe_y,uaxe_z
   integer :: i,j,k,nxz,nxy,nxyz,nboucle
   integer :: longueur
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
!**************************STRET*******************************************
   real(8),dimension(ny) :: ppy,pp2y,pp3y,pp4y,fux1,d1,d,hprime
   real(8),dimension(my) :: yp,yeta
   real(8),dimension(2) :: ja,jb
   real(8),dimension(nx) :: d5
   real(8),dimension(ny/2,m1+m2+1) :: a,a1,a2
   real(8),dimension(ny/2,m1) :: al,al1
   real(8),dimension(ny/2) :: indx,indx1,e,c
   real(8) :: yinf,beta,den,xnum,alpha,xcx,den1,den2,den3,den4,xnum1,cst
   real(8) :: pi,fourpi
!**************************************************************************
!
   character(len=80) path_network,nom_network,path_files,&
        nom_file,nom_script,nom_film,&
        nom_fichier
!
   character(len=80) lect_vort,avs_vort,lect_q,avs_q,lect_pdf,lect_moy,lect_exp
   character(len=80) lect,fichier
   character(len=3)  suffixe
!
   nxyz=nx*ny*nz
   nxy =nx*ny
   nxz =nx*nz
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
   pi=acos(-1.)
   fourpi=4.*pi
!     yinf cest la borne inf du domaine que l'on veut obtenir
   yinf=-yly/2.
      
      beta=2.5
!!!!!      beta=2.
      den=2.*beta*yinf
      xnum=-yinf-sqrt(pi*pi*beta*beta+yinf*yinf)
      alpha=abs(xnum/den)
      xcx=1./beta/alpha
!      print *,'alpha',alpha,beta,xcx,yly,yinf,den,xnum
 !     pause
      
!
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
      print *, 'min',yp(nym/2+1)-yp(nym/2),'=?',dx
!      pause'yy'

      do j=1,ny
         ppy(j)=1./(alpha/pi+(1./pi/beta)*sin(pi*yeta(j))* &
            sin(pi*yeta(j)))
         pp2y(j)=ppy(j)*ppy(j)
         pp4y(j)=pp2y(j)*(-2./beta*cos(pi*yeta(j))*sin(pi*yeta(j)))
         hprime(j)=1./ppy(j)
!         print *,ppy(j),pp2y(j)
      enddo
!      pause
!
   
   call schemas(ffx,fcx,fbx,ffy,fcy,fby,ffz,fcz,fbz,sfx,scx,sbx,sfy,&
     scy,sby,sfz,scz,sbz,fsx,fwx,fsy,fwy,fsz,fwz,ssx,swx,ssy,swy,ssz,swz,&
     ffxp,fsxp,fwxp,ffyp,fsyp,fwyp,ffzp,fszp,fwzp,sfxp,ssxp,swxp,sfyp,ssyp,&
     swyp,sfzp,sszp,swzp,nx,nxm,ny,nym,nz,cfx6,ccx6,cbx6,cfxp6,ccxp6,cbxp6,&
     csxp6,cwxp6,csx6,cwx6,cifx6,cicx6,cibx6,cifxp6,cicxp6,cibxp6,cisxp6,&
     ciwxp6,cisx6,ciwx6,cfi6,cci6,cbi6,cfip6,ccip6,cbip6,csip6,cwip6,csi6,&
     cwi6,cifi6,cici6,cibi6,cifip6,cicip6,cibip6,cisip6,ciwip6,cisi6,ciwi6,&
     cfy6,ccy6,cby6,cfyp6,ccyp6,cbyp6,csyp6,cwyp6,csy6,&
     cwy6,cify6,cicy6,ciby6,cifyp6,cicyp6,cibyp6,cisyp6,ciwyp6,cisy6,ciwy6,&
     cfi6y,cci6y,cbi6y,cfip6y,ccip6y,cbip6y,csip6y,cwip6y,csi6y,cwi6y,&
     cifi6y,cici6y,cibi6y,cifip6y,cicip6y,cibip6y,cisip6y,ciwip6y,cisi6y,&
     ciwi6y,ppy,pp2y)

   do i=1,nx
      uaxe_y(i)=0.
      uaxe_z(i)=0.
   enddo
!
   do nboucle=1,1000,1
      call numcar(nboucle,suffixe)
      nom_file='test0'//suffixe
      longueur=index(nom_file,' ')-1
      open(10,file=nom_file(1:longueur),&
           form='unformatted',&
           status='unknown')
      read(10) ux,uy
      close(10)
! 
      call  calcwz(ux,uy,uz,gx,gy,gz,&
           tx,ty,tz,di,sx,sy,sz,&
           ffx,fsx,fwx,ffy,fsy,fwy,ffz,fsz,fwz,&
           ffxp,fsxp,fwxp,&
           ffyp,fsyp,fwyp,&
           ffzp,fszp,fwzp,&
           nx,ny,nz)
!
      call numcar(nboucle,suffixe)
      nom_file='test'//suffixe
      longueur=index(nom_file,' ')-1
      open(20,file=nom_file(1:longueur)//'.avs',&
           form='unformatted',status='unknown')
           write(20) tz
      close(20)
   enddo
!
end program stats_cm
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
   USE derivex_m 
   USE derivey_m 
   USE derivez_m 
   USE param1_m 
   USE param2_m 
   USE limit_m 
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
   integer :: iecoule,iles,ifiltre,ivirtuel,iscalaire,ientree,idebut,ifin
   integer :: nschema,nxboite,icontrole,iavs
   real(8) :: a,u1,u2,rayon,re,sc
!
   character(len=80) path_network,nom_network,path_files,&
        nom_file,nom_script,nom_film,&
        nom_fichier
!
   character(len=80) lect_vort,avs_vort,lect_q,avs_q,lect_pdf,lect_moy,lect_exp
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
!     
!********************************************************************
!
subroutine calcwz(ux,uy,uz,gx,gy,gz,&
     tx,ty,tz,di,sx,sy,sz,&
     ffx,fsx,fwx,ffy,fsy,fwy,ffz,fsz,fwz,&
     ffxp,fsxp,fwxp,&
     ffyp,fsyp,fwyp,&
     ffzp,fszp,fwzp,&
     nx,ny,nz)
!
!********************************************************************
!
   implicit none
!
   integer :: nx,ny,nz,i,j,k
!
   real(8),dimension(nx,ny,nz) :: ux,uy,uz,gx,gy,gz,tx,ty,tz,di
   real(8),dimension(ny,nz) :: sx
   real(8),dimension(nx,nz) :: sy
   real(8),dimension(nx,ny) :: sz
   real(8),dimension(nx) :: ffx,fcx,fbx,fsx,fwx
   real(8),dimension(nx) :: ffxp,fsxp,fwxp
   real(8),dimension(ny) :: ffy,fcy,fby,fsy,fwy
   real(8),dimension(ny) :: ffyp,fsyp,fwyp
   real(8),dimension(nz) :: ffz,fcz,fbz,fsz,fwz
   real(8),dimension(nz) :: ffzp,fszp,fwzp
   real(8) :: wzmin,wzmax
!
   wzmin= 1.e10
   wzmax=-1.e10
!
   call derx (tx,uy,gx,sx,ffxp,fsxp,fwxp,nx,ny,nz,1)
   call dery (ty,ux,gx,di,sy,ffyp,fsyp,fwyp,nx,ny,nz,1)
!
   do k=1,nz
   do j=1,ny
   do i=1,nx
      tz(i,j,k)=(tx(i,j,k)-ty(i,j,k))
   enddo
!   print *,ffyp(j),fsyp(j),fwyp(j)
   enddo
!   pause
   enddo
!
   return
end subroutine calcwz
!
!********************************************************************
!
subroutine calcwn(ux,uy,uz,gx,gy,gz,&
     tx,ty,tz,di,sx,sy,sz,&
     ffx,fsx,fwx,ffy,fsy,fwy,ffz,fsz,fwz,&
     ffxp,fsxp,fwxp,&
     ffyp,fsyp,fwyp,&
     ffzp,fszp,fwzp,&
     nx,ny,nz)
!
!********************************************************************
!
   implicit none
!
   integer :: nx,ny,nz,i,j,k
   real(8) :: wzmin,wzmax
!
   real(8),dimension(nx,ny,nz) :: ux,uy,uz,gx,gy,gz,tx,ty,tz,di
   real(8),dimension(ny,nz) :: sx
   real(8),dimension(nx,nz) :: sy
   real(8),dimension(nx,ny) :: sz
   real(8),dimension(nx) :: ffx,fcx,fbx,fsx,fwx
   real(8),dimension(nx) :: ffxp,fsxp,fwxp
   real(8),dimension(ny) :: ffy,fcy,fby,fsy,fwy
   real(8),dimension(ny) :: ffyp,fsyp,fwyp
   real(8),dimension(nz) :: ffz,fcz,fbz,fsz,fwz
   real(8),dimension(nz) :: ffzp,fszp,fwzp
!
   wzmin= 1.e10
   wzmax=-1.e10
!
   call derx (tx,uy,gx,sx,ffxp,fsxp,fwxp,nx,ny,nz,1)
   call dery (ty,ux,gx,di,sy,ffyp,fsyp,fwyp,nx,ny,nz,1)
!
   do k=1,nz
   do j=1,ny
   do i=1,nx
      tz(i,j,k)=tx(i,j,k)-ty(i,j,k)
      if (tz(i,j,k).lt.wzmin) wzmin=tz(i,j,k)
      if (tz(i,j,k).gt.wzmax) wzmax=tz(i,j,k)
   enddo
   enddo
   enddo
!
   call derz(gx,ux,gz,sz,ffzp,fszp,fwzp,nx,ny,nz,1)
   call derx(gy,uz,gz,sx,ffxp,fsxp,fwxp,nx,ny,nz,1)
!
   do k=1,nz
   do j=1,ny
   do i=1,nx
      ty(i,j,k)=gx(i,j,k)-gy(i,j,k)
   enddo
   enddo
   enddo
!
   call dery (gy,uz,gx,di,sy,ffyp,fsyp,fwyp,nx,ny,nz,1)
   call derz (gz,uy,gx,sz,ffzp,fszp,fwzp,nx,ny,nz,1)
!
   do k=1,nz
   do j=1,ny
   do i=1,nx
      tx(i,j,k)=gy(i,j,k)-gz(i,j,k)
   enddo
   enddo
   enddo
!
   do k=1,nz
   do j=1,ny
   do i=1,nx
      tz(i,j,k)=sqrt(tx(i,j,k)*tx(i,j,k)+&
           ty(i,j,k)*ty(i,j,k)+&
           tz(i,j,k)*tz(i,j,k))
   enddo
   enddo
   enddo 
!
   return
end subroutine calcwn
!
!********************************************************************
!
subroutine calcq(ux,uy,uz,gx,gy,gz,&
     wx,wy,wz,di,sx,sy,sz,&
     ffx,fsx,fwx,ffy,fsy,fwy,ffz,fsz,fwz,&
     ffxp,fsxp,fwxp,&
     ffyp,fsyp,fwyp,&
     ffzp,fszp,fwzp,&
     nx,ny,nz)
!
!********************************************************************
!
   implicit none
!
   integer :: nx,ny,nz,i,j,k
   real(8) :: wzmin,wzmax
!
   real(8),dimension(nx,ny,nz) :: ux,uy,uz,gx,gy,gz,tx,ty,tz,di,wx,wy,wz
   real(8),dimension(ny,nz) :: sx
   real(8),dimension(nx,nz) :: sy
   real(8),dimension(nx,ny) :: sz
   real(8),dimension(nx) :: ffx,fcx,fbx,fsx,fwx
   real(8),dimension(nx) :: ffxp,fsxp,fwxp
   real(8),dimension(ny) :: ffy,fcy,fby,fsy,fwy
   real(8),dimension(ny) :: ffyp,fsyp,fwyp
   real(8),dimension(nz) :: ffz,fcz,fbz,fsz,fwz
   real(8),dimension(nz) :: ffzp,fszp,fwzp
!
   wzmin= 1.e10
   wzmax=-1.e10
!
   call derx (wx,ux,gx,sx,ffxp,fsxp,fwxp,nx,ny,nz,1)
   call dery (wy,uy,gy,di,sy,ffyp,fsyp,fwyp,nx,ny,nz,1)
   call derz (wz,uz,gz,sz,ffzp,fszp,fwzp,nx,ny,nz,1)
!
   call dery (gx,ux,gz,di,sy,ffyp,fsyp,fwyp,nx,ny,nz,1)
   call derx (gy,uy,gz,sx,ffxp,fsxp,fwxp,nx,ny,nz,1)
!
   do k=1,nz
   do j=1,ny
   do i=1,nx
      wx(i,j,k)=-(wx(i,j,k)*wx(i,j,k)+&
           wy(i,j,k)*wy(i,j,k)+&
           wz(i,j,k)*wz(i,j,k)+&
           gx(i,j,k)*gy(i,j,k))
   enddo
   enddo
   enddo
!
   call derz (gx,ux,gy,sz,ffzp,fszp,fwzp,nx,ny,nz,1)
   call derx (gz,uz,gy,sx,ffxp,fsxp,fwxp,nx,ny,nz,1)
!
   call dery (wz,uz,gy,di,sy,ffyp,fsyp,fwyp,nx,ny,nz,1)
   call derz (wy,uy,gy,sz,ffzp,fszp,fwzp,nx,ny,nz,1)
!
   do k=1,nz
   do j=1,ny
   do i=1,nx
      wx(i,j,k)=wx(i,j,k)-&
           (gx(i,j,k)*gz(i,j,k)+&
           wz(i,j,k)*wy(i,j,k))
   enddo
   enddo
   enddo
!
   wzmin=1.e8
   wzmax=-1.e8
   do k=1,nz
   do j=1,ny
   do i=1,nx
      if (wx(i,j,k).lt.wzmin/1.) wzmin=wx(i,j,k)
      if (wx(i,j,k).gt.wzmax/1.) wzmax=wx(i,j,k)
   enddo
   enddo
   enddo 
   print *,'Qmax,Qmin =',wzmax,wzmin
!
   return
end subroutine calcq
!
!*****************************************************************
!
subroutine numcar(num,car)
!
!*****************************************************************
!
   implicit none
!
   character(len=3) car
   integer :: num
!
   if (num.ge.100) then
      write (car,1) num
1     format (i3)
   else
      if (num.ge.10) then
         write (car,2) num
2        format ('0',i2)
      else
         write (car,3) num
3        format ('00',i1)
      endif
   endif
!      
   return
end subroutine numcar
!
!*********************************************************************
!
subroutine annule (gx,nxyz)
!
!*********************************************************************
!
  implicit none
!
  real(8),dimension(nxyz) :: gx
  integer :: ijk,nxyz
!
  do ijk=1,nxyz
     gx(ijk)=0.
  enddo
!
  return
end subroutine annule
!
!********************************************************************
!********************************************************************
!
subroutine schemas(ffx,fcx,fbx,ffy,fcy,fby,ffz,fcz,fbz,sfx,scx,sbx,sfy,&
     scy,sby,sfz,scz,sbz,fsx,fwx,fsy,fwy,fsz,fwz,ssx,swx,ssy,swy,ssz,swz,&
     ffxp,fsxp,fwxp,ffyp,fsyp,fwyp,ffzp,fszp,fwzp,sfxp,ssxp,swxp,sfyp,ssyp,&
     swyp,sfzp,sszp,swzp,nx,nxm,ny,nym,nz,cfx6,ccx6,cbx6,cfxp6,ccxp6,cbxp6,&
     csxp6,cwxp6,csx6,cwx6,cifx6,cicx6,cibx6,cifxp6,cicxp6,cibxp6,cisxp6,&
     ciwxp6,cisx6,ciwx6,cfi6,cci6,cbi6,cfip6,ccip6,cbip6,csip6,cwip6,csi6,&
     cwi6,cifi6,cici6,cibi6,cifip6,cicip6,cibip6,cisip6,ciwip6,cisi6,ciwi6,&
     cfy6,ccy6,cby6,cfyp6,ccyp6,cbyp6,csyp6,cwyp6,csy6,&
     cwy6,cify6,cicy6,ciby6,cifyp6,cicyp6,cibyp6,cisyp6,ciwyp6,cisy6,ciwy6,&
     cfi6y,cci6y,cbi6y,cfip6y,ccip6y,cbip6y,csip6y,cwip6y,csi6y,cwi6y,&
     cifi6y,cici6y,cibi6y,cifip6y,cicip6y,cibip6y,cisip6y,ciwip6y,cisi6y,&
     ciwi6y,ppy,pp2y)
!
!********************************************************************
!
  USE paramx_m
  USE paramy_m
  USE paramz_m
  USE derivex_m
  USE derivey_m
  USE derivez_m
  USE dericex6_m 
  USE interpol6_m 
  USE dericey6_m 
  USE interpoly6_m 
!
  implicit none
!
  real(8),dimension(nx) :: ffx,fcx,fbx,sfx,scx,sbx,fsx,fwx,ssx,swx 
  real(8),dimension(nx) :: ffxp,fsxp,fwxp,sfxp,ssxp,swxp
  real(8),dimension(ny) :: ffy,fcy,fby,sfy,scy,sby,fsy,fwy,ssy,swy 
  real(8),dimension(ny) :: ffyp,fsyp,fwyp,sfyp,ssyp,swyp
  real(8),dimension(nz) :: ffz,fcz,fbz,sfz,scz,sbz,fsz,fwz,ssz,swz 
  real(8),dimension(nz) :: ffzp,fszp,fwzp,sfzp,sszp,swzp
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
!**************************STRET***************************
   real(8),dimension(ny) :: ppy,pp2y
!**********************************************************  
!
      integer  :: nx,nxm,ny,nym,nz
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, j, k 
!-----------------------------------------------
! 
!
      alfa1x= 2.
      af1x  =-(5./2.  )/dx
      bf1x  = (   2.  )/dx
      cf1x  = (1./2.  )/dx
      df1x  = 0.
      alfa2x= 1./4.
      af2x  = (3./4.  )/dx
!
      alfanx= 2.
      afnx  =-(5./2.  )/dx
      bfnx  = (   2.  )/dx
      cfnx  = (1./2.  )/dx
      dfnx  = 0.
      alfamx= 1./4.
      afmx  = (3./4.  )/dx
!
      alfaix= 1./3.
      afix  = (7./9.  )/dx
      bfix  = (1./36. )/dx
!
      alsa1x= 11.
      as1x  = (13.    )/dx2
      bs1x  =-(27.    )/dx2
      cs1x  = (15.    )/dx2
      ds1x  =-(1.     )/dx2
      alsa2x= 1./10.
      as2x  = (6./5.  )/dx2
!
      alsanx= 11.
      asnx  = (13.    )/dx2
      bsnx  =-(27.    )/dx2
      csnx  = (15.    )/dx2
      dsnx  =-(1.     )/dx2
      alsamx= 1./10.
      asmx  = (6./5.  )/dx2
!
      alsaix= 2./11.
      asix  = (12./11.)/dx2
      bsix  = (3./44. )/dx2
!
      alfa1y= 2.
      af1y  =-(5./2.  )/dy
      bf1y  = (   2.  )/dy
      cf1y  = (1./2.  )/dy
      df1y  = 0.
      alfa2y= 1./4.
      af2y  = (3./4.  )/dy
!
      alfany= 2.
      afny  =-(5./2.  )/dy
      bfny  = (   2.  )/dy
      cfny  = (1./2.  )/dy
      dfny  = 0.
      alfamy= 1./4.
      afmy  = (3./4.  )/dy
!
      alfajy= 1./3.
      afjy  = (7./9.  )/dy
      bfjy  = (1./36. )/dy
!
      alsa1y= 11.
      as1y  = (13.    )/dy2
      bs1y  =-(27.    )/dy2
      cs1y  = (15.    )/dy2
      ds1y  =-(1.     )/dy2
      alsa2y= 1./10.
      as2y  = (6./5.  )/dy2
!
      alsany= 11.
      asny  = (13.    )/dy2
      bsny  =-(27.    )/dy2
      csny  = (15.    )/dy2
      dsny  =-(1.     )/dy2
      alsamy= 1./10.
      asmy  = (6./5.  )/dy2
!
      alsajy= 2./11.
      asjy  = (12./11.)/dy2
      bsjy  = (3./44. )/dy2
!
!*************COMPACT 6 en X***************************
   alcaix6=9./62. 
   acix6=(63./62.)/dx 
   bcix6=(17./62.)/3./dx 
!
   cfx6(1)=alcaix6 
   cfx6(2)=alcaix6 
   cfx6(nxm-2)=alcaix6 
   cfx6(nxm-1)=alcaix6 
   cfx6(nxm)=0. 
   ccx6(1)=1. + alcaix6 
   ccx6(2)=1. 
   ccx6(nxm-2)=1. 
   ccx6(nxm-1)=1. 
   ccx6(nxm)=1. + alcaix6 
   cbx6(1)=alcaix6 
   cbx6(2)=alcaix6 
   cbx6(nxm-2)=alcaix6 
   cbx6(nxm-1)=alcaix6 
   cbx6(nxm)=0. 
   do i=3,nxm-3 
      cfx6(i)=alcaix6 
      ccx6(i)=1. 
      cbx6(i)=alcaix6 
   enddo
!********************************
   cfi6(1)=alcaix6 + alcaix6 
   cfi6(2)=alcaix6 
   cfi6(nx-2)=alcaix6 
   cfi6(nx-1)=alcaix6 
   cfi6(nx)=0. 
   cci6(1)=1. 
   cci6(2)=1. 
   cci6(nx-2)=1. 
   cci6(nx-1)=1. 
   cci6(nx)=1. 
   cbi6(1)=alcaix6 
   cbi6(2)=alcaix6 
   cbi6(nx-2)=alcaix6 
   cbi6(nx-1)=alcaix6 + alcaix6 
   cbi6(nx)=0. 
   do i=3,nx-3 
      cfi6(i)=alcaix6 
      cci6(i)=1. 
      cbi6(i)=alcaix6 
   enddo
!**************INTERPOLATION ORDRE 6 en X**************
   ailcaix6=3./10. 
   aicix6=3./4. 
   bicix6=1./20. 
!
   cifx6(1)=ailcaix6 
   cifx6(2)=ailcaix6 
   cifx6(nxm-2)=ailcaix6 
   cifx6(nxm-1)=ailcaix6 
   cifx6(nxm)=0. 
   cicx6(1)=1. + ailcaix6 
   cicx6(2)=1. 
   cicx6(nxm-2)=1. 
   cicx6(nxm-1)=1. 
   cicx6(nxm)=1. + ailcaix6 
   cibx6(1)=ailcaix6 
   cibx6(2)=ailcaix6 
   cibx6(nxm-2)=ailcaix6 
   cibx6(nxm-1)=ailcaix6 
   cibx6(nxm)=0. 
   do i=3,nxm-3 
      cifx6(i)=ailcaix6 
      cicx6(i)=1. 
      cibx6(i)=ailcaix6 
   enddo
!*************************
   cifi6(1)=ailcaix6 + ailcaix6 
   cifi6(2)=ailcaix6 
   cifi6(nx-2)=ailcaix6 
   cifi6(nx-1)=ailcaix6 
   cifi6(nx)=0. 
   cici6(1)=1. 
   cici6(2)=1. 
   cici6(nx-2)=1. 
   cici6(nx-1)=1. 
   cici6(nx)=1. 
   cibi6(1)=ailcaix6 
   cibi6(2)=ailcaix6 
   cibi6(nx-2)=ailcaix6 
   cibi6(nx-1)=ailcaix6 + ailcaix6 
   cibi6(nx)=0. 
   do i=3,nx-3 
      cifi6(i)=ailcaix6 
      cici6(i)=1. 
      cibi6(i)=ailcaix6 
   enddo
!*************************************************
!*************COMPACT 6 en Y**********************
   alcaiy6=9./62. 
   aciy6=(63./62.)/dy 
   bciy6=(17./62.)/3./dy 
!
   if (ncly==1) then 
      cfy6(1)=alcaiy6 
      cfy6(2)=alcaiy6 
      cfy6(nym-2)=alcaiy6 
      cfy6(nym-1)=alcaiy6 
      cfy6(nym)=0. 
      ccy6(1)=1. + alcaiy6 
      ccy6(2)=1. 
      ccy6(nym-2)=1. 
      ccy6(nym-1)=1. 
      ccy6(nym)=1. + alcaiy6 
      cby6(1)=alcaiy6 
      cby6(2)=alcaiy6 
      cby6(nym-2)=alcaiy6 
      cby6(nym-1)=alcaiy6 
      cby6(nym)=0. 
      do j=3,nym-3 
         cfy6(j)=alcaiy6 
         ccy6(j)=1. 
         cby6(j)=alcaiy6 
      enddo
!********************************
      cfi6y(1)=alcaiy6 + alcaiy6 
      cfi6y(2)=alcaiy6 
      cfi6y(ny-2)=alcaiy6 
      cfi6y(ny-1)=alcaiy6 
      cfi6y(ny)=0. 
      cci6y(1)=1. 
      cci6y(2)=1. 
      cci6y(ny-2)=1. 
      cci6y(ny-1)=1. 
      cci6y(ny)=1. 
      cbi6y(1)=alcaiy6 
      cbi6y(2)=alcaiy6 
      cbi6y(ny-2)=alcaiy6 
      cbi6y(ny-1)=alcaiy6 + alcaiy6 
      cbi6y(ny)=0. 
      do j=3,ny-3 
         cfi6y(j)=alcaiy6 
         cci6y(j)=1. 
         cbi6y(j)=alcaiy6 
      enddo
   endif
!**************INTERPOLATION ORDRE 6 en Y**************
   ailcaiy6=3./10. 
   aiciy6=3./4. 
   biciy6=1./20. 
!
   if (ncly==1) then 
      cify6(1)=ailcaiy6 
      cify6(2)=ailcaiy6 
      cify6(nym-2)=ailcaiy6 
      cify6(nym-1)=ailcaiy6 
      cify6(nym)=0. 
      cicy6(1)=1. + ailcaiy6 
      cicy6(2)=1. 
      cicy6(nym-2)=1. 
      cicy6(nym-1)=1. 
      cicy6(nym)=1. + ailcaiy6 
      ciby6(1)=ailcaiy6 
      ciby6(2)=ailcaiy6 
      ciby6(nym-2)=ailcaiy6 
      ciby6(nym-1)=ailcaiy6 
      ciby6(nym)=0. 
      do j=3,nym-3 
         cify6(j)=ailcaiy6 
         cicy6(j)=1. 
         ciby6(j)=ailcaiy6 
      enddo
!*************************
      cifi6y(1)=ailcaiy6 + ailcaiy6 
      cifi6y(2)=ailcaiy6 
      cifi6y(ny-2)=ailcaiy6 
      cifi6y(ny-1)=ailcaiy6 
      cifi6y(ny)=0. 
      cici6y(1)=1. 
      cici6y(2)=1. 
      cici6y(ny-2)=1. 
      cici6y(ny-1)=1. 
      cici6y(ny)=1. 
      cibi6y(1)=ailcaiy6 
      cibi6y(2)=ailcaiy6 
      cibi6y(ny-2)=ailcaiy6 
      cibi6y(ny-1)=ailcaiy6 + ailcaiy6 
      cibi6y(ny)=0. 
      do j=3,ny-3 
         cifi6y(j)=ailcaiy6 
         cici6y(j)=1. 
         cibi6y(j)=ailcaiy6 
      enddo
   endif
!*************************************************
!
      if (nz.gt.1) then
         alfa1z= 2.
         af1z  =-(5./2.  )/dz
         bf1z  = (   2.  )/dz
         cf1z  = (1./2.  )/dz
         df1z  = 0.
         alfa2z= 1./4.
         af2z  = (3./4.  )/dz
!
         alfanz= 2.
         afnz  =-(5./2.  )/dz
         bfnz  = (   2.  )/dz
         cfnz  = (1./2.  )/dz
         dfnz  = 0.
         alfamz= 1./4.
         afmz  = (3./4.  )/dz
!
         alfakz= 1./3.
         afkz  = (7./9.  )/dz
         bfkz  = (1./36. )/dz
!
         alsa1z= 11.
         as1z  = (13.    )/dz2
         bs1z  =-(27.    )/dz2
         cs1z  = (15.    )/dz2
         ds1z  =-(1.     )/dz2
         alsa2z= 1./10.
         as2z  = (6./5.  )/dz2
!
         alsanz= 11.
         asnz  = (13.    )/dz2
         bsnz  =-(27.    )/dz2
         csnz  = (15.    )/dz2
         dsnz  =-(1.     )/dz2
         alsamz= 1./10.
         asmz  = (6./5.  )/dz2
!
         alsakz= 2./11.
         askz  = (12./11.)/dz2
         bskz  = (3./44. )/dz2
      endif
!
      if (nclx.eq.0) then
         ffx(1)   =alfaix
         ffx(2)   =alfaix
         ffx(nx-2)=alfaix
         ffx(nx-1)=alfaix
         ffx(nx)  =0.
         fcx(1)   =2.
         fcx(2)   =1.
         fcx(nx-2)=1.
         fcx(nx-1)=1.
         fcx(nx  )=1.+alfaix*alfaix
         fbx(1)   =alfaix
         fbx(2)   =alfaix
         fbx(nx-2)=alfaix
         fbx(nx-1)=alfaix
         fbx(nx  )=0.
         do i=3,nx-3
            ffx(i)=alfaix
            fcx(i)=1.
            fbx(i)=alfaix
         enddo
      endif
!
      if (nclx.eq.1) then
         ffx(1)   =alfaix+alfaix
         ffx(2)   =alfaix
         ffx(nx-2)=alfaix
         ffx(nx-1)=alfaix
         ffx(nx)  =0.
         fcx(1)   =1.
         fcx(2)   =1.
         fcx(nx-2)=1.
         fcx(nx-1)=1.
         fcx(nx  )=1.
         fbx(1)   =alfaix 
         fbx(2)   =alfaix
         fbx(nx-2)=alfaix
         fbx(nx-1)=alfaix+alfaix
         fbx(nx  )=0.
         do i=3,nx-3
            ffx(i)=alfaix
            fcx(i)=1.
            fbx(i)=alfaix
         enddo
      endif
!
      if (nclx.eq.2) then
         ffx(1)   =alfa1x
         ffx(2)   =alfa2x
         ffx(nx-2)=alfaix
         ffx(nx-1)=alfamx
         ffx(nx)  =0.
         fcx(1)   =1.
         fcx(2)   =1.
         fcx(nx-2)=1.
         fcx(nx-1)=1.
         fcx(nx  )=1.
         fbx(1)   =alfa2x 
         fbx(2)   =alfaix
         fbx(nx-2)=alfamx
         fbx(nx-1)=alfanx
         fbx(nx  )=0.
         do i=3,nx-3
            ffx(i)=alfaix
            fcx(i)=1.
            fbx(i)=alfaix
         enddo
      endif
!
      if (ncly.eq.0) then
         ffy(1)   =alfajy
         ffy(2)   =alfajy
         ffy(ny-2)=alfajy
         ffy(ny-1)=alfajy
         ffy(ny)  =0.
         fcy(1)   =2.
         fcy(2)   =1.
         fcy(ny-2)=1.
         fcy(ny-1)=1.
         fcy(ny  )=1.+alfajy*alfajy
         fby(1)   =alfajy
         fby(2)   =alfajy
         fby(ny-2)=alfajy
         fby(ny-1)=alfajy
         fby(ny  )=0.
         do j=3,ny-3
            ffy(j)=alfajy
            fcy(j)=1.
            fby(j)=alfajy
         enddo
      endif
!
      if (ncly.eq.1) then
         ffy(1)   =alfajy+alfajy
         ffy(2)   =alfajy
         ffy(ny-2)=alfajy
         ffy(ny-1)=alfajy
         ffy(ny)  =0.
         fcy(1)   =1.
         fcy(2)   =1.
         fcy(ny-2)=1.
         fcy(ny-1)=1.
         fcy(ny  )=1.
         fby(1)   =alfajy 
         fby(2)   =alfajy
         fby(ny-2)=alfajy
         fby(ny-1)=alfajy+alfajy
         fby(ny  )=0.
         do j=3,ny-3
            ffy(j)=alfajy
            fcy(j)=1.
            fby(j)=alfajy
         enddo
      endif
!
      if (ncly.eq.2) then
         ffy(1)   =alfa1y
         ffy(2)   =alfa2y
         ffy(ny-2)=alfajy
         ffy(ny-1)=alfamy
         ffy(ny)  =0.
         fcy(1)   =1.
         fcy(2)   =1.
         fcy(ny-2)=1.
         fcy(ny-1)=1.
         fcy(ny  )=1.
         fby(1)   =alfa2y 
         fby(2)   =alfajy
         fby(ny-2)=alfamy
         fby(ny-1)=alfany
         fby(ny  )=0.
         do j=3,ny-3
            ffy(j)=alfajy
            fcy(j)=1.
            fby(j)=alfajy
         enddo
      endif
!
      if (nz.gt.1) then
         if (nclz.eq.0) then
            ffz(1)   =alfakz
            ffz(2)   =alfakz
            ffz(nz-2)=alfakz
            ffz(nz-1)=alfakz
            ffz(nz)  =0.
            fcz(1)   =2.
            fcz(2)   =1.
            fcz(nz-2)=1.
            fcz(nz-1)=1.
            fcz(nz  )=1.+alfakz*alfakz
            fbz(1)   =alfakz
            fbz(2)   =alfakz
            fbz(nz-2)=alfakz
            fbz(nz-1)=alfakz
            fbz(nz  )=0.
            do k=3,nz-3
               ffz(k)=alfakz
               fcz(k)=1.
               fbz(k)=alfakz
            enddo
         endif
!
         if (nclz.eq.1) then
            ffz(1)   =alfakz+alfakz
            ffz(2)   =alfakz
            ffz(nz-2)=alfakz
            ffz(nz-1)=alfakz
            ffz(nz)  =0.
            fcz(1)   =1.
            fcz(2)   =1.
            fcz(nz-2)=1.
            fcz(nz-1)=1.
            fcz(nz  )=1.
            fbz(1)   =alfakz 
            fbz(2)   =alfakz
            fbz(nz-2)=alfakz
            fbz(nz-1)=alfakz+alfakz
            fbz(nz  )=0.
            do k=3,nz-3
               ffz(k)=alfakz
               fcz(k)=1.
               fbz(k)=alfakz
            enddo
         endif
!
         if (nclz.eq.2) then
            ffz(1)   =alfa1z
            ffz(2)   =alfa2z
            ffz(nz-2)=alfakz
            ffz(nz-1)=alfamz
            ffz(nz)  =0.
            fcz(1)   =1.
            fcz(2)   =1.
            fcz(nz-2)=1.
            fcz(nz-1)=1.
            fcz(nz  )=1.
            fbz(1)   =alfa2z 
            fbz(2)   =alfakz
            fbz(nz-2)=alfamz
            fbz(nz-1)=alfanz
            fbz(nz  )=0.
            do k=3,nz-3
               ffz(k)=alfakz
               fcz(k)=1.
               fbz(k)=alfakz
            enddo
         endif
      endif
!
      if (nclx.eq.0) then
         sfx(1)   =alsaix
         sfx(2)   =alsaix
         sfx(nx-2)=alsaix
         sfx(nx-1)=alsaix
         sfx(nx)  =0.
         scx(1)   =2.
         scx(2)   =1.
         scx(nx-2)=1.
         scx(nx-1)=1.
         scx(nx  )=1.+alsaix*alsaix
         sbx(1)   =alsaix
         sbx(2)   =alsaix
         sbx(nx-2)=alsaix
         sbx(nx-1)=alsaix
         sbx(nx  )=0.
         do i=3,nx-3
            sfx(i)=alsaix
            scx(i)=1.
            sbx(i)=alsaix
         enddo
      endif
!
      if (nclx.eq.1) then
         sfx(1)   =alsaix+alsaix
         sfx(2)   =alsaix
         sfx(nx-2)=alsaix
         sfx(nx-1)=alsaix
         sfx(nx)  =0.
         scx(1)   =1.
         scx(2)   =1.
         scx(nx-2)=1.
         scx(nx-1)=1.
         scx(nx  )=1.
         sbx(1)   =alsaix
         sbx(2)   =alsaix
         sbx(nx-2)=alsaix
         sbx(nx-1)=alsaix+alsaix
         sbx(nx  )=0.
         do i=3,nx-3
            sfx(i)=alsaix
            scx(i)=1.
            sbx(i)=alsaix
         enddo
      endif
!
      if (nclx.eq.2) then
         sfx(1)   =alsa1x
         sfx(2)   =alsa2x
         sfx(nx-2)=alsaix
         sfx(nx-1)=alsamx
         sfx(nx)  =0.
         scx(1)   =1.
         scx(2)   =1.
         scx(nx-2)=1.
         scx(nx-1)=1.
         scx(nx  )=1.
         sbx(1)   =alsa2x 
         sbx(2)   =alsaix
         sbx(nx-2)=alsamx
         sbx(nx-1)=alsanx
         sbx(nx  )=0.
         do i=3,nx-3
            sfx(i)=alsaix
            scx(i)=1.
            sbx(i)=alsaix
         enddo
      endif
!
      if (ncly.eq.0) then
         sfy(1)   =alsajy
         sfy(2)   =alsajy
         sfy(ny-2)=alsajy
         sfy(ny-1)=alsajy
         sfy(ny)  =0.
         scy(1)   =2.
         scy(2)   =1.
         scy(ny-2)=1.
         scy(ny-1)=1.
         scy(ny  )=1.+alsajy*alsajy
         sby(1)   =alsajy
         sby(2)   =alsajy
         sby(ny-2)=alsajy
         sby(ny-1)=alsajy
         sby(ny  )=0.
         do j=3,ny-3
            sfy(j)=alsajy
            scy(j)=1.
            sby(j)=alsajy
         enddo
      endif
!
      if (ncly.eq.1) then
         sfy(1)   =alsajy+alsajy
         sfy(2)   =alsajy
         sfy(ny-2)=alsajy
         sfy(ny-1)=alsajy
         sfy(ny)  =0.
         scy(1)   =1.
         scy(2)   =1.
         scy(ny-2)=1.
         scy(ny-1)=1.
         scy(ny  )=1.
         sby(1)   =alsajy 
         sby(2)   =alsajy
         sby(ny-2)=alsajy
         sby(ny-1)=alsajy+alsajy
         sby(ny  )=0.
         do j=3,ny-3
            sfy(j)=alsajy
            scy(j)=1.
            sby(j)=alsajy
         enddo
      endif
!
      if (ncly.eq.2) then
         sfy(1)   =alsa1y
         sfy(2)   =alsa2y
         sfy(ny-2)=alsajy
         sfy(ny-1)=alsamy
         sfy(ny)  =0.
         scy(1)   =1.
         scy(2)   =1.
         scy(ny-2)=1.
         scy(ny-1)=1.
         scy(ny  )=1.
         sby(1)   =alsa2y 
         sby(2)   =alsajy
         sby(ny-2)=alsamy
         sby(ny-1)=alsany
         sby(ny  )=0.
         do j=3,ny-3
            sfy(j)=alsajy
            scy(j)=1.
            sby(j)=alsajy
         enddo
      endif
!
      if (nz.gt.1) then
         if (nclz.eq.0) then
            sfz(1)   =alsakz
            sfz(2)   =alsakz
            sfz(nz-2)=alsakz
            sfz(nz-1)=alsakz
            sfz(nz)  =0.
            scz(1)   =2.
            scz(2)   =1.
            scz(nz-2)=1.
            scz(nz-1)=1.
            scz(nz  )=1.+alsakz*alsakz
            sbz(1)   =alsakz
            sbz(2)   =alsakz
            sbz(nz-2)=alsakz
            sbz(nz-1)=alsakz
            sbz(nz  )=0.
            do k=3,nz-3
               sfz(k)=alsakz
               scz(k)=1.
               sbz(k)=alsakz
            enddo
         endif
!
         if (nclz.eq.1) then
            sfz(1)   =alsakz+alsakz
            sfz(2)   =alsakz
            sfz(nz-2)=alsakz
            sfz(nz-1)=alsakz
            sfz(nz)  =0.
            scz(1)   =1.
            scz(2)   =1.
            scz(nz-2)=1.
            scz(nz-1)=1.
            scz(nz  )=1.
            sbz(1)   =alsakz 
            sbz(2)   =alsakz
            sbz(nz-2)=alsakz
            sbz(nz-1)=alsakz+alsakz
            sbz(nz  )=0.
            do k=3,nz-3
               sfz(k)=alsakz
               scz(k)=1.
               sbz(k)=alsakz
            enddo
         endif
!
         if (nclz.eq.2) then
            sfz(1)   =alsa1z
            sfz(2)   =alsa2z
            sfz(nz-2)=alsakz
            sfz(nz-1)=alsamz
            sfz(nz)  =0.
            scz(1)   =1.
            scz(2)   =1.
            scz(nz-2)=1.
            scz(nz-1)=1.
            scz(nz  )=1.
            sbz(1)   =alsa2z 
            sbz(2)   =alsakz
            sbz(nz-2)=alsamz
            sbz(nz-1)=alsanz
            sbz(nz  )=0.
            do k=3,nz-3
               sfz(k)=alsakz
               scz(k)=1.
               sbz(k)=alsakz
            enddo
         endif
      endif
!
      do i=1,nx
         ffxp(i)=ffx(i)
         sfxp(i)=sfx(i)
      enddo

      do j=1,ny
         ffyp(j)=ffy(j)
         sfyp(j)=sfy(j)
      enddo
      do k=1,nz
         ffzp(k)=ffz(k)
         sfzp(k)=sfz(k)
      enddo
!
   do i=1,nxm 
!***************************
      cfxp6(i)=cfx6(i) 
      cifxp6(i)=cifx6(i) 
!***************************
   enddo
   do i=1,nx 
!***************************
      cifip6(i)=cifi6(i) 
      cfip6(i)=cfi6(i) 
!***************************
   enddo
   do j=1,nym 
!***************************
      cfyp6(j)=cfy6(j) 
      cifyp6(j)=cify6(j) 
!***************************
   enddo
   do j=1,ny 
!***************************
      cifip6y(j)=cifi6y(j) 
      cfip6y(j)=cfi6y(j) 
!***************************
   enddo
!
      if (nclx.eq.1) then
         ffxp(1)=0.
         sfx (1)=0.
      endif
!***************************
   cfxp6(1)=0. 
   cfip6(1)=0. 
!***************************
   if (ncly==1) then  
!***************************
      cfyp6(1)=0. 
      cfip6y(1)=0. 
!***************************
   endif
!
      if (ncly.eq.1) then
         ffyp(1)=0.
         sfy (1)=0.
      endif

      if (nclz.eq.1) then
         ffzp(1)=0.
         sfz (1)=0.
      endif
!
!******************STRET STRET STRET STRET STRET*******************
      do j=1,ny-1
         ffy(j)=ffy(j)*ppy(j+1)
         ffyp(j)=ffyp(j)*ppy(j+1)
         fcy(j)=fcy(j)*ppy(j)
         fby(j)=fby(j)*ppy(j)
         sfy(j)=sfy(j)*pp2y(j+1)
         sfyp(j)=sfyp(j)*pp2y(j+1)
         scy(j)=scy(j)*pp2y(j)
         sby(j)=sby(j)*pp2y(j)
      enddo
      fcy(ny)=fcy(ny)*ppy(ny)
      scy(ny)=scy(ny)*pp2y(ny)
!******************************************************************
      call prepare (fbx,fcx,ffx,fsx,fwx,nx)
!*********************************************
   call prepare (cbx6,ccx6,cfx6,csx6,cwx6,nxm) 
   call prepare (cibx6,cicx6,cifx6,cisx6,ciwx6,nxm) 
   call prepare (cbi6,cci6,cfi6,csi6,cwi6,nx) 
   call prepare (cibi6,cici6,cifi6,cisi6,ciwi6,nx) 
!
   call prepare (cby6,ccy6,cfy6,csy6,cwy6,nym) 
   call prepare (ciby6,cicy6,cify6,cisy6,ciwy6,nym) 
   call prepare (cbi6y,cci6y,cfi6y,csi6y,cwi6y,ny) 
   call prepare (cibi6y,cici6y,cifi6y,cisi6y,ciwi6y,ny) 
!*********************************************
      call prepare (fby,fcy,ffy,fsy,fwy,ny)
      if (nz.gt.1) call prepare (fbz,fcz,ffz,fsz,fwz,nz)
!
      call prepare (fbx,fcx,ffxp,fsxp,fwxp,nx)
      call prepare (fby,fcy,ffyp,fsyp,fwyp,ny)
      if (nz.gt.1) call prepare (fbz,fcz,ffzp,fszp,fwzp,nz)
!
      if (nclx.eq.1) then
         fbx(nx-1)=0.
         call prepare (fbx,fcx,ffxp,fsxp,fwxp,nx)
      endif
!
!****************************
   cbx6(nxm-1)=0. 
   cibx6(nxm)=0 
   cbi6(nx-1)=0. 
   cibi6(nx)=0 
   call prepare (cbx6,ccx6,cfxp6,csxp6,cwxp6,nxm) 
   call prepare (cibx6,cicx6,cifxp6,cisxp6,ciwxp6,nxm) 
   call prepare (cbi6,cci6,cfip6,csip6,cwip6,nx) 
   call prepare (cibi6,cici6,cifip6,cisip6,ciwip6,nx)
!****************************
!
      if (ncly.eq.1) then
!****************************
      cby6(nym-1)=0. 
      ciby6(nym)=0 
      cbi6y(ny-1)=0. 
      cibi6y(ny)=0 
!****************************
         fby(ny-1)=0.
         call prepare (fby,fcy,ffyp,fsyp,fwyp,ny)
!************************************************
      call prepare (cby6,ccy6,cfyp6,csyp6,cwyp6,nym) 
      call prepare (ciby6,cicy6,cifyp6,cisyp6,ciwyp6,nym) 
      call prepare (cbi6y,cci6y,cfip6y,csip6y,cwip6y,ny) 
      call prepare (cibi6y,cici6y,cifip6y,cisip6y,ciwip6y,ny) 
!************************************************
      endif
!
      if (nz.gt.1) then
         if (nclz.eq.1) then
            fbz(nz-1)=0.
            call prepare (fbz,fcz,ffzp,fszp,fwzp,nz)
         endif
      endif
!
      call prepare (sbx,scx,sfx,ssx,swx,nx)
      call prepare (sby,scy,sfy,ssy,swy,ny)
      if (nz.gt.1) call prepare (sbz,scz,sfz,ssz,swz,nz)
!
      call prepare (sbx,scx,sfxp,ssxp,swxp,nx)
      call prepare (sby,scy,sfyp,ssyp,swyp,ny)
      if (nz.gt.1) call prepare (sbz,scz,sfzp,sszp,swzp,nz)
!
      if (nclx.eq.1) then
         sbx(nx-1)=0.
         call prepare (sbx,scx,sfx ,ssx ,swx ,nx)
      endif
!
      if (ncly.eq.1) then
         sby(ny-1)=0.
         call prepare (sby,scy,sfy ,ssy ,swy ,ny)
      endif
!
      if (nclz.eq.1) then
         sbz(nz-1)=0.
         call prepare (sbz,scz,sfz ,ssz ,swz ,nz)
      endif
!
      return
      end
!
!*******************************************************************
!
      subroutine prepare (b,c,f,s,w,n)
! 
!*******************************************************************
!
      real(8),dimension(n) :: b,c,f,s,w
!
      do i=1,n
         w(i)=c(i)
      enddo
!
      do i=2,n
         s(i)=b(i-1)/w(i-1)
         w(i)=w(i)-f(i-1)*s(i)
      enddo
!
      do i=1,n
         w(i)=1./w(i)    
      enddo
!
      return
      end
!********************************************************************
!
subroutine derx(tx,ux,rx,sx,ffx,fsx,fwx,nx,ny,nz,npaire) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
  USE paramx_m
  USE derivex_m
!...Translated by PSUITE Trans90                  4.3ZH  9:45:39   1/26/ 4  
!...Switches: -yv INDDO=0 -nbejkno
  implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  integer , intent(in) :: nx 
  integer , intent(in) :: ny 
  integer , intent(in) :: nz 
  integer , intent(in) :: npaire 
  real , intent(inout) :: tx(nx,ny,nz) 
  real(8) , intent(in) :: ux(nx,ny,nz) 
  real(8) , intent(inout) :: rx(nx,ny,nz) 
  real(8) , intent(inout) :: sx(ny,nz) 
  real(8) , intent(in) :: ffx(nx) 
  real(8) , intent(in) :: fsx(nx) 
  real(8) , intent(in) :: fwx(nx) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer :: k, j, i 
!-----------------------------------------------
!
!********************************************************************
!
   if (nclx==0) then 
      do k=1,nz 
      do j=1,ny 
         tx(1,j,k)=afix*(ux(2,j,k)-ux(nx,j,k))& 
                  +bfix*(ux(3,j,k)-ux(nx-1,j,k)) 
         rx(1,j,k)=-1. 
         tx(2,j,k)=afix*(ux(3,j,k)-ux(1,j,k))&
                  +bfix*(ux(4,j,k)-ux(nx,j,k)) 
         rx(2,j,k)=0. 
      enddo
      enddo
      do i=3,nx-2 
      do k=1,nz 
      do j=1,ny 
         tx(i,j,k)=afix*(ux(i+1,j,k)-ux(i-1,j,k))&
                  +bfix*(ux(i+2,j,k)-ux(i-2,j,k)) 
          rx(i,j,k)=0. 
       enddo
       enddo
       enddo
       do k=1,nz 
       do j=1,ny 
          tx(nx-1,j,k)=afix*(ux(nx,j,k)-ux(nx-2,j,k))&
                      +bfix*(ux(1,j,k)-ux(nx-3,j,k)) 
          rx(nx-1,j,k)=0. 
          tx(nx,j,k)=afix*(ux(1,j,k)-ux(nx-1,j,k))&
                    +bfix*(ux(2,j,k)-ux(nx-2,j,k)) 
          rx(nx,j,k)=alfaix           
       enddo
       enddo 
       do i=2, nx 
       do k=1, nz 
       do j=1, ny 
          tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*fsx(i) 
          rx(i,j,k)=rx(i,j,k)-rx(i-1,j,k)*fsx(i) 
       enddo 
       enddo 
       enddo 
       do k=1,nz 
       do j=1,ny 
          tx(nx,j,k)=tx(nx,j,k)*fwx(nx) 
          rx(nx,j,k)=rx(nx,j,k)*fwx(nx) 
       enddo
       enddo 
       do i=nx-1,1,-1 
       do k=1,nz 
       do j=1,ny 
          tx(i,j,k)=(tx(i,j,k)-ffx(i)*tx(i+1,j,k))*fwx(i) 
          rx(i,j,k)=(rx(i,j,k)-ffx(i)*rx(i+1,j,k))*fwx(i) 
       enddo
       enddo
       enddo
       do k=1,nz 
       do j=1,ny 
          sx(j,k)=(tx(1,j,k)-alfaix*tx(nx,j,k))&
                 /(1.+rx(1,j,k)-alfaix*rx(nx,j,k)) 
       enddo
       enddo
       do k=1,nz 
       do j=1,ny 
       do i=1,nx 
          tx(i,j,k)=tx(i,j,k)-sx(j,k)*rx(i,j,k) 
       enddo 
       enddo 
       enddo 
   endif
!
   if (nclx==1) then 
      if (npaire==1) then 
         do k=1,nz 
         do j=1,ny 
            tx(1,j,k)=0. 
            tx(2,j,k)=afix*(ux(3,j,k)-ux(1,j,k))&
                     +bfix*(ux(4,j,k)-ux(2,j,k)) 
         enddo
         enddo
         do i=3,nx-2 
         do k=1,nz 
         do j=1,ny 
            tx(i,j,k)=afix*(ux(i+1,j,k)-ux(i-1,j,k))&
                     +bfix*(ux(i+2,j,k)-ux(i-2,j,k)) 
         enddo
         enddo
         enddo
         do k=1,nz 
         do j=1,ny 
            tx(nx-1,j,k)=afix*(ux(nx,j,k)-ux(nx-2,j,k))&
                        +bfix*(ux(nx-1,j,k)-ux(nx-3,j,k)) 
            tx(nx,j,k)=0. 
         enddo 
         enddo 
         do i=2,nx 
         do k=1,nz 
         do j=1,ny 
            tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*fsx(i) 
         enddo
         enddo
         enddo
         do k=1,nz 
         do j=1,ny 
            tx(nx,j,k)=tx(nx,j,k)*fwx(nx) 
         enddo
         enddo
         do i=nx-1,1,-1 
         do k=1,nz 
         do j=1,ny 
            tx(i,j,k)=(tx(i,j,k)-ffx(i)*tx(i+1,j,k))*fwx(i) 
         enddo
         enddo
         enddo
      endif
      if (npaire==0) then 
         do k=1,nz 
         do j=1,ny 
            tx(1,j,k)=afix*(ux(2,j,k)+ux(2,j,k))&
                     +bfix*(ux(3,j,k)+ux(3,j,k)) 
            tx(2,j,k)=afix*(ux(3,j,k)-ux(1,j,k))&
                     +bfix*(ux(4,j,k)+ux(2,j,k)) 
         enddo
         enddo
         do i=3,nx-2 
         do k=1,nz 
         do j=1,ny 
            tx(i,j,k)=afix*(ux(i+1,j,k)-ux(i-1,j,k))&
                     +bfix*(ux(i+2,j,k)-ux(i-2,j,k)) 
         enddo
         enddo
         enddo
         do k=1,nz 
         do j=1,ny 
            tx(nx-1,j,k)=afix*(ux(nx,j,k)-ux(nx-2,j,k))&
                        +bfix*((-ux(nx-1,j,k))-ux(nx-3,j,k)) 
            tx(nx,j,k)=afix*((-ux(nx-1,j,k))-ux(nx-1,j,k))&
                      +bfix*((-ux(nx-2,j,k))-ux(nx-2,j,k)) 
         enddo
         enddo
         do i=2,nx 
         do k=1,nz 
         do j=1,ny 
            tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*fsx(i) 
         enddo
         enddo
         enddo
         do k=1,nz 
         do j=1,ny 
            tx(nx,j,k)=tx(nx,j,k)*fwx(nx) 
         enddo
         enddo
         do i=nx-1,1,-1 
         do k=1,nz 
         do j=1,ny 
            tx(i,j,k)=(tx(i,j,k)-ffx(i)*tx(i+1,j,k))*fwx(i) 
         enddo
         enddo
         enddo
      endif
   endif
!
   if (nclx==2) then 
      do k=1,nz 
      do j=1,ny 
         tx(1,j,k)=af1x*ux(1,j,k)+bf1x*ux(2,j,k)+cf1x*ux(3,j,k) 
         tx(2,j,k)=af2x*(ux(3,j,k)-ux(1,j,k)) 
      enddo 
      enddo 
      do i=3,nx-2 
      do k=1,nz 
      do j=1,ny 
         tx(i,j,k)=afix*(ux(i+1,j,k)-ux(i-1,j,k))&
                  +bfix*(ux(i+2,j,k)-ux(i-2,j,k)) 
      enddo 
      enddo 
      enddo 
      do k=1,nz 
      do j=1,ny 
         tx(nx-1,j,k)=afmx*(ux(nx,j,k)-ux(nx-2,j,k)) 
         tx(nx,j,k)=(-afnx*ux(nx,j,k))-bfnx*ux(nx-1,j,k)-cfnx*ux(nx-2,j,k) 
      enddo
      enddo
      do i=2,nx 
      do k=1,nz 
      do j=1,ny 
         tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*fsx(i) 
      enddo
      enddo
      enddo
      do k=1,nz 
      do j=1,ny 
         tx(nx,j,k)=tx(nx,j,k)*fwx(nx) 
      enddo 
      enddo 
      do i=nx-1,1,-1 
      do k=1,nz 
      do j=1,ny 
         tx(i,j,k)=(tx(i,j,k)-ffx(i)*tx(i+1,j,k))*fwx(i) 
      enddo 
      enddo 
      enddo 
   endif 
!
end subroutine derx 
!
!********************************************************************
!
subroutine dery(ty,uy,ry,di,sy,ffy,fsy,fwy,nx,ny,nz,npaire) 
!...Translated by PSUITE Trans90                  4.3ZH  9:45:39   1/26/ 4  
!...Switches: -yv INDDO=0 -nbejkno
  USE paramy_m

  implicit none


!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  integer  :: nx 
  integer  :: ny 
  integer  :: nz 
  integer  :: npaire 
  real(8)  :: ty(nx,ny,nz) 
  real(8) , intent(in) :: uy(nx,ny,nz) 
  real(8)  :: ry(nx,nz,ny) 
  real(8)  :: di(nx,nz,ny) 
  real(8)  :: sy(nx,nz) 
  real(8)  :: ffy(ny) 
  real(8)  :: fsy(ny) 
  real(8)  :: fwy(ny) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer :: k, j, i 
!-----------------------------------------------
!
!********************************************************************
!
   do k=1,nz 
   do j=1,ny 
   do i=1,nx 
      di(i,k,j)=uy(i,j,k) 
   enddo
   enddo 
   enddo 
!
   call dery1 (ry,di,ty,sy,ffy,fsy,fwy,nx,ny,nz,npaire) 
!
   do k=1,nz 
   do j=1,ny 
   do i=1,nx 
       ty(i,j,k)=ry(i,k,j)*yly 
   enddo 
   enddo
   enddo 
!
   return  
end subroutine dery
!
!********************************************************************
!
subroutine dery1(ty,uy,ry,sy,ffy,fsy,fwy,nx,ny,nz,npaire) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
  USE paramy_m 
  USE derivey_m 
!...Translated by PSUITE Trans90                  4.3ZH  9:45:39   1/26/ 4  
!...Switches: -yv INDDO=0 -nbejkno
  implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  integer , intent(in) :: nx 
  integer , intent(in) :: ny 
  integer , intent(in) :: nz 
  integer , intent(in) :: npaire 
  real(8) , intent(inout) :: ty(nx,nz,ny) 
  real(8) , intent(in) :: uy(nx,nz,ny) 
  real(8) , intent(inout) :: ry(nx,nz,ny) 
  real(8) , intent(inout) :: sy(nx,nz) 
  real(8) , intent(in) :: ffy(ny) 
  real(8) , intent(in) :: fsy(ny) 
  real(8) , intent(in) :: fwy(ny) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer :: k, i, j 
!-----------------------------------------------
!
!********************************************************************
!
!      include 'precision.inc'
!
!
!
   if (ncly==0) then 
      do k=1,nz 
      do i=1,nx 
         ty(i,k,1)=afjy*(uy(i,k,2)-uy(i,k,ny))&
                  +bfjy*(uy(i,k,3)-uy(i,k,ny-1)) 
         ry(i,k,1)=-1. 
         ty(i,k,2)=afjy*(uy(i,k,3)-uy(i,k,1))&
                  +bfjy*(uy(i,k,4)-uy(i,k,ny)) 
         ry(i,k,2)=0. 
      enddo
      enddo 
      do j=3,ny-2 
      do k=1,nz 
      do i=1,nx 
         ty(i,k,j)=afjy*(uy(i,k,j+1)-uy(i,k,j-1))&
                  +bfjy*(uy(i,k,j+2)-uy(i,k,j-2)) 
         ry(i,k,j)=0. 
      enddo
      enddo 
      enddo
      do k=1,nz 
      do i=1,nx 
         ty(i,k,ny-1)=afjy*(uy(i,k,ny)-uy(i,k,ny-2))&
                     +bfjy*(uy(i,k,1)-uy(i,k,ny-3)) 
         ry(i,k,ny-1)=0. 
         ty(i,k,ny)=afjy*(uy(i,k,1)-uy(i,k,ny-1))&
                   +bfjy*(uy(i,k,2)-uy(i,k,ny-2)) 
         ry(i,k,ny)=alfajy 
      enddo 
      enddo 
      do j=2,ny 
      do k=1,nz 
      do i=1,nx 
         ty(i,k,j)=ty(i,k,j)-ty(i,k,j-1)*fsy(j) 
         ry(i,k,j)=ry(i,k,j)-ry(i,k,j-1)*fsy(j) 
      enddo 
      enddo 
      enddo
      do k=1,nz 
      do i=1,nx 
         ty(i,k,ny)=ty(i,k,ny)*fwy(ny) 
         ry(i,k,ny)=ry(i,k,ny)*fwy(ny) 
      enddo 
      enddo 
      do j=ny-1,1,-1 
      do k=1,nz 
      do i=1,nx 
         ty(i,k,j)=(ty(i,k,j)-ffy(j)*ty(i,k,j+1))*fwy(j) 
         ry(i,k,j)=(ry(i,k,j)-ffy(j)*ry(i,k,j+1))*fwy(j) 
      enddo 
      enddo 
      enddo
      do k=1,nz 
      do i=1,nx 
         sy(i,k)=(ty(i,k,1)-alfajy*ty(i,k,ny))&
                /(1.+ry(i,k,1)-alfajy*ry(i,k,ny)) 
      enddo 
      enddo 
      do j=1,ny 
      do k=1,nz 
      do i=1,nx 
         ty(i,k,j)=ty(i,k,j)-sy(i,k)*ry(i,k,j) 
      enddo 
      enddo
      enddo
   endif
!
   if (ncly==1) then 
      if (npaire==1) then 
         do k=1,nz 
         do i=1,nx 
            ty(i,k,1)=0. 
            ty(i,k,2)=afjy*(uy(i,k,3)-uy(i,k,1))&
                     +bfjy*(uy(i,k,4)-uy(i,k,2)) 
         enddo 
         enddo 
         do k=1,nz 
         do i=1,nx 
         do j=3,ny-2 
            ty(i,k,j)=afjy*(uy(i,k,j+1)-uy(i,k,j-1))&
                     +bfjy*(uy(i,k,j+2)-uy(i,k,j-2)) 
         enddo 
         enddo 
         enddo 
         do k=1,nz 
         do i=1,nx 
            ty(i,k,ny-1)=afjy*(uy(i,k,ny)-uy(i,k,ny-2))&
                        +bfjy*(uy(i,k,ny-1)-uy(i,k,ny-3)) 
            ty(i,k,ny)=0. 
         enddo 
         enddo 
         do j=2,ny 
         do k=1,nz 
         do i=1,nx 
            ty(i,k,j)=ty(i,k,j)-ty(i,k,j-1)*fsy(j) 
         enddo 
         enddo 
         enddo 
         do k=1,nz 
         do i=1,nx 
            ty(i,k,ny)=ty(i,k,ny)*fwy(ny) 
         enddo 
         enddo 
         do j=ny-1,1,-1 
         do k=1,nz 
         do i=1,nx 
            ty(i,k,j)=(ty(i,k,j)-ffy(j)*ty(i,k,j+1))*fwy(j) 
         enddo 
         enddo 
         enddo 
      endif
      if (npaire==0) then 
         do k=1,nz 
         do i=1,nx 
            ty(i,k,1)=afjy*(uy(i,k,2)+uy(i,k,2))&
                     +bfjy*(uy(i,k,3)+uy(i,k,3)) 
            ty(i,k,2)=afjy*(uy(i,k,3)-uy(i,k,1))&
                     +bfjy*(uy(i,k,4)+uy(i,k,2)) 
         enddo 
         enddo 
         do k=1,nz 
         do i=1,nx 
         do j=3,ny-2 
            ty(i,k,j)=afjy*(uy(i,k,j+1)-uy(i,k,j-1))&
                     +bfjy*(uy(i,k,j+2)-uy(i,k,j-2)) 
         enddo 
         enddo 
         enddo 
         do k=1,nz 
         do i=1,nx 
            ty(i,k,ny-1)=afjy*(uy(i,k,ny)-uy(i,k,ny-2))&
                        +bfjy*((-uy(i,k,ny-1))-uy(i,k,ny-3)) 
            ty(i,k,ny)=afjy*((-uy(i,k,ny-1))-uy(i,k,ny-1))&
                      +bfjy*((-uy(i,k,ny-2))-uy(i,k,ny-2)) 
         enddo 
         enddo 
         do j=2,ny 
         do k=1,nz 
         do i=1,nx 
            ty(i,k,j)=ty(i,k,j)-ty(i,k,j-1)*fsy(j) 
         enddo 
         enddo 
         enddo 
         do k=1,nz 
         do i=1,nx 
            ty(i,k,ny)=ty(i,k,ny)*fwy(ny) 
         enddo
         enddo 
         do j=ny-1,1,-1 
         do k=1,nz 
         do i=1,nx 
            ty(i,k,j)=(ty(i,k,j)-ffy(j)*ty(i,k,j+1))*fwy(j) 
         enddo 
         enddo 
         enddo 
      endif
   endif
!
   if (ncly==2) then 
      do k=1,nz 
      do i=1,nx 
         ty(i,k,1)=af1y*uy(i,k,1)+bf1y*uy(i,k,2)+cf1y*uy(i,k,3) 
         ty(i,k,2)=af2y*(uy(i,k,3)-uy(i,k,1)) 
      enddo 
      enddo 
      do j=3,ny-2 
      do k=1,nz 
      do i=1,nx 
         ty(i,k,j)=afjy*(uy(i,k,j+1)-uy(i,k,j-1))&
                  +bfjy*(uy(i,k,j+2)-uy(i,k,j-2)) 
      enddo 
      enddo 
      enddo 
      do k=1,nz 
      do i=1,nx 
         ty(i,k,ny-1)=afmy*(uy(i,k,ny)-uy(i,k,ny-2)) 
         ty(i,k,ny)=(-afny*uy(i,k,ny))-bfny*uy(i,k,ny-1)-cfny*uy(i,k,ny-2) 
      enddo 
      enddo 
      do j=2,ny 
      do k=1,nz 
      do i=1,nx 
         ty(i,k,j)=ty(i,k,j)-ty(i,k,j-1)*fsy(j) 
      enddo 
      enddo 
      enddo 
      do k=1,nz 
      do i=1,nx 
         ty(i,k,ny)=ty(i,k,ny)*fwy(ny) 
      enddo 
      enddo 
      do j=ny-1,1,-1 
      do k=1,nz 
      do i=1,nx 
         ty(i,k,j)=(ty(i,k,j)-ffy(j)*ty(i,k,j+1))*fwy(j) 
      enddo 
      enddo 
      enddo
   endif
!
end subroutine dery1
!
!********************************************************************
!
subroutine derz(tz,uz,rz,sz,ffz,fsz,fwz,nx,ny,nz,npaire) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
  USE paramz_m 
  USE derivez_m 
!...Translated by PSUITE Trans90                  4.3ZH  9:45:39   1/26/ 4  
!...Switches: -yv INDDO=0 -nbejkno
  implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  integer , intent(in) :: nx 
  integer , intent(in) :: ny 
  integer , intent(in) :: nz 
  integer , intent(in) :: npaire 
  real(8) , intent(inout) :: tz(nx,ny,nz) 
  real(8) , intent(in) :: uz(nx,ny,nz) 
  real(8) , intent(inout) :: rz(nx,ny,nz) 
  real(8) , intent(inout) :: sz(nx,ny) 
  real(8) , intent(in) :: ffz(nz) 
  real(8) , intent(in) :: fsz(nz) 
  real(8) , intent(in) :: fwz(nz) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer :: j, i, k 
!-----------------------------------------------
!
!********************************************************************
!
!
   if (nclz==0) then
      do j=1,ny
      do i=1,nx
         tz(i,j,1)=afkz*(uz(i,j,2)-uz(i,j,nz  ))&
                  +bfkz*(uz(i,j,3)-uz(i,j,nz-1))
         rz(i,j,1)=-1.
         tz(i,j,2)=afkz*(uz(i,j,3)-uz(i,j,1 ))&
                  +bfkz*(uz(i,j,4)-uz(i,j,nz))
         rz(i,j,2)=0.
      enddo
      enddo
      do k=3,nz-2
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=afkz*(uz(i,j,k+1)-uz(i,j,k-1))&
                  +bfkz*(uz(i,j,k+2)-uz(i,j,k-2))
         rz(i,j,k)=0.
      enddo
      enddo
      enddo
      do j=1,ny
      do i=1,nx
         tz(i,j,nz-1)=afkz*(uz(i,j,nz)-uz(i,j,nz-2))&
                     +bfkz*(uz(i,j,1 )-uz(i,j,nz-3))
         rz(i,j,nz-1)=0.
         tz(i,j,nz  )=afkz*(uz(i,j,1)-uz(i,j,nz-1))&
                     +bfkz*(uz(i,j,2)-uz(i,j,nz-2))
         rz(i,j,nz  )=alfakz
      enddo
      enddo
      do k=2,nz
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=tz(i,j,k)-tz(i,j,k-1)*fsz(k)
         rz(i,j,k)=rz(i,j,k)-rz(i,j,k-1)*fsz(k)
      enddo
      enddo
      enddo
      do j=1,ny
      do i=1,nx
         tz(i,j,nz)=tz(i,j,nz)*fwz(nz)
         rz(i,j,nz)=rz(i,j,nz)*fwz(nz)
      enddo
      enddo
      do k=nz-1,1,-1
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=(tz(i,j,k)-ffz(k)*tz(i,j,k+1))*fwz(k)
         rz(i,j,k)=(rz(i,j,k)-ffz(k)*rz(i,j,k+1))*fwz(k)
      enddo
      enddo
      enddo
      do j=1,ny
      do i=1,nx
         sz(i,j)=(   tz(i,j,1)-alfakz*tz(i,j,nz))/&
                 (1.+rz(i,j,1)-alfakz*rz(i,j,nz))
      enddo
      enddo
      do k=1,nz
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=tz(i,j,k)-sz(i,j)*rz(i,j,k)
      enddo
      enddo
      enddo
   endif
!
   if (nclz==1) then
      if (npaire==1) then
         do j=1,ny
         do i=1,nx
            tz(i,j,1)=0.
            tz(i,j,2)=afkz*(uz(i,j,3)-uz(i,j,1))&
                     +bfkz*(uz(i,j,4)-uz(i,j,2))
         enddo
         enddo
         do k=3,nz-2
         do j=1,ny
         do i=1,nx
               tz(i,j,k)=afkz*(uz(i,j,k+1)-uz(i,j,k-1))&
                        +bfkz*(uz(i,j,k+2)-uz(i,j,k-2))
         enddo
         enddo
         enddo
         do j=1,ny
         do i=1,nx
            tz(i,j,nz-1)=afkz*(uz(i,j,nz  )-uz(i,j,nz-2))&
                        +bfkz*(uz(i,j,nz-1)-uz(i,j,nz-3))
            tz(i,j,nz  )=0.
         enddo
         enddo
         do k=2,nz
         do j=1,ny
         do i=1,nx
            tz(i,j,k)=tz(i,j,k)-tz(i,j,k-1)*fsz(k)
         enddo
         enddo
         enddo
         do j=1,ny
         do i=1,nx
            tz(i,j,nz)=tz(i,j,nz)*fwz(nz)
         enddo
         enddo
         do k=nz-1,1,-1
         do j=1,ny
         do i=1,nx
            tz(i,j,k)=(tz(i,j,k)-ffz(k)*tz(i,j,k+1))*fwz(k)
         enddo
         enddo
         enddo
      endif
      if (npaire==0) then
         do j=1,ny
         do i=1,nx
            tz(i,j,1)=afkz*(uz(i,j,2)+uz(i,j,2))&
                     +bfkz*(uz(i,j,3)+uz(i,j,3))
            tz(i,j,2)=afkz*(uz(i,j,3)-uz(i,j,1))&
                     +bfkz*(uz(i,j,4)+uz(i,j,2))
         enddo
         enddo
         do k=3,nz-2
         do j=1,ny
         do i=1,nx
            tz(i,j,k)=afkz*(uz(i,j,k+1)-uz(i,j,k-1))&
                     +bfkz*(uz(i,j,k+2)-uz(i,j,k-2))
         enddo
         enddo
         enddo
         do j=1,ny
         do i=1,nx
            tz(i,j,nz-1)=afkz*( uz(i,j,nz  )-uz(i,j,nz-2))&
                        +bfkz*(-uz(i,j,nz-1)-uz(i,j,nz-3))
            tz(i,j,nz  )=afkz*(-uz(i,j,nz-1)-uz(i,j,nz-1))&
                        +bfkz*(-uz(i,j,nz-2)-uz(i,j,nz-2))
         enddo
         enddo
         do k=2,nz
         do j=1,ny
         do i=1,nx
            tz(i,j,k)=tz(i,j,k)-tz(i,j,k-1)*fsz(k)
         enddo
         enddo
         enddo
         do j=1,ny
         do i=1,nx
            tz(i,j,nz)=tz(i,j,nz)*fwz(nz)
         enddo
         enddo
         do k=nz-1,1,-1
         do j=1,ny
         do i=1,nx
            tz(i,j,k)=(tz(i,j,k)-ffz(k)*tz(i,j,k+1))*fwz(k)
         enddo
         enddo
         enddo
      endif
   endif
!
   if (nclz==2) then
      do j=1,ny
      do i=1,nx
         tz(i,j,1)=af1z*uz(i,j,1)+bf1z*uz(i,j,2)&
                  +cf1z*uz(i,j,3)
         tz(i,j,2)=af2z*(uz(i,j,3)-uz(i,j,1))
      enddo
      enddo
      do k=3,nz-2
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=afkz*(uz(i,j,k+1)-uz(i,j,k-1))&
                  +bfkz*(uz(i,j,k+2)-uz(i,j,k-2))
      enddo
      enddo
      enddo
      do j=1,ny
      do i=1,nx
         tz(i,j,nz-1)= afmz*(uz(i,j,nz)-uz(i,j,nz-2))
         tz(i,j,nz  )=-afnz*uz(i,j,nz)-bfnz*uz(i,j,nz-1)&
                      -cfnz*uz(i,j,nz-2)
      enddo
      enddo
      do k=2,nz
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=tz(i,j,k)-tz(i,j,k-1)*fsz(k)
      enddo
      enddo
      enddo
      do j=1,ny
      do i=1,nx
         tz(i,j,nz)=tz(i,j,nz)*fwz(nz)
      enddo
      enddo
      do k=nz-1,1,-1
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=(tz(i,j,k)-ffz(k)*tz(i,j,k+1))*fwz(k)
      enddo
      enddo
      enddo
   endif
!
   return
end subroutine derz


!
!********************************************************************
!
subroutine derxx(tx,ux,rx,sx,sfx,ssx,swx,nx,ny,nz,npaire) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
  USE paramx_m 
  USE derivex_m 
!...Translated by PSUITE Trans90                  4.3ZH  9:45:39   1/26/ 4  
!...Switches: -yv INDDO=0 -nbejkno
  implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  integer , intent(in) :: nx 
  integer , intent(in) :: ny 
  integer , intent(in) :: nz 
  integer , intent(in) :: npaire 
  real(8) , intent(inout) :: tx(nx,ny,nz) 
  real(8) , intent(in) :: ux(nx,ny,nz) 
  real(8) , intent(inout) :: rx(nx,ny,nz) 
  real(8) , intent(inout) :: sx(ny,nz) 
  real(8) , intent(in) :: sfx(nx) 
  real(8) , intent(in) :: ssx(nx) 
  real(8) , intent(in) :: swx(nx) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer :: k, j, i 
!-----------------------------------------------
!
!********************************************************************
!
   if (nclx==0) then
      do k=1,nz
      do j=1,ny
         tx(1,j,k)=asix*(ux(2,j,k)-ux(1   ,j,k)&
                        -ux(1,j,k)+ux(nx  ,j,k))&
                  +bsix*(ux(3,j,k)-ux(1   ,j,k)&
                        -ux(1,j,k)+ux(nx-1,j,k))
         rx(1,j,k)=-1.
         tx(2,j,k)=asix*(ux(3,j,k)-ux(2 ,j,k)&
                        -ux(2,j,k)+ux(1 ,j,k))&
                  +bsix*(ux(4,j,k)-ux(2 ,j,k)&
                        -ux(2,j,k)+ux(nx,j,k))
         rx(2,j,k)=0.
      enddo
      enddo
      do i=3,nx-2
      do k=1,nz
      do j=1,ny
         tx(i,j,k)=asix*(ux(i+1,j,k)-ux(i  ,j,k)&
                        -ux(i  ,j,k)+ux(i-1,j,k))&
                  +bsix*(ux(i+2,j,k)-ux(i  ,j,k)&
                        -ux(i  ,j,k)+ux(i-2,j,k))
         rx(i,j,k)=0.
      enddo
      enddo
      enddo
      do k=1,nz
      do j=1,ny
         tx(nx-1,j,k)=asix*(ux(nx  ,j,k)-ux(nx-1,j,k)&
                           -ux(nx-1,j,k)+ux(nx-2,j,k))&
                     +bsix*(ux(1   ,j,k)-ux(nx-1,j,k)&
                           -ux(nx-1,j,k)+ux(nx-3,j,k))
         rx(nx-1,j,k)=0.
         tx(nx  ,j,k)=asix*(ux(1 ,j,k)-ux(nx  ,j,k)&
                            -ux(nx,j,k)+ux(nx-1,j,k))&
                     +bsix*(ux(2 ,j,k)-ux(nx  ,j,k)&
                           -ux(nx,j,k)+ux(nx-2,j,k))
         rx(nx  ,j,k)=alsaix
      enddo
      enddo
      do i=2,nx
      do k=1,nz
      do j=1,ny
         tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*ssx(i)
         rx(i,j,k)=rx(i,j,k)-rx(i-1,j,k)*ssx(i)
      enddo
      enddo
      enddo
      do k=1,nz
      do j=1,ny
         tx(nx,j,k)=tx(nx,j,k)*swx(nx)
         rx(nx,j,k)=rx(nx,j,k)*swx(nx)
      enddo
      enddo
      do i=nx-1,1,-1
      do k=1,nz
      do j=1,ny
         tx(i,j,k)=(tx(i,j,k)-sfx(i)*tx(i+1,j,k))*swx(i)
         rx(i,j,k)=(rx(i,j,k)-sfx(i)*rx(i+1,j,k))*swx(i)
      enddo
      enddo
      enddo
      do k=1,nz
      do j=1,ny
         sx(j,k)=(   tx(1,j,k)-alsaix*tx(nx,j,k))/&
                 (1.+rx(1,j,k)-alsaix*rx(nx,j,k))
      enddo
      enddo
      do k=1,nz
      do j=1,ny
      do i=1,nx
         tx(i,j,k)=tx(i,j,k)-sx(j,k)*rx(i,j,k)
      enddo
      enddo
      enddo
   endif
!
   if (nclx==1) then
      if (npaire==1) then
         do k=1,nz
         do j=1,ny
            tx(1,j,k)=asix*(ux(2,j,k)-ux(1,j,k)&
                           -ux(1,j,k)+ux(2,j,k))&
                     +bsix*(ux(3,j,k)-ux(1,j,k)&
                           -ux(1,j,k)+ux(3,j,k))
            tx(2,j,k)=asix*(ux(3,j,k)-ux(2,j,k)&
                           -ux(2,j,k)+ux(1,j,k))&
                     +bsix*(ux(4,j,k)-ux(2,j,k)&
                           -ux(2,j,k)+ux(2,j,k))
         enddo
         enddo
         do i=3,nx-2
         do k=1,nz
         do j=1,ny
            tx(i,j,k)=asix*(ux(i+1,j,k)-ux(i  ,j,k)&
                           -ux(i  ,j,k)+ux(i-1,j,k))&
                     +bsix*(ux(i+2,j,k)-ux(i  ,j,k)&
                           -ux(i  ,j,k)+ux(i-2,j,k))
         enddo
         enddo
         enddo
         do k=1,nz
         do j=1,ny
            tx(nx-1,j,k)=asix*(ux(nx  ,j,k)-ux(nx-1,j,k)&
                              -ux(nx-1,j,k)+ux(nx-2,j,k))&
                        +bsix*(ux(nx-1,j,k)-ux(nx-1,j,k)&
                              -ux(nx-1,j,k)+ux(nx-3,j,k))
            tx(nx  ,j,k)=asix*(ux(nx-1,j,k)-ux(nx  ,j,k)&
                               -ux(nx  ,j,k)+ux(nx-1,j,k))&
                         +bsix*(ux(nx-2,j,k)-ux(nx  ,j,k)&
                               -ux(nx  ,j,k)+ux(nx-2,j,k))
         enddo
         enddo
         do i=2,nx
         do k=1,nz
         do j=1,ny
            tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*ssx(i)
         enddo
         enddo
         enddo
         do k=1,nz
         do j=1,ny
            tx(nx,j,k)=tx(nx,j,k)*swx(nx)
         enddo
         enddo
         do i=nx-1,1,-1
         do k=1,nz
         do j=1,ny
            tx(i,j,k)=(tx(i,j,k)-sfx(i)*tx(i+1,j,k))*swx(i)
         enddo
         enddo
         enddo
      endif
      if (npaire==0) then
         do k=1,nz
         do j=1,ny
            tx(1,j,k)=0.
            tx(2,j,k)=asix*(ux(3,j,k)-ux(2,j,k)&
                           -ux(2,j,k)+ux(1,j,k))&
                     +bsix*(ux(4,j,k)-ux(2,j,k)&
                           -ux(2,j,k)-ux(2,j,k))
         enddo
         enddo
         do i=3,nx-2
         do k=1,nz
         do j=1,ny
            tx(i,j,k)=asix*(ux(i+1,j,k)-ux(i  ,j,k)&
                           -ux(i  ,j,k)+ux(i-1,j,k))&
                     +bsix*(ux(i+2,j,k)-ux(i  ,j,k)&
                           -ux(i  ,j,k)+ux(i-2,j,k))
         enddo
         enddo
         enddo
         do k=1,nz
         do j=1,ny
            tx(nx-1,j,k)=asix*(ux(nx  ,j,k)-ux(nx-1,j,k)&
                              -ux(nx-1,j,k)+ux(nx-2,j,k))&
                        +bsix*(-ux(nx-1,j,k)-ux(nx-1,j,k)&
                               -ux(nx-1,j,k)+ux(nx-3,j,k))
            tx(nx  ,j,k)=0.
         enddo
         enddo
         do i=2,nx
         do k=1,nz
         do j=1,ny
            tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*ssx(i)
         enddo
         enddo
         enddo
         do k=1,nz
         do j=1,ny
            tx(nx,j,k)=tx(nx,j,k)*swx(nx)
         enddo
         enddo
         do i=nx-1,1,-1
         do k=1,nz
         do j=1,ny
            tx(i,j,k)=(tx(i,j,k)-sfx(i)*tx(i+1,j,k))*swx(i)
         enddo
         enddo
         enddo
      endif
   endif
!
   if (nclx==2) then
      do k=1,nz
      do j=1,ny
         tx(1,j,k)=as1x*ux(1,j,k)+bs1x*ux(2,j,k)&
                  +cs1x*ux(3,j,k)+ds1x*ux(4,j,k)
         tx(2,j,k)=as2x*(ux(3,j,k)-ux(2,j,k)&
                        -ux(2,j,k)+ux(1,j,k))
      enddo
      enddo
!
      do i=3,nx-2
      do k=1,nz
      do j=1,ny
         tx(i,j,k)=asix*(ux(i+1,j,k)+ux(i-1,j,k))&
                  +bsix*(ux(i+2,j,k)+ux(i-2,j,k))&
                  -(asix*2.+bsix*2.)*ux(i,j,k)
      enddo
      enddo
      enddo
!
      do k=1,nz
      do j=1,ny
         tx(nx-1,j,k)=asmx*(ux(nx  ,j,k)-ux(nx-1,j,k)&
                           -ux(nx-1,j,k)+ux(nx-2,j,k))
         tx(nx  ,j,k)=asnx*ux(nx  ,j,k)+bsnx*ux(nx-1,j,k)&
                     +csnx*ux(nx-2,j,k)+dsnx*ux(nx-3,j,k)
      enddo
      enddo
      do i=2,nx
      do k=1,nz
      do j=1,ny
         tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*ssx(i)
      enddo
      enddo
      enddo
      do k=1,nz
      do j=1,ny
         tx(nx,j,k)=tx(nx,j,k)*swx(nx)
      enddo
      enddo
      do i=nx-1,1,-1
      do k=1,nz
      do j=1,ny
         tx(i,j,k)=(tx(i,j,k)-sfx(i)*tx(i+1,j,k))*swx(i)
      enddo
      enddo
      enddo
   endif
!
   return  
end subroutine derxx
!
!********************************************************************
!
subroutine deryy(ty,uy,ry,di,sy,sfy,ssy,swy,nx,ny,nz,npaire) 
!...Translated by PSUITE Trans90                  4.3ZH  9:45:39   1/26/ 4  
!...Switches: -yv INDDO=0 -nbejkno
  implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  integer  :: nx 
  integer  :: ny 
  integer  :: nz 
  integer  :: npaire 
  real(8)  :: ty(nx,ny,nz) 
  real(8) , intent(in) :: uy(nx,ny,nz) 
  real(8)  :: ry(nx,nz,ny) 
  real(8)  :: di(nx,nz,ny) 
  real(8)  :: sy(nx,nz) 
  real(8)  :: sfy(ny) 
  real(8)  :: ssy(ny) 
  real(8)  :: swy(ny) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer :: k, j, i 
!-----------------------------------------------
!
!********************************************************************
!
   do k=1,nz 
   do j=1,ny 
   do i=1,nx 
      di(i,k,j)=uy(i,j,k) 
   enddo 
   enddo 
   enddo 
!
   call deryy1 (ry,di,ty,sy,sfy,ssy,swy,nx,ny,nz,npaire) 
!
   do k=1,nz 
   do j=1,ny 
   do i=1,nx 
      ty(i,j,k)=ry(i,k,j) 
   enddo 
   enddo 
   enddo 
!
   return  
end subroutine deryy


!
!********************************************************************
!
subroutine deryy1(ty,uy,ry,sy,sfy,ssy,swy,nx,ny,nz,npaire) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
  USE paramy_m 
  USE derivey_m 
!...Translated by PSUITE Trans90                  4.3ZH  9:45:39   1/26/ 4  
!...Switches: -yv INDDO=0 -nbejkno
  implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  integer , intent(in) :: nx 
  integer , intent(in) :: ny 
  integer , intent(in) :: nz 
  integer , intent(in) :: npaire 
  real(8) , intent(inout) :: ty(nx,nz,ny) 
  real(8) , intent(in) :: uy(nx,nz,ny) 
  real(8) , intent(inout) :: ry(nx,nz,ny) 
  real(8) , intent(inout) :: sy(nx,nz) 
  real(8) , intent(in) :: sfy(ny) 
  real(8) , intent(in) :: ssy(ny) 
  real(8) , intent(in) :: swy(ny) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer :: k, i, j 
!-----------------------------------------------
!
!********************************************************************
!
   if (ncly==0) then
      do k=1,nz
      do i=1,nx
         ty(i,k,1)=asjy*(uy(i,k,2)-uy(i,k,1   )&
                        -uy(i,k,1)+uy(i,k,ny  ))&
                  +bsjy*(uy(i,k,3)-uy(i,k,1   )&
                        -uy(i,k,1)+uy(i,k,ny-1))
         ry(i,k,1)=-1.
         ty(i,k,2)=asjy*(uy(i,k,3)-uy(i,k,2 )&
                        -uy(i,k,2)+uy(i,k,1 ))&
                  +bsjy*(uy(i,k,4)-uy(i,k,2 )&
                        -uy(i,k,2)+uy(i,k,ny))
         ry(i,k,2)=0.
      enddo
      enddo
      do j=3,ny-2
      do k=1,nz
      do i=1,nx
         ty(i,k,j)=asjy*(uy(i,k,j+1)-uy(i,k,j  )&
                        -uy(i,k,j  )+uy(i,k,j-1))&
                   +bsjy*(uy(i,k,j+2)-uy(i,k,j  )&
                              -uy(i,k,j  )+uy(i,k,j-2))
         ry(i,k,j)=0.
      enddo
      enddo
      enddo
      do k=1,nz
      do i=1,nx
         ty(i,k,ny-1)=asjy*(uy(i,k,ny  )-uy(i,k,ny-1)&
                           -uy(i,k,ny-1)+uy(i,k,ny-2))&
                     +bsjy*(uy(i,k,1   )-uy(i,k,ny-1)&
                           -uy(i,k,ny-1)+uy(i,k,ny-3))
         ry(i,k,ny-1)=0.
         ty(i,k,ny  )=asjy*(uy(i,k,1 )-uy(i,k,ny  )&
                           -uy(i,k,ny)+uy(i,k,ny-1))&
                     +bsjy*(uy(i,k,2 )-uy(i,k,ny  )&
                           -uy(i,k,ny)+uy(i,k,ny-2))
         ry(i,k,ny  )=alsajy
      enddo
      enddo
      do j=2,ny
      do k=1,nz
      do i=1,nx
         ty(i,k,j)=ty(i,k,j)-ty(i,k,j-1)*ssy(j)
         ry(i,k,j)=ry(i,k,j)-ry(i,k,j-1)*ssy(j)
      enddo
      enddo
      enddo
      do k=1,nz
      do i=1,nx
         ty(i,k,ny)=ty(i,k,ny)*swy(ny)
         ry(i,k,ny)=ry(i,k,ny)*swy(ny)
      enddo
      enddo
      do j=ny-1,1,-1
      do k=1,nz
      do i=1,nx
         ty(i,k,j)=(ty(i,k,j)-sfy(j)*ty(i,k,j+1))*swy(j)
         ry(i,k,j)=(ry(i,k,j)-sfy(j)*ry(i,k,j+1))*swy(j)
      enddo
      enddo
      enddo
      do k=1,nz
      do i=1,nx
         sy(i,k)=(   ty(i,k,1)-alsajy*ty(i,k,ny))/&
                 (1.+ry(i,k,1)-alsajy*ry(i,k,ny))
      enddo
      enddo
      do j=1,ny
      do k=1,nz
      do i=1,nx
         ty(i,k,j)=ty(i,k,j)-sy(i,k)*ry(i,k,j)
      enddo
      enddo
      enddo
   endif
!
   if (ncly==1) then
      if (npaire==1) then
         do k=1,nz
         do i=1,nx
            ty(i,k,1)=asjy*(uy(i,k,2)-uy(i,k,1)&
                           -uy(i,k,1)+uy(i,k,2))&
                     +bsjy*(uy(i,k,3)-uy(i,k,1)&
                           -uy(i,k,1)+uy(i,k,3))
            ty(i,k,2)=asjy*(uy(i,k,3)-uy(i,k,2)&
                           -uy(i,k,2)+uy(i,k,1))&
                     +bsjy*(uy(i,k,4)-uy(i,k,2)&
                           -uy(i,k,2)+uy(i,k,2))
         enddo
         enddo
         do k=1,nz
         do i=1,nx
         do j=3,ny-2
            ty(i,k,j)=asjy*(uy(i,k,j+1)-uy(i,k,j  )&
                           -uy(i,k,j  )+uy(i,k,j-1))&
                     +bsjy*(uy(i,k,j+2)-uy(i,k,j  )&
                           -uy(i,k,j  )+uy(i,k,j-2))
         enddo
         enddo
         enddo
         do k=1,nz
         do i=1,nx
            ty(i,k,ny-1)=asjy*(uy(i,k,ny  )-uy(i,k,ny-1)&
                              -uy(i,k,ny-1)+uy(i,k,ny-2))&
                        +bsjy*(uy(i,k,ny-1)-uy(i,k,ny-1)&
                              -uy(i,k,ny-1)+uy(i,k,ny-3))
            ty(i,k,ny  )=asjy*(uy(i,k,ny-1)-uy(i,k,ny  )&
                              -uy(i,k,ny  )+uy(i,k,ny-1))&
                        +bsjy*(uy(i,k,ny-2)-uy(i,k,ny  )&
                              -uy(i,k,ny  )+uy(i,k,ny-2))
         enddo
         enddo
         do j=2,ny
         do k=1,nz
         do i=1,nx
            ty(i,k,j)=ty(i,k,j)-ty(i,k,j-1)*ssy(j)
         enddo
         enddo
         enddo
         do k=1,nz
         do i=1,nx
            ty(i,k,ny)=ty(i,k,ny)*swy(ny)
         enddo
         enddo
         do j=ny-1,1,-1
         do k=1,nz
         do i=1,nx
            ty(i,k,j)=(ty(i,k,j)-sfy(j)*ty(i,k,j+1))*swy(j)
         enddo
         enddo
         enddo
      endif
      if (npaire==0) then
         do k=1,nz
         do i=1,nx
            ty(i,k,1)=0.
            ty(i,k,2)=asjy*(uy(i,k,3)-uy(i,k,2)&
                           -uy(i,k,2)+uy(i,k,1))&
                     +bsjy*(uy(i,k,4)-uy(i,k,2)&
                           -uy(i,k,2)-uy(i,k,2))
         enddo
         enddo
         do k=1,nz
         do i=1,nx
         do j=3,ny-2
            ty(i,k,j)=asjy*(uy(i,k,j+1)-uy(i,k,j  )&
                           -uy(i,k,j  )+uy(i,k,j-1))&
                     +bsjy*(uy(i,k,j+2)-uy(i,k,j  )&
                           -uy(i,k,j  )+uy(i,k,j-2))
         enddo
         enddo
         enddo
         do k=1,nz
         do i=1,nx
            ty(i,k,ny-1)=asjy*(uy(i,k,ny  )-uy(i,k,ny-1)&
                              -uy(i,k,ny-1)+uy(i,k,ny-2))&
                        +bsjy*(-uy(i,k,ny-1)-uy(i,k,ny-1)&
                               -uy(i,k,ny-1)+uy(i,k,ny-3))
            ty(i,k,ny  )=0.
         enddo
         enddo
         do j=2,ny
         do k=1,nz
         do i=1,nx
            ty(i,k,j)=ty(i,k,j)-ty(i,k,j-1)*ssy(j)
         enddo
         enddo
         enddo
         do k=1,nz
         do i=1,nx
            ty(i,k,ny)=ty(i,k,ny)*swy(ny)
         enddo
         enddo
         do j=ny-1,1,-1
         do k=1,nz
         do i=1,nx
            ty(i,k,j)=(ty(i,k,j)-sfy(j)*ty(i,k,j+1))*swy(j)
         enddo
         enddo
         enddo
      endif
   endif
!
   if (ncly==2) then
      do k=1,nz
      do i=1,nx
         ty(i,k,1)=as1y*uy(i,k,1)+bs1y*uy(i,k,2)&
                  +cs1y*uy(i,k,3)+ds1y*uy(i,k,4)
         ty(i,k,2)=as2y*(uy(i,k,3)-uy(i,k,2)&
                        -uy(i,k,2)+uy(i,k,1))
      enddo
      enddo
      do k=1,nz
      do i=1,nx
      do j=3,ny-2
         ty(i,k,j)=asjy*(uy(i,k,j+1)-uy(i,k,j  )&
                        -uy(i,k,j  )+uy(i,k,j-1))&
                  +bsjy*(uy(i,k,j+2)-uy(i,k,j  )&
                        -uy(i,k,j  )+uy(i,k,j-2))
      enddo
      enddo
      enddo
      do k=1,nz
      do i=1,nx
         ty(i,k,ny-1)=asmy*(uy(i,k,ny  )-uy(i,k,ny-1)&
                           -uy(i,k,ny-1)+uy(i,k,ny-2))
         ty(i,k,ny  )=asny*uy(i,k,ny  )+bsny*uy(i,k,ny-1)&
                     +csny*uy(i,k,ny-2)+dsny*uy(i,k,ny-3)
      enddo
      enddo
      do j=2,ny
      do k=1,nz
      do i=1,nx
         ty(i,k,j)=ty(i,k,j)-ty(i,k,j-1)*ssy(j)
      enddo
      enddo
      enddo
      do k=1,nz
      do i=1,nx
         ty(i,k,ny)=ty(i,k,ny)*swy(ny)
      enddo
      enddo
      do j=ny-1,1,-1
      do k=1,nz
      do i=1,nx
         ty(i,k,j)=(ty(i,k,j)-sfy(j)*ty(i,k,j+1))*swy(j)
      enddo
      enddo
      enddo
   endif
!
   return  
end subroutine deryy1
!
!********************************************************************
!
subroutine derzz(tz,uz,rz,sz,sfz,ssz,swz,nx,ny,nz,npaire) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
  USE paramz_m 
  USE derivez_m 
!...Translated by PSUITE Trans90                  4.3ZH  9:45:39   1/26/ 4  
!...Switches: -yv INDDO=0 -nbejkno
  implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  integer , intent(in) :: nx 
  integer , intent(in) :: ny 
  integer , intent(in) :: nz 
  integer , intent(in) :: npaire 
  real(8) , intent(inout) :: tz(nx,ny,nz) 
  real(8) , intent(in) :: uz(nx,ny,nz) 
  real(8) , intent(inout) :: rz(nx,ny,nz) 
  real(8) , intent(inout) :: sz(nx,ny) 
  real(8) , intent(in) :: sfz(nz) 
  real(8) , intent(in) :: ssz(nz) 
  real(8) , intent(in) :: swz(nz) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer :: j, i, k 
!-----------------------------------------------
!
!********************************************************************
!
   if (nclz==0) then
      do j=1,ny
      do i=1,nx
         tz(i,j,1)=askz*(uz(i,j,2)-uz(i,j,1   )&
                        -uz(i,j,1)+uz(i,j,nz  ))&
                  +bskz*(uz(i,j,3)-uz(i,j,1   )&
                        -uz(i,j,1)+uz(i,j,nz-1))
         rz(i,j,1)=-1.
         tz(i,j,2)=askz*(uz(i,j,3)-uz(i,j,2 )&
                        -uz(i,j,2)+uz(i,j,1 ))&
                  +bskz*(uz(i,j,4)-uz(i,j,2 )&
                        -uz(i,j,2)+uz(i,j,nz))
         rz(i,j,2)=0.
      enddo
      enddo
      do k=3,nz-2
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=askz*(uz(i,j,k+1)-uz(i,j,k  )&
                        -uz(i,j,k  )+uz(i,j,k-1))&
                  +bskz*(uz(i,j,k+2)-uz(i,j,k  )&
                        -uz(i,j,k  )+uz(i,j,k-2))
          rz(i,j,k)=0.
      enddo
      enddo
      enddo
      do j=1,ny
      do i=1,nx
         tz(i,j,nz-1)=askz*(uz(i,j,nz  )-uz(i,j,nz-1)&
                           -uz(i,j,nz-1)+uz(i,j,nz-2))&
                     +bskz*(uz(i,j,1   )-uz(i,j,nz-1)&
                           -uz(i,j,nz-1)+uz(i,j,nz-3))
         rz(i,j,nz-1)=0.
         tz(i,j,nz  )=askz*(uz(i,j,1 )-uz(i,j,nz  )&
                           -uz(i,j,nz)+uz(i,j,nz-1))&
                     +bskz*(uz(i,j,2 )-uz(i,j,nz  )&
                           -uz(i,j,nz)+uz(i,j,nz-2))
         rz(i,j,nz  )=alsakz
      enddo
      enddo
      do k=2,nz
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=tz(i,j,k)-tz(i,j,k-1)*ssz(k)
         rz(i,j,k)=rz(i,j,k)-rz(i,j,k-1)*ssz(k)
      enddo
      enddo
      enddo
      do j=1,ny
      do i=1,nx
         tz(i,j,nz)=tz(i,j,nz)*swz(nz)
         rz(i,j,nz)=rz(i,j,nz)*swz(nz)
      enddo
      enddo
      do k=nz-1,1,-1
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=(tz(i,j,k)-sfz(k)*tz(i,j,k+1))*swz(k)
         rz(i,j,k)=(rz(i,j,k)-sfz(k)*rz(i,j,k+1))*swz(k)
      enddo
      enddo
      enddo
      do j=1,ny
      do i=1,nx
         sz(i,j)=(   tz(i,j,1)-alsakz*tz(i,j,nz))/&
                 (1.+rz(i,j,1)-alsakz*rz(i,j,nz))
      enddo
      enddo
      do k=1,nz
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=tz(i,j,k)-sz(i,j)*rz(i,j,k)
      enddo
      enddo
      enddo
   endif
!
   if (nclz==1) then
      if (npaire==1) then
         do j=1,ny
         do i=1,nx
            tz(i,j,1)=askz*(uz(i,j,2)-uz(i,j,1)&
                           -uz(i,j,1)+uz(i,j,2))&
                     +bskz*(uz(i,j,3)-uz(i,j,1)&
                           -uz(i,j,1)+uz(i,j,3))
            tz(i,j,2)=askz*(uz(i,j,3)-uz(i,j,2)&
                           -uz(i,j,2)+uz(i,j,1))&
                     +bskz*(uz(i,j,4)-uz(i,j,2)&
                           -uz(i,j,2)+uz(i,j,2))
         enddo
         enddo
         do k=3,nz-2
         do j=1,ny
         do i=1,nx
            tz(i,j,k)=askz*(uz(i,j,k+1)-uz(i,j,k  )&
                           -uz(i,j,k  )+uz(i,j,k-1))&
                     +bskz*(uz(i,j,k+2)-uz(i,j,k  )&
                           -uz(i,j,k  )+uz(i,j,k-2))
         enddo
         enddo
         enddo
         do j=1,ny
         do i=1,nx
            tz(i,j,nz-1)=askz*(uz(i,j,nz  )-uz(i,j,nz-1)&
                              -uz(i,j,nz-1)+uz(i,j,nz-2))&
                        +bskz*(uz(i,j,nz-1)-uz(i,j,nz-1)&
                              -uz(i,j,nz-1)+uz(i,j,nz-3))
            tz(i,j,nz  )=askz*(uz(i,j,nz-1)-uz(i,j,nz  )&
                              -uz(i,j,nz  )+uz(i,j,nz-1))&
                        +bskz*(uz(i,j,nz-2)-uz(i,j,nz  )&
                              -uz(i,j,nz  )+uz(i,j,nz-2))
         enddo
         enddo
         do k=2,nz
         do j=1,ny
         do i=1,nx
            tz(i,j,k)=tz(i,j,k)-tz(i,j,k-1)*ssz(k)
         enddo
         enddo
         enddo
         do j=1,ny
         do i=1,nx
            tz(i,j,nz)=tz(i,j,nz)*swz(nz)
         enddo
         enddo
         do k=nz-1,1,-1
         do j=1,ny
         do i=1,nx
            tz(i,j,k)=(tz(i,j,k)-sfz(k)*tz(i,j,k+1))*swz(k)
         enddo
         enddo
         enddo
      endif
      if (npaire==0) then
         do j=1,ny
         do i=1,nx
            tz(i,j,1)=0.
            tz(i,j,2)=askz*(uz(i,j,3)-uz(i,j,2)&
                           -uz(i,j,2)+uz(i,j,1))&
                     +bskz*(uz(i,j,4)-uz(i,j,2)&
                           -uz(i,j,2)-uz(i,j,2))
         enddo
         enddo
         do k=3,nz-2
         do j=1,ny
         do i=1,nx
            tz(i,j,k)=askz*(uz(i,j,k+1)-uz(i,j,k  )&
                           -uz(i,j,k  )+uz(i,j,k-1))&
                     +bskz*(uz(i,j,k+2)-uz(i,j,k  )&
                           -uz(i,j,k  )+uz(i,j,k-2))
         enddo
         enddo
         enddo
         do j=1,ny
         do i=1,nx
            tz(i,j,nz-1)=askz*(uz(i,j,nz  )-uz(i,j,nz-1)&
                              -uz(i,j,nz-1)+uz(i,j,nz-2))&
                        +bskz*(-uz(i,j,nz-1)-uz(i,j,nz-1)&
                              -uz(i,j,nz-1)+uz(i,j,nz-3))
            tz(i,j,nz  )=0.
         enddo
         enddo
         do k=2,nz
         do j=1,ny
         do i=1,nx
            tz(i,j,k)=tz(i,j,k)-tz(i,j,k-1)*ssz(k)
         enddo
         enddo
         enddo
         do j=1,ny
         do i=1,nx
            tz(i,j,nz)=tz(i,j,nz)*swz(nz)
         enddo
         enddo
         do k=nz-1,1,-1
         do j=1,ny
         do i=1,nx
            tz(i,j,k)=(tz(i,j,k)-sfz(k)*tz(i,j,k+1))*swz(k)
         enddo
         enddo
         enddo
      endif
   endif
!
   if (nclz==2) then
      do j=1,ny
      do i=1,nx
         tz(i,j,1)=as1z*uz(i,j,1)+bs1z*uz(i,j,2)&
                  +cs1z*uz(i,j,3)+ds1z*uz(i,j,4)
         tz(i,j,2)=as2z*(uz(i,j,3)-uz(i,j,2)&
                        -uz(i,j,2)+uz(i,j,1))
      enddo
      enddo
      do k=3,nz-2
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=askz*(uz(i,j,k+1)-uz(i,j,k  )&
                        -uz(i,j,k  )+uz(i,j,k-1))&
                  +bskz*(uz(i,j,k+2)-uz(i,j,k  )&
                        -uz(i,j,k  )+uz(i,j,k-2))
      enddo
      enddo
      enddo
      do j=1,ny
      do i=1,nx
         tz(i,j,nz-1)=asmz*(uz(i,j,nz  )-uz(i,j,nz-1)&
                           -uz(i,j,nz-1)+uz(i,j,nz-2))
         tz(i,j,nz  )=asnz*uz(i,j,nz  )+bsnz*uz(i,j,nz-1)&
                     +csnz*uz(i,j,nz-2)+dsnz*uz(i,j,nz-3)
      enddo
      enddo
      do k=2,nz
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=tz(i,j,k)-tz(i,j,k-1)*ssz(k)
      enddo
      enddo
      enddo
      do j=1,ny
      do i=1,nx
         tz(i,j,nz)=tz(i,j,nz)*swz(nz)
      enddo
      enddo
      do k=nz-1,1,-1
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=(tz(i,j,k)-sfz(k)*tz(i,j,k+1))*swz(k)
      enddo
      enddo
      enddo
   endif
!
   return  
end subroutine derzz
!
