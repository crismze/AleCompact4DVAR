! *****************************
program extraction  
! *****************************
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
   integer,parameter :: nx=257,ny=257,nz=1
   integer,parameter :: m1=2,m2=2;
   integer,parameter :: nxm=nx-1,nym=ny-1,nzm=nz
   integer,parameter :: mx=nx+1,my=ny+1,mz=nz
   integer,parameter :: nyp=257,nzp=1,mw=4*mx*my*mz,nvec=256
   integer,parameter :: mxf=2*nx,myf=2*ny,mzf=2*nz,mb=2*mx
   integer,parameter :: mxy=mx*my,myz=my*mz
   integer,parameter :: nwork=max((3*nym+3)*(3*nxm+3),2*ny*nx)
   integer,parameter :: ntrigsX=4*nx,ntrigsY=4*ny,nfax=19
!
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
   integer :: i,j,k,nxz,nxy,nxyz,nboucle
   integer :: longueur
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
   real(8) :: yinf,beta,den,xnum,alpha,xcx,den1,den2,den3,den4,xnum1,cst
   real(8) :: pi,twopi,fourpi,r,x11,x,val
   integer :: i5,i6,i55
   integer :: j1,j2,j3,j4,jj
!**************************************************************************
!
   character(len=80) lect,fichier
   character(len=4)  suffixe
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
     twopi=2.*pi
     fourpi=4.*pi
!     yinf cest la borne inf du domaine que l'on veut obtenir
     yinf=-yly/2.
!
      beta=2.
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

!
! ***********************************************************************
!      
      open(10,file='cmlcn010',form='unformatted')
      read(10) ux,uy
      close(10)
!
      do i=1,nx
      do j=1,ny
         sc1(i,j,1)=ux(i,j,1)
      enddo
      enddo
   u1=0.
1002 format(120F12.7)
1003 format(7F12.7)
   open (10,file='profiluxpois.dat',form='formatted')
   do j=1,ny
   write (10,1002) yp(j),ux(1,j,1)-u1,((1.)*dx)+ux(1,j,1)-u1,((2.)*dx)+ux(2,j,1)-u1,(3.*dx)+ux(3,j,1)-u1,(4*dx)+ux(4,j,1)-u1,&
                   ((5.)*dx)+ux(5,j,1)-u1,((6.)*dx)+ux(6,j,1)-u1,((7.)*dx)+ux(7,j,1)-u1,((8.)*dx)+ux(8,j,1)-u1,&
                   ((9.)*dx)+ux(9,j,1)-u1,((10.)*dx)+ux(10,j,1)-u1,((11.)*dx)+ux(11,j,1)-u1,((12.)*dx)+ux(12,j,1)-u1, &
                   ((12.)*dx)+ux(12,j,1)-u1,((13.)*dx)+ux(13,j,1)-u1,((14.)*dx)+ux(14,j,1)-u1,((15.)*dx)+ux(15,j,1)-u1, &
                   ((16.)*dx)+ux(16,j,1)-u1,((17.)*dx)+ux(17,j,1)-u1,((18.)*dx)+ux(18,j,1)-u1,((19.)*dx)+ux(19,j,1)-u1, &
                   ((20.)*dx)+ux(20,j,1)-u1,((21.)*dx)+ux(21,j,1)-u1,((22.)*dx)+ux(22,j,1)-u1,((23.)*dx)+ux(23,j,1)-u1, &
                   ((24.)*dx)+ux(24,j,1)-u1,((25.)*dx)+ux(25,j,1)-u1,((26.)*dx)+ux(26,j,1)-u1,((27.)*dx)+ux(27,j,1)-u1, &
                   ((28.)*dx)+ux(28,j,1)-u1,((29.)*dx)+ux(29,j,1)-u1,((30.)*dx)+ux(30,j,1)-u1,((31.)*dx)+ux(31,j,1)-u1, &
                   ((32.)*dx)+ux(32,j,1)-u1,((33.)*dx)+ux(33,j,1)-u1,((34.)*dx)+ux(34,j,1)-u1,((35.)*dx)+ux(35,j,1)-u1, &
                   ((36.)*dx)+ux(36,j,1)-u1,((37.)*dx)+ux(37,j,1)-u1,((38.)*dx)+ux(38,j,1)-u1,((39.)*dx)+ux(39,j,1)-u1, &
                   ((40.)*dx)+ux(40,j,1)-u1,((41.)*dx)+ux(41,j,1)-u1,((42.)*dx)+ux(42,j,1)-u1,((43.)*dx)+ux(43,j,1)-u1, &
                   ((44.)*dx)+ux(44,j,1)-u1,((45.)*dx)+ux(45,j,1)-u1,((46.)*dx)+ux(46,j,1)-u1,((47.)*dx)+ux(47,j,1)-u1, &
                   ((48.)*dx)+ux(48,j,1)-u1,((49.)*dx)+ux(49,j,1)-u1,((50.)*dx)+ux(50,j,1)-u1,((51.)*dx)+ux(51,j,1)-u1, &
                   ((52.)*dx)+ux(52,j,1)-u1,((53.)*dx)+ux(53,j,1)-u1,((54.)*dx)+ux(54,j,1)-u1,((55.)*dx)+ux(55,j,1)-u1, &
                   ((56.)*dx)+ux(56,j,1)-u1,((57.)*dx)+ux(57,j,1)-u1,((58.)*dx)+ux(58,j,1)-u1,((59.)*dx)+ux(59,j,1)-u1, &
                   ((60.)*dx)+ux(60,j,1)-u1,((61.)*dx)+ux(61,j,1)-u1,((62.)*dx)+ux(62,j,1)-u1,((63.)*dx)+ux(53,j,1)-u1, &
                   ((64.)*dx)+ux(64,j,1)-u1,((65.)*dx)+ux(65,j,1)-u1,((66.)*dx)+ux(66,j,1)-u1,((67.)*dx)+ux(67,j,1)-u1, &
                   ((68.)*dx)+ux(68,j,1)-u1,((69.)*dx)+ux(69,j,1)-u1,((70.)*dx)+ux(70,j,1)-u1,((71.)*dx)+ux(71,j,1)-u1, &
                   ((72.)*dx)+ux(72,j,1)-u1,((73.)*dx)+ux(73,j,1)-u1,((74.)*dx)+ux(74,j,1)-u1,((75.)*dx)+ux(75,j,1)-u1, &
                   ((76.)*dx)+ux(76,j,1)-u1,((77.)*dx)+ux(77,j,1)-u1,((78.)*dx)+ux(78,j,1)-u1,((79.)*dx)+ux(79,j,1)-u1, &
                   ((80.)*dx)+ux(80,j,1)-u1,((81.)*dx)+ux(81,j,1)-u1,((82.)*dx)+ux(82,j,1)-u1,((83.)*dx)+ux(83,j,1)-u1, &
                   ((84.)*dx)+ux(84,j,1)-u1,((85.)*dx)+ux(85,j,1)-u1,((86.)*dx)+ux(86,j,1)-u1,((87.)*dx)+ux(87,j,1)-u1, &
                   ((88.)*dx)+ux(88,j,1)-u1,((89.)*dx)+ux(89,j,1)-u1,((90.)*dx)+ux(90,j,1)-u1,((91.)*dx)+ux(91,j,1)-u1
   enddo
   close(10)
   open (15,file='profiluxsortie.dat',form='formatted')
   do j=1,ny
   write (15,1003) yp(j),((nx-6.)*dx)+ux(nx-6,j,1),((nx-5.)*dx)+ux(nx-5,j,1),((nx-4.)*dx)+ux(nx-4,j,1),((nx-3.)*dx)+ux(nx-3,j,1),&
                   ((nx-1.)*dx)+ux(nx-1,j,1),((nx)*dx)+ux(nx,j,1)
   enddo
   close(15)
!
end program extraction
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
!
!
   character(len=80) lect_vort,avs_vort,lect_q,avs_q,lect_pdf,lect_moy,lect_exp
   character(len=80) lect,fichier,a
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

