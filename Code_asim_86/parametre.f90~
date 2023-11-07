!
!********************************************************************
!
subroutine parametre(nx,ny,nz,nxr,nyr,rdxy) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
  USE paramx_m 
  USE paramy_m 
  USE paramz_m 
  USE param1_m 
  USE param2_m 
  USE delta_m 
  USE barreau_m 
  USE ecoulement_m 
  USE controle_m 
  USE sauvegarde_m 
  USE module_avs_m 
!...Translated by PSUITE Trans90                  4.3ZH 16:31:24   1/27/ 4  
!...Switches: -yv INDDO=0 -nbejkno
  implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  integer , intent(in) :: nx,nxr 
  integer , intent(in) :: ny,nyr 
  integer , intent(in) :: nz 
  integer , intent(in) :: rdxy 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer :: longueur 
!  real(8) :: re, sc, pi, theta, dtv, cfl 
  real(8) :: sc, pi, theta, dtv, cfl 
  character :: a*80 
!-----------------------------------------------
!
!********************************************************************
!
1000 format(a,8x) 
1001 format(i6) 
1002 format(f12.8) 
1003 format(a,80x) 
      open(10,file='general.prm',status='unknown',form='formatted') 
      read (10,1000) a 
      read (10,1000) a 
      read (10,1000) a 
      read (10,1000) a 
      read (10,1001) nclx 
      read (10,1000) a 
      read (10,1001) ncly 
      read (10,1000) a 
      read (10,1001) nclz 
      read (10,1000) a 
      read (10,1002) u1 
      read (10,1000) a 
      read (10,1002) u2 
      read (10,1000) a 
      read (10,1002) rayon 
      read (10,1000) a 
      read (10,1002) bruit 
      read (10,1000) a 
      read (10,1002) xlx 
      read (10,1000) a 
      read (10,1002) yly 
      read (10,1000) a 
      read (10,1002) zlz 
      read (10,1000) a 
      read (10,1002) re 
      read (10,1000) a 
      read (10,1002) sc 
      read (10,1000) a 
      read (10,1001) iecoule 
      read (10,1000) a 
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
   print *,'nx,ny,nz=',nx,ny,nz 
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
   print *,'ivirtuel=',ivirtuel 
!
   open(10,file='fichiers.prm',status='unknown',form='formatted') 
      read (10,1000) a 
      read (10,1000) a 
      read (10,1000) a 
      read (10,1000) a 
      read (10,1001) isave 
      read (10,1000) a 
      read (10,1001) ilit 
      read (10,1000) a 
      read (10,1001) idebmod 
      read (10,1000) a 
      read (10,1001) imodulo 
      read (10,1000) a 
      read (10,1001) idemarre 
      read (10,1000) a 
      read (10,1001) icommence 
      read (10,1000) a 
      read (10,1001) irecord 
      read (10,1000) a 
      read (10,1001) istat 
      read (10,1000) a 
      read (10,1003) nchamp 
      read (10,1000) a 
      read (10,1003) nchamp1 
      read (10,1000) a 
      read (10,1003) filecharge 
      read (10,1000) a 
      read (10,1003) filesauve 
      read (10,1000) a 
      read (10,1003) filebruit 
   close(10) 
!
   print *,'idemarre=',idemarre 
!
   open(10,file='assimilation.prm',status='unknown',form='formatted') 
      read (10,1000) a 
      read (10,1000) a 
      read (10,1000) a 
      read (10,1000) a 
      read (10,1001) iobs 
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
   pi=acos(-1.) 
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
   if (iavs==1) then 
      open(10,file='avs.prm',status='unknown',form='formatted') 
         read (10,1000) a 
         read (10,1000) a 
         read (10,1000) a 
         read (10,1000) a 
         read (10,1001) idebavs 
         read (10,1000) a 
         read (10,1001) imodavs 
         read (10,1000) a 
         read (10,1003) path_network 
         read (10,1000) a 
         read (10,1003) nom_network 
         read (10,1000) a 
         read (10,1003) path_files 
         read (10,1000) a 
         read (10,1003) nom_file 
         read (10,1000) a 
         read (10,1003) nom_script 
         read (10,1000) a 
         read (10,1003) nom_film 
         read (10,1000) a 
         read (10,1000) a 
         read (10,1000) a 
         read (10,1000) a 
         read (10,1002) val_iso 
         read (10,1000) a 
         read (10,1001) ndsize 
         read (10,1000) a 
         read (10,1000) a 
         read (10,1001) nx_min 
         read (10,1000) a 
         read (10,1001) nx_max 
         read (10,1000) a 
         read (10,1001) ny_min 
         read (10,1000) a 
         read (10,1001) ny_max 
         read (10,1000) a 
         read (10,1001) nz_min 
         read (10,1000) a 
         read (10,1001) nz_max 
      close(10) 
   endif 
!
   xnu=1./re 
   xkappa=xnu/sc 
!
   if (nclx==0) dx=xlx/nx 
   if (nclx==1 .or. nclx==2) dx=xlx/(nx-1.) 
!
   if (ncly==0) dy=yly/ny 
   if (ncly==1 .or. ncly==2) dy=yly/(ny-1.) 
!
   if (nz>1) then 
      if (nclz==0) dz=zlz/nz 
      if (nclz==1 .or. nclz==2) dz=zlz/(nz-1.) 
   else 
      dz=zlz 
   endif
!
   dx2=dx*dx 
   dy2=dy*dy 
   dz2=dz*dz 
!
!  resolucion grosera
   dxg=rdxy*dx
   dyg=rdxy*dy
!
!  taille dominio reducido
   xlxr=(nxr-1)*dx
   ylyr=(nyr-1)*dy
!
   if (ncly==0) delty=yly/ny 
   if (ncly==1) delty=yly/(ny-1) 
!
   if (nclz==0) deltz=zlz/nz 
   if (nclz==1) deltz=zlz/(nz-1) 
!
   y1=yly/2. 
   z1=zlz/2. 
!
   if (iecoule==5) then 
      if (nz>1) then 
         delty=0. 
         y1=0. 
      else 
         deltz=0. 
         z1=0. 
      endif
   endif
!
   print *,'dx,dy,dz = ',dx,dy,dz 
   print *,'dxg,dyg = ',dxg,dyg 
   print *,'xlxr,ylyr = ',xlxr,ylyr 
!
!   if (nschema==1) dt=0.5*dx/(pi*u1)
   if (nschema==1) dt=0.01
!
   if (nschema==2 .or. nschema==3) dt=0.0625
!dt=0.5*sqrt(3.)*dx/(1.989*u1) 
!
   print *,'dt (AB,RK3 ou RK4)=',dt 
!
   if (ivirtuel==1) then 
      dtv=(((-betav)-sqrt(betav*betav-2*alphav*xkk))/alphav)*sqrt(3.) 
      print *,'dtv=',dtv 
      if (dtv<dt) dt=dtv
      if (dt<dtv) dt=dt
   endif
!
   cfl=u1*dt/dx 
   print *,'CFL = ',cfl 
   print *,'dt=',dt 
!
!   dt=dt/2.
   print *,'****************' 
   print *,'dt modifie HARD=',dt 
   print *,'****************' 
!
   if (ientree==2) then 
      longueur=index(filebruit,' ')-1 
      open(40,file=filebruit(1:longueur),form='unformatted',status='unknown') 
   endif
   if (nxboite/=0 .or. ientree==3) then 
      longueur=index(filesauve,' ')-1 
      open(50,file=filesauve(1:longueur),form='unformatted',status='unknown') 
   endif
!
   return  
end subroutine parametre
