!***************************************************************************
!
subroutine mesure (ux,uy,uz,py1,pz1,sp,gx,gy,gz,sg,hx,hy,hz,sh,&
     tx,ty,tz,di,tr,us,ps,&
     sx,sy,sz,fenetre1,error1,&
     ffx,fsx,fwx,ffy,fsy,fwy,ffz,fsz,fwz,&
     ffxp,fsxp,fwxp,&
     ffyp,fsyp,fwyp,&
     ffzp,fszp,fwzp,&
     txy1,txy2,txy3,txy4,&
     txz1,txz2,txz3,&
     tyz1,tyz2,tyz3,&
     itime,nx,ny,nz,nxm,nym,nzm,mx,my,mz,num,num_fich,ppm)
!
!***************************************************************************

!
   USE sauvegarde_m
   USE module_avs_m
   USE ecoulement_m
!
   implicit none
!

!     
   real(8),dimension(nx,ny,nz) :: ux,uy,uz,hx,hy,hz,gx,gy,gz
   real(8),dimension(nx,ny,nz) :: tx,ty,tz,sp,sg,sh,di,py,py1,pz1
   real(8),dimension(nxm,nym,nzm) :: ppm
   real(8),dimension(nx,4,1000) :: tyz1,tyz2,tyz3
   real(8),dimension(ny,4,1000) :: txz1,txz2,txz3
   real(8),dimension(nz,4,1000) :: txy1,txy2,txy3,txy4
   real(8),dimension(mx,my,mz) :: tr,us,ps
   real(8),dimension(ny,nz) :: sx
   real(8),dimension(nx,nz) :: sy
   real(8),dimension(nx,ny) :: sz
   real(8),dimension(nx) ::  ffx,fsx,fwx,ffxp,fsxp,fwxp
   real(8),dimension(ny) :: ffy,fsy,fwy,ffyp,fsyp,fwyp
   real(8),dimension(nz) :: ffz,fsz,fwz,ffzp,fszp,fwzp
!
   integer fenetre1,error1
!
  integer  :: itime,i,j,k,nxz,longueur1 
  integer  :: nx,ny,nz,mx,my,mz,nxm,nym,nzm 
  integer  :: num 
  integer  :: num_fich 
!
   character(len=80) nom_file2
   character(len=4) car
!      
!      print *,itime,imodulo,idebmod
!      pause'1'
   if (mod(itime,imodulo).eq.0 .and. itime.ge.idebmod) then 
!         pause'11'
      call ecrit (ux,uy,uz,py1,sp,gx,gy,gz,sg,hx,hy,hz,sh,&
          nx,ny,nz)
      isave=isave+1
!          pause'22'
   endif
!      pause'2'
      if (mod(itime,1000).eq.0) then
         call sauve (ux,uy,uz,py1,pz1,gx,gy,gz,hx,hy,hz,sp&
        ,ppm,nx,ny,nz,nxm,nym,nzm,itime)
      endif         
!
   if (iavs.eq.1) then
      if (itime.ge.idebavs .and. mod(itime,imodavs).eq.0) then
         call calcwn(ux,uy,uz,gx,gy,gz,&
              tx,ty,tz,di,sx,sy,sz,&
              ffx,fsx,fwx,ffy,fsy,fwy,ffz,fsz,fwz,&
              ffxp,fsxp,fwxp,&
              ffyp,fsyp,fwyp,&
              ffzp,fszp,fwzp,&
              nx,ny,nz)
         call numcar(num_fich,car)
         longueur1=index(nchamp,' ')-1
         nom_file2=nchamp(1:longueur1)//'_wn_'//car
         call sauve_avs (tz,nom_file2,nx,ny,nz)
         nxz=nx*nz
         num_fich=num_fich+1
      endif
   endif
!
   return
end subroutine mesure
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
   real(8),dimension(nx,ny,nz) :: ux,uy,uz,tx,ty,tz,gx,gy,gz,di
   real(8),dimension(ny,nz) :: sx
   real(8),dimension(nx,nz) :: sy
   real(8),dimension(nx,ny) :: sz
   real(8),dimension(nx) :: ffx,fsx,fwx,ffxp,fsxp,fwxp
   real(8),dimension(ny) :: ffy,fsy,fwy,ffyp,fsyp,fwyp
   real(8),dimension(nz) :: ffz,fsz,fwz,ffzp,fszp,fwzp
!
  integer  :: nx 
  integer  :: ny 
  integer  :: nz 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer :: k, j, i 
  real :: wzmin, wzmax 
!-----------------------------------------------
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
      tz(i,j,k)=abs(tx(i,j,k)-ty(i,j,k))
!      if (tz(i,j,k).lt.wzmin) wzmin=tz(i,j,k)
!      if (tz(i,j,k).gt.wzmax) wzmax=tz(i,j,k)
   enddo
   enddo
   enddo
!
   return
end subroutine calcwz
!
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
   real(8),dimension(nx,ny,nz) :: ux,uy,uz,tx,ty,tz,gx,gy,gz,di
   real(8),dimension(ny,nz) :: sx
   real(8),dimension(nx,nz) :: sy
   real(8),dimension(nx,ny) :: sz
   real(8),dimension(nx) :: ffx,fsx,fwx,ffxp,fsxp,fwxp
   real(8),dimension(ny) :: ffy,fsy,fwy,ffyp,fsyp,fwyp
   real(8),dimension(nz) :: ffz,fsz,fwz,ffzp,fszp,fwzp
!
  integer  :: nx 
  integer  :: ny 
  integer  :: nz 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer :: k, j, i 
  real :: wzmin, wzmax 
!-----------------------------------------------
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
subroutine ecrit (ux,uy,uz,py1,sp,gx,gy,gz,sg,hx,hy,hz,sh,&
     nx,ny,nz)
!
!********************************************************************
!
   USE ecoulement_m
   USE sauvegarde_m
!
   implicit none
!
   integer :: nx 
   integer :: ny 
   integer :: nz 
   integer :: k,j,i  
!
   real(8),dimension(nx,ny,nz) :: ux,uy,uz,hx,hy,hz,gx,gy,gz,sp,sg,sh
   real(8),dimension(nx,ny,nz) :: py1
!

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer :: longueur 
  real :: wzmin, wzmax 
!-----------------------------------------------
!
   character(len=4) suffix
   character(len=20) nfichier
!
   call numcar (isave,suffix)
   longueur=index(nchamp,' ')-1
   nfichier=nchamp(1:longueur)//suffix
   longueur=index(nfichier,' ')-1
!
   wzmin= 1.e10
   wzmax=-1.e10
!
!a    open(10,file=nfichier(1:longueur),&
!a         form='unformatted',status='unknown')
   open(10,file=nfichier(1:longueur),form='formatted')
   if (nz.gt.1) then
      if (iscalaire.eq.0) then
         write(10) ux,uy,uz
      else
         write(10) ux,uy,uz,sp
      endif
   else
      if (iscalaire.eq.0) then
!a         write(10) ux,uy
         do j=1,ny
         do i=1,nx
            write(10,*) ux(i,j,1),uy(i,j,1)
         enddo
         enddo 
      else
         write(10) ux,uy,sp
      endif
   endif
   close(10)
!      
   if ((nschema.eq.1).and.((istat.eq.2).or.(istat.eq.4))) then
      open(20,file=nfichier(1:longueur)//'g',&
           form='unformatted',status='unknown')
      open(30,file=nfichier(1:longueur)//'h',&
           form='unformatted',status='unknown')
      if (nz.gt.1) then
         if (iscalaire.eq.0) then
            write(20) gx,gy,gz
            write(30) hx,hy,hz
         else
            write(20) gx,gy,gz,sg
            write(30) hx,hy,hz,sh
         endif
      else
         if (iscalaire.eq.0) then
            write(20) gx,gy
            write(30) hx,hy
         else
            write(20) gx,gy,sg
            write(30) hx,hy,sh
         endif
      endif
      close(20)
      close(30)
   endif
!            
   return
end subroutine ecrit

!********************************************************************
!
subroutine lit (ux,uy,uz,py1,pz1,sp,gx,gy,gz,sg,hx,hy,hz,sh,&
     bxo,byo,bzo,bs1,bxr,byr,bzr,bsr,bxp,byp,bzp,bsp,&
     nx,ny,nz,nxm,nym,nzm,nyp,nzp,ppm)
!
!********************************************************************
!
   USE ecoulement_m
   USE sauvegarde_m
!
   implicit none
!
   real(8),dimension(nx,ny,nz) :: ux,uy,uz,hx,hy,hz,gx,gy,gz,py1,pz1
   real(8),dimension(nxm,nym,nzm) :: ppm
   real(8),dimension(nx,ny,nz) :: sp,sg,sh
   real(8),dimension(ny,nz) :: bx1,by1,bz1,bxr,byr,bzr,bxp,byp,bzp
   real(8),dimension(ny,nz) :: bs1,bsp,bsr
   real(8),dimension(nyp,nzp) :: bxo,byo,bzo
!
  integer , intent(in) :: nx 
  integer , intent(in) :: ny 
  integer , intent(in) :: nz 
  integer , intent(in) :: nyp 
  integer , intent(in) :: nzp,nxm,nym,nzm 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer :: longueur,i,j,k,iturb
!-----------------------------------------------
!
   longueur=index(filecharge,' ')-1
!a    open(10,file=filecharge(1:longueur),&
!a         form='unformatted',status='unknown')
   open(10,file=filecharge(1:longueur),form='formatted')
   if (nz.gt.1) then
      if (iscalaire.eq.0) then
         read(10) ux,uy,uz,py1,pz1,ppm,gx,gy,gz
      else
         read(10) ux,uy,uz,sp
      endif
   else
      if (iscalaire.eq.0) then
!a         read(10) ux,uy,py1,ppm,gx,gy
         do j=1,ny
         do i=1,nx
            read(10,*) ux(i,j,1),uy(i,j,1)
         enddo
         enddo 
      else
         read(10) ux,uy,sp
      endif
   endif
!

!
   if (nschema.eq.1) then
      open(20,file=filecharge(1:longueur)//'g',&
           form='unformatted',status='unknown')
      open(30,file=filecharge(1:longueur)//'h',&
           form='unformatted',status='unknown')
      if (nz.gt.1) then
         if (iscalaire.eq.0) then
            read(20) gx,gy,gz
            read(30) hx,hy,hz
         else
            read(20) gx,gy,gz,sg
            read(30) hx,hy,hz,sh
         endif
!
      else
         if (iscalaire.eq.0) then
!a             read(20) gx,gy
!a             read(30) hx,hy
         else
            read(20) gx,gy,sg
            read(30) hx,hy,sh
         endif
      endif
   endif
!
!   if (nz.gt.1) then
!      if (ientree.eq.3) then
!         if (iscalaire.eq.0) then
!            do i=1,idemarre
!               read(50) bx1,by1,bz1
!            enddo
!            read(50) bxr,byr,bzr
!            read(50) bxp,byp,bzp
!         else
!            do i=1,idemarre
!               read(50) bx1,by1,bz1,bs1
!            enddo
!            read(50) bxr,byr,bzr,bsr
!            if (nschema.eq.2) read(50) bxp,byp,bzp,bsp
!        endif
!     endif
!      if (ientree.eq.2) then
!         do i=1,icommence
!            read(40) bxo,byo,bzo
!         enddo
!      endif
!   else
!      if (ientree.eq.3) then
!         if (iscalaire.eq.0) then
!            do i=1,idemarre
!               read(50) bx1,by1
!            enddo
!            read(50) bxr,byr
!            read(50) bxp,byp
!         else
!            do i=1,idemarre
!               read(50) bx1,by1,bs1
!            enddo
!            read(50) bxr,byr,bsr
!            if (nschema.eq.2) read(50) bxp,byp,bsp
!        endif
!      endif
!         if (ientree.eq.2) then
!            do i=1,icommence
!               read(40) bxo,byo
!               iturb=iturb+1
!            enddo
!!!!!!! ici
!            read(40) bxr,byr
!            iturb=iturb+1
!            read(40) bxp,byp
!            iturb=iturb+1
!         endif
!  endif
!         
   return
end subroutine lit
!
!********************************************************************
!
subroutine sauve_champ (ux,uy,uz,sp,cx,cy,cz,cs,nx,ny,nz)
!
!********************************************************************
!
   USE ecoulement_m
   USE sauvegarde_m
!
   implicit none
!
   
!
   real(8),dimension(nx,ny,nz) :: ux,uy,uz,sp
   real(8),dimension(ny,nz) :: cx,cy,cz,cs
!
  integer , intent(in) :: nx 
  integer , intent(in) :: ny 
  integer , intent(in) :: nz 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer :: k, j 
!-----------------------------------------------   
!
   if (nz.gt.1) then
      if (iscalaire.eq.0) then
         do k=1,nz
         do j=1,ny
            cx(j,k)=ux(nxboite,j,k)
            cy(j,k)=uy(nxboite,j,k)
            cz(j,k)=uz(nxboite,j,k)
         enddo
         enddo
         write(50) cx,cy,cz
      else
         do k=1,nz
         do j=1,ny
            cx(j,k)=ux(nxboite,j,k)
            cy(j,k)=uy(nxboite,j,k)
            cz(j,k)=uz(nxboite,j,k)
            cs(j,k)=sp(nxboite,j,k)
         enddo
         enddo
         write(50) cx,cy,cz,cs
      endif
   else
      if (iscalaire.eq.0) then
         do j=1,ny
            cx(j,1)=ux(nxboite,j,1)
            cy(j,1)=uy(nxboite,j,1)
         enddo
         write(50) cx,cy
      else
         do j=1,ny
            cx(j,1)=ux(nxboite,j,1)
            cy(j,1)=uy(nxboite,j,1)
            cs(j,1)=sp(nxboite,j,1)
         enddo
         write(50) cx,cy,cs
      endif
   endif
!
   return
end subroutine sauve_champ
!
!********************************************************************
!
subroutine sauve_stats(ux,uy,uz,sp,nx,ny,nz)
!
!********************************************************************
!
   USE paramx_m
!
   implicit none
!
  integer , intent(in) :: nx 
  integer , intent(in) :: ny 
  integer , intent(in) :: nz 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer :: i 
  real :: x 
!-----------------------------------------------
   real(8),dimension(nx,ny,nz) :: ux,uy,uz,sp
!
   x=0.
   do i=1,nx
      write(80,*) x,ux(i,ny/2,nz/2)
      x=x+dx
   enddo
!
   return
end subroutine sauve_stats
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
  integer  :: nx 
  integer  :: ny 
  integer  :: nz 
 !-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer :: k, j, i 
  real :: wzmin, wzmax 
!----------------------------------------------- 
   real(8),dimension(nx,ny,nz) :: ux,uy,uz,gx,gy,gz,wx,wy,wz,di
   real(8),dimension(ny,nz) :: sx
   real(8),dimension(nx,nz) :: sy
   real(8),dimension(nx,ny) :: sz
   real(8),dimension(nx) :: ffx,fsx,fwx,ffxp,fsxp,fwxp
   real(8),dimension(ny) :: ffy,fsy,fwy,ffyp,fsyp,fwyp
   real(8),dimension(nz) :: ffz,fsz,fwz,ffzp,fszp,fwzp
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
!*************************************************************
!
subroutine ouvre_fichier
!
!*************************************************************
!

   open(80,file='usuraxe.dat',form='formatted',&
       status='unknown')
!
   open (15,file='spectre.dat',status='unknown',&
        form='formatted')
!
   return
end subroutine ouvre_fichier
!
!*************************************************************
!
subroutine ferme_fichier
!
!*************************************************************
!
  implicit none
!
   close(80)
   close(15)
!
   return
end subroutine ferme_fichier
!
!*************************************************************
!
   subroutine moyt(uxmt,uymt,uzmt,uxux,uyuy,uzuz,&
        uxuy,uxuz,uyuz,ux,uy,uz,nx,ny,nz)
!
!*************************************************************
!
     implicit none
!
   integer :: nx,ny,nz,i,nxyz
   real(8),dimension(nx,ny,nz) :: uxmt,uymt,uzmt,uxux,uxuy,uxuz,uyuy,uyuz,uzuz
   real(8),dimension(nx,ny,nz) :: ux,uy,uz
!
      nxyz=nx*ny*nz
!
      do i=1,nxyz
         uxmt(i,1,1)=uxmt(i,1,1)+ux(i,1,1)
         uymt(i,1,1)=uymt(i,1,1)+uy(i,1,1)
         uzmt(i,1,1)=uzmt(i,1,1)+uz(i,1,1)
         uxux(i,1,1)=uxux(i,1,1)+ux(i,1,1)*ux(i,1,1)
         uyuy(i,1,1)=uyuy(i,1,1)+uy(i,1,1)*uy(i,1,1)
         uzuz(i,1,1)=uzuz(i,1,1)+uz(i,1,1)*uz(i,1,1)
         uxuy(i,1,1)=uxuy(i,1,1)+ux(i,1,1)*uy(i,1,1)
         uxuz(i,1,1)=uxuz(i,1,1)+ux(i,1,1)*uz(i,1,1)
         uyuz(i,1,1)=uyuz(i,1,1)+uy(i,1,1)*uz(i,1,1)
      enddo
!
      return
    end subroutine moyt
!
!********************************************************************
!
   subroutine sauve (ux,uy,uz,py1,pz1,gx,gy,gz,hx,hy,hz,sp&
        ,ppm,nx,ny,nz,nxm,nym,nzm,itime)
!
!*******************************************************************
!
     USE ecoulement_m
!
     implicit none
!
     real(8),dimension(nx,ny,nz) :: py1,ux,uy,uz,gx,gy,gz,hx,hy,hz,sp,sg,sh,pz1
     real(8),dimension(nxm,nym,nzm) :: ppm
     integer :: nx,ny,nz,itime,nxm,nym,nzm

     if (itime.eq.1) then
        open(10,file='debut_swirl.dat',&
             form='unformatted',status='unknown')
     else
        open(10,file='sauve',&
             form='unformatted',status='unknown')
     endif
     if (nz.gt.1) then
        if (iscalaire.eq.0) then
           write(10) ux,uy,uz,py1,pz1,ppm,gx,gy,gz
        else
           write(10) ux,uy,uz,sp
        endif
     else
        if (iscalaire.eq.0) then
           write(10) ux,uy,py1,ppm,gx,gy
        else
           write(10) ux,uy,sp
        endif
     endif
     close(10)
!
     if (nschema.eq.1) then
        open(20,file='sauve'//'g',&
             form='unformatted',status='unknown')
        open(30,file='sauve'//'h',&
             form='unformatted',status='unknown')
        if (nz.gt.1) then
           if (iscalaire.eq.0) then
              write(20) gx,gy,gz!
              write(30) hx,hy,hz
           else
              write(20) gx,gy,gz,sg
              write(30) hx,hy,hz,sh
           endif
        else
           if (iscalaire.eq.0) then
              write(20) gx,gy
              write(30) hx,hy
           else
              write(20) gx,gy,sg
              write(30) hx,hy,sh
           endif
        endif
        close(20)
        close(30)
     endif
     return
   end subroutine sauve
!
