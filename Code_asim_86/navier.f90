!********************************************************************
!
subroutine datosinc3d (ffx,fsx,fwx,ffy,fsy,fwy,ffxp,fsxp,fwxp,ffyp,fsyp,fwyp,&
             ffx_r,fsx_r,fwx_r,ffy_r,fsy_r,fwy_r,&
             ffxp_r,fsxp_r,fwxp_r,ffyp_r,fsyp_r,fwyp_r,&
             ffx_rg,fsx_rg,fwx_rg,ffy_rg,fsy_rg,fwy_rg,&
             ffxp_rg,fsxp_rg,fwxp_rg,ffyp_rg,fsyp_rg,fwyp_rg,&
             sfx,ssx,swx,sfy,ssy,swy,sfxp,ssxp,swxp,sfyp,ssyp,swyp,&
             cfx6,csx6,cwx6,cfxp6,csxp6,cwxp6,&
             cifx6,cisx6,ciwx6,cifxp6,cisxp6,ciwxp6,&
             cfy6,csy6,cwy6,cfyp6,csyp6,cwyp6,&
             cify6,cisy6,ciwy6,cifyp6,cisyp6,ciwyp6,&   
             cfip6,csip6,cwip6,cfip6y,csip6y,cwip6y,&
             cifip6,cisip6,ciwip6,cifip6y,cisip6y,ciwip6y,&
             fifx,ficx,fibx,fiffx,fibbx,fiz1x,fiz2x,&
             filax,filaxp,&
             fifxp,ficxp,fibxp,fiffxp,fibbxp,&
             fify,ficy,fiby,fiffy,fibby,fiz1y,fiz2y,&
             filay,filayp,&
             fifyp,ficyp,fibyp,fiffyp,fibbyp,&
             a_tot,al_tot,indx_tot,a2_tot,al1_tot,indx1_tot,d5,&
             pp4y,adt,bdt,gdt,&
             nx,ny,nz,nxm,nym,mx,my,mz,m1,m2,zero,nxr,nyr,myr,rdxy,nxrg,nyrg,myrg)
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
   USE aleatoire_m
   USE pwtcom_m
! 
   implicit none
!
   real(8),dimension(nx) :: ffx,fcx,fbx,sfx,scx,sbx,fsx,fwx,ssx,swx
   real(8),dimension(nx) :: ffxp,fsxp,fwxp,sfxp,ssxp,swxp
   real(8),dimension(ny) :: ffy,fcy,fby,sfy,scy,sby,fsy,fwy,ssy,swy
   real(8),dimension(ny) :: ffyp,fsyp,fwyp,sfyp,ssyp,swyp
   real(8),dimension(nz) :: ffz,fcz,fbz,sfz,scz,sbz,fsz,fwz,ssz,swz
   real(8),dimension(nz) :: ffzp,fszp,fwzp,sfzp,sszp,swzp
   real(8),dimension(nxr) :: ffx_r,fcx_r,fbx_r,fsx_r,fwx_r 
   real(8),dimension(nxr) :: ffxp_r,fsxp_r,fwxp_r
   real(8),dimension(nyr) :: ffy_r,fcy_r,fby_r,fsy_r,fwy_r 
   real(8),dimension(nyr) :: ffyp_r,fsyp_r,fwyp_r
   real(8),dimension(nxrg) :: ffx_rg,fcx_rg,fbx_rg,fsx_rg,fwx_rg 
   real(8),dimension(nxrg) :: ffxp_rg,fsxp_rg,fwxp_rg
   real(8),dimension(nyrg) :: ffy_rg,fcy_rg,fby_rg,fsy_rg,fwy_rg 
   real(8),dimension(nyrg) :: ffyp_rg,fsyp_rg,fwyp_rg
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
!  ********************* Filtrage *********************
   real(8),dimension(nx) :: fifx,ficx,fibx,fiffx,fibbx,fiz1x,fiz2x
   real(8),dimension(nx,2) :: filax,filaxp
   real(8),dimension(nx) :: fifxp,ficxp,fibxp,fiffxp,fibbxp
   real(8),dimension(ny) :: fify,ficy,fiby,fiffy,fibby,fiz1y,fiz2y
   real(8),dimension(ny,2) :: filay,filayp
   real(8),dimension(ny) :: fifyp,ficyp,fibyp,fiffyp,fibbyp
   real(8),dimension(nz) :: fifz,ficz,fibz,fiffz,fibbz,fiz1z,fiz2z
   real(8),dimension(nz,2) :: filaz,filazp
   real(8),dimension(nz) :: fifzp,ficzp,fibzp,fiffzp,fibbzp
!  ********************* STRET *********************
   real(8),dimension(ny) :: ppy,pp2y,pp3y,pp4y,fux1,d1,d,hprime
   real(8),dimension(nyr) :: ppy_r
   real(8),dimension(nyrg) :: ppy_rg
   real(8),dimension(nym) :: ppyi,pp2yi,pp3yi,pp4yi,hprimei
   real(8),dimension(my) :: yp,ypi,yeta,yetai
   real(8),dimension(myr) :: yeta_r
   real(8),dimension(myrg) :: yeta_rg
   real(8),dimension(2) :: ja,jb
   real(8),dimension(nx) :: d5
   real(8),dimension(ny/2,m1+m2+1) :: a,a1,a2
   real(8),dimension(ny/2,m1) :: al,al1
   real(8),dimension(ny/2) :: indx,indx1,e,c
   real(8),dimension(ny/2,m1+m2+1,nxm) :: a_tot,a2_tot
   real(8),dimension(ny/2,m1,nxm) :: al_tot,al1_tot
   real(8),dimension(ny/2,nxm) :: indx_tot,indx1_tot
   real(8) :: pi,yinf,beta,den,xnum,alpha,xcx,den1,den2,den3,den4,xnum1,cst
!  *************************************************
   real(8) :: adt,bdt,cdt,gdt,dum,zero
   integer :: itime,k,j,i
   integer :: nx,ny,nz,nxm,nym,mx,my,mz,m1,m2,nxr,nyr,myr,rdxy,nxrg,nyrg,myrg
   real(8),external :: epmach
!
!
!  Declaration de format
1000 format(4f12.8)
!
   call parametre(nx,ny,nz,nxr,nyr,rdxy)
!
!  determination du zero machine
   dum=0.
   zero=epmach(dum)*100
   print *,'zero = ',zero
!
!  ***************************************************************
!  **************** fonction de stretching x=h(s) ****************
!  ***************************************************************
!
      pi=acos(-1.)
      yinf=-yly/2.  ! borne inferior du domaine physique dans la direction inhomogene
      beta=12.5
      den=2.*beta*yinf
      xnum=-yinf-sqrt(pi*pi*beta*beta+yinf*yinf)
      alpha=abs(xnum/den)  ! ec. (3.71) tesis sylvain
      xcx=1./beta/alpha
!
!     yeta(j) ----> s : variable de l'espace de calcul, point principaux
!     yp(j) ----> x : variable physique, point principaux
!     ppy(j) ----> h' : dh/ds, point principaux
! 
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
            if (yeta(j).lt.0.5) yp(j)=xnum1-cst-yinf  ! ec. (3.67) tesis sylvain
            if (yeta(j).eq.0.5) yp(j)=0.-yinf
            if (yeta(j).gt.0.5) yp(j)=xnum1+cst-yinf
         enddo
      endif
!
      if (alpha.eq.0.) then
         yp(1)=-1.e10
         do j=2,ny
            yeta(j)=(j-1.)*(1./ny)
            yp(j)=-beta*cos(pi*yeta(j))/sin(yeta(j)*pi)  ! ec. (3.66) tesis sylvain
         enddo
      endif
!
      do j=1,ny
         ppy(j)=1./(alpha/pi+(1./pi/beta)*sin(pi*yeta(j))* &
            sin(pi*yeta(j)))  ! ec. (3.63) tesis sylvain
         pp2y(j)=ppy(j)*ppy(j)
         pp4y(j)=pp2y(j)*(-2./beta*cos(pi*yeta(j))*sin(pi*yeta(j)))
         hprime(j)=1./ppy(j)
      enddo
!
      open (15,file='stret.dat',form='formatted')
      do j=1,ny
         write(15,1000) yp(j)
      enddo
      close(15)
!
!     yetai(j) ----> s : variable de l'espace de calcul, point secondaires (grille decalee)
!     ypi(j) ----> x : variable physique, point secondaires
!     ppyi(j) ----> h' : dh/ds, point secondaires
!
      if (alpha.ne.0.) then 
         do j=1,ny
            if (ncly.eq.0) yetai(j)=(j-0.5)*(1./ny)
            if (ncly.eq.1) yetai(j)=(j-0.5)*(1./(ny-1.))
            den1=sqrt(alpha*beta+1.)
            xnum=den1/sqrt(alpha/pi)/sqrt(beta)/sqrt(pi)     
            den=2.*sqrt(alpha/pi)*sqrt(beta)*pi*sqrt(pi)
            den3=((sin(pi*yetai(j)))*(sin(pi*yetai(j)))/beta/pi)+alpha/pi
            den4=2.*alpha*beta-cos(2.*pi*yetai(j))+1.
            xnum1=(atan(xnum*tan(pi*yetai(j))))*den4/den1/den3/den
            cst=sqrt(beta)*pi/(2.*sqrt(alpha)*sqrt(alpha*beta+1.))
            if (yetai(j).lt.0.5) ypi(j)=xnum1-cst-yinf
            if (yetai(j).eq.0.5) ypi(j)=0.-yinf
            if (yetai(j).gt.0.5) ypi(j)=xnum1+cst-yinf
         enddo
      endif
!
      if (alpha.eq.0.) then
         ypi(1)=-1.e10
         do j=2,ny
            yetai(j)=(j-1.)*(1./ny)
            ypi(j)=-beta*cos(pi*yetai(j))/sin(yetai(j)*pi)
         enddo
      endif      
!
      do j=1,nym
         ppyi(j)=1./(alpha/pi+(1./pi/beta)*sin(pi*yetai(j))* &
            sin(pi*yetai(j)))
         pp2yi(j)=ppyi(j)*ppyi(j)
         pp4yi(j)=pp2yi(j)*(-2./beta*cos(pi*yetai(j))*sin(pi*yetai(j)))
         hprimei(j)=1./ppyi(j)
      enddo
!
!  ***************************************************************
!  ********************* determina ppy_r(j) *********************
!  ***************************************************************
!
      yinf=-ylyr/2.  ! borne inferior du domaine physique dans la direction inhomogene
      den=2.*beta*yinf
      xnum=-yinf-sqrt(pi*pi*beta*beta+yinf*yinf)
      alpha=abs(xnum/den)  ! ec. (3.71) tesis sylvain
      xcx=1./beta/alpha
!
!     yeta_r(j) ----> s : variable de l'espace de calcul, point principaux
!     ppy_r(j) ----> h' : dh/ds, point principaux
! 
      if (alpha.ne.0.) then 
         do j=1,nyr
            if (ncly.eq.0) yeta_r(j)=(j-1.)*(1./nyr)
            if (ncly.eq.1) yeta_r(j)=(j-1.)*(1./(nyr-1.))
         enddo
      endif
!
      if (alpha.eq.0.) then
         do j=2,nyr
            yeta_r(j)=(j-1.)*(1./nyr)
         enddo
      endif
!
      do j=1,nyr
         ppy_r(j)=1./(alpha/pi+(1./pi/beta)*sin(pi*yeta_r(j))* &
            sin(pi*yeta_r(j)))  ! ec. (3.63) tesis sylvain
      enddo
!
!  ***************************************************************
!  ********************* determina ppy_rg(j) *********************
!  ***************************************************************
!
      yinf=-ylyr/2.  ! borne inferior du domaine physique dans la direction inhomogene
      den=2.*beta*yinf
      xnum=-yinf-sqrt(pi*pi*beta*beta+yinf*yinf)
      alpha=abs(xnum/den)  ! ec. (3.71) tesis sylvain
      xcx=1./beta/alpha
!
!     yeta_rg(j) ----> s : variable de l'espace de calcul, point principaux
!     ppy_rg(j) ----> h' : dh/ds, point principaux
! 
      if (alpha.ne.0.) then 
         do j=1,nyrg
            if (ncly.eq.0) yeta_rg(j)=(j-1.)*(1./nyrg)
            if (ncly.eq.1) yeta_rg(j)=(j-1.)*(1./(nyrg-1.))
         enddo
      endif
!
      if (alpha.eq.0.) then
         do j=2,nyrg
            yeta_rg(j)=(j-1.)*(1./nyrg)
         enddo
      endif
!
      do j=1,nyrg
         ppy_rg(j)=1./(alpha/pi+(1./pi/beta)*sin(pi*yeta_rg(j))* &
            sin(pi*yeta_rg(j)))  ! ec. (3.63) tesis sylvain
      enddo
!
!  ***************************************************************
!  ***************************************************************
!  ***************************************************************
!
   call filtres(fiffx,fifx,ficx,fibx,fibbx,filax,&
        fiffxp,fifxp,ficxp,fibxp,fibbxp,filaxp,&
        fiffy,fify,ficy,fiby,fibby,filay,&
        fiffyp,fifyp,ficyp,fibyp,fibbyp,filayp,&
        fiffz,fifz,ficz,fibz,fibbz,filaz,&
        fiffzp,fifzp,ficzp,fibzp,fibbzp,filazp,&
        fiz1x,fiz2x,fiz1y,fiz2y,fiz1z,fiz2z,&
        nx,ny,nz)
!
!  coef. esquemas para dominio reconstruido
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
        ciwi6y,ppy,pp2y,ppyi,pp2yi,yp,ypi)
!
!  coef. esquemas para dominio reducido fino
   call schemas_r(ffx_r,fcx_r,fbx_r,ffy_r,fcy_r,fby_r,&
        fsx_r,fwx_r,fsy_r,fwy_r,&
        ffxp_r,fsxp_r,fwxp_r,ffyp_r,fsyp_r,fwyp_r,&
        nxr,nyr,ppy_r)
!
!  coef. esquemas para dominio reducido grosero
   call schemas_rg(ffx_rg,fcx_rg,fbx_rg,ffy_rg,fcy_rg,fby_rg,&
        fsx_rg,fwx_rg,fsy_rg,fwy_rg,&
        ffxp_rg,fsxp_rg,fwxp_rg,ffyp_rg,fsyp_rg,fwyp_rg,&
        nxrg,nyrg,ppy_rg)
!
!  Construye una descomposicion LU de la matriz de coef. B=A^2 
   call construye_matriz_B(a_tot,al_tot,indx_tot,a2_tot,al1_tot,indx1_tot,d5,&
        alpha,beta,m1,m2,nx,nxm,ny,nym,nz,mx,my,mz)
!
   adt=1.5*dt
   bdt=-0.5*dt
   cdt=0.*dt
   gdt=adt+bdt+cdt
!
   return
end subroutine datosinc3d
!
!********************************************************************
!
!  Original code Y=F(X0)
!
!  Entrada: ux, uy, gx, gy, bx1, by1 
!  Salida: uxn, uyn, g8x, g8y
subroutine inc3dorigin (ux,uy,gx,gy,&
             uxn,uyn,g8x,g8y,&
             bx1,by1,itime,&
             ffx,fsx,fwx,ffy,fsy,fwy,ffxp,fsxp,fwxp,ffyp,fsyp,fwyp,&
             sfx,ssx,swx,sfy,ssy,swy,sfxp,ssxp,swxp,sfyp,ssyp,swyp,&
             cfx6,csx6,cwx6,cfxp6,csxp6,cwxp6,&
             cifx6,cisx6,ciwx6,cifxp6,cisxp6,ciwxp6,&
             cfy6,csy6,cwy6,cfyp6,csyp6,cwyp6,&
             cify6,cisy6,ciwy6,cifyp6,cisyp6,ciwyp6,&   
             cfip6,csip6,cwip6,cfip6y,csip6y,cwip6y,&
             cifip6,cisip6,ciwip6,cifip6y,cisip6y,ciwip6y,&
             a_tot,al_tot,indx_tot,a2_tot,al1_tot,indx1_tot,d5,&
             pp4y,adt,bdt,gdt,&
             nx,ny,nz,nxm,nym,mx,my,mz,m1,m2,ditime)
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
   USE aleatoire_m
   USE pwtcom_m
! 
   implicit none
!
   real(8),dimension(nx,ny,nz) :: ux,uy,gx,gy
   real(8),dimension(nx,ny,nz) :: uxn,uyn,u8x,u8y,g8x,g8y
   real(8),dimension(mx,my) :: ps2d1,ps2d2,ps2d3,us2d1,us2d2,us2d3
   real(8),dimension(ny,ditime) :: bx1,by1
   real(8),dimension(nx) :: ffx,fsx,fwx,ffxp,fsxp,fwxp
   real(8),dimension(ny) :: ffy,fsy,fwy,ffyp,fsyp,fwyp
   real(8),dimension(nx) :: sfx,ssx,swx,sfxp,ssxp,swxp
   real(8),dimension(ny) :: sfy,ssy,swy,sfyp,ssyp,swyp
   real(8),dimension(nxm) :: cfx6,csx6,cwx6,cifx6,cisx6,ciwx6 
   real(8),dimension(nxm) :: cifxp6,cisxp6,ciwxp6,cfxp6,csxp6,cwxp6
   real(8),dimension(nym) :: cfy6,csy6,cwy6,cify6,cisy6,ciwy6 
   real(8),dimension(nym) :: cifyp6,cisyp6,ciwyp6,cfyp6,csyp6,cwyp6
   real(8),dimension(nx) :: cfip6,csip6,cwip6
   real(8),dimension(nx) :: cifip6,cisip6,ciwip6
   real(8),dimension(ny) :: cfip6y,csip6y,cwip6y
   real(8),dimension(ny) :: cifip6y,cisip6y,ciwip6y
!  ********************* STRET *********************
   real(8),dimension(ny) :: pp4y
   real(8),dimension(nx) :: d5
   real(8),dimension(ny/2,m1+m2+1,nxm) :: a_tot,a2_tot
   real(8),dimension(ny/2,m1,nxm) :: al_tot,al1_tot
   real(8),dimension(ny/2,nxm) :: indx_tot,indx1_tot
!  *************************************************
   real(8) :: adt,bdt,gdt,t
   integer :: itime
   integer :: nx,ny,nz,nxm,nym,mx,my,mz,m1,m2,ditime
!
!
   print *,'itime = ',itime
!
!  Entrada: ux, uy, gx, gy. Salida: u8x, u8y, g8x, g8y, ps2d1
   call funcion1 (ux,uy,gx,gy,u8x,u8y,g8x,g8y,ps2d1,&
	  ffx,fsx,fwx,ffy,fsy,fwy,&
	  ffxp,fsxp,fwxp,&
	  ffyp,fsyp,fwyp,&
	  sfx,ssx,swx,sfy,ssy,swy,&
	  sfxp,ssxp,swxp,&
	  sfyp,ssyp,swyp,&
	  cfx6,csx6,cwx6,&
	  cfxp6,csxp6,cwxp6,&
	  cifx6,cisx6,ciwx6,cifxp6,cisxp6,ciwxp6,&
	  cfy6,csy6,cwy6,cfyp6,csyp6,cwyp6,&
	  cify6,cisy6,ciwy6,cifyp6,cisyp6,ciwyp6,&    
	  pp4y,bx1,by1,itime,ditime,nx,ny,nz,nxm,nym,mx,my,mz,adt,bdt,gdt)
! 
!  Entrada: ps2d1. Salida: ps2d2
   call funcion2(ps2d1,ps2d2,mx,my,nxm,nym)
! 
!  Entrada: ps2d2. Salida: ps2d3
   call funcion3(ps2d2,ps2d3,mx,my,nx,ny,nxm,nym)
! 
!  Entrada: ps2d3. Salida: us2d1
   call funcion4(ps2d3,us2d1,&
	  a_tot,al_tot,indx_tot,a2_tot,al1_tot,indx1_tot,d5,&
	  m1,m2,nx,nxm,ny,nym,mx,my)
! 
!  Entrada: us2d1. Salida: us2d2
   call funcion5(us2d1,us2d2,mx,my,nx,ny,nxm,nym)
! 
!  Entrada: us2d2. Salida: us2d3
   call funcion6(us2d2,us2d3,mx,my,nxm,nym)
!
!  Entrada: us2d3, u8x, u8y. Salida: uxn, uyn
   call funcion7(us2d3,u8x,u8y,uxn,uyn,&
	  cfip6,csip6,cwip6,cfip6y,csip6y,cwip6y,&
	  cifip6,cisip6,ciwip6,cifip6y,cisip6y,ciwip6y,&
	  nx,nxm,mx,ny,nym,my,nz)
!
!  Entrada: uxn, uyn. Imprime DIV U
   call divU (uxn,uyn,&
	  cfx6,csx6,cwx6,&
	  cfxp6,csxp6,cwxp6,&
	  cifx6,cisx6,ciwx6,cifxp6,cisxp6,ciwxp6,&
	  cfy6,csy6,cwy6,cfyp6,csyp6,cwyp6,&
	  cify6,cisy6,ciwy6,cifyp6,cisyp6,ciwyp6,&
	  nx,ny,nz,nxm,nym)
!                
   return
end subroutine inc3dorigin
!
!********************************************************************
!
!  Tangent code Yd=F_D(X0,X0d)
!
!  Entrada: ux, uy, gx, gy, bx1, by1, uxd, uyd, gxd, gyd, bx1d, by1d
!  Salida: uxnd, uynd, g8xd, g8yd
subroutine inc3dtangent (ux,uy,gx,gy,&
             uxd,uyd,gxd,gyd,&
             uxnd,uynd,g8xd,g8yd,&
             bx1,by1,bx1d,by1d,itime,&
             ffx,fsx,fwx,ffy,fsy,fwy,ffxp,fsxp,fwxp,ffyp,fsyp,fwyp,&
             sfx,ssx,swx,sfy,ssy,swy,sfxp,ssxp,swxp,sfyp,ssyp,swyp,&
             cfx6,csx6,cwx6,cfxp6,csxp6,cwxp6,&
             cifx6,cisx6,ciwx6,cifxp6,cisxp6,ciwxp6,&
             cfy6,csy6,cwy6,cfyp6,csyp6,cwyp6,&
             cify6,cisy6,ciwy6,cifyp6,cisyp6,ciwyp6,&   
             cfip6,csip6,cwip6,cfip6y,csip6y,cwip6y,&
             cifip6,cisip6,ciwip6,cifip6y,cisip6y,ciwip6y,&
             a_tot,al_tot,indx_tot,a2_tot,al1_tot,indx1_tot,d5,&
             pp4y,adt,bdt,gdt,&
             nx,ny,nz,nxm,nym,mx,my,mz,m1,m2,ditime)
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
   USE aleatoire_m
   USE pwtcom_m
! 
   implicit none
!
   real(8),dimension(nx,ny,nz) :: ux,uy,gx,gy
   real(8),dimension(nx,ny,nz) :: uxn,uyn,u8x,u8y,g8x,g8y
   real(8),dimension(mx,my) :: ps2d1
   real(8),dimension(nx,ny,nz) :: uxd,uyd,gxd,gyd
   real(8),dimension(nx,ny,nz) :: uxnd,uynd,u8xd,u8yd,g8xd,g8yd
   real(8),dimension(mx,my) :: ps2d1d,ps2d2d,ps2d3d,us2d1d,us2d2d,us2d3d
   real(8),dimension(ny,ditime) :: bx1,by1,bx1d,by1d
   real(8),dimension(nx) :: ffx,fsx,fwx,ffxp,fsxp,fwxp
   real(8),dimension(ny) :: ffy,fsy,fwy,ffyp,fsyp,fwyp
   real(8),dimension(nx) :: sfx,ssx,swx,sfxp,ssxp,swxp
   real(8),dimension(ny) :: sfy,ssy,swy,sfyp,ssyp,swyp
   real(8),dimension(nxm) :: cfx6,csx6,cwx6,cifx6,cisx6,ciwx6 
   real(8),dimension(nxm) :: cifxp6,cisxp6,ciwxp6,cfxp6,csxp6,cwxp6
   real(8),dimension(nym) :: cfy6,csy6,cwy6,cify6,cisy6,ciwy6 
   real(8),dimension(nym) :: cifyp6,cisyp6,ciwyp6,cfyp6,csyp6,cwyp6
   real(8),dimension(nx) :: cfip6,csip6,cwip6
   real(8),dimension(nx) :: cifip6,cisip6,ciwip6
   real(8),dimension(ny) :: cfip6y,csip6y,cwip6y
   real(8),dimension(ny) :: cifip6y,cisip6y,ciwip6y
!  ********************* STRET *********************
   real(8),dimension(ny) :: pp4y
   real(8),dimension(nx) :: d5
   real(8),dimension(ny/2,m1+m2+1,nxm) :: a_tot,a2_tot
   real(8),dimension(ny/2,m1,nxm) :: al_tot,al1_tot
   real(8),dimension(ny/2,nxm) :: indx_tot,indx1_tot
!  *************************************************
   real(8) :: adt,bdt,gdt,t
   integer :: itime
   integer :: nx,ny,nz,nxm,nym,mx,my,mz,m1,m2,ditime
!
!
   print *,'itime = ',itime
!
!  Entrada: ux, uy, gx, gy, bx1, by1, uxd, uyd, gxd, gyd, bx1d, by1d
!  Salida: u8xd, u8yd, g8xd, g8yd, ps2d1d
   call FUNCION1_D(ux,uxd,uy,uyd,gx,gxd,gy,gyd,bx1,bx1d,by1,by1d,&
             u8x,u8xd,u8y,u8yd,g8x,g8xd,g8y,g8yd,ps2d1,ps2d1d,&
             ffx,fsx,fwx,ffy,fsy,fwy,&
             ffxp,fsxp,fwxp,&
             ffyp,fsyp,fwyp,&
             sfx,ssx,swx,sfy,ssy,swy,&
             sfxp,ssxp,swxp,&
             sfyp,ssyp,swyp,&
             cfx6,csx6,cwx6,&
             cfxp6,csxp6,cwxp6,&
             cifx6,cisx6,ciwx6,cifxp6,cisxp6,ciwxp6,&
             cfy6,csy6,cwy6,cfyp6,csyp6,cwyp6,&
             cify6,cisy6,ciwy6,cifyp6,cisyp6,ciwyp6,&    
             pp4y,itime,ditime,nx,ny,nz,nxm,nym,mx,my,mz,adt,bdt,gdt)
!  
!  Entrada: ps2d1d. Salida: ps2d2d
   call funcion2(ps2d1d,ps2d2d,mx,my,nxm,nym)
!  
!  Entrada: ps2d2d. Salida: ps2d3d
   call FUNCION3_D(ps2d2d,ps2d3d,mx,my,nx,ny,nxm,nym)
!  
!  Entrada: ps2d3d. Salida: us2d1d
   call FUNCION4_D(ps2d3d,us2d1d,a_tot,al_tot,indx_tot,&
             a2_tot,al1_tot,indx1_tot,d5,m1,m2,nx,nxm,ny,nym,mx,my)
!  
!  Entrada: us2d1d. Salida: us2d2d
   call FUNCION5_D(us2d1d,us2d2d,mx,my,nx,ny,nxm,nym)
! 
!  Entrada: us2d2d. Salida: us2d3d
   call funcion6(us2d2d,us2d3d,mx,my,nxm,nym)
!
!  Entrada: us2d3d, u8xd, u8yd. Salida: uxnd, uynd
   call FUNCION7_D(us2d3d,u8xd,u8yd,uxnd,uynd,&
             cfip6,csip6,cwip6,cfip6y,csip6y,cwip6y,cifip6,cisip6,&
             ciwip6,cifip6y,cisip6y,ciwip6y,nx,nxm,mx,ny,nym,my,nz)
!                
   return
end subroutine inc3dtangent
!
!********************************************************************
!
!  Adjoint code X0b=F_B(X0,Yb)
!
!  Entrada: ux, uy, gx, gy, bx1, by1, uxnb, uynb, g8xb, g8yb
!  Salida: uxb, uyb, gxb, gyb, bx1b, by1b
subroutine inc3dadjoint (ux,uy,gx,gy,&
             uxnb,uynb,g8xb,g8yb,&
             uxb,uyb,gxb,gyb,&
             bx1,by1,bx1b,by1b,itime,&
             ffx,fsx,fwx,ffy,fsy,fwy,ffxp,fsxp,fwxp,ffyp,fsyp,fwyp,&
             sfx,ssx,swx,sfy,ssy,swy,sfxp,ssxp,swxp,sfyp,ssyp,swyp,&
             cfx6,csx6,cwx6,cfxp6,csxp6,cwxp6,&
             cifx6,cisx6,ciwx6,cifxp6,cisxp6,ciwxp6,&
             cfy6,csy6,cwy6,cfyp6,csyp6,cwyp6,&
             cify6,cisy6,ciwy6,cifyp6,cisyp6,ciwyp6,&   
             cfip6,csip6,cwip6,cfip6y,csip6y,cwip6y,&
             cifip6,cisip6,ciwip6,cifip6y,cisip6y,ciwip6y,&
             a_tot,al_tot,indx_tot,a2_tot,al1_tot,indx1_tot,d5,&
             pp4y,adt,bdt,gdt,&
             nx,ny,nz,nxm,nym,mx,my,mz,m1,m2,ditime)
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
   USE aleatoire_m
   USE pwtcom_m
! 
   implicit none
!
   real(8),dimension(nx,ny,nz) :: ux,uy,gx,gy
   real(8),dimension(nx,ny,nz) :: uxn,uyn,u8x,u8y,g8x,g8y
   real(8),dimension(mx,my) :: ps2d1
   real(8),dimension(nx,ny,nz) :: uxb,uyb,gxb,gyb
   real(8),dimension(nx,ny,nz) :: uxnb,uynb,u8xb,u8yb,g8xb,g8yb
   real(8),dimension(mx,my) :: ps2d1b,ps2d2b,ps2d3b,us2d1b,us2d2b,us2d3b
   real(8),dimension(ny,ditime) :: bx1,by1,bx1b,by1b
   real(8),dimension(nx) :: ffx,fsx,fwx,ffxp,fsxp,fwxp
   real(8),dimension(ny) :: ffy,fsy,fwy,ffyp,fsyp,fwyp
   real(8),dimension(nx) :: sfx,ssx,swx,sfxp,ssxp,swxp
   real(8),dimension(ny) :: sfy,ssy,swy,sfyp,ssyp,swyp
   real(8),dimension(nxm) :: cfx6,csx6,cwx6,cifx6,cisx6,ciwx6 
   real(8),dimension(nxm) :: cifxp6,cisxp6,ciwxp6,cfxp6,csxp6,cwxp6
   real(8),dimension(nym) :: cfy6,csy6,cwy6,cify6,cisy6,ciwy6 
   real(8),dimension(nym) :: cifyp6,cisyp6,ciwyp6,cfyp6,csyp6,cwyp6
   real(8),dimension(nx) :: cfip6,csip6,cwip6
   real(8),dimension(nx) :: cifip6,cisip6,ciwip6
   real(8),dimension(ny) :: cfip6y,csip6y,cwip6y
   real(8),dimension(ny) :: cifip6y,cisip6y,ciwip6y
!  ********************* STRET *********************
   real(8),dimension(ny) :: pp4y
   real(8),dimension(nx) :: d5
   real(8),dimension(ny/2,m1+m2+1,nxm) :: a_tot,a2_tot
   real(8),dimension(ny/2,m1,nxm) :: al_tot,al1_tot
   real(8),dimension(ny/2,nxm) :: indx_tot,indx1_tot
!  *************************************************
   real(8) :: adt,bdt,gdt,t,normAdj
   integer :: itime,j,i
   integer :: nx,ny,nz,nxm,nym,mx,my,mz,m1,m2,ditime
!
!
   print *,'itime = ',itime
!
!  Entrada: uxnb, uynb. Salida: us2d3b, u8xb, u8yb
   call FUNCION7_B(us2d3b,u8xb,u8yb,uxnb,uynb,&
             cfip6,csip6,cwip6,cfip6y,csip6y,cwip6y,cifip6,cisip6,&
             ciwip6,cifip6y,cisip6y,ciwip6y,nx,nxm,mx,ny,nym,my,nz)
! 
!  Entrada: us2d3b. Salida: us2d2b
   call funcion6inv(us2d3b,us2d2b,mx,my,nxm,nym)
!  
!  Entrada: us2d2b. Salida: us2d1b
   call FUNCION5_B(us2d1b,us2d2b,mx,my,nx,ny,nxm,nym)
!  
!  Entrada: us2d1b. Salida: ps2d3b
   call FUNCION4_B(ps2d3b,us2d1b,a_tot,al_tot,indx_tot,&
             a2_tot,al1_tot,indx1_tot,d5,m1,m2,nx,nxm,ny,nym,mx,my)
! 
!  Entrada: ps2d3b. Salida: ps2d2b
   call FUNCION3_B(ps2d2b,ps2d3b,mx,my,nx,ny,nxm,nym)
! 
!  Entrada: ps2d2b. Salida: ps2d1b
   call funcion2inv(ps2d2b,ps2d1b,mx,my,nxm,nym)
!
!  Entrada: ux, uy, gx, gy, bx1, by1, u8xb, u8yb, g8xb, g8yb, ps2d1b
!  Salida: uxb, uyb, gxb, gyb, bx1b, by1b
   call FUNCION1_B(ux,uxb,uy,uyb,gx,gxb,gy,gyb,bx1,bx1b,by1,by1b,&
             u8x,u8xb,u8y,u8yb,g8x,g8xb,g8y,g8yb,ps2d1,ps2d1b,&
             ffx,fsx,fwx,ffy,fsy,fwy,&
             ffxp,fsxp,fwxp,&
             ffyp,fsyp,fwyp,&
             sfx,ssx,swx,sfy,ssy,swy,&
             sfxp,ssxp,swxp,&
             sfyp,ssyp,swyp,&
             cfx6,csx6,cwx6,&
             cfxp6,csxp6,cwxp6,&
             cifx6,cisx6,ciwx6,cifxp6,cisxp6,ciwxp6,&
             cfy6,csy6,cwy6,cfyp6,csyp6,cwyp6,&
             cify6,cisy6,ciwy6,cifyp6,cisyp6,ciwyp6,&    
             pp4y,itime,ditime,nx,ny,nz,nxm,nym,mx,my,mz,adt,bdt,gdt)
!                
   return
end subroutine inc3dadjoint
!
!********************************************************************
!
!  Original code Y=F(X0) modificado para calcular campos de presion
!
!  Entrada: ux, uy, gx, gy, bx1, by1 
!  Salida: uxn, uyn, ppn, g8x, g8y
subroutine inc3dorigin_tray (ux,uy,gx,gy,&
             uxn,uyn,ppn,g8x,g8y,&
             bx1,by1,itime,&
             ffx,fsx,fwx,ffy,fsy,fwy,ffxp,fsxp,fwxp,ffyp,fsyp,fwyp,&
             sfx,ssx,swx,sfy,ssy,swy,sfxp,ssxp,swxp,sfyp,ssyp,swyp,&
             cfx6,csx6,cwx6,cfxp6,csxp6,cwxp6,&
             cifx6,cisx6,ciwx6,cifxp6,cisxp6,ciwxp6,&
             cfy6,csy6,cwy6,cfyp6,csyp6,cwyp6,&
             cify6,cisy6,ciwy6,cifyp6,cisyp6,ciwyp6,&   
             cfip6,csip6,cwip6,cfip6y,csip6y,cwip6y,&
             cifip6,cisip6,ciwip6,cifip6y,cisip6y,ciwip6y,&
             a_tot,al_tot,indx_tot,a2_tot,al1_tot,indx1_tot,d5,&
             pp4y,adt,bdt,gdt,&
             nx,ny,nz,nxm,nym,mx,my,mz,m1,m2,ditime)
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
   USE aleatoire_m
   USE pwtcom_m
! 
   implicit none
!
   real(8),dimension(nx,ny,nz) :: ux,uy,gx,gy
   real(8),dimension(nx,ny,nz) :: uxn,uyn,ppn,u8x,u8y,g8x,g8y
   real(8),dimension(mx,my) :: ps2d1,ps2d2,ps2d3,us2d1,us2d2,us2d3
   real(8),dimension(ny,ditime) :: bx1,by1
   real(8),dimension(nx) :: ffx,fsx,fwx,ffxp,fsxp,fwxp
   real(8),dimension(ny) :: ffy,fsy,fwy,ffyp,fsyp,fwyp
   real(8),dimension(nx) :: sfx,ssx,swx,sfxp,ssxp,swxp
   real(8),dimension(ny) :: sfy,ssy,swy,sfyp,ssyp,swyp
   real(8),dimension(nxm) :: cfx6,csx6,cwx6,cifx6,cisx6,ciwx6 
   real(8),dimension(nxm) :: cifxp6,cisxp6,ciwxp6,cfxp6,csxp6,cwxp6
   real(8),dimension(nym) :: cfy6,csy6,cwy6,cify6,cisy6,ciwy6 
   real(8),dimension(nym) :: cifyp6,cisyp6,ciwyp6,cfyp6,csyp6,cwyp6
   real(8),dimension(nx) :: cfip6,csip6,cwip6
   real(8),dimension(nx) :: cifip6,cisip6,ciwip6
   real(8),dimension(ny) :: cfip6y,csip6y,cwip6y
   real(8),dimension(ny) :: cifip6y,cisip6y,ciwip6y
!  ********************* STRET *********************
   real(8),dimension(ny) :: pp4y
   real(8),dimension(nx) :: d5
   real(8),dimension(ny/2,m1+m2+1,nxm) :: a_tot,a2_tot
   real(8),dimension(ny/2,m1,nxm) :: al_tot,al1_tot
   real(8),dimension(ny/2,nxm) :: indx_tot,indx1_tot
!  *************************************************
   real(8) :: adt,bdt,gdt,t
   integer :: itime
   integer :: nx,ny,nz,nxm,nym,mx,my,mz,m1,m2,ditime
!
!
   print *,'itime = ',itime
!
!  Entrada: ux, uy, gx, gy. Salida: u8x, u8y, g8x, g8y, ps2d1
   call funcion1 (ux,uy,gx,gy,u8x,u8y,g8x,g8y,ps2d1,&
	  ffx,fsx,fwx,ffy,fsy,fwy,&
	  ffxp,fsxp,fwxp,&
	  ffyp,fsyp,fwyp,&
	  sfx,ssx,swx,sfy,ssy,swy,&
	  sfxp,ssxp,swxp,&
	  sfyp,ssyp,swyp,&
	  cfx6,csx6,cwx6,&
	  cfxp6,csxp6,cwxp6,&
	  cifx6,cisx6,ciwx6,cifxp6,cisxp6,ciwxp6,&
	  cfy6,csy6,cwy6,cfyp6,csyp6,cwyp6,&
	  cify6,cisy6,ciwy6,cifyp6,cisyp6,ciwyp6,&    
	  pp4y,bx1,by1,itime,ditime,nx,ny,nz,nxm,nym,mx,my,mz,adt,bdt,gdt)
! 
!  Entrada: ps2d1. Salida: ps2d2
   call funcion2(ps2d1,ps2d2,mx,my,nxm,nym)
! 
!  Entrada: ps2d2. Salida: ps2d3
   call funcion3(ps2d2,ps2d3,mx,my,nx,ny,nxm,nym)
! 
!  Entrada: ps2d3. Salida: us2d1
   call funcion4(ps2d3,us2d1,&
	  a_tot,al_tot,indx_tot,a2_tot,al1_tot,indx1_tot,d5,&
	  m1,m2,nx,nxm,ny,nym,mx,my)
! 
!  Entrada: us2d1. Salida: us2d2
   call funcion5(us2d1,us2d2,mx,my,nx,ny,nxm,nym)
! 
!  Entrada: us2d2. Salida: us2d3
   call funcion6(us2d2,us2d3,mx,my,nxm,nym)
!
!  Entrada: us2d3, u8x, u8y. Salida: uxn, uyn, ppn
   call funcion7_tray(us2d3,u8x,u8y,uxn,uyn,ppn,&
	  cfip6,csip6,cwip6,cfip6y,csip6y,cwip6y,&
	  cifip6,cisip6,ciwip6,cifip6y,cisip6y,ciwip6y,&
	  nx,nxm,mx,ny,nym,my,nz,gdt)
!
!  Entrada: uxn, uyn. Imprime DIV U
   call divU (uxn,uyn,&
	  cfx6,csx6,cwx6,&
	  cfxp6,csxp6,cwxp6,&
	  cifx6,cisx6,ciwx6,cifxp6,cisxp6,ciwxp6,&
	  cfy6,csy6,cwy6,cfyp6,csyp6,cwyp6,&
	  cify6,cisy6,ciwy6,cifyp6,cisyp6,ciwyp6,&
	  nx,ny,nz,nxm,nym)
!                
   return
end subroutine inc3dorigin_tray
!
!********************************************************************
!
subroutine initial (ux,uy,uz,sp,gx,gy,gz,sg,hx,hy,hz,sh,&
     vx,vy,vz,bx1,by1,bz1,bs1,bxo,byo,bzo,&
     bxr,byr,bzr,bsr,bxp,byp,bzp,bsp,&
     bxs,bys,bzs,nx,ny,nz,nyp,nzp,yp,iimin,i5,i55,jj,j1,j2,j3,j4,j5,angle3)
!
!********************************************************************
!
   USE paramx_m
   USE paramy_m
   USE paramz_m
   USE sauvegarde_m
   USE ecoulement_m
   USE delta_m
   USE param2_m
   USE aleatoire_m
   USE turbi_m    
   USE controle_m    
   USE barreau_m    
!
   implicit none
!
   real(8),dimension(nx,ny,nz) :: ux,uy,uz,gx,gy,gz,vx,vy,vz,hx,hy,hz
   real(8),dimension(nx,ny,nz) :: sp,sg,sh
   real(8),dimension(ny,nz) :: bx1,by1,bz1,bxr,byr,bzr
   real(8),dimension(ny,nz) :: bxp,byp,bzp,bxs,bys,bzs
   real(8),dimension(ny,nz) :: bs1,bsp,bsr
   real(8),dimension(nyp,nzp) :: bxo,byo,bzo   
   integer  :: nx,iimin,i5,i55 
   integer  :: ny,j1,j2,j3,j4,j5,jj
   integer  :: nz,nzdt,kp,n,iread
   integer,intent(in) :: nyp 
   integer,intent(in) :: nzp 
   real(8),dimension(ny+1) :: yp
!
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: k,j,i 
      real(8) :: seed,xinit,sigma,y,z,r,um,r1,r2,r3,t,h,angle3,pi
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      real(8),external :: alea,rand2
!-----------------------------------------------      
   pi=acos(-1.)
   sigma=0.2
! 
! ra c'est le rayon de ma plaque      
      iread=0
      print *,'iread entree de initial',iread
      if (ientree.eq.1) then
         if (nz.gt.1) then
            do k=1,nz
            do j=1,ny
               y=(j-1)*delty-y1
               z=(k-1)*deltz-z1
               r=sqrt(y*y+z*z)
               um=exp(-0.2*(r-rayon)*(r-rayon))
               do i=1,nx
                  r1=rand2(idum)-0.5
                  r2=rand2(idum)-0.5
                  r3=rand2(idum)-0.5
                  ux(i,j,k)=um*r1*bruit
                  uy(i,j,k)=um*r2*bruit
                  uz(i,j,k)=um*r3*bruit
               enddo
            enddo
            enddo
         else
!************************STRET STRET*******YP*********************
            do j=1,ny
               r=abs(yp(j)-y1)
               um=exp(-sigma*r*r)
               do i=1,nx
                  r1=rand2(idum)-0.5
                  r2=rand2(idum)-0.5
                  r3=rand2(idum)-0.5
                  ux(i,j,1)=um*r1*bruit
                  uy(i,j,1)=um*r2*bruit
                  uz(i,j,1)=um*r3*bruit
               enddo
            enddo
         endif
      else 
         if (nz.gt.1) then
            do i=1,icommence-1
               read(40) bxo,byo,bzo
            enddo
         else
            do i=1,icommence-1
               read(40) bxo,byo
            enddo
         endif
         do k=1,nz
         do j=1,ny
            bx1(j,k)=0.
            by1(j,k)=0.
            bz1(j,k)=0.
            bs1(j,k)=0.
         enddo
         enddo
         if (nz.gt.1) then
            iread=0
            do i=1,nx
               read(40) bxo,byo,bzo
               iread=iread+1
               do n=1,9
                  read(40) bxo,byo,bzo
                  iread=iread+1
               enddo
! Modif entree
!
               nzdt=(nz-1)/2
               do k=1,nzdt-(nzp+1)/2
               do j=1,ny
                  bx1(j,k)=bxo(j,1)
                  by1(j,k)=byo(j,1)
                  bz1(j,k)=bzo(j,1)
               enddo
               enddo
               do k=nzdt+(nzp+1)/2,nz
               do j=1,ny
                  bx1(j,k)=bxo(j,nzp)
                  by1(j,k)=byo(j,nzp)
                  bz1(j,k)=bzo(j,nzp)
               enddo
               enddo
               kp=1
               do k=nzdt-(nzp-1)/2,nzdt+(nzp-1)/2
                  do j=1,ny
                     bx1(j,k)=bxo(j,kp)
                     by1(j,k)=byo(j,kp)
                     bz1(j,k)=bzo(j,kp)
                  enddo
                  kp=kp+1
               enddo 
! Fin Modif entree
!
               do k=1,nz
               do j=1,ny
                  y=yp(j)-y1
                  z=(k-1)*deltz-z1
                  um=exp(-0.2*z*z)
                  ux(nx-i+1,j,k)=bruit*um*bx1(j,k)
                  uy(nx-i+1,j,k)=bruit*um*by1(j,k)
                  uz(nx-i+1,j,k)=bruit*um*bz1(j,k)
               enddo
               enddo
            enddo
         else
            do i=1,nx
               read(40) bxo,byo
               do j=1,ny/2-nyp/2
                  bx1(j,1)=bxo(1,1)
                  by1(j,1)=byo(1,1)
                  bx1(ny-j+1,1)=bx1(j,1)
                  by1(ny-j+1,1)=by1(j,1)
               enddo
               do j=1,nyp
                  bx1(j+ny/2-nyp/2,1)=bxo(j,1)
                  by1(j+ny/2-nyp/2,1)=byo(j,1)
               enddo 
               do j=1,ny
                  y=yp(j)-y1
                  r=sqrt(y*y)
                  um=exp(-0.2*(r-rayon)*(r-rayon))
                  ux(nx-i+1,j,1)=bruit*um*bx1(j,1)
                  uy(nx-i+1,j,1)=bruit*um*by1(j,1)
               enddo
            enddo
         endif
      endif
!
      print *,'iread',iread
!
      t=0.
      call ecoule (bx1,by1,bz1,bs1,nx,ny,nz,t,yp,jj,j1,j2,j3,j4,j5,iimin,i5,i55,h,angle3)
!
      if (iscalaire.eq.0) then
         if (nschema.eq.1) then
            do k=1,nz
             do j=1,ny
             do i=1,nx
                ux(i,j,k)=ux(i,j,k)+bx1(j,k)
                uy(i,j,k)=uy(i,j,k)+by1(j,k)
                uz(i,j,k)=uz(i,j,k)+bz1(j,k)
                gx(i,j,k)=ux(i,j,k)
                gy(i,j,k)=uy(i,j,k)
                gz(i,j,k)=uz(i,j,k)
                hx(i,j,k)=ux(i,j,k)
                hy(i,j,k)=uy(i,j,k)
                hz(i,j,k)=uz(i,j,k)
             enddo
             enddo
             if (ivirtuel.eq.1) then 
                do j=j5,j4
                do i=2,nx
                   ux(i,j,k)=0.
                   uy(i,j,k)=0.
                   uz(i,j,k)=uz(i,j,k)+bz1(j,k)
                   gx(i,j,k)=ux(i,j,k)
                   gy(i,j,k)=uy(i,j,k)
                   gz(i,j,k)=uz(i,j,k)
                   hx(i,j,k)=ux(i,j,k)
                   hy(i,j,k)=uy(i,j,k)
                   hz(i,j,k)=uz(i,j,k)
                enddo
                enddo
             endif
            enddo
         else
            do k=1,nz
            do j=1,ny
            do i=1,nx
               ux(i,j,k)=ux(i,j,k)+bx1(j,k)
               uy(i,j,k)=uy(i,j,k)+by1(j,k)
               uz(i,j,k)=uz(i,j,k)+bz1(j,k)
            enddo
            if (ientree.eq.1) then
               bxp(j,k)=0.
               byp(j,k)=0.
               bzp(j,k)=0.
               bxs(j,k)=0.
               bys(j,k)=0.
               bzs(j,k)=0.
            else
               bxp(j,k)=bx1(j,k)
               byp(j,k)=by1(j,k)
               bzp(j,k)=bz1(j,k)
               bxr(j,k)=bx1(j,k)
               byr(j,k)=by1(j,k)
               bzr(j,k)=bz1(j,k)
            endif
            enddo
            enddo
         endif
      else
         if (nschema.eq.1) then
            do k=1,nz
            do j=1,ny
            do i=1,nx
               ux(i,j,k)=ux(i,j,k)+bx1(j,k)
               uy(i,j,k)=uy(i,j,k)+by1(j,k)
               uz(i,j,k)=uz(i,j,k)+bz1(j,k)
               sp(i,j,k)=bs1(j,k)
               gx(i,j,k)=ux(i,j,k)
               gy(i,j,k)=uy(i,j,k)
               gz(i,j,k)=uz(i,j,k)
               sg(i,j,k)=sp(i,j,k)
               hx(i,j,k)=ux(i,j,k)
               hy(i,j,k)=uy(i,j,k)
               hz(i,j,k)=uz(i,j,k)
               sh(i,j,k)=sp(i,j,k)
            enddo
            enddo
            enddo
         else
            do k=1,nz
            do j=1,ny
            do i=1,nx
               ux(i,j,k)=ux(i,j,k)+bx1(j,k)
               uy(i,j,k)=uy(i,j,k)+by1(j,k)
               uz(i,j,k)=uz(i,j,k)+bz1(j,k)
               sp(i,j,k)=bs1(j,k)
               gx(i,j,k)=ux(i,j,k)
               gy(i,j,k)=uy(i,j,k)
               gz(i,j,k)=uz(i,j,k)
               sg(i,j,k)=sp(i,j,k)
            enddo
            if (ientree.eq.1) then
               bxp(j,k)=0.
               byp(j,k)=0.
               bzp(j,k)=0.
               bsp(j,k)=0.
               bxs(j,k)=0.
               bys(j,k)=0.
               bzs(j,k)=0.
               bsr(j,k)=0.
            else
               bxp(j,k)=bx1(j,k)
               byp(j,k)=by1(j,k)
               bzp(j,k)=bz1(j,k)
               bsp(j,k)=bs1(j,k)
               bxr(j,k)=bx1(j,k)
               byr(j,k)=by1(j,k)
               bzr(j,k)=bz1(j,k)
               bsr(j,k)=bs1(j,k)
            endif
            enddo
            enddo
         endif
      endif
!
   return
end subroutine initial
!
!********************************************************************
!
subroutine visqueux (vx,vy,vz,ux,uy,uz,xnut,sp,sv,&
     tx,ty,tz,tr,di,sx,sy,sz,&
     sfx,ssx,swx,sfy,ssy,swy,sfz,ssz,swz,&
     sfxp,ssxp,swxp,&
     sfyp,ssyp,swyp,&
     sfzp,sszp,swzp,&
     ffy,fsy,fwy,ffyp,fsyp,fwyp,pp4y,ty1,nx,ny,nz,yp)
!********************************************************************
!
   USE param1_m
   USE ecoulement_m
   USE paramy_m
   USE paramx_m
!
   implicit none
!
   real(8),dimension(nx,ny,nz) :: ux,uy,uz,vx,vy,vz,tx,ty,tz,tr,sp,sv
   real(8),dimension(nx,ny,nz) :: di,xnut,ty1
   real(8),dimension(ny,nz) :: sx
   real(8),dimension(nx,nz) :: sy
   real(8),dimension(nx,ny) :: sz
   real(8),dimension(nx) :: sfx,ssx,swx,sfxp,ssxp,swxp
   real(8),dimension(ny) :: sfy,ssy,swy,sfyp,ssyp,swyp
   real(8),dimension(nz) :: sfz,ssz,swz,sfzp,sszp,swzp 
   real(8),dimension(nz) :: ffy,fsy,fwy,ffyp,fsyp,fwyp
!
!************************STRET********************************
   real(8),dimension(ny) :: pp4y
   real(8),dimension(ny+1) :: yp
   real(8) :: yl2y,pi,x
!*************************************************************
   integer  :: nx 
   integer  :: ny 
   integer  :: nz
!
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: nxyz,ijk,j,i,k 
!----------------------------------------------- 
   nxyz=nx*ny*nz
!
   yl2y=yly*yly
!
   call derxx (tx,ux,tr,sx,sfx ,ssx ,swx ,nx,ny,nz,0)
   call deryy (ty,ux,tr,di,sy,sfyp,ssyp,swyp,nx,ny,nz,1)
!  **************** correccion por stretching ****************
   call dery (ty1,ux,tr,di,sy,ffyp,fsyp,fwyp,nx,ny,nz,1)
   do k=1,nz
   do j=1,ny
   do i=1,nx
      ty(i,j,k)=ty(i,j,k)*yl2y-pp4y(j)*ty1(i,j,k)/yl2y  ! ec. (3.58) tesis sylvain
   enddo
   enddo
   enddo
!  ***********************************************************
   do j=1,ny
   do i=1,nx
      vx(i,j,1)=xnu*(tx(i,j,1)+ty(i,j,1))-vx(i,j,1)
   enddo
   enddo
!
!
   call derxx (tx,uy,tr,sx,sfxp,ssxp,swxp,nx,ny,nz,1)
   call deryy (ty,uy,tr,di,sy,sfy,ssy,swy,nx,ny,nz,0)
!  **************** correccion por stretching ****************
   call dery (ty1,uy,tr,di,sy,ffy,fsy,fwy,nx,ny,nz,0)
   do k=1,nz
   do j=1,ny
   do i=1,nx
      ty(i,j,k)=ty(i,j,k)*yl2y-pp4y(j)*ty1(i,j,k)/yl2y
   enddo
   enddo
   enddo
!  ***********************************************************
   do j=1,ny
   do i=1,nx
      vy(i,j,1)=xnu*(tx(i,j,1)+ty(i,j,1))-vy(i,j,1)
   enddo
   enddo
!
   return
end subroutine visqueux
!
!********************************************************************
!
!  Entrada: ux, uy, gx, gy. Salida: u8x, u8y, g8x, g8y, ps2d1
subroutine funcion1 (ux,uy,gx,gy,u8x,u8y,g8x,g8y,ps2d1,&
             ffx,fsx,fwx,ffy,fsy,fwy,&
             ffxp,fsxp,fwxp,&
             ffyp,fsyp,fwyp,&
             sfx,ssx,swx,sfy,ssy,swy,&
             sfxp,ssxp,swxp,&
             sfyp,ssyp,swyp,&
             cfx6,csx6,cwx6,&
             cfxp6,csxp6,cwxp6,&
             cifx6,cisx6,ciwx6,cifxp6,cisxp6,ciwxp6,&
             cfy6,csy6,cwy6,cfyp6,csyp6,cwyp6,&
             cify6,cisy6,ciwy6,cifyp6,cisyp6,ciwyp6,&    
             pp4y,bx1,by1,itime,ditime,nx,ny,nz,nxm,nym,mx,my,mz,adt,bdt,gdt)
!********************************************************************
!
   USE paramx_m
   USE paramy_m
   USE paramz_m
   USE param1_m
   USE ecoulement_m
   USE delta_m
   USE barreau_m
   USE dericex6_m
   USE interpol6_m
   USE dericey6_m
   USE interpoly6_m
!
   implicit none
!
   real(8),dimension(nx,ny,nz) :: ux,uy,vx,vy,rx,ry
   real(8),dimension(nx,ny,nz) :: u8x,u8y,gx,gy,g8x,g8y
   real(8),dimension(nx,ny,nz) :: cx,cy,cz,tx,ty,tz
   real(8),dimension(nx,ny,nz) :: tz1,ty1,ty2,tpx,tpy,tpy1,tpy2
   real(8),dimension(nxm,nym,nz) :: wtx,wty
   real(8),dimension(nx,nym,nz) :: wty1
   real(8),dimension(nxm,ny,nz) :: wtx1
   real(8),dimension(nxm,nym,nz) :: ppm
   real(8),dimension(mx,my,mz) :: ps,ps1
   real(8),dimension(mx,my) :: ps2d1
   real(8),dimension(ny,ditime) :: bx1,by1
   real(8),dimension(ny,nz) :: bxn,byn
   real(8),dimension(nx) :: ffx,fsx,fwx,ffxp,fsxp,fwxp
   real(8),dimension(ny) :: ffy,fsy,fwy,ffyp,fsyp,fwyp
   real(8),dimension(nx) :: sfx,ssx,swx,sfxp,ssxp,swxp
   real(8),dimension(ny) :: sfy,ssy,swy,sfyp,ssyp,swyp
   real(8),dimension(nxm) :: cfx6,csx6,cwx6,cifx6,cisx6,ciwx6 
   real(8),dimension(nxm) :: cifxp6,cisxp6,ciwxp6,cfxp6,csxp6,cwxp6
   real(8),dimension(nym) :: cfy6,csy6,cwy6,cify6,cisy6,ciwy6 
   real(8),dimension(nym) :: cifyp6,cisyp6,ciwyp6,cfyp6,csyp6,cwyp6
!*************** STRET *******************
   real(8),dimension(ny) :: pp4y
   real(8) :: yl2y
!*****************************************
   integer,intent(in) :: itime,ditime,nx,ny,nz,nxm,nym,mx,my,mz 
   real(8),intent(in) :: adt,bdt,gdt 
!
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
   integer :: i,j,k 
   real(8) :: udxv,vphase,cxv,vphase_sum
!----------------------------------------------- 
!
!  ******************** urotu *********************
   call derx (tz,uy,ffxp,fsxp,fwxp,nx,ny,nz,1)
   call dery (cz,ux,ffyp,fsyp,fwyp,nx,ny,nz,1)
!
!  ******************* visqueux *******************
   call derxx (tx,ux,sfx ,ssx ,swx ,nx,ny,nz,0)
   call deryy (ty,ux,sfyp,ssyp,swyp,nx,ny,nz,1)
   call dery (ty1,ux,ffyp,fsyp,fwyp,nx,ny,nz,1)
   call derxx (tpx,uy,sfxp,ssxp,swxp,nx,ny,nz,1)
   call deryy (tpy,uy,sfy,ssy,swy,nx,ny,nz,0)
   call dery (tpy1,uy,ffy,fsy,fwy,nx,ny,nz,0)
!
   yl2y=yly*yly
!
   do j=1,ny
   do i=1,nx
!     ******************** urotu *********************
      tz1(i,j,1)= tz(i,j,1)-cz(i,j,1)
      cx(i,j,1)=-tz1(i,j,1)*uy(i,j,1)
      cy(i,j,1)= tz1(i,j,1)*ux(i,j,1)
!
!     ******************* visqueux *******************
      ty2(i,j,1)=ty(i,j,1)*yl2y-pp4y(j)*ty1(i,j,1)/yl2y  
      tpy2(i,j,1)=tpy(i,j,1)*yl2y-pp4y(j)*tpy1(i,j,1)/yl2y
      vx(i,j,1)=xnu*(tx(i,j,1)+ty2(i,j,1))
      vy(i,j,1)=xnu*(tpx(i,j,1)+tpy2(i,j,1))
!
!     ******************** F(u^n) ********************
      rx(i,j,1)=vx(i,j,1)-cx(i,j,1)
      ry(i,j,1)=vy(i,j,1)-cy(i,j,1)
!
!     ******************** adams ********************
      u8x(i,j,1)=adt*rx(i,j,1)+bdt*gx(i,j,1)+ux(i,j,1) 
      u8y(i,j,1)=adt*ry(i,j,1)+bdt*gy(i,j,1)+uy(i,j,1) 
      g8x(i,j,1)=rx(i,j,1)
      g8y(i,j,1)=ry(i,j,1)
   enddo
   enddo
!
   udxv=1./dx
!
   vphase_sum=0.
   do j=1,ny
      vphase_sum=vphase_sum+ux(nx,j,1)
   enddo
   vphase=vphase_sum/ny
!
   cxv=vphase*gdt*udxv
!
   do j=1,ny
!     ******************** sortie ********************
      bxn(j,1)=ux(nx,j,1)-cxv*(ux(nx,j,1)-ux(nx-1,j,1))
      byn(j,1)=uy(nx,j,1)-cxv*(uy(nx,j,1)-uy(nx-1,j,1))
!     **************** cond. de borde ****************
      u8x(nx,j,1)=bxn(j,1)
      u8y(nx,j,1)=byn(j,1)
      u8x(1,j,1)=bx1(j,itime)
      u8y(1,j,1)=by1(j,itime)
   enddo
!
!  ******************* divbric16bis ********************
   call intery6(wty1,u8x,cifyp6,cisyp6,ciwyp6,nx,ny,nym,nz,1)
   call decx6(wtx,wty1,cfx6,csx6,cwx6,nx,nxm,ny,nym,nz,0)  
   call inter6(wtx1,u8y,cifxp6,cisxp6,ciwxp6,nx,nxm,ny,nz,1)
   call decy6(wty,wtx1,cfy6,csy6,cwy6,nx,nxm,ny,nym,nz,0)
!     
   do j=1,nym
   do i=1,nxm
      ppm(i,j,1)=wtx(i,j,1)+wty(i,j,1)
   enddo
   enddo
!
!  ******************* temporellebis16STR ********************
!  ************** D(nxm,nym,nz) ---> D(mx,my) ****************
   do k=1,mz
   do j=1,my
   do i=1,mx
      ps(i,j,k)=0.
      ps1(i,j,k)=0.
   enddo
   enddo
   enddo
!
   do j=1,my
   do i=1,mx
      ps2d1(i,j)=0.
   enddo
   enddo
!
   do k=1,nz
   do j=1,nym
   do i=1,nxm/2
      ps(i,j,k)=ppm(2*(i-1)+1,j,k)
   enddo
   enddo
   enddo
!
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
!
   do k=1,nz
   do i=1,nxm
   do j=nym/2+1,nym
      ps1(i,j,k)=ps(i,2*nym-2*j+2,k)
   enddo
   enddo
   enddo
!
   do i=1,nxm
   do j=1,nym
      ps2d1(i,j)=ps1(i,j,1)
   enddo
   enddo     
!
   return
end subroutine funcion1
!
!********************************************************************
!
!  Entrada: ps2din. Salida: ps2d
subroutine funcion2(ps2din,ps2d,mx,my,nxm,nym)
!********************************************************************
!
   implicit none
!
   real(8),dimension(mx,my) :: ps2d,ps2din
   integer :: i,j,mx,my,nxm,nym
!***************** FFT *******************
   real(8) :: scale
   real(8),dimension(100+2*(nxm+nym)) :: table
   real(8),dimension(512*max(nxm,nym)) :: work
!*****************************************
!
!  ******************* temporellebis16STR ********************
!  **************** D(mx,my) ---> D^(mx,my) ******************
!
!  Guardo la entrada (GF)
   do i=1,mx
   do j=1,my
      ps2d(i,j)=ps2din(i,j)
   enddo
   enddo     
!
   scale=1./sqrt(real(nxm*nym))
!
   call scfft2d(0,nxm,nym,scale,ps2d,mx,ps2d,mx/2,table,work,0)
   call scfft2d(-1,nxm,nym,scale,ps2d,mx,ps2d,mx/2,table,work,0)
!
!  Construyo matriz que conserve la norma (CN) 
   do i=3,nxm
   do j=1,my
      ps2d(i,j) = sqrt(2.)*ps2d(i,j);
   enddo
   enddo
!
end subroutine funcion2
!
!********************************************************************
!
!  Entrada: ps2d. Salida: ps2din
subroutine funcion2inv(ps2d,ps2din,mx,my,nxm,nym)
!********************************************************************
!
   implicit none
!
   real(8),dimension(mx,my) :: ps2d,ps2din
   integer :: i,j,mx,my,nxm,nym
!***************** FFT *******************
   real(8) :: scale
   real(8),dimension(100+2*(nxm+nym)) :: table
   real(8),dimension(512*max(nxm,nym)) :: work
!*****************************************
!
!  (CN^-1)
   do i=3,nxm
   do j=1,my
      ps2d(i,j) = (1./sqrt(2.))*ps2d(i,j);
   enddo
   enddo
!
   scale=1./sqrt(real(nxm*nym))
!
   call scfft2d(0,nxm,nym,scale,ps2d,mx,ps2d,mx/2,table,work,0)
   call csfft2d(1,nxm,nym,scale,ps2d,mx/2,ps2d,mx,table,work,0)
!
!  (GF^-1)
   ps2din(:,:)=0.
   do i=1,nxm
   do j=1,nym
      ps2din(i,j)=ps2d(i,j)
   enddo
   enddo   
!
end subroutine funcion2inv
!
!********************************************************************
!
!  Entrada: ps2d2. Salida: ps2d3
subroutine funcion3(ps2d2,ps2d3,mx,my,nx,ny,nxm,nym)
!********************************************************************
!
   implicit none
!
   real(8),dimension(mx,my) :: ps2d2,ps2d3
   real(8),dimension(mx,my) :: ur,ue,ur1,ue1
   real(8) :: xx1,xx2
   real(8) :: pi
   integer :: i,j,mx,my,nx,ny,nxm,nym
!
!  ******************* temporellebis16STR ********************
!  *************** D^(mx,my) ---> D^d(mx,my) *****************
!
   pi=acos(-1.)
!
!  (CN^-1). Asi recupero la salida de scfft2d corresp. a funcion2
   do i=3,nxm
   do j=1,my
      ps2d2(i,j) = (1./sqrt(2.))*ps2d2(i,j);
   enddo
   enddo
!
   do j=1,my
   do i=1,mx
      ps2d3(i,j)=ps2d2(i,j)
   enddo
   enddo
!
   do j=1,my
   do i=1,mx
      ur(i,j)=0.
      ue(i,j)=0.
      ur1(i,j)=0.
      ue1(i,j)=0.
   enddo
   enddo
!
   do j=1,nym
   do i=1,nxm/2
      ur(i,j)=ps2d2(2*i-1,j)
   enddo
   enddo
!      
   do j=1,nym
      ur(nxm/2+1,j)=ps2d2(nxm+1,j)
   enddo
!
   do j=1,nym
   do i=2,nxm/2
      ur(i+nxm/2,j)=ur(nxm/2-i+2,j)
   enddo
   enddo
!
   do j=1,nym
   do i=1,nxm/2
      ue(i,j)=ps2d2(2*i,j)
   enddo
   enddo
!
   do j=1,nym
      ue(nxm/2+1,j)=ps2d2(nx+1,j)
   enddo
!
   do j=1,nym
   do i=2,nxm/2
      ue(i+nxm/2,j)=-ue(nxm/2-i+2,j)
   enddo
   enddo
!
   do j=1,ny
   do i=1,nxm/2+1
      xx1=ur(i,j)*cos((j-1)/2.*pi/(nym))/2.+ue(i,j)*sin((j-1)/2.*pi/(nym))/2.
      xx2=ur(i,nym-j+2)*cos((j-1)/2.*pi/(nym))/2.-ue(i,nym-j+2)*sin((j-1)/2.*pi/(nym))/2.           
      ur1(i,j)=xx1+xx2      
   enddo
   enddo
!
   do j=1,ny
   do i=2,nxm/2
      ur1(i+nxm/2,j)=ur1(nxm/2-i+2,j)
   enddo
   enddo
!
   do j=1,ny
   do i=1,nxm/2+1
      xx1=-ur(i,j)*sin((j-1)/2.*pi/(nym))/2.+ue(i,j)*cos((j-1)/2.*pi/(nym))/2.
      xx2=ur(i,nym-j+2)*sin((j-1)/2.*pi/(nym))/2.+ue(i,nym-j+2)*cos((j-1)/2.*pi/(nym))/2.           
      ue1(i,j)=xx1+xx2 
   enddo
   enddo
!
   do j=1,ny
   do i=2,nxm/2
      ue1(i+nxm/2,j)=-ue1(nxm/2-i+2,j)
   enddo
   enddo
!
   do i=1,nxm
      ur1(i,1)=ur1(i,1)*2.
      ue1(i,1)=ue1(i,1)*2.
   enddo
!  
   do j=1,ny
   do i=1,nxm
      xx1=(ur1(i,j)+ur1(nxm-i+2,j))*cos((i-1)/2.*pi/(nxm))/2.
      xx2=(ue1(i,j)-ue1(nxm-i+2,j))*sin((i-1)/2.*pi/(nxm))/2.
      ps2d3(i,j)=xx1+xx2
   enddo
   enddo
!
   do j=1,ny
      ps2d3(1,j)=2.*ps2d3(1,j)
   enddo
!
end subroutine funcion3
!
!********************************************************************
!
!  Entrada: ps2d. Salida: us2d
subroutine funcion4(ps2d,us2d,&
             a_tot,al_tot,indx_tot,a2_tot,al1_tot,indx1_tot,d5,&
             m1,m2,nx,nxm,ny,nym,mx,my)
!********************************************************************
!
   implicit none
!
   real(8),dimension(mx,my) :: ps2d,us2d
   real(8),dimension(ny/2,m1+m2+1) :: a,a1,a2
   real(8),dimension(ny/2,m1) :: al,al1
   real(8),dimension(ny/2) :: indx,indx1,e,c,e1,c1
   real(8),dimension(ny/2,m1+m2+1,nxm) :: a_tot,a2_tot
   real(8),dimension(ny/2,m1,nxm) :: al_tot,al1_tot
   real(8),dimension(ny/2,nxm) :: indx_tot,indx1_tot
   real(8),dimension(nx) :: d5
   integer :: i,j,k,nx,mx,nxm,ny,my,nym,m1,m2
!
!  ******************* temporellebis16STR ********************
!  *************** D^d(mx,my) ---> p^d(mx,my) ****************
!
   do j=1,my
   do i=1,mx
      us2d(i,j)=0.
   enddo
   enddo
!
   do i=1,nxm ! Resuelve el sistema B p^d = D^d
!
!     Recupera resultados de la descomposicion LU de la matriz de coef. B=A^2
      a = a_tot(:,:,i)
      al = al_tot(:,:,i)
      indx = indx_tot(:,i)
      a2 = a2_tot(:,:,i)
      al1 = al1_tot(:,:,i)
      indx1 = indx1_tot(:,i)
!
      do j=1,ny/2
         e(j)=0.
         c(j)=0.
         e1(j)=0.
         c1(j)=0.
      enddo
!
!     Separamos los a+ib provenientes de D^
      do j=1,ny/2
         e(j)=ps2d(i,2*j-1)
         c(j)=ps2d(i,2*j)
      enddo
!     
!     Solves the band diagonal linear equations A·x = e.
!     Given the arrays a, al, and indx as returned from bandec, and given a right-hand-side vector e.
!     The solution vector x se escribe en e1.
      call banbks1(a,ny,my,m1,m2,al,indx,e,e1)
! 
      call banbks1(a2,ny,my,m1,m2,al1,indx1,c,c1)
!
!     On reconstruit les a+ib provenientes de p^
      do j=1,ny-1,2
         us2d(i,j)=e1((j+1)/2)
      enddo
      do j=2,ny,2
         us2d(i,j)=c1(j/2)
      enddo 
!
      if (d5(i).eq.0) then
         us2d(i,1)=0.
         us2d(i,ny)=0.
      endif
!
   enddo ! Resuelve el sistema B p^d = D^d
!
end subroutine funcion4
!
!********************************************************************
!
!  Entrada: us2d1. Salida: us2d2
subroutine funcion5(us2d1,us2d2,mx,my,nx,ny,nxm,nym)
!********************************************************************
!
   implicit none
!
   real(8),dimension(mx,my) :: us2d1,us2d2
   real(8),dimension(mx,my) :: ur1,ue1,w1,w11
   real(8) :: pi,xx1,xx2,xx3,xx4
   integer :: i,j,mx,my,nx,ny,nxm,nym
!
!  ******************* temporellebis16STR ********************
!  *************** p^d(mx,my) ---> p^(mx,my) *****************
!
   pi=acos(-1.)
!
   do j=1,my
   do i=1,mx
      ur1(i,j)=0.
      ue1(i,j)=0.
   enddo
   enddo
!
   do j=1,ny
   do i=1,nxm/2+1
      w1(i,j)=us2d1(i,j)
   enddo
   enddo
!
   do j=1,ny
   do i=1,nxm/2+1
      w11(i,j)=us2d1(nxm-i+2,j)
   enddo
   enddo
!
   do j=1,ny
   do i=1,nxm/2+1
      ur1(i,j)=w1(i,j)*cos((i-1)/2.*pi/(nxm))+w11(i,j)*cos((nxm-i+1)/2.*pi/(nxm))
      ue1(i,j)=w1(i,j)*sin((i-1)/2.*pi/(nxm))-w11(i,j)*sin((nxm-i+1)/2.*pi/(nxm))
   enddo
   enddo
!
   do j=1,ny
   do i=2,nxm/2
      ur1(i+nxm/2,j)=ur1(nxm/2-i+2,j)
   enddo
   enddo
!
   do j=1,ny
   do i=2,nxm/2
      ue1(i+nxm/2,j)=-ue1(nxm/2-i+2,j)
   enddo
   enddo
!
   do j=1,ny
   do i=1,nxm/2+1
      xx1=ur1(i,j)*cos((j-1)/2.*pi/(nym))-ue1(i,j)*cos((nym-j+1)/2.*pi/(nym))
      xx2=ur1(i,nym-j+2)*sin((j-1)/2.*pi/(nym))+ue1(i,nym-j+2)*cos((j-1)/2.*pi/(nym)) 
      xx3=-ur1(i,nym-j+2)*sin((nym-j+1)/2.*pi/(nym))+ue1(i,nym-j+2)*sin((j-1)/2.*pi/(nym))
      xx4=ur1(i,j)*sin((j-1)/2.*pi/(nym))+ue1(i,j)*cos((j-1)/2.*pi/(nym))
      us2d2(2*i-1,j)=xx1+xx2
      us2d2(2*i,j)=xx3+xx4
   enddo
   enddo
!
   do i=1,mx
     us2d2(i,ny)=0.
   enddo
!
   do i=1,mx
      us2d2(i,ny)=us2d2(i,(ny+1)/2)
   enddo
!
!  (CN). Asi recupero la entrada de csfft2d corresp. a funcion6
   do i=3,nxm
   do j=1,my
      us2d2(i,j) = sqrt(2.)*us2d2(i,j);
   enddo
   enddo
!
end subroutine funcion5
!
!********************************************************************
!
!  Entrada: ps2d. Salida: ps2din
subroutine funcion6(ps2d,ps2din,mx,my,nxm,nym)
!********************************************************************
!
   implicit none
!
   real(8),dimension(mx,my) :: ps2d,ps2din
   integer :: i,j,mx,my,nxm,nym
!***************** FFT *******************
   real(8) :: scale
   real(8),dimension(100+2*(nxm+nym)) :: table
   real(8),dimension(512*max(nxm,nym)) :: work
!*****************************************
!
!  ****************** temporellebis16STR *******************
!  *************** p^(mx,my) ---> p(mx,my) *****************
!
!  (CN^-1)
   do i=3,nxm
   do j=1,my
      ps2d(i,j) = (1./sqrt(2.))*ps2d(i,j);
   enddo
   enddo
!
   scale=1./sqrt(real(nxm*nym))
!
   call scfft2d(0,nxm,nym,scale,ps2d,mx,ps2d,mx/2,table,work,0)
   call csfft2d(1,nxm,nym,scale,ps2d,mx/2,ps2d,mx,table,work,0)
!
!  (GF^-1)
   ps2din(:,:)=0.
   do i=1,nxm
   do j=1,nym
      ps2din(i,j)=ps2d(i,j)
   enddo
   enddo   
!
end subroutine funcion6
!
!********************************************************************
!
!  Entrada: ps2din. Salida: ps2d
subroutine funcion6inv(ps2din,ps2d,mx,my,nxm,nym)
!********************************************************************
!
   implicit none
!
   real(8),dimension(mx,my) :: ps2d,ps2din
   integer :: i,j,mx,my,nxm,nym
!***************** FFT *******************
   real(8) :: scale
   real(8),dimension(100+2*(nxm+nym)) :: table
   real(8),dimension(512*max(nxm,nym)) :: work
!*****************************************
!
!  Guardo la entrada (GF)
   do i=1,mx
   do j=1,my
      ps2d(i,j)=ps2din(i,j)
   enddo
   enddo     
!
   scale=1./sqrt(real(nxm*nym))
!
   call scfft2d(0,nxm,nym,scale,ps2d,mx,ps2d,mx/2,table,work,0)
   call scfft2d(-1,nxm,nym,scale,ps2d,mx,ps2d,mx/2,table,work,0)
!
!  Construyo matriz que conserve la norma (CN) 
   do i=3,nxm
   do j=1,my
      ps2d(i,j) = sqrt(2.)*ps2d(i,j);
   enddo
   enddo
!
end subroutine funcion6inv
!
!********************************************************************
!
!  Entrada: us2d3, u8x, u8y. Salida: uxn, uyn
subroutine funcion7(us2d3,u8x,u8y,uxn,uyn,&
             cfip6,csip6,cwip6,cfip6y,csip6y,cwip6y,&
             cifip6,cisip6,ciwip6,cifip6y,cisip6y,ciwip6y,&
             nx,nxm,mx,ny,nym,my,nz)
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
   real(8),dimension(mx,my) :: us2d3,ux2d,pp2d
   real(8),dimension(nx) :: cfip6,csip6,cwip6
   real(8),dimension(nx) :: cifip6,cisip6,ciwip6
   real(8),dimension(ny) :: cfip6y,csip6y,cwip6y
   real(8),dimension(ny) :: cifip6y,cisip6y,ciwip6y
   integer :: i,j,k,nx,mx,nxm,ny,my,nym,nz
!
!  ******************* temporellebis16STR ********************
!  ************** p(mx,my) ---> p(nxm,nym,nz) ****************
   do j=1,ny
   do i=1,nxm/2
      pp2d(2*i-1,j)=us2d3(i,j)
   enddo
   enddo
!
   do j=1,ny
   do i=1,nxm/2
      pp2d(2*i,j)=us2d3(nxm-i+1,j) 
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
end subroutine funcion7
!
!********************************************************************
!
!  Entrada: us2d3, u8x, u8y. Salida: uxn, uyn, ppn
subroutine funcion7_tray(us2d3,u8x,u8y,uxn,uyn,ppn,&
             cfip6,csip6,cwip6,cfip6y,csip6y,cwip6y,&
             cifip6,cisip6,ciwip6,cifip6y,cisip6y,ciwip6y,&
             nx,nxm,mx,ny,nym,my,nz,gdt)
!********************************************************************
!
   USE paramx_m
   USE paramy_m
!
   implicit none
!
   real(8),dimension(nx,ny,nz) :: uxn,uyn,u8x,u8y,px,py,pp,ppn
   real(8),dimension(nxm,nym,nz) :: ppm
   real(8),dimension(nxm,ny,nz) :: wk1
   real(8),dimension(nx,nym,nz) :: wk2
   real(8),dimension(mx,my) :: us2d3,ux2d,pp2d
   real(8),dimension(nx) :: cfip6,csip6,cwip6
   real(8),dimension(nx) :: cifip6,cisip6,ciwip6
   real(8),dimension(ny) :: cfip6y,csip6y,cwip6y
   real(8),dimension(ny) :: cifip6y,cisip6y,ciwip6y
   real(8) :: gdt,udt
   integer :: i,j,k,nx,mx,nxm,ny,my,nym,nz
!
!  ******************* temporellebis16STR ********************
!  ************** p(mx,my) ---> p(nxm,nym,nz) ****************
   do j=1,ny
   do i=1,nxm/2
      pp2d(2*i-1,j)=us2d3(i,j)
   enddo
   enddo
!
   do j=1,ny
   do i=1,nxm/2
      pp2d(2*i,j)=us2d3(nxm-i+1,j) 
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
!  ***********************************************************
!  *********** Total pressure --> Static pressure ************
!  ***********************************************************
!
!  ************* p(nxm,nym,nz) ---> p(nx,ny,nz) **************
!  interpolation en y --> (p)^I_y
   call interiy6(wk1,ppm,cifip6y,cisip6y,ciwip6y,nxm,nym,ny,nz,1)
!
!  interpolation en x --> (p)^I_xy
   call interi6(pp,wk1,cifip6,cisip6,ciwip6,nxm,nx,ny,nz,1)
!
   udt=1./gdt
!
   do j=1,ny
   do i=1,nx
!       ppn(i,j,1)=udt*pp(i,j,1)-0.5*(uxn(i,j,1)*uxn(i,j,1)+uyn(i,j,1)*uyn(i,j,1))
      ppn(i,j,1)=udt*pp(i,j,1)
   enddo
   enddo
!
   return
end subroutine funcion7_tray
!
!********************************************************************
!
subroutine gradp(pp,px,py,py1,pz,tr,di,sx,sy,sz,&
     ffx,fsx,fwx,ffy,fsy,fwy,ffz,fsz,fwz,&
     ffxp,fsxp,fwxp,&
     ffyp,fsyp,fwyp,&
     ffzp,fszp,fwzp,&
     nx,ny,nz,gdt)
!
!********************************************************************
!
   USE paramx_m
!
   implicit none
!
   real(8),dimension(nx,ny,nz) :: px,py,pz,pp,tr,di,py1
   real(8),dimension(ny,nz) :: sx
   real(8),dimension(nx,nz) :: sy
   real(8),dimension(nx,ny) :: sz
   real(8),dimension(nx) :: ffx,fsx,fwx,ffxp,fsxp,fwxp
   real(8),dimension(ny) :: ffy,fsy,fwy,ffyp,fsyp,fwyp
   real(8),dimension(nz) :: ffz,fsz,fwz,ffzp,fszp,fwzp       
!
   integer  :: nx,ny,nz,i,j,k 
   real(8) :: gdt
!
   call derxb (px,pp,tr,sx,ffxp,fsxp,fwxp,nx,ny,nz,1)
   call dery (py,pp,tr,di,sy,ffyp,fsyp,fwyp,nx,ny,nz,1)
   if (nz.gt.1) then
      call derz (pz,pp,tr,sz,ffzp,fszp,fwzp,nx,ny,nz,1)
   endif
!
   do k=1,nz
   do j=1,ny
   do i=1,nx
      py1(i,j,k)=py(i,j,k)/gdt 
   enddo
   enddo
   enddo
!
   return
end subroutine gradp
!
!********************************************************************
!
subroutine urotu(cx,cy,cz,ux,uy,uz,sp,sv,&
     xnut,xnutx,xnuty,xnutz,&
     tx,ty,tz,tr,di,sx,sy,sz,&
     ffx,fsx,fwx,ffy,fsy,fwy,ffz,fsz,fwz,&
     ffxp,fsxp,fwxp,&
     ffyp,fsyp,fwyp,&
     ffzp,fszp,fwzp,&
     nx,ny,nz)
! 
!********************************************************************
!
!
   USE ecoulement_m
!
   implicit none
!
   integer  :: nx 
   integer  :: ny 
   integer  :: nz 
!
   real(8),dimension(nx,ny,nz) :: ux,uy,uz,cx,cy,cz,tx,ty,tz
   real(8),dimension(nx,ny,nz) :: tr,sp,sv,di,xnut,xnutx,xnuty,xnutz
   real(8),dimension(ny,nz) :: sx
   real(8),dimension(nx,nz) :: sy
   real(8),dimension(nx,ny) :: sz
   real(8),dimension(nx) :: ffx,fsx,fwx,ffxp,fsxp,fwxp
   real(8),dimension(ny) :: ffy,fsy,fwy,ffyp,fsyp,fwyp
   real(8),dimension(nz) :: ffz,fsz,fwz,ffzp,fszp,fwzp     
!
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: nxyz,ijk,j,i,k 
!-----------------------------------------------   
!
   call derx (tz,uy,tr,sx,ffxp,fsxp,fwxp,nx,ny,nz,1)
   call dery (cz,ux,tr,di,sy,ffyp,fsyp,fwyp,nx,ny,nz,1)
   do j=1,ny
   do i=1,nx
      tz(i,j,1)= tz(i,j,1)-cz(i,j,1)
      cx(i,j,1)=-tz(i,j,1)*uy(i,j,1)
      cy(i,j,1)= tz(i,j,1)*ux(i,j,1)
   enddo
   enddo
!
   return
end subroutine urotu
!
!
!********************************************************************
!
subroutine rk4 (ux,uy,uz,sp,gx,gy,gz,sg,hx,hy,hz,sh,&
     vx,vy,vz,sv,nxyz,adt,bdt,itr)
! 
!********************************************************************
!
!
!
   USE ecoulement_m
!
   implicit none
!
   real(8),dimension(nxyz) :: ux,uy,uz,vx,vy,vz,gx,gy,gz,hx,hy,hz
   real(8),dimension(nxyz) :: sp,sg,sv,sh
!
   integer,intent(in) :: nxyz 
   integer,intent(in) :: itr 
   real(8),intent(in) :: adt 
   real(8),intent(in) :: bdt 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
   integer :: ijk 
!-----------------------------------------------
!
   if (iscalaire.eq.0) then
      if (itr.eq.1) then
         do ijk=1,nxyz
            ux(ijk)=adt*vx(ijk)+hx(ijk)
            uy(ijk)=adt*vy(ijk)+hy(ijk)
            uz(ijk)=adt*vz(ijk)+hz(ijk)
            gx(ijk)=bdt*vx(ijk)+hx(ijk)
            gy(ijk)=bdt*vy(ijk)+hy(ijk)
            gz(ijk)=bdt*vz(ijk)+hz(ijk)
         enddo
      else
         if (itr.eq.4) then
            do ijk=1,nxyz
               ux(ijk)=bdt*vx(ijk)+gx(ijk)
               uy(ijk)=bdt*vy(ijk)+gy(ijk)
               uz(ijk)=bdt*vz(ijk)+gz(ijk)
            enddo
         else
            do ijk=1,nxyz
               ux(ijk)=adt*vx(ijk)+hx(ijk)
               uy(ijk)=adt*vy(ijk)+hy(ijk)
               uz(ijk)=adt*vz(ijk)+hz(ijk)
               gx(ijk)=bdt*vx(ijk)+gx(ijk)
               gy(ijk)=bdt*vy(ijk)+gy(ijk)
               gz(ijk)=bdt*vz(ijk)+gz(ijk)
            enddo
         endif
      endif
   else
      if (itr.eq.1) then
         do ijk=1,nxyz
            ux(ijk)=adt*vx(ijk)+hx(ijk)
            uy(ijk)=adt*vy(ijk)+hy(ijk)
            uz(ijk)=adt*vz(ijk)+hz(ijk)
            sp(ijk)=adt*sv(ijk)+sh(ijk)
            gx(ijk)=bdt*vx(ijk)+hx(ijk)
            gy(ijk)=bdt*vy(ijk)+hy(ijk)
            gz(ijk)=bdt*vz(ijk)+hz(ijk)
            sg(ijk)=bdt*sv(ijk)+sh(ijk)
         enddo
      else
         if (itr.eq.4) then
            do ijk=1,nxyz
               ux(ijk)=adt*vx(ijk)+gx(ijk)
               uy(ijk)=adt*vy(ijk)+gy(ijk)
               uz(ijk)=adt*vz(ijk)+gz(ijk)
               sp(ijk)=adt*sv(ijk)+sg(ijk)
            enddo
         else
            do ijk=1,nxyz
               ux(ijk)=adt*vx(ijk)+hx(ijk)
               uy(ijk)=adt*vy(ijk)+hy(ijk)
               uz(ijk)=adt*vz(ijk)+hz(ijk)
               sp(ijk)=adt*sv(ijk)+sh(ijk)
               gx(ijk)=bdt*vx(ijk)+gx(ijk)
               gy(ijk)=bdt*vy(ijk)+gy(ijk)
               gz(ijk)=bdt*vz(ijk)+gz(ijk)
               sg(ijk)=bdt*sv(ijk)+sg(ijk)
            enddo
         endif
      endif
   endif
!
   return
end subroutine rk4
!
!********************************************************************
!
subroutine rkutta (ux,uy,uz,sp,gx,gy,gz,sg,&
     vx,vy,vz,sv,nxyz,adt,bdt,nz)
! 
!******************************************************************
!
!
!
   USE ecoulement_m
!
   implicit none
!
   real(8),dimension(nxyz) :: ux,uy,uz,vx,vy,vz,gx,gy,gz
   real(8),dimension(nxyz) :: sp,sg,sv
!
   integer,intent(in) :: nxyz 
   integer,intent(in) :: nz 
   real(8),intent(in) :: adt 
   real(8),intent(in) :: bdt 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
   integer :: ijk 
!-----------------------------------------------
!
   if (iscalaire.eq.0) then
      if (nz.gt.1) then
         if (bdt.ne.0.) then
            do ijk=1,nxyz
               ux(ijk)=adt*vx(ijk)+bdt*gx(ijk)+ux(ijk) 
               uy(ijk)=adt*vy(ijk)+bdt*gy(ijk)+uy(ijk)
               uz(ijk)=adt*vz(ijk)+bdt*gz(ijk)+uz(ijk)
               gx(ijk)=vx(ijk)
               gy(ijk)=vy(ijk)
               gz(ijk)=vz(ijk)
            enddo
         else
            do ijk=1,nxyz
               ux(ijk)=adt*vx(ijk)+ux(ijk) 
               uy(ijk)=adt*vy(ijk)+uy(ijk)
               uz(ijk)=adt*vz(ijk)+uz(ijk)
               gx(ijk)=vx(ijk)
               gy(ijk)=vy(ijk)
               gz(ijk)=vz(ijk)
            enddo
         endif
      else
         if (bdt.ne.0.) then
            do ijk=1,nxyz
               ux(ijk)=adt*vx(ijk)+bdt*gx(ijk)+ux(ijk) 
               uy(ijk)=adt*vy(ijk)+bdt*gy(ijk)+uy(ijk)
               gx(ijk)=vx(ijk)
               gy(ijk)=vy(ijk)
            enddo
         else
            do ijk=1,nxyz
               ux(ijk)=adt*vx(ijk)+ux(ijk) 
               uy(ijk)=adt*vy(ijk)+uy(ijk)
               gx(ijk)=vx(ijk)
               gy(ijk)=vy(ijk)
            enddo
         endif
      endif
   else
      if (nz.gt.1) then
         if (bdt.ne.0.) then
            do ijk=1,nxyz
               ux(ijk)=adt*vx(ijk)+bdt*gx(ijk)+ux(ijk) 
               uy(ijk)=adt*vy(ijk)+bdt*gy(ijk)+uy(ijk)
               uz(ijk)=adt*vz(ijk)+bdt*gz(ijk)+uz(ijk)
               sp(ijk)=adt*sv(ijk)+bdt*sg(ijk)+sp(ijk)
               gx(ijk)=vx(ijk)
               gy(ijk)=vy(ijk)
               gz(ijk)=vz(ijk)
               sg(ijk)=sv(ijk)
            enddo
         else
            do ijk=1,nxyz
               ux(ijk)=adt*vx(ijk)+ux(ijk) 
               uy(ijk)=adt*vy(ijk)+uy(ijk)
               uz(ijk)=adt*vz(ijk)+uz(ijk)
               sp(ijk)=adt*sv(ijk)+sp(ijk)
               gx(ijk)=vx(ijk)
               gy(ijk)=vy(ijk)
               gz(ijk)=vz(ijk)
               sg(ijk)=sv(ijk)
            enddo
         endif
      else
         if (bdt.ne.0.) then
            do ijk=1,nxyz
               ux(ijk)=adt*vx(ijk)+bdt*gx(ijk)+ux(ijk) 
               uy(ijk)=adt*vy(ijk)+bdt*gy(ijk)+uy(ijk)
               sp(ijk)=adt*sv(ijk)+bdt*sg(ijk)+sp(ijk)
               gx(ijk)=vx(ijk)
               gy(ijk)=vy(ijk)
               sg(ijk)=sv(ijk)
            enddo
         else
            do ijk=1,nxyz
               ux(ijk)=adt*vx(ijk)+ux(ijk) 
               uy(ijk)=adt*vy(ijk)+uy(ijk)
               sp(ijk)=adt*sv(ijk)+sp(ijk)
               gx(ijk)=vx(ijk)
               gy(ijk)=vy(ijk)
               sg(ijk)=sv(ijk)
            enddo
         endif
      endif
   endif
!
   return
end subroutine rkutta
!
!*******************************************************************
!
!  Entrada: ux, uy, vx, vy, gx, gy. Salida: u8x, u8y, gx, gy
subroutine adams (ux,uy,u8x,u8y,gx,gy,vx,vy,nxyz,adt,bdt)
!
!*******************************************************************
!
   USE ecoulement_m
!
   implicit none
!
   real(8),dimension(nxyz) :: ux,uy,u8x,u8y,vx,vy,gx,gy
!
   integer,intent(in) :: nxyz 
   real(8),intent(in) :: adt 
   real(8),intent(in) :: bdt 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
   integer :: ijk 
!-----------------------------------------------
!
   do ijk=1,nxyz
      u8x(ijk)=adt*vx(ijk)+bdt*gx(ijk)+ux(ijk) 
      u8y(ijk)=adt*vy(ijk)+bdt*gy(ijk)+uy(ijk)
      gx(ijk)=vx(ijk)
      gy(ijk)=vy(ijk)
   enddo
!
   return
end subroutine adams
!
!********************************************************************
!
subroutine filtrage (cx,cy,cz,ux,uy,uz,sp,sv,&
     tx,ty,tz,tr,di,sx,sy,sz,vx,vy,vz,&
     fiffx,fifx,ficx,fibx,fibbx,&
     fiffxp,fifxp,ficxp,fibxp,fibbxp,&
     fiffy,fify,ficy,fiby,fibby,&
     fiffyp,fifyp,ficyp,fibyp,fibbyp,&
     fiffz,fifz,ficz,fibz,fibbz,&
     fiffzp,fifzp,ficzp,fibzp,fibbzp,&
     fiz1x,fiz2x,filax,filaxp,&
     fiz1y,fiz2y,filay,filayp,&
     fiz1z,fiz2z,filaz,filazp,&
     nx,ny,nz)
! 
!********************************************************************
!
  implicit none
!
  integer  :: nx 
  integer  :: ny 
  integer  :: nz   
!
   real(8),dimension(nx) :: fifx,ficx,fibx,fiffx,fibbx,fiz1x,fiz2x
   real(8),dimension(nx,2) ::filax,filaxp
   real(8),dimension(nx) :: fifxp,ficxp,fibxp,fiffxp,fibbxp
   real(8),dimension(ny) :: fify,ficy,fiby,fiffy,fibby,fiz1y,fiz2y
   real(8),dimension(ny,2) ::filay,filayp
   real(8),dimension(ny) :: fifyp,ficyp,fibyp,fiffyp,fibbyp
   real(8),dimension(nz) :: fifz,ficz,fibz,fiffz,fibbz,fiz1z,fiz2z
   real(8),dimension(nz,2) ::filaz,filazp
   real(8),dimension(nz) :: fifzp,ficzp,fibzp,fiffzp,fibbzp
   real(8),dimension(ny,nz) :: sx,vx 
   real(8),dimension(nx,nz) :: sy,vy
   real(8),dimension(nx,ny) :: sz,vz
   real(8),dimension(nx,ny,nz) :: ux,uy,uz,cx,cy,cz,tx,ty,tz,tr,sp,sv,di
!

!
   call filx(cx,ux,tr,sx,vx,fiffx,fifx,ficx,fibx,fibbx,&
        filax,fiz1x,fiz2x,nx,ny,nz,0)
   call filx(cy,uy,tr,sx,vx,fiffx,fifx,ficx,fibx,fibbx,&
        filax,fiz1x,fiz2x,nx,ny,nz,0)
!
   call fily(tx,cx,tr,sy,vy,fiffyp,fifyp,ficyp,fibyp,fibbyp,&
        filayp,fiz1y,fiz2y,nx,ny,nz,1)
   call fily(ty,cy,tr,sy,vy,fiffy,fify,ficy,fiby,fibby,&
        filay,fiz1y,fiz2y,nx,ny,nz,0)
!
   if (nz.gt.1) then
      call filx(cz,uz,tr,sx,vx,fiffx,fifx,ficx,fibx,fibbx,&
           filax,fiz1x,fiz2x,nx,ny,nz,0)
      call fily(tz,cz,tr,sy,vy,fiffy,fify,ficy,fiby,fibby,&
           filay,fiz1y,fiz2y,nx,ny,nz,0)
      call filz(ux,tx,tr,sz,vz,fiffzp,fifzp,ficzp,fibzp,fibbzp,&
           filazp,fiz1z,fiz2z,nx,ny,nz,1)
      call filz(uy,ty,tr,sz,vz,fiffzp,fifzp,ficzp,fibzp,fibbzp,&
           filazp,fiz1z,fiz2z,nx,ny,nz,1)
      call filz(uz,tz,tr,sz,vz,fiffz,fifz,ficz,fibz,fibbz,&
           filaz,fiz1z,fiz2z,nx,ny,nz,0)
   endif
!
   return
end subroutine filtrage
!
!********************************************************************
!
subroutine filtrage1 (ux,uy,tx,ty,&
     fiffx,fifx,ficx,fibx,fibbx,&
     fiffxp,fifxp,ficxp,fibxp,fibbxp,&
     fiffy,fify,ficy,fiby,fibby,&
     fiffyp,fifyp,ficyp,fibyp,fibbyp,&
     fiz1x,fiz2x,filax,filaxp,&
     fiz1y,fiz2y,filay,filayp,&
     nx,ny,nz)
! 
!********************************************************************
!
  implicit none
!
  integer :: nx,ny,nz   
!
  real(8),dimension(nx) :: fifx,ficx,fibx,fiffx,fibbx,fiz1x,fiz2x
  real(8),dimension(nx,2) ::filax,filaxp
  real(8),dimension(nx) :: fifxp,ficxp,fibxp,fiffxp,fibbxp
  real(8),dimension(ny) :: fify,ficy,fiby,fiffy,fibby,fiz1y,fiz2y
  real(8),dimension(ny,2) ::filay,filayp
  real(8),dimension(ny) :: fifyp,ficyp,fibyp,fiffyp,fibbyp
  real(8),dimension(ny,nz) :: sx,vx 
  real(8),dimension(nx,nz) :: sy,vy
  real(8),dimension(nx,ny,nz) :: ux,uy,cx,cy,tx,ty,tr
!
  call filx(cx,ux,tr,sx,vx,fiffx,fifx,ficx,fibx,fibbx,&
         filax,fiz1x,fiz2x,nx,ny,nz,0)
  call filx(cy,uy,tr,sx,vx,fiffx,fifx,ficx,fibx,fibbx,&
         filax,fiz1x,fiz2x,nx,ny,nz,0)
!
  call fily(tx,cx,tr,sy,vy,fiffyp,fifyp,ficyp,fibyp,fibbyp,&
         filayp,fiz1y,fiz2y,nx,ny,nz,1)
  call fily(ty,cy,tr,sy,vy,fiffy,fify,ficy,fiby,fibby,&
         filay,fiz1y,fiz2y,nx,ny,nz,0)
!
  return
end subroutine filtrage1
!
!********************************************************************
!
subroutine filtrage2 (uy,ty,&
     fiffx,fifx,ficx,fibx,fibbx,&
     fiffxp,fifxp,ficxp,fibxp,fibbxp,&
     fiffy,fify,ficy,fiby,fibby,&
     fiffyp,fifyp,ficyp,fibyp,fibbyp,&
     fiz1x,fiz2x,filax,filaxp,&
     fiz1y,fiz2y,filay,filayp,&
     nx,ny,nz)
! 
!********************************************************************
!
  implicit none
!
  integer :: nx,ny,nz   
!
  real(8),dimension(nx) :: fifx,ficx,fibx,fiffx,fibbx,fiz1x,fiz2x
  real(8),dimension(nx,2) ::filax,filaxp
  real(8),dimension(nx) :: fifxp,ficxp,fibxp,fiffxp,fibbxp
  real(8),dimension(ny) :: fify,ficy,fiby,fiffy,fibby,fiz1y,fiz2y
  real(8),dimension(ny,2) ::filay,filayp
  real(8),dimension(ny) :: fifyp,ficyp,fibyp,fiffyp,fibbyp
  real(8),dimension(ny,nz) :: sx,vx 
  real(8),dimension(nx,nz) :: sy,vy
  real(8),dimension(nx,ny,nz) :: uy,cy,ty,tr
!
  call filx(cy,uy,tr,sx,vx,fiffx,fifx,ficx,fibx,fibbx,&
         filax,fiz1x,fiz2x,nx,ny,nz,0)
!
  call fily(ty,cy,tr,sy,vy,fiffy,fify,ficy,fiby,fibby,&
         filay,fiz1y,fiz2y,nx,ny,nz,0)
!
  return
end subroutine filtrage2
! 
!********************************************************************
!
subroutine corgp (ux,uy,uz,px,py,pz,&
     bx1,by1,bz1,bxn,byn,bzn,&
     nx,ny,nz)
! 
!********************************************************************
!
!
!
   USE paramx_m
!
   implicit none
!
   real(8),dimension(nx,ny,nz) :: ux,uy,uz,px,py,pz
   real(8),dimension(ny,nz) :: bx1,by1,bz1,bxn,byn,bzn
   integer,intent(in) :: nx 
   integer,intent(in) :: ny 
   integer,intent(in) :: nz 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
   integer :: k,j,i 
!-----------------------------------------------
!
   if (nz.gt.1) then
      do k=1,nz
      do j=1,ny
      do i=1,nx
         ux(i,j,k)=-px(i,j,k)+ux(i,j,k)
         uy(i,j,k)=-py(i,j,k)+uy(i,j,k) 
         uz(i,j,k)=-pz(i,j,k)+uz(i,j,k) 
      enddo
      enddo
      enddo
   else
      do j=1,ny
      do i=1,nx
         ux(i,j,1)=-px(i,j,1)+ux(i,j,1)
         uy(i,j,1)=-py(i,j,1)+uy(i,j,1) 
      enddo
      enddo
   endif
!
   return
end subroutine corgp
!
!********************************************************************
!
subroutine div (pp,ux,uy,uz,tx,ty,tz,tr,di,sx,sy,sz,t,r,&
     ffx,fsx,fwx,ffy,fsy,fwy,ffz,fsz,fwz,&
     ffxp,fsxp,fwxp,&
     ffyp,fsyp,fwyp,&
     ffzp,fszp,fwzp,&
     nx,ny,nz,nvec,itime)
! 
!********************************************************************
!
!
  implicit none
!
   real(8),dimension(nx,ny,nz) :: ux,uy,uz,tx,ty,tz,tr,pp,di
   real(8),dimension(nvec) :: t,r 
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
   integer,intent(in) :: nvec ,itime
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: imax,jmax,kmax,i,nxyz,ijk,l,j 
      real(8) :: umax,tmax,tmoy 
!-----------------------------------------------   
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
     

   call derxb (tx,ux,tr,sx,ffx,fsx,fwx,nx,ny,nz,0)
!       if (itime.eq.1) then
!       open (31,file='courbe.dat',form='formatted',status='unknown')
!           do i=1,nx
!              x=(i-1)*dx
!              write(31,*) x,tx(i,10,1)             
!           enddo
!           close(31)
!           call system('gnuplot courbe.gnu')
!           stop
!           endif
   call dery (ty,uy,tr,di,sy,ffy,fsy,fwy,nx,ny,nz,0)
     
!
   if (nz.gt.1) then
      call derz (tz,uz,tr,sz,ffz ,fsz ,fwz ,nx,ny,nz,0)
      nxyz=nx*ny*nz
      do ijk=1,nxyz
         pp(ijk,1,1)=tx(ijk,1,1)+ty(ijk,1,1)+tz(ijk,1,1)
      enddo
      do ijk=1,nxyz,nvec
         imax=min(ijk+nvec-1,nxyz)
         l=imax-ijk+1
         do i=1,l
            t(i)=t(i)+abs(pp(ijk+i-1,1,1))
            r(i)=(r(i)+pp(ijk+i-1,1,1))/2.&
                 +abs(r(i)-pp(ijk+i-1,1,1))/2.
         enddo
      enddo
      do i=1,nvec
         tmoy=tmoy+t(i)
         tmax=(tmax+r(i))/2.+abs(tmax-r(i))/2.
      enddo
   else
      do j=1,ny
      do i=1,nx
         pp(i,j,1)=tx(i,j,1)+ty(i,j,1)
         tmoy=tmoy+abs(pp(i,j,1))
         tmax=(tmax+pp(i,j,1))/2.+abs(tmax-pp(i,j,1))/2.
      enddo
      enddo
      kmax=1
   endif
!
   tmoy=tmoy/(nx*ny*nz)
   print *,'divu = ',tmax,tmoy
!
   return
end subroutine div
!
!*********************************************************
!
subroutine sortie (ux,uy,bxn,byn,nx,ny,nz,gdt)
!
!*********************************************************
!
!
   USE paramx_m
   USE paramy_m
   USE paramz_m
   USE ecoulement_m
   USE delta_m
!
   implicit none
!
   real(8),dimension(nx,ny,nz) :: ux,uy
   real(8),dimension(ny,nz) :: bxn,byn
!
   integer,intent(in) :: nx,ny,nz 
   real(8),intent(in) :: gdt 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
   integer :: j 
   real(8) :: udx,vphase,cx 
!-----------------------------------------------
!
   udx=1./dx
   vphase=0.5*((2.*u2)-1.)
   cx=vphase*gdt*udx
!
   do j=1,ny
      bxn(j,1)=ux(nx,j,1)-cx*(ux(nx,j,1)-ux(nx-1,j,1))
      byn(j,1)=uy(nx,j,1)-cx*(uy(nx,j,1)-uy(nx-1,j,1))
   enddo
!
   return
end subroutine sortie
!
!**********************************************************************
!
subroutine ecoule (bx1,by1,bz1,bs1,nx,ny,nz,t,yp,jj,j1,j2,j3,j4,j5,iimin,i5,i55,h,angle3)
!
!**********************************************************************
!
!
!
   USE paramx_m
   USE paramy_m
   USE paramz_m
   USE ecoulement_m
   USE delta_m
   USE controle_m
   USE barreau_m
!
   implicit none
!
   real(8),dimension(ny,nz) :: bx1,by1,bz1,bs1
   real(8),dimension(ny+1) :: yp 
!
   integer :: nx 
   integer :: ny 
   integer :: nz
   real(8) ::t    
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
   integer :: i,j,k
   integer :: ii,i21,i31,i5,i6,i7,i55,iimin
   integer :: j1,j2,j21,j3,j4,j5,jj,jj1,jj2 
   real(8) :: pi,twopi,y,omega,ct1,st1,ct2,st2,sot1,sot2,yl1,zl1
   real(8) :: yl2,yc1,yp1,yc2,yp2,r1,r2,um1,um2,ut,z,r,zl2,ut1,ut2
   real(8) :: zc1,zc2,uh,ud,um,zp1,s0,so,yplaca1,yplaca2,xlb1,xlb2,ym
   real(8) :: val,vmax,b,c,d,e,zlzp,zlzph,udszlzp,udsczlzp
   real(8) :: h,yy,angle3,coef,A,f,dt,alpha,dui
!-----------------------------------------------
!
   pi=acos(-1.)
   twopi=2.*acos(-1.)

!
   if (iecoule.eq.1) then
      do j=1,ny
         y=abs((j-1)*delty-y1)
         if (y.eq.0.) then
            bx1(j,1)=u1
         else
            bx1(j,1)=u2+((u1-u2)/2.)*&
                 (1.-tanh(rayon/4.*&
                 (y/rayon-rayon/y)))
         endif
         by1(j,1)=0.
      enddo
!
   endif
   if (iecoule.eq.2) then
      if (ientree.ne.3) then
         do j=1,ny
            y=abs((j-1)*delty-y1)
            if (y.eq.0.) then
               bx1(j,1)=u1
            else
               bx1(j,1)=u2+((u1-u2)/2.)*&
                    (1.-tanh(rayon/4.*&
                    (y/rayon-rayon/y)))
            endif
         enddo
      endif
      omega=twopi/12.
      ct1=cos(theta1)
      st1=sin(theta1)
      ct2=cos(theta2)
      st2=sin(theta2)
      sot1=1.+0.2*sin(omega*t)
      sot1=1.
      sot2=1.-0.2*sin(omega*t)
      sot2=1.
      yl1=yly/2.-distance*rayon
      zl1=zlz/2.
      yl2=yly/2.+distance*rayon
!      
      do j=1,ny
         yc1=(j-1)*delty-yl1
         yp1=ct1*yc1
         yc2=(j-1)*delty-yl2
         yp2=ct2*yc2
         r1=abs(yp1)
         r2=abs(yp2)
         if (r1.eq.0.) then
            um1=u1-u2
         else
            um1=((u1-u2)/2.)*&
                 (1.-tanh(rlat/4.*&
                 (r1/rlat-rlat/r1)))
         endif
         if (r2.eq.0.) then
            um2=u1-u2
         else
            um2=((u1-u2)/2.)*&
                 (1.-tanh(rlat/4.*&
                 (r2/rlat-rlat/r2)))
         endif
         bx1(j,1)=u2+um1*ct1
         by1(j,1)=um1*st1+um2*st2
      enddo
   endif
   if (iecoule.eq.3) then
      uampl=u1/2.
      ufreq=1/20.
      ut=u1+uampl*sin(t*ufreq)
!
      do k=1,nz
      do j=1,ny
         y=(j-1)*delty-y1
         z=(k-1)*deltz-z1
         r=sqrt(y*y+z*z)
         if (r.eq.0.) then
            bx1(j,k)=u1
         else
            bx1(j,k)=u2+((u1-u2)/2.)*&
                 (1.-tanh(rayon/4.*&
                 (r/rayon-rayon/r)))
         endif
         by1(j,k)=0.
         bz1(j,k)=0.
      enddo
      enddo
   endif
   if (iecoule.eq.4) then
      print *,'t=',t
      if (ientree.ne.3) then
         do k=1,nz
         do j=1,ny
            y=(j-1)*delty-y1
            z=(k-1)*deltz-z1
            r=sqrt(y*y+z*z)
            if (r.eq.0.) then
               bx1(j,k)=u1
            else
               bx1(j,k)=u2+((u1-u2)/2.)*&
                    (1.-tanh(rayon/4.*&
                    (r/rayon-rayon/r)))
            endif
            by1(j,k)=0.
            bz1(j,k)=0.
         enddo
         enddo
      else
         do k=1,nz
         do j=1,ny
            bx1(j,k)=0.
            by1(j,k)=0.
            bz1(j,k)=0.
         enddo
         enddo
      endif
      theta1=(0.5*pi-phi)
      theta2=-(0.5*pi-phi)
      ct1=cos(theta1)
      st1=sin(theta1)
      ct2=cos(theta2)
      st2=sin(theta2)
      sot1=1.
      sot2=1.
      yl1=y1-distance*rayon
      zl1=z1
      yl2=y1+distance*rayon
      zl2=z1
      call control (ut1,ut2,t)
!
      do k=1,nz
      do j=1,ny
         yc1=(j-1)*delty-yl1
         zc1=(k-1)*deltz-zl1
         yp1=ct1*yc1
         yc2=(j-1)*delty-yl2
         zc2=(k-1)*deltz-zl2
         yp2=ct2*yc2
         r1=sqrt(yp1*yp1+zc1*zc1)
         r2=sqrt(yp2*yp2+zc2*zc2)
         if (r1.eq.0.) then
            um1=ut1-u2
         else
            um1=((ut1-u2)/2.)*&
                 (1.-tanh(rlat/4.*&
                 (r1/rlat-rlat/r1)))
         endif
         if (r2.eq.0.) then
            um2=ut2-u2
         else
            um2=((ut2-u2)/2.)*&
                 (1.-tanh(rlat/4.*&
                 (r2/rlat-rlat/r2)))
         endif
         bx1(j,k)=bx1(j,k)+um1*ct1+um2*ct2
         by1(j,k)=um1*st1+um2*st2
         bz1(j,k)=0.
      enddo
      enddo
   endif
   if (iecoule.eq.5) then
      uh=0.5*(u1+u2)
      ud=0.5*(u1-u2)
      if (nz.gt.1) then
         do k=1,nz
         do j=1,ny
            z=(k-1)*deltz-z1
            bx1(j,k)=uh+ud*tanh(2.*z)
            by1(j,k)=0.
            bz1(j,k)=0.
            bs1(j,k)=0.5+0.5*tanh(2.*z)
         enddo
         enddo
      else
         do j=1,ny
            y=(j-1)*delty-y1
            bx1(j,1)=uh+ud*tanh(2.*y)
            by1(j,1)=0.
            bz1(j,1)=0.
         enddo
      endif
      if (icontrole.eq.1) then
         theta1=(0.5*pi-phi)
         ct1=cos(theta1)
         st1=sin(theta1)
         sot1=1.
         yl1=yly/2.
         zl1=zlz/2.-distance*3.
         call control (ut1,ut2,t)
         do k=1,nz
         do j=1,ny
            yc1=(j-1)*yly/ny-yl1
            zc1=(k-1)*deltz-zl1
            zp1=ct1*zc1
            r1=sqrt(yc1*yc1+zp1*zp1)
            if (r1.eq.0.) then
               um1=ut1-u2
            else
               um1=((ut1-u2)/2.)*&
                    (1.-tanh(rlat*&
                    (r1/rlat-rlat/r1)))
            endif
            bx1(j,k)=bx1(j,k)+um1*ct1
            by1(j,k)=0.
            bz1(j,k)=um1*st1
         enddo
         enddo
         if (iscalaire.eq.1) then
            do k=1,nz
            do j=1,ny
               bs1(j,k)=um1
            enddo
            enddo
         endif
      endif
   endif
!
   if (iscalaire.eq.1 .and. iecoule.lt.5) then
      if (nz.gt.1) then
         do k=1,nz
         do j=1,ny
            y=(j-1)*delty-y1
            z=(k-1)*deltz-z1
            r=sqrt(y*y+z*z)
            if (r.eq.0.) then
               bs1(j,k)=u1-u2
            else
               bs1(j,k)=((u1-u2)/2.)*&
                    (1.-tanh(rayon/4.*&
                    (r/rayon-rayon/r)))
            endif
         enddo
         enddo
      else
         do j=1,ny
            y=abs((j-1)*delty-y1)
            if (y.eq.0.) then
               bs1(j,1)=u1-u2
            else
               bs1(j,1)=((u1-u2)/2.)*&
                    (1.-tanh(rayon/4.*&
                    (y/rayon-rayon/y)))
            endif
         enddo
      endif
   endif
!************************************************************
   if (iecoule.eq.6) then
      s0=1.
      if (nz.gt.1) then
         do k=1,nz
            z=(k-1)*deltz-z1
            do j=1,ny
               y=yp(j)-y1
               bx1(j,k)=1-s0*exp(-log(2.)*y**2)
               by1(j,k)=0.
               bz1(j,k)=0.
            enddo
         enddo
      else
         do j=1,ny
            y=yp(j)-y1
!               bx1(j,1)=1-s0*exp(-log(2.)*y**2)
            bx1(j,1)=1.-s0/2.*(1.-tanh((abs(y)-1)/0.33))
            by1(j,1)=0.
            bz1(j,1)=0.
         enddo
      endif
   endif
!************************************************************
!************************************************************
   if (iecoule.eq.7) then   !ecoulement constant
!************************************************************
!  ra c'est le rayon de ma plaque      
   j21=1
   do j=ny,jj,-1
      r=abs(yp(j)-cey)
      if (r.gt.(h/cos(phi))) j21=j+1
   enddo

!**********************************************************
!*******EPAISSEUR DES COUCHES LIMITES**********************
!**********************************************************
! poiseuille 
!   do j=j1,j2
!      bx1(j,1)=0.
!      by1(j,1)=0.
!   enddo
!
   val=0.
   do j=j2+1,j3-1
      if (val.lt.(yp(j)-yp(j2+1))*(yp(j3-1)-yp(j))) then 
      val=(yp(j)-yp(j2+1))*(yp(j3-1)-yp(j))
      endif
   enddo

! Ecoulement constant
!      um=0.5*(u1+u2)
!      do j=1,ny
!         bx1(j,1)=um
!         by1(j,1)=0.
!      enddo
!***********TRAVAIL SUR LA VITESSE U U U U U U************************
!12/19 c'est delta99 de perret
!poly d'ordre 6 soluce 0.772588621523 ( soluce de u/ue-0.99=0 )
!delta 100 egal a 0.81 en gros
!**********************************************************************
!on a la relation theta=985/9009 * delta99 ( theta=0.109 * delta99 )
!on veut theta=1.23 ce qui donne delta99=11.2498172589
!**********************************************************************      
      jj2=j4
      do j=j4,ny
         r=abs(yp(j)-yp(j4))
         if (r.le.(5./19.)) jj2=jj2+1
      enddo
      do j=j4,jj2
         ym=(yp(j)-yp(j4))/(yp(jj2)-yp(j4))
         ym=ym*0.772588621523
         bx1(j,1)=2.*ym*u1-5.*ym*ym*ym*ym*u1+6.*ym*ym*ym*ym*ym*u1-2.*ym*ym*ym*ym*ym*ym*u1
         if ((j>j2+1).and.(bx1(j,1) < bx1(j-1,1))) bx1(j,1)=u1
      enddo
      do j=jj2+1,ny
         bx1(j,1)=u1
      enddo
!***********************
!      do j=j1,j2
!         bx1(j,1)=0. Pas dans le cas d un poiseuille
!      enddo
!***********************
      do j=j5,j21
         bx1(j,1)=0.
         by1(j,1)=0.
      enddo
!
      do j=j21+1,j3-1
         um=((yp(j)-yp(j2+1))*(yp(j3-1)-yp(j))/val)
         bx1(j,1)=um*u2*cos(angle3)
         by1(j,1)=um*u2*sin(angle3)
      enddo
!
      do j=j3,j4
         bx1(j,1)=0.
         by1(j,1)=0.
      enddo
!
!**********************************************************************
!on a la relation theta=985/9009 * delta99 ( theta=0.109 * delta99 )
!on veut theta=1.06 ce qui donne delta99=9.69496446701
!**********************************************************************   
      jj1=1
      do j=2,j5
         r=abs(yp(j)-yp(j5))
         if (r.gt.((19./19.))) jj1=j+1
      enddo
!
      do j=j5,jj1,-1
         ym=(yp(j)-yp(j5))/(yp(jj1)-yp(j5))
         ym=ym*0.772588621523
         bx1(j,1)=2*ym*u2-5.*ym*ym*ym*ym*u2+6*ym*ym*ym*ym*ym*u2-2.*ym*ym*ym*ym*ym*ym*u2 
        if ((j<j1-1).and.(bx1(j,1) < bx1(j+1,1))) bx1(j,1)=u2
      enddo
      do j=1,jj1-1
         bx1(j,1)=u2
      enddo
!      do j=j5,jj1-1
!         write (11,*) yp(j),bx1(j,1),by1(j,1)
!      enddo
!      
!   stop 'Navier iecoule.eq.7'
!*********************************************************************
!***********TRAVAIL SUR LA VITESSE V V V V V V************************
!!!!!!!!!!!!PAS DE TRAVAIL!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*********************************************************************      
!      
   endif
   if (iecoule.eq.8) then   !ecoulement sillage  CAS2D uniquement
      zlzp=99.          !
      zlzph=0.5*zlzp        !
      uh=0.5*(u1+u2)
      ud=0.5*(u1-u2)
      udszlzp=ud/zlzp 
      coef=2.        !Coefficient
      udsczlzp=ud/(coef*zlzp)
!      A=uh/(yp(ny/2)-yp(ny/2-1))
      dui=1.
      A=40.   
      f=0.161
      alpha=0.00025
      do j=1,ny
         ym=yp(j)-(yly*0.5)
!         bx1(j,1)=uh+ud*tanh(ym*2.)
      bx1(j,1)=uh+udsczlzp*log(cosh(coef*ym)/cosh(coef*(ym-zlzp)))&
                    + 0.5*(tanh(2.*ym)-1.)*dui  ! ec. (3.5) tesis toutou
         by1(j,1)=A*exp(-0.2*ym**2.)*(alpha*sin(2.*pi*f*t) &
                   + alpha/2.*sin(2.*pi*(f/2.)*t))  ! ec. (3.6) tesis toutou
         bz1(j,1)=0.
!         write (15,*) ym,bx1(j,1),by1(j,1),(bx1(j,1)+by1(j,1)) 
      enddo
   endif
!   stop
   return
end subroutine ecoule
!
!****************************************************************
!
      subroutine blasius(bx1,xlb,Re,um,yplaca,nx,ny,nz,yp,j1)
!
!!****************************************************************        
        USE paramy_m

      implicit none

      integer neq                ! Numero de equacoes
      integer nx,ny,nz           ! Dimensoes das matrizes
      integer i,j,n              ! Contadores de laco
      integer j1                 ! Coordenada da placa

      parameter (neq=3)


      real(8) h1,h2,h               ! Espacamento vertical de eta
      real(8) yplaca                ! Localizacao da placa plana
      real(8) xlb                   ! Distancia do bordo da placa
      real(8) Re                    ! Numero de Reynolds
      real(8) um                    ! Velocidade de referencia

      real(8) y,dyb,xa
      
      real(8) eta(ny)               ! Coordenada Vertical 
      real(8) y1(ny)                ! Segunda derivada
      real(8) y2(ny)                ! Primeira derivada
      real(8) y3(ny)                ! Funcao

      real(8) k(4,neq)              ! Coeficientes do metodo RK3
      real(8) bx1(ny,nz)
      real(8),dimension(ny+1) :: yp,yp1
      integer :: jj1,jj2

      do j=1,ny
        eta(j) = 0.
        y1(j)  = 0.
        y2(j)  = 0.
        y3(j)  = 0.
!        bx1(j,1) = 0.
      enddo
!###########################################################################
! Calculo dos parametros da camada limite
!###########################################################################
!
!      print *, 'Re_D     = ', Re
!      print *, 'xlb      = ', xlb
!      print *, 'Re_x     = ', xlb*Re
!      print *, 'delta    = ', 5.*xlb/sqrt(xlb*Re)
!
!############################################################################
! Definindo as Condicoes Inicais da Solucao de Blasius
!############################################################################

      !!j1 = 91!nint(yplaca/dyb) + 1

      eta(j1) = 0.               ! Coordenada vertical 
      y1(j1)  = 0.4696           ! Cond. Inicial - Segunda derivada
      y2(j1)  = 0.               ! Cond. Inicial - Primeira derivada
      y3(j1)  = 0.               ! Cond. Inicial - Funcao

!##########################################################################
! Calcula o incremento de eta
!##########################################################################
      dyb = yp(ny/2+1)-yp(ny/2)
      h1 = 0.
      h2 = (1*dyb/xlb)*sqrt(0.5*xlb*Re)
      h = (h2-h1)

!      print *, 'd_eta = ',h
     
      do j=j1,ny-1
        do n=1,neq
          call funcao(k(1,n),eta(j),y1(j),y2(j),y3(j),n)
        enddo
      
        do n=1,neq
          call funcao(k(2,n),eta(j)+0.5*h,y1(j)+0.5*h*k(1,1),&
                         y2(j)+0.5*h*k(1,2),y3(j)+0.5*h*k(1,3),n)
        enddo
      
        do n=1,neq
          call funcao(k(3,n),eta(j)+0.5*h,y1(j)+0.5*h*k(2,1),&
                         y2(j)+0.5*h*k(2,2),y3(j)+0.5*h*k(2,3),n)
        enddo

        do n=1,neq
          call funcao(k(4,n),eta(j)+h,y1(j)+h*k(3,1),&
                         y2(j)+h*k(3,2),y3(j)+h*k(3,3),n)
        enddo

        y1(j+1) = y1(j) + (h*(k(1,1)+2.*k(2,1)+2.*k(3,1)+k(4,1)))/6.
        y2(j+1) = y2(j) + (h*(k(1,2)+2.*k(2,2)+2.*k(3,2)+k(4,2)))/6.
        y3(j+1) = y3(j) + (h*(k(1,3)+2.*k(2,3)+2.*k(3,3)+k(4,3)))/6.
        eta(j+1) = eta(j) + h
        
        bx1(j,1) = y2(j)*um
        
       if(eta(j).gt.11) then
         bx1(j,1) = um
       endif
      enddo

      bx1(ny,1) = bx1(ny-1,1)

!      open(10,file='perfil.dat')
!      do j=j1,ny
!        y = eta(j)*xlb/sqrt(xlb*Re*0.5)
!        write(10,'(f13.7 f13.7 f13.7 f13.7 f13.7)')
!     1                      y,eta(j),y3(j),y2(j),y1(j)
!      enddo
!      close(10)
!******************************************
!      bx1(:,:)=0.
!     do j=j1,ny
!        yp1(j)=yp(j)-yp(j1)
!     enddo
!     do j=j1,ny
!       yp(j)=yp1(j)
!     enddo
!      do j=j1,ny
!         yp(j)=yp(j)/yly
!         bx1(j,1)=2.*yp(j)-2.*yp(j)*yp(j)*yp(j)+9.*yp(j)*yp(j)*yp(j)*yp(j)
!      enddo
!*******************************************
!      do j=1,ny
!         print *,j,bx1(j,1)
!      enddo
!      pause'1'

    end subroutine blasius
!***********************************************************************
!***********************************************************************
!
   subroutine blasiusreverse(bx1,xlb,Re,um,yplaca,nx,ny,nz,yp,j1)
!
!***********************************************************************
     Use paramy_m
     implicit none
!
     integer :: neq                ! Numeros des equations
     integer :: nx,ny,nz,i,j,n        
     integer :: j1                 ! coordonnee de la plaque
!
     parameter (neq=3)
!
      real(8) :: h1,h2,h               ! Espacamento vertical de eta
      real(8) :: yplaca                ! Localizacao da placa plana
      real(8) :: xlb                   ! Distancia do bordo da placa
      real(8) :: Re                    ! Nombre de Reynolds
      real(8) :: um                    ! Vitesse de reference
!
      real(8) :: y,dyb
!      
      real(8) :: eta(ny)               ! Coordenada Vertical 
      real(8) :: y1(ny)                ! derivee seconde
      real(8) :: y2(ny)                ! derivee premiere
      real(8) :: y3(ny)                ! Funcao

      real(8) :: k(4,neq)              ! Coeficientes do metodo RK3
      real(8) :: bx1(ny,nz)
      real(8),dimension(ny+1) :: yp

      do j=1,ny
        eta(j) = 0.
        y1(j)  = 0.
        y2(j)  = 0.
        y3(j)  = 0.
!        bx1(j,1) = 0.
      enddo
!###########################################################################
! Calcul des parametres de la couche limite
!###########################################################################
!
!      print *, 'Re_D     = ', Re
!      print *, 'xlb      = ', xlb
!      print *, 'Re_x     = ', xlb*Re
!      print *, 'delta    = ', 5.*xlb/sqrt(xlb*Re)
!
!############################################################################
! Definition des conditions initiales de la soluction de Blasius
!############################################################################

     

      eta(j1) = 0.               ! Coordonnee verticale 
      y1(j1)  = 0.4696           ! Cond. Initiale - Derivee premiere
      y2(j1)  = 0.               ! Cond. Initiale - Derivee seconde
      y3(j1)  = 0.               ! Cond. Initiale - Funcao

!##########################################################################
! Calcula o incremento de eta
!##########################################################################
      dyb=yp(ny/2+1)-yp(ny/2)
      h1 = 0.
      h2 = (1*dyb/xlb)*sqrt(0.5*xlb*Re)
      h = (h2-h1)

!      print *, 'd_eta = ',h
     
      do j=j1,2,-1
        do n=1,neq
          call funcao(k(1,n),eta(j),y1(j),y2(j),y3(j),n)
        enddo
      
        do n=1,neq
          call funcao(k(2,n),eta(j)+0.5*h,y1(j)+0.5*h*k(1,1),&
                         y2(j)+0.5*h*k(1,2),y3(j)+0.5*h*k(1,3),n)
        enddo
      
        do n=1,neq
          call funcao(k(3,n),eta(j)+0.5*h,y1(j)+0.5*h*k(2,1),&
                         y2(j)+0.5*h*k(2,2),y3(j)+0.5*h*k(2,3),n)
        enddo

        do n=1,neq
          call funcao(k(4,n),eta(j)+h,y1(j)+h*k(3,1),&
                         y2(j)+h*k(3,2),y3(j)+h*k(3,3),n)
        enddo

        y1(j-1) = y1(j) + (h*(k(1,1)+2.*k(2,1)+2.*k(3,1)+k(4,1)))/6.
        y2(j-1) = y2(j) + (h*(k(1,2)+2.*k(2,2)+2.*k(3,2)+k(4,2)))/6.
        y3(j-1) = y3(j) + (h*(k(1,3)+2.*k(2,3)+2.*k(3,3)+k(4,3)))/6.
        eta(j-1) = eta(j) + h
        
        bx1(j,1) = y2(j)*um
        
       if(eta(j).gt.11) then
         bx1(j,1) = um
       endif
      enddo

      bx1(1,1) = bx1(2,1)

   end subroutine blasiusreverse
!**************************************************************************
!********************************************************************
!
subroutine divergence (pp,ux,uy,uz,tx,ty,tz,tr,di,sx,sy,sz,t,r,&
     ffx,fsx,fwx,ffy,fsy,fwy,ffz,fsz,fwz,&
     ffxp,fsxp,fwxp,&
     ffyp,fsyp,fwyp,&
     ffzp,fszp,fwzp,&
     nx,ny,nz,nvec)
! 
!********************************************************************
!
   real(8),dimension(nx,ny,nz) :: ux,uy,uz,tx,ty,tz,tr,pp,di
   real(8),dimension(nvec) :: t,r 
   real(8),dimension(ny,nz) :: sx
   real(8),dimension(nx,nz) :: sy
   real(8),dimension(nx,ny) :: sz
   real(8),dimension(nx) :: ffx,fsx,fwx,ffxp,fsxp,fwxp
   real(8),dimension(ny) :: ffy,fsy,fwy,ffyp,fsyp,fwyp
   real(8),dimension(nz) :: ffz,fsz,fwz,ffzp,fszp,fwzp
!
   call derx (tx,ux,tr,sx,ffx,fsx,fwx,nx,ny,nz,0)
     
   call dery (ty,uy,tr,di,sy,ffy,fsy,fwy,nx,ny,nz,0)
!
   if (nz.gt.1) then
      call derz (tz,uz,tr,sz,ffz ,fsz ,fwz ,nx,ny,nz,0)
      nxyz=nx*ny*nz
      do ijk=1,nxyz
         pp(ijk,1,1)=tx(ijk,1,1)+ty(ijk,1,1)+tz(ijk,1,1)
      enddo
   else
      do j=1,ny
      do i=1,nx
         pp(i,j,1)=tx(i,j,1)+ty(i,j,1)           
      enddo
      enddo
        
   endif
!
   return
end subroutine divergence
! 
!********************************************************************
!
subroutine div1 (pp,ux,uy,uz,tx,ty,tz,tr,di,sx,sy,sz,t,r,&
     ffx,fsx,fwx,ffy,fsy,fwy,ffz,fsz,fwz,&
     ffxp,fsxp,fwxp,&
     ffyp,fsyp,fwyp,&
     ffzp,fszp,fwzp,nx,ny,nz,nxm,nym,mx,my,&
     nwork,ntrigsX,ntrigsY,nvec,itime)
!
!********************************************************************
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
  USE paramx_m 
  USE paramy_m 
  USE paramz_m 
!
  implicit none
!
   real(8),dimension(nx,ny,nz) :: ux,uy,uz,tx,ty,tz,tr,pp,di
   real(8),dimension(nvec) :: t,r 
   real(8),dimension(ny,nz) :: sx
   real(8),dimension(nx,nz) :: sy
   real(8),dimension(nx,ny) :: sz
   real(8),dimension(nx) :: ffx,fsx,fwx,ffxp,fsxp,fwxp
   real(8),dimension(ny) :: ffy,fsy,fwy,ffyp,fsyp,fwyp
   real(8),dimension(nz) :: ffz,fsz,fwz,ffzp,fszp,fwzp
   real(8),dimension(mx,my) :: us2d,ps2d
!
   integer  :: nx,nxm,nym,mx,my,itime 
   integer  :: ny,nwork,ntrigsX,ntrigsY
   integer  :: nz,k
   integer  :: nvec
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: imax,jmax,kmax,i,nxyz,ijk,l,j 
      real(8) :: umax,tmax,tmoy 
!-----------------------------------------------   
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
   us2d(:,:)=0.
   ps2d(:,:)=0.

   call derxb (tx,ux,tr,sx,ffx,fsx,fwx,nx,ny,nz,0)
   call dery (ty,uy,tr,di,sy,ffy,fsy,fwy,nx,ny,nz,0)
     
!
   if (nz.gt.1) then
      call derz (tz,uz,tr,sz,ffz ,fsz ,fwz ,nx,ny,nz,0)
      nxyz=nx*ny*nz
      do ijk=1,nxyz
         pp(ijk,1,1)=tx(ijk,1,1)+ty(ijk,1,1)+tz(ijk,1,1)
      enddo
      do ijk=1,nxyz,nvec
         imax=min(ijk+nvec-1,nxyz)
         l=imax-ijk+1
         do i=1,l
            t(i)=t(i)+abs(pp(ijk+i-1,1,1))
            r(i)=(r(i)+pp(ijk+i-1,1,1))/2.&
                 +abs(r(i)-pp(ijk+i-1,1,1))/2.
         enddo
      enddo
      do i=1,nvec
         tmoy=tmoy+t(i)
        tmax=(tmax+r(i))/2.+abs(tmax-r(i))/2.
      enddo
   else
      do j=1,ny
      do i=1,nx
         pp(i,j,1)=tx(i,j,1)+ty(i,j,1)
         us2d(i,j)=pp(i,j,1)
      enddo
      enddo

!      print*,ny,mx,my,nxm,nym,nwork,ntrigsX,ntrigsY,nclx,ncly
      call SLFFT2D(us2d,nx,ny,mx,my,nxm,nym,nwork,&
           ntrigsX,ntrigsY,nclx,ncly,-1)
!
      us2d(1,1)=0.
      us2d(nx,1)=0.
      us2d(nx,ny)=0.
      do i=1,mx
         us2d(i,my-1)=0.
         us2d(i,my)=0.
      enddo
!      do j=1,my
!         us2d(nx+1,j)=0.
!      enddo
!
!      do i=1,mx
!      do j=1,my
!        print *,us2d(i,j)
!      enddo
!      pause
!      enddo
!      pause
      call SLFFT2D(us2d,nx,ny,mx,my,nxm,nym,nwork,&
           ntrigsX,ntrigsY,nclx,ncly,1)
!
     do j=1,ny
     do i=1,nx
        pp(i,j,1)=us2d(i,j)
     enddo
     enddo
     do k=1,nz
     do j=1,ny
     do i=1,nx
        pp(i,j,k)=pp(i,j,1)
     enddo
     enddo
     enddo
!
      do j=1,ny
      do i=1,nx
         tmoy=tmoy+abs(pp(i,j,1))
         tmax=(tmax+pp(i,j,1))/2.+abs(tmax-pp(i,j,1))/2.
      enddo
      enddo
      kmax=1
   endif
!
   tmoy=tmoy/(nx*ny*nz)
   print *,'divu = ',tmax,tmoy
!
   return
end subroutine div1
!
!*******************************************************************
!
subroutine les (cx,cy,cz,ux,uy,uz,sp,sv,&
     tx,ty,tz,tr,di,sx,sy,sz,vx,vy,vz,&
     fiffx,fifx,ficx,fibx,fibbx,&
     fiffxp,fifxp,ficxp,fibxp,fibbxp,&
     fiffy,fify,ficy,fiby,fibby,&
     fiffyp,fifyp,ficyp,fibyp,fibbyp,&
     fiffz,fifz,ficz,fibz,fibbz,&
     fiffzp,fifzp,ficzp,fibzp,fibbzp,&
     fiz1x,fiz2x,filax,filaxp,&
     fiz1y,fiz2y,filay,filayp,&
     fiz1z,fiz2z,filaz,filazp,&
     nx,ny,nz)
!
!*******************************************************************
!
  USE paramx_m
  USE paramy_m
  USE paramz_m
  USE parfix_m 
  USE parfiy_m 
  USE parfiz_m 
!
   implicit none
!
   integer :: nx,ny,nz,i,j
   real(8) :: pi,x,y
!
   real(8),dimension(nx,ny,nz) :: ux,uy,uz,cx,cy,cz,tx,ty,tz
   real(8),dimension(nx,ny,nz) :: tr,sp,sv,di
   real(8),dimension(ny,nz) :: sx,vx
   real(8),dimension(nx,nz) :: sy,vy
   real(8),dimension(nx,ny) :: sz,vz
   real(8),dimension(nx) :: fifx,ficx,fibx,fiffx,fibbx,fiz1x,fiz2x
   real(8),dimension(nx) :: fifxp,ficxp,fibxp,fiffxp,fibbxp
   real(8),dimension(nx,2) ::filax,filaxp
   real(8),dimension(nx) :: fify,ficy,fiby,fiffy,fibby,fiz1y,fiz2y 
   real(8),dimension(ny) :: fifyp,ficyp,fibyp,fiffyp,fibbyp
   real(8),dimension(ny,2) ::filay,filayp 
   real(8),dimension(nx) :: fifz,ficz,fibz,fiffz,fibbz,fiz1z,fiz2z 
   real(8),dimension(ny) :: fifzp,ficzp,fibzp,fiffzp,fibbzp
   real(8),dimension(ny,2) ::filaz,filazp 
!
!   print *,nz
   if (nz.gt.1) then
      call filx(cx,ux,tr,sx,vx,fiffx,fifx,ficx,fibx,fibbx,&
           filax,fiz1x,fiz2x,nx,ny,nz,0)
      call filx(cy,uy,tr,sx,vx,fiffx,fifx,ficx,fibx,fibbx,&
           filax,fiz1x,fiz2x,nx,ny,nz,1)
      call filx(cz,uz,tr,sx,vx,fiffx,fifx,ficx,fibx,fibbx,&
           filax,fiz1x,fiz2x,nx,ny,nz,1)

      if (ncly.eq.0) then
         call fily(tx,cx,tr,sy,vy,fiffy,fify,ficy,fiby,fibby,&
              filay,fiz1y,fiz2y,nx,ny,nz,1)
         call fily(ty,cy,tr,sy,vy,fiffy,fify,ficy,fiby,fibby,&
              filay,fiz1y,fiz2y,nx,ny,nz,0)
         call fily(tz,cz,tr,sy,vy,fiffy,fify,ficy,fiby,fibby,&
              filay,fiz1y,fiz2y,nx,ny,nz,1)
      endif
      if (ncly.eq.1) then 
         call fily(tx,cx,tr,sy,vy,fiffyp,fifyp,ficyp,fibyp,fibbyp,&
              filayp,fiz1y,fiz2y,nx,ny,nz,1)
         call fily(ty,cy,tr,sy,vy,fiffy,fify,ficy,fiby,fibby,&
              filay,fiz1y,fiz2y,nx,ny,nz,0)
         call fily(tz,cz,tr,sy,vy,fiffyp,fifyp,ficyp,fibyp,fibbyp,&
             filayp,fiz1y,fiz2y,nx,ny,nz,1)
      endif

      if (nclz.eq.0) then
         call filz(ux,tx,tr,sz,vz,fiffz,fifz,ficz,fibz,fibbz,&
              filaz,fiz1z,fiz2z,nx,ny,nz,1)
         call filz(uy,ty,tr,sz,vz,fiffz,fifz,ficz,fibz,fibbz,&
              filaz,fiz1z,fiz2z,nx,ny,nz,1)
         call filz(uz,tz,tr,sz,vz,fiffz,fifz,ficz,fibz,fibbz,&
              filaz,fiz1z,fiz2z,nx,ny,nz,0)
      endif
!
      if (nclz.eq.1) then
         call filz(ux,tx,tr,sz,vz,fiffzp,fifzp,ficzp,fibzp,fibbzp,&
              filazp,fiz1z,fiz2z,nx,ny,nz,1)
         call filz(uy,ty,tr,sz,vz,fiffzp,fifzp,ficzp,fibzp,fibbzp,&
              filazp,fiz1z,fiz2z,nx,ny,nz,1)
         call filz(uz,tz,tr,sz,vz,fiffz,fifz,ficz,fibz,fibbz,&
              filaz,fiz1z,fiz2z,nx,ny,nz,0)
      endif
   else

      call filx(cx,ux,tr,sx,vx,fiffx,fifx,ficx,fibx,fibbx,&
           filax,fiz1x,fiz2x,nx,ny,nz,0)
      call filx(cy,uy,tr,sx,vx,fiffx,fifx,ficx,fibx,fibbx,&
           filax,fiz1x,fiz2x,nx,ny,nz,1)

!
      if (ncly.eq.0) then   
         call fily(ux,cx,tr,sy,vy,fiffy,fify,ficy,fiby,fibby,&
              filay,fiz1y,fiz2y,nx,ny,nz,1)
         call fily(uy,cy,tr,sy,vy,fiffy,fify,ficy,fiby,fibby,&
              filay,fiz1y,fiz2y,nx,ny,nz,0)
      endif
!    
      if (ncly.eq.1) then
         call fily(cx,ux,tr,sy,vy,fiffyp,fifyp,ficyp,fibyp,fibbyp,&
              filayp,fiz1y,fiz2y,nx,ny,nz,1)
         call fily(cy,uy,tr,sy,vy,fiffy,fify,ficy,fiby,fibby,&
              filay,fiz1y,fiz2y,nx,ny,nz,0)
      endif
   endif
!
   return
end subroutine les
! 
!********************************************************************
!
subroutine skew(cx,cy,cz,ux,uy,uz,sp,sv,&
     xnut,xnutx,xnuty,xnutz,hx,hy,hz,&
     tx,ty,tz,tr,di,sx,sy,sz,&
     ffx,fsx,fwx,ffy,fsy,fwy,ffz,fsz,fwz,&
     ffxp,fsxp,fwxp,&
     ffyp,fsyp,fwyp,&
     ffzp,fszp,fwzp,&
     nx,ny,nz,itime)
! 
!********************************************************************
!
   USE ecoulement_m
!
   implicit none
!
   integer :: nx,ny,nz,i,j,k,nxyz,ijk,itime
!
   real(8),dimension(nx,ny,nz) :: ux,uy,uz,cx,cy,cz,hx,hy,hz,tx,ty,tz
   real(8),dimension(nx,ny,nz) :: tr,sp,sv,di,xnut,xnutx,xnuty,xnutz
   real(8),dimension(ny,nz) :: sx
   real(8),dimension(nx,nz) :: sy
   real(8),dimension(nx,ny) :: sz
   real(8),dimension(nx) :: ffx,fsx,fwx,ffxp,fsxp,fwxp
   real(8),dimension(ny) :: ffy,fsy,fwy,ffyp,fsyp,fwyp
   real(8),dimension(ny) :: ffz,fsz,fwz,ffzp,fszp,fwzp      
!
   nxyz=nx*ny*nz
!


   if (iscalaire.eq.1) then
      if (nz.gt.1) then
         do ijk=1,nxyz
            tx(ijk,1,1)=sp(ijk,1,1)*ux(ijk,1,1)
            ty(ijk,1,1)=sp(ijk,1,1)*uy(ijk,1,1)
            tz(ijk,1,1)=sp(ijk,1,1)*uz(ijk,1,1)
         enddo
         call derx (cx,tx,tr,sx,ffxp,fsxp,fwxp,nx,ny,nz,1)
         call dery (cy,ty,tr,di,sy,ffyp,fsyp,fwyp,nx,ny,nz,1)
         call derz (cz,tz,tr,sz,ffzp,fszp,fwzp,nx,ny,nz,1)
         do ijk=1,nxyz
            sv(ijk,1,1)=cx(ijk,1,1)+cy(ijk,1,1)+cz(ijk,1,1)
         enddo
      else
         do j=1,ny
         do i=1,nx
            tx(i,j,1)=sp(i,j,1)*ux(i,j,1)
            ty(i,j,1)=sp(i,j,1)*uy(i,j,1)
         enddo
         enddo
         call derx (cx,tx,tr,sx,ffxp,fsxp,fwxp,nx,ny,nz,1)
         call dery (cy,ty,tr,di,sy,ffyp,fsyp,fwyp,nx,ny,nz,1)
         do j=1,ny
         do i=1,nx
            sv(i,j,1)=cx(i,j,1)+cy(i,j,1)
         enddo
         enddo
      endif
   endif
   if (nz.gt.1) then
      do k=1,nz 
      do j=1,ny
      do i=1,nx
         tx(i,j,k)=ux(i,j,k)*ux(i,j,k)
         ty(i,j,k)=uy(i,j,k)*uy(i,j,k)
         tz(i,j,k)=uz(i,j,k)*uz(i,j,k)
         xnutx(i,j,k)=uy(i,j,k)*uz(i,j,k)
         xnuty(i,j,k)=ux(i,j,k)*uz(i,j,k)
         xnutz(i,j,k)=ux(i,j,k)*uy(i,j,k)
      enddo
      enddo
      enddo
!
      call derx (cx,tx,tr,sx,ffxp,fsxp,fwxp,nx,ny,nz,1)
      call dery (cy,ty,tr,di,sy,ffyp,fsyp,fwyp,nx,ny,nz,1)
      call derz (cz,tz,tr,sz,ffzp,fszp,fwzp,nx,ny,nz,1)
      call derx (tx,xnutz,tr,sx,ffxp,fsxp,fwxp,nx,ny,nz,1)
      call dery (ty,xnutz,tr,di,sy,ffyp,fsyp,fwyp,nx,ny,nz,1)
      call derx (tz,xnuty,tr,sx,ffxp,fsxp,fwxp,nx,ny,nz,1)
      call derz (xnut,xnuty,tr,sz,ffzp,fszp,fwzp,nx,ny,nz,1)
      call dery (hx,xnutx,tr,di,sy,ffyp,fsyp,fwyp,nx,ny,nz,1)
      call derz (hy,xnutx,tr,sz,ffzp,fszp,fwzp,nx,ny,nz,1)
!
      do k=1,nz 
      do j=1,ny
      do i=1,nx
         cx(i,j,k)=cx(i,j,k)+ty(i,j,k)+xnut(i,j,k)
         cy(i,j,k)=cy(i,j,k)+tx(i,j,k)+hy(i,j,k)
         cz(i,j,k)=cz(i,j,k)+tz(i,j,k)+hx(i,j,k)
!            cx(i,j,k)=ty(i,j,k)+xnut(i,j,k)
!            cy(i,j,k)=tx(i,j,k)+sv(i,j,k)
!            cz(i,j,k)=tz(i,j,k)+sp(i,j,k)
      enddo
      enddo
      enddo
!
      call derx (tx,ux,tr,sx,ffxp,fsxp,fwxp,nx,ny,nz,1)
      call dery (ty,uy,tr,di,sy,ffyp,fsyp,fwyp,nx,ny,nz,1)
      call derz (tz,uz,tr,sz,ffzp,fszp,fwzp,nx,ny,nz,1)
!
      call derx (xnutz,uy,tr,sx,ffxp,fsxp,fwxp,nx,ny,nz,1)
      call dery (xnutx,uz,tr,di,sy,ffyp,fsyp,fwyp,nx,ny,nz,1)
      call derz (xnuty,ux,tr,sz,ffzp,fszp,fwzp,nx,ny,nz,1)
!
      call derx (xnut,uz,tr,sx,ffxp,fsxp,fwxp,nx,ny,nz,1)
      call dery (hx,ux,tr,di,sy,ffyp,fsyp,fwyp,nx,ny,nz,1)
      call derz (hy,uy,tr,sz,ffzp,fszp,fwzp,nx,ny,nz,1)
!
      do k=1,nz 
      do j=1,ny
      do i=1,nx
         cx(i,j,k)=0.5*(cx(i,j,k)+ux(i,j,k)*tx(i,j,k)+&
              uy(i,j,k)*hx(i,j,k)+&
              uz(i,j,k)*xnuty(i,j,k))
         cy(i,j,k)=0.5*(cy(i,j,k)+uy(i,j,k)*ty(i,j,k)+&
              ux(i,j,k)*xnutz(i,j,k)+&
              uz(i,j,k)*hy(i,j,k))
         cz(i,j,k)=0.5*(cz(i,j,k)+uz(i,j,k)*tz(i,j,k)+&
              uy(i,j,k)*xnutx(i,j,k)+&
              ux(i,j,k)*xnut(i,j,k))
      enddo
      enddo
      enddo
   else

      do k=1,nz
      do j=1,ny
      do i=1,nx
         tx(i,j,k)=ux(i,j,k)*ux(i,j,k)
         ty(i,j,k)=uy(i,j,k)*uy(i,j,k)
         xnutz(i,j,k)=ux(i,j,k)*uy(i,j,k)
      enddo
      enddo
      enddo
!
      call derx (cx,tx,tr,sx,ffxp,fsxp,fwxp,nx,ny,nz,1)
      call dery (cy,ty,tr,di,sy,ffyp,fsyp,fwyp,nx,ny,nz,1)
      call derx (tx,xnutz,tr,sx,ffxp,fsxp,fwxp,nx,ny,nz,1)
      call dery (ty,xnutz,tr,di,sy,ffy,fsy,fwy,nx,ny,nz,0)
!


      do k=1,nz
      do j=1,ny
      do i=1,nx
         cx(i,j,k)=cx(i,j,k)+ty(i,j,k)
         cy(i,j,k)=cy(i,j,k)+tx(i,j,k)
      enddo
      enddo
      enddo
!
      call derx (tx,ux,tr,sx,ffxp,fsxp,fwxp,nx,ny,nz,1)
      call dery (ty,uy,tr,di,sy,ffy,fsy,fwy,nx,ny,nz,0)
      call derx (xnutz,uy,tr,sx,ffxp,fsxp,fwxp,nx,ny,nz,1)
      call dery (hx,ux,tr,di,sy,ffyp,fsyp,fwyp,nx,ny,nz,1)
!
      do k=1,nz
      do j=1,ny,10
      do i=1,nx,4
         cx(i,j,k)=0.5*(cx(i,j,k)+ux(i,j,k)*tx(i,j,k)+&
              uy(i,j,k)*hx(i,j,k))
         cy(i,j,k)=0.5*(cy(i,j,k)+uy(i,j,k)*ty(i,j,k)+&
              ux(i,j,k)*xnutz(i,j,k))
!         print *,cx(i,j,1),cy(i,j,1)
      enddo
!      pause
      enddo
      enddo      
   endif
!
   return
end subroutine skew
!********************************************************************
!
subroutine corgpbis (ux,uy,px,py,nx,ny,nz,nlock)
! 
!********************************************************************
!
   USE paramx_m
!
   implicit none
!
   integer :: nx,ny,nz,i,j,k,nlock
   real(8),dimension(nx,ny,nz) :: ux,uy,px,py
! 
   if (nlock.eq.1) then
      do i=1,nx
      do j=1,ny
         uy(i,j,1)=-py(i,j,1)+uy(i,j,1) 
         ux(i,j,1)=-px(i,j,1)+ux(i,j,1)
      enddo
      enddo
   endif
!
   if (nlock.eq.2) then
      do i=1,nx
      do j=1,ny
         uy(i,j,1)=py(i,j,1)+uy(i,j,1) 
         ux(i,j,1)=px(i,j,1)+ux(i,j,1)
      enddo
      enddo
   endif
!
   return
end subroutine corgpbis
! 
!**********************************************************************
!**************************************************************************
!
      subroutine  funcao(outvar,x,y1,y2,y3,nf)
!
!**************************************************************************
!
      implicit none
!
      integer nf
!      
      real(8) outvar,x
      real(8) y1,y2,y3,y4

      if (nf.eq.1) then
        outvar = -y3*y1
      elseif (nf.eq.2) then
        outvar = y1
      else
        outvar = y2
      endif
      
      return

      end
      
