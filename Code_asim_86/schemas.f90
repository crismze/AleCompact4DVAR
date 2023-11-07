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
     ciwi6y,ppy,pp2y,ppyi,pp2yi,yp,ypi)
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
  real(8),dimension(ny) :: cisip6y,ciwip6y,cisi6y,ciwi6y,yp,ypi   
!
!**************************STRET***************************
  real(8),dimension(ny) :: ppy,pp2y
  real(8),dimension(nym) :: ppyi,pp2yi
!**********************************************************  
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
!
!*********** grille decalee, matriz Cx, annexe A tesis sylvain ***********
   alcaix6=9./62.  ! coef. ec. (3.32) tesis sylvain 
   acix6=(63./62.)/dx 
   bcix6=(17./62.)/3./dx 
!
   cfx6(1)=alcaix6 
   cfx6(2)=alcaix6 
   cfx6(nxm-2)=alcaix6 
   cfx6(nxm-1)=alcaix6 
   cfx6(nxm)=0. 
   if (nclx==0) ccx6(1)=2.  
   if (nclx==1) ccx6(1)=1. + alcaix6 
   if (nclx==2) ccx6(1)=1. + alcaix6 
   ccx6(2)=1. 
   ccx6(nxm-2)=1. 
   ccx6(nxm-1)=1. 
   if (nclx==0) ccx6(nxm)=1. + alcaix6*alcaix6  
   if (nclx==1) ccx6(nxm)=1. + alcaix6 
   if (nclx==2) ccx6(nxm)=1. + alcaix6
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
!
!*********** grille decalee, matriz Cx, page 55 tesis sylvain ***********
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
!
!*********** grille decalee, matriz CIx, annexe A tesis sylvain ***********
   ailcaix6=3./10.  ! coef. ec. (3.33) tesis sylvain 
   aicix6=3./4. 
   bicix6=1./20. 
!
   cifx6(1)=ailcaix6 
   cifx6(2)=ailcaix6 
   cifx6(nxm-2)=ailcaix6 
   cifx6(nxm-1)=ailcaix6 
   cifx6(nxm)=0. 
   if (nclx==0) cicx6(1)=2.  
   if (nclx==1) cicx6(1)=1. + ailcaix6 
   if (nclx==2) cicx6(1)=1. + ailcaix6 
   cicx6(2)=1. 
   cicx6(nxm-2)=1. 
   cicx6(nxm-1)=1. 
   if (nclx==0) cicx6(nxm)=1. + ailcaix6*ailcaix6
   if (nclx==1) cicx6(nxm)=1. + ailcaix6 
   if (nclx==2) cicx6(nxm)=1. + ailcaix6
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
!
!*********** grille decalee, matriz CIx, page 55 tesis sylvain ***********
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
!
!*********** grille decalee, matriz Cy, annexe A tesis sylvain ***********
   alcaiy6=9./62. 
   aciy6=(63./62.)/dy 
   bciy6=(17./62.)/3./dy 
!
   cfy6(1)=alcaiy6 
   cfy6(2)=alcaiy6 
   cfy6(nym-2)=alcaiy6 
   cfy6(nym-1)=alcaiy6 
   cfy6(nym)=0. 
   if (ncly==0) ccy6(1)=2.  
   if (ncly==1) ccy6(1)=1. + alcaiy6 
   if (ncly==2) ccy6(1)=1. + alcaiy6
   ccy6(2)=1. 
   ccy6(nym-2)=1. 
   ccy6(nym-1)=1. 
   if (ncly==0) ccy6(nym)=1. + alcaiy6*alcaiy6  
   if (ncly==1) ccy6(nym)=1. + alcaiy6 
   if (ncly==2) ccy6(nym)=1. + alcaiy6
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
!
!*********** grille decalee, matriz Cy, page 55 tesis sylvain ***********
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
!
!*********** grille decalee, matriz CIy, annexe A tesis sylvain ***********
   ailcaiy6=3./10. 
   aiciy6=3./4. 
   biciy6=1./20. 
!
   cify6(1)=ailcaiy6 
   cify6(2)=ailcaiy6 
   cify6(nym-2)=ailcaiy6 
   cify6(nym-1)=ailcaiy6 
   cify6(nym)=0. 
   if (ncly==0) cicy6(1)=2.  
   if (ncly==1) cicy6(1)=1. + ailcaiy6 
   if (ncly==2) cicy6(1)=1. + ailcaiy6 
   cicy6(2)=1. 
   cicy6(nym-2)=1. 
   cicy6(nym-1)=1. 
   if (ncly==0) cicy6(nym)=1. + ailcaiy6*ailcaiy6
   if (ncly==1) cicy6(nym)=1. + ailcaiy6 
   if (ncly==2) cicy6(nym)=1. + ailcaiy6
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
!
!*********** grille decalee, matriz CIy, page 55 tesis sylvain ***********
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
!
!
!*****************************************
!*****************************************
      do i=1,nx
         ffxp(i)=ffx(i)
         sfxp(i)=sfx(i)
      enddo
      do i=1,nxm   
         cfxp6(i)=cfx6(i)
         cifxp6(i)=cifx6(i)
      enddo
      do i=1,nx
         cifip6(i)=cifi6(i)
         cfip6(i)=cfi6(i)
      enddo
!*****************************************
      do j=1,ny
         ffyp(j)=ffy(j)
         sfyp(j)=sfy(j)
      enddo
      do j=1,nym
         cfyp6(j)=cfy6(j)
         cifyp6(j)=cify6(j)
      enddo
      do j=1,ny
         cifip6y(j)=cifi6y(j)
         cfip6y(j)=cfi6y(j)
      enddo
!*****************************************
      do k=1,nz
         ffzp(k)=ffz(k)
         sfzp(k)=sfz(k)
      enddo
!
      if (nclx.eq.1) then
         ffxp(1)=0.
         sfx (1)=0.
      endif
      cfxp6(1)=0.
      cfip6(1)=0.
!*****************************************
      if (ncly.eq.1) then
         ffyp(1)=0.
         sfy (1)=0.
      endif
!*****************************************
      cfyp6(1)=0.
      cfip6y(1)=0.
!*****************************************
      if (nclz.eq.1) then
         ffzp(1)=0.
         sfz (1)=0.
      endif
!*****************************************
!*****************************************
!
!
!  ***********************************************************
!  **************** correccion por stretching ****************
!  ***********************************************************

!*********** matrices Ay1, Ay1' *************
      do j=1,ny-1
         ffy(j)=ffy(j)*ppy(j+1)  ! ec. (3.54) tesis sylvain
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
!
!*********** grille decalee, matriz Cy1, annexe A tesis sylvain *************
!      do j=1,nym !version sylvain
!         cfy6(j)=cfy6(j)*ppyi(j+1)
!         cfyp6(j)=cfyp6(j)*ppyi(j+1)
!         ccy6(j)=ccy6(j)*ppyi(j)
!         cby6(j)=cby6(j)*ppyi(j)
!      enddo 
      do j=1,nym-1 !version philippe
         cfy6(j)=cfy6(j)*ppyi(j+1)  ! ec. (3.56) tesis sylvain
         cfyp6(j)=cfyp6(j)*ppyi(j+1)
      enddo 
      do j=1,nym !version philippe
         ccy6(j)=ccy6(j)*ppyi(j)
         cby6(j)=cby6(j)*ppyi(j)
      enddo 
!
!*********** grille decalee, matriz Cy1, page 55 tesis sylvain *************
      do j=1,nym
         cfi6y(j)=cfi6y(j)*ppy(j+1)
         cfip6y(j)=cfip6y(j)*ppy(j+1)
         cci6y(j)=cci6y(j)*ppy(j)
         cbi6y(j)=cbi6y(j)*ppy(j)
      enddo
!      cci6y(nym)=cci6y(ny)*ppy(ny)
!
!  ***********************************************************
!  ***********************************************************
!  ***********************************************************
!
!
!*************************************************************
!*************************************************************
      call prepare (fbx,fcx,ffx ,fsx ,fwx ,nx)
      call prepare (fbx,fcx,ffxp,fsxp,fwxp,nx)
      call prepare (cbx6,ccx6,cfx6 ,csx6 ,cwx6 ,nxm)
      call prepare (cbx6,ccx6,cfxp6,csxp6,cwxp6,nxm)
      call prepare (cibx6,cicx6,cifx6 ,cisx6 ,ciwx6 ,nxm)
      call prepare (cibx6,cicx6,cifxp6,cisxp6,ciwxp6,nxm)
      call prepare (cbi6,cci6,cfi6 ,csi6 ,cwi6 ,nx)
      call prepare (cbi6,cci6,cfip6,csip6,cwip6,nx)
      call prepare (cibi6,cici6,cifi6 ,cisi6 ,ciwi6 ,nx)
      call prepare (cibi6,cici6,cifip6,cisip6,ciwip6,nx)
!*************************************************************
      call prepare (fby,fcy,ffy ,fsy ,fwy ,ny)
      call prepare (fby,fcy,ffyp,fsyp,fwyp,ny)
      call prepare (cby6,ccy6,cfy6 ,csy6 ,cwy6 ,nym)
      call prepare (cby6,ccy6,cfyp6,csyp6,cwyp6,nym)
      call prepare (ciby6,cicy6,cify6 ,cisy6 ,ciwy6 ,nym)
      call prepare (ciby6,cicy6,cifyp6,cisyp6,ciwyp6,nym)
      call prepare (cbi6y,cci6y,cfi6y ,csi6y ,cwi6y ,ny)
      call prepare (cbi6y,cci6y,cfip6y,csip6y,cwip6y,ny)
      call prepare (cibi6y,cici6y,cifi6y ,cisi6y ,ciwi6y ,ny)
      call prepare (cibi6y,cici6y,cifip6y,cisip6y,ciwip6y,ny)
!*************************************************************
      if (nz.gt.1) call prepare (fbz,fcz,ffz,fsz,fwz,nz)
!*************************************************************
      if (nclx.eq.1) then
         fbx(nx-1)=0.
         cbx6(nxm-1)=0.
         cibx6(nxm)=0
         cbi6(nx-1)=0.
         cibi6(nx)=0
         call prepare (fbx,fcx,ffxp,fsxp,fwxp,nx)
         call prepare (cbx6,ccx6,cfxp6,csxp6,cwxp6,nxm)
         call prepare (cibx6,cicx6,cifxp6,cisxp6,ciwxp6,nxm)
         call prepare (cbi6,cci6,cfip6,csip6,cwip6,nx)
         call prepare (cibi6,cici6,cifip6,cisip6,ciwip6,nx)
      endif
!
      if (nclx.eq.2) then
         cbx6(nxm-1)=0.
         cibx6(nxm)=0
         cbi6(nx-1)=0.
         cibi6(nx)=0
         call prepare (cbx6,ccx6,cfxp6,csxp6,cwxp6,nxm)
         call prepare (cibx6,cicx6,cifxp6,cisxp6,ciwxp6,nxm)
         call prepare (cbi6,cci6,cfip6,csip6,cwip6,nx)
         call prepare (cibi6,cici6,cifip6,cisip6,ciwip6,nx)
      endif
!
      if (ncly.eq.1) then
         cby6(nym-1)=0.
         ciby6(nym)=0
         cbi6y(ny-1)=0.
         cibi6y(ny)=0 
         fby(ny-1)=0.
         call prepare (fby,fcy,ffyp,fsyp,fwyp,ny)
         call prepare (cby6,ccy6,cfyp6,csyp6,cwyp6,nym)
         call prepare (ciby6,cicy6,cifyp6,cisyp6,ciwyp6,nym)
         call prepare (cbi6y,cci6y,cfip6y,csip6y,cwip6y,ny)
         call prepare (cibi6y,cici6y,cifip6y,cisip6y,ciwip6y,ny)
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
!
!********************************************************************
!
subroutine schemas_r(ffx,fcx,fbx,ffy,fcy,fby,&
        fsx,fwx,fsy,fwy,&
        ffxp,fsxp,fwxp,ffyp,fsyp,fwyp,&
        nx,ny,ppy)
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
  real(8),dimension(nx) :: ffx,fcx,fbx,fsx,fwx 
  real(8),dimension(nx) :: ffxp,fsxp,fwxp
  real(8),dimension(ny) :: ffy,fcy,fby,fsy,fwy 
  real(8),dimension(ny) :: ffyp,fsyp,fwyp
!
!**************************STRET***************************
  real(8),dimension(ny) :: ppy
!**********************************************************  
  integer  :: nx,ny
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer :: i, j 
!-----------------------------------------------
! 
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
!
!*****************************************
!*****************************************
      do i=1,nx
         ffxp(i)=ffx(i)
      enddo
!*****************************************
      do j=1,ny
         ffyp(j)=ffy(j)
      enddo
!*****************************************
      if (ncly.eq.1) then
         ffyp(1)=0.
      endif
!*****************************************
!*****************************************
!
!
!  ***********************************************************
!  **************** correccion por stretching ****************
!  ***********************************************************

!*********** matrices Ay1, Ay1' *************
      do j=1,ny-1
         ffy(j)=ffy(j)*ppy(j+1)  ! ec. (3.54) tesis sylvain
         ffyp(j)=ffyp(j)*ppy(j+1)
         fcy(j)=fcy(j)*ppy(j)
         fby(j)=fby(j)*ppy(j)
      enddo
      fcy(ny)=fcy(ny)*ppy(ny)
!
!  ***********************************************************
!  ***********************************************************
!  ***********************************************************
!
!
!*************************************************************
!*************************************************************
      call prepare (fbx,fcx,ffx ,fsx ,fwx ,nx)
      call prepare (fbx,fcx,ffxp,fsxp,fwxp,nx)
!*************************************************************
      call prepare (fby,fcy,ffy ,fsy ,fwy ,ny)
      call prepare (fby,fcy,ffyp,fsyp,fwyp,ny)
!*************************************************************
      if (ncly.eq.1) then
         fby(ny-1)=0.
         call prepare (fby,fcy,ffyp,fsyp,fwyp,ny)
      endif
!
      return
      end
!
!********************************************************************
!
subroutine schemas_rg(ffx,fcx,fbx,ffy,fcy,fby,&
        fsx,fwx,fsy,fwy,&
        ffxp,fsxp,fwxp,ffyp,fsyp,fwyp,&
        nx,ny,ppy)
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
  real(8),dimension(nx) :: ffx,fcx,fbx,fsx,fwx 
  real(8),dimension(nx) :: ffxp,fsxp,fwxp
  real(8),dimension(ny) :: ffy,fcy,fby,fsy,fwy 
  real(8),dimension(ny) :: ffyp,fsyp,fwyp
!
!**************************STRET***************************
  real(8),dimension(ny) :: ppy
!**********************************************************  
  integer  :: nx,ny
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer :: i, j 
!-----------------------------------------------
! 
!
      alfa1x_rg= 2.
      af1x_rg  =-(5./2.  )/dxg
      bf1x_rg  = (   2.  )/dxg
      cf1x_rg  = (1./2.  )/dxg
      df1x_rg  = 0.
      alfa2x_rg= 1./4.
      af2x_rg  = (3./4.  )/dxg
!
      alfanx_rg= 2.
      afnx_rg  =-(5./2.  )/dxg
      bfnx_rg  = (   2.  )/dxg
      cfnx_rg  = (1./2.  )/dxg
      dfnx_rg  = 0.
      alfamx_rg= 1./4.
      afmx_rg  = (3./4.  )/dxg
!
      alfaix_rg= 1./3.
      afix_rg  = (7./9.  )/dxg
      bfix_rg  = (1./36. )/dxg
!
      alfa1y_rg= 2.
      af1y_rg  =-(5./2.  )/dyg
      bf1y_rg  = (   2.  )/dyg
      cf1y_rg  = (1./2.  )/dyg
      df1y_rg  = 0.
      alfa2y_rg= 1./4.
      af2y_rg  = (3./4.  )/dyg
!
      alfany_rg= 2.
      afny_rg  =-(5./2.  )/dyg
      bfny_rg  = (   2.  )/dyg
      cfny_rg  = (1./2.  )/dyg
      dfny_rg  = 0.
      alfamy_rg= 1./4.
      afmy_rg  = (3./4.  )/dyg
!
      alfajy_rg= 1./3.
      afjy_rg  = (7./9.  )/dyg
      bfjy_rg  = (1./36. )/dyg
!
      if (nclx.eq.2) then
         ffx(1)   =alfa1x_rg
         ffx(2)   =alfa2x_rg
         ffx(nx-2)=alfaix_rg
         ffx(nx-1)=alfamx_rg
         ffx(nx)  =0.
         fcx(1)   =1.
         fcx(2)   =1.
         fcx(nx-2)=1.
         fcx(nx-1)=1.
         fcx(nx  )=1.
         fbx(1)   =alfa2x_rg 
         fbx(2)   =alfaix_rg
         fbx(nx-2)=alfamx_rg
         fbx(nx-1)=alfanx_rg
         fbx(nx  )=0.
         do i=3,nx-3
            ffx(i)=alfaix_rg
            fcx(i)=1.
            fbx(i)=alfaix_rg
         enddo
      endif
!
      if (ncly.eq.1) then
         ffy(1)   =alfajy_rg+alfajy_rg
         ffy(2)   =alfajy_rg
         ffy(ny-2)=alfajy_rg
         ffy(ny-1)=alfajy_rg
         ffy(ny)  =0.
         fcy(1)   =1.
         fcy(2)   =1.
         fcy(ny-2)=1.
         fcy(ny-1)=1.
         fcy(ny  )=1.
         fby(1)   =alfajy_rg 
         fby(2)   =alfajy_rg
         fby(ny-2)=alfajy_rg
         fby(ny-1)=alfajy_rg+alfajy_rg
         fby(ny  )=0.
         do j=3,ny-3
            ffy(j)=alfajy_rg
            fcy(j)=1.
            fby(j)=alfajy_rg
         enddo
      endif
!
!
!*****************************************
!*****************************************
      do i=1,nx
         ffxp(i)=ffx(i)
      enddo
!*****************************************
      do j=1,ny
         ffyp(j)=ffy(j)
      enddo
!*****************************************
      if (ncly.eq.1) then
         ffyp(1)=0.
      endif
!*****************************************
!*****************************************
!
!
!  ***********************************************************
!  **************** correccion por stretching ****************
!  ***********************************************************

!*********** matrices Ay1, Ay1' *************
      do j=1,ny-1
         ffy(j)=ffy(j)*ppy(j+1)  ! ec. (3.54) tesis sylvain
         ffyp(j)=ffyp(j)*ppy(j+1)
         fcy(j)=fcy(j)*ppy(j)
         fby(j)=fby(j)*ppy(j)
      enddo
      fcy(ny)=fcy(ny)*ppy(ny)
!
!  ***********************************************************
!  ***********************************************************
!  ***********************************************************
!
!
!*************************************************************
!*************************************************************
      call prepare (fbx,fcx,ffx ,fsx ,fwx ,nx)
      call prepare (fbx,fcx,ffxp,fsxp,fwxp,nx)
!*************************************************************
      call prepare (fby,fcy,ffy ,fsy ,fwy ,ny)
      call prepare (fby,fcy,ffyp,fsyp,fwyp,ny)
!*************************************************************
      if (ncly.eq.1) then
         fby(ny-1)=0.
         call prepare (fby,fcy,ffyp,fsyp,fwyp,ny)
      endif
!
      return
      end
!
