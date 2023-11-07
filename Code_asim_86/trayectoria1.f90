!********************************************************************
! Guarda secuencia de campos de las trayectorias para armar video
PROGRAM trayectoria1
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
   integer,parameter :: nxrg=67,nyrg=47 ! malla reducida grosera (obs.)
   integer,parameter :: nxr=199,nyr=139 ! malla reducida fina (modelo)
   integer,parameter :: rdxy=3 ! relacion entre la resolucion fina (1/dx) y la resolucion grosera (1/dxg)
   integer,parameter :: jyc=50,jyt=100 ! delimita zonas de veloc. constante y de transicion en malla reconstruida
   integer,parameter :: nx=nxr,ny=nyr+2*jyt,nz=1 ! malla reconstruida (modelo)
   integer,parameter :: i_redg1=5,i_redg2=nxrg-4 ! bordes izquierdo y derecho del campo asimilado en malla reducida grosera
   integer,parameter :: j_redg1=5,j_redg2=nyrg-4 ! bordes inferior y superior del campo asimilado en malla reducida grosera
   integer,parameter :: nxrg_asim=i_redg2-i_redg1+1,nyrg_asim=j_redg2-j_redg1+1 ! dominio espacial de asimilacion, malla grosera
   integer,parameter :: i_red1=rdxy*i_redg1-rdxy+1,i_red2=rdxy*i_redg2-rdxy+1 ! bordes izquierdo y derecho del campo asimilado en malla reducida fina
   integer,parameter :: j_red1=rdxy*j_redg1-rdxy+1,j_red2=rdxy*j_redg2-rdxy+1 ! bordes inferior y superior del campo asimilado en malla reducida fina
   integer,parameter :: nxr_asim=i_red2-i_red1+1,nyr_asim=j_red2-j_red1+1 ! dominio espacial de asimilacion, malla fina
   integer,parameter :: i_rec1=i_red1,i_rec2=i_red2 ! bordes izquierdo y derecho del campo asimilado en malla reconstruida
   integer,parameter :: j_rec1=j_red1+jyt-1,j_rec2=j_red2+jyt-1 ! bordes inferior y superior del campo asimilado en malla reconstruida
   integer,parameter :: nxm=nx-1,nym=ny-1
   integer,parameter :: mx=nx+1,my=ny+1,mz=nz+2
   integer,parameter :: myr=nyr+1
   integer,parameter :: myrg=nyrg+1
   integer,parameter :: m1=2,m2=2
   integer,parameter :: fenetre=10 ! taille de la sequence des obs. para asimilacion
   integer,parameter :: fenetre_obs=10 ! taille de la sequence des obs. para armar video
   integer,parameter :: fenetrem=fenetre-1
   integer,parameter :: imodulob=128 ! intervalle temporel entre obs.
   integer,parameter :: imodulo_video=2 ! intervalle temporel entre campos para armar video
   integer,parameter :: ditime=fenetrem*imodulob
   integer,parameter :: ditime_asim=ditime-imodulob ! dominio temporal de asimilacion
   integer,parameter :: imoduloe=16 ! intervalle spatial entre obs.
   integer,parameter :: ditimev=fenetrem*imoduloe+1
   integer,parameter :: r_g=rdxy ! radio de la zona considerada en la reduccion de 1/dx a 1/dxg
   integer,parameter :: n_nodos_max=(2*r_g+1)**2 ! numero maximo de nodos contenidos en la zona de reduccion espacial
   integer,parameter :: rq_g=rdxy ! radio de la zona considerada en el suavizado de la cond. inicial
   integer,parameter :: n_nodos_q_max=(2*rq_g+1)**2 ! numero maximo de nodos contenidos en la zona de suavizado de la cond. inicial
   integer,parameter :: d_g=imodulob ! distancia temporal considerada en el suavizado de la cond. de entrada
   integer,parameter :: n_nodos_t_max=2*d_g+1 ! numero maximo de nodos contenidos en el intervalo de suavizado de la cond. de entrada
!  ********************* LBFGS *********************
   integer,parameter :: nxyr=nxr*nyr 
   integer,parameter :: nytr=nyr*ditime 
   integer,parameter :: ndim=4*nxyr+2*nytr ! number of variables
   integer,parameter :: msave=15 ! number of corrections used in the BFGS update
   integer,parameter :: nwork=ndim*(2*msave+1)+2*msave ! length of array W used as workspace
!  *************************************************
!
   real(8),dimension(fenetrem,nx,ny,nz) :: ux,uy,gx,gy,wz
   real(8),dimension(fenetrem,nx,ny,nz) :: uxb,uyb,gxb,gyb,wzb
   real(8),dimension(imodulob,nx,ny,nz) :: ux_rec,uy_rec,gx_rec,gy_rec
   real(8),dimension(fenetre_obs,nxrg,nyrg,nz) :: ux_obsg,uy_obsg,wz_obsg
   real(8),dimension(nxr,nyr,nz) :: ux0r_obs,uy0r_obs
   real(8),dimension(nxr,nyr,nz) :: ux0r,uy0r,gx0r,gy0r,wz0r
   real(8),dimension(nxr,nyr,nz) :: ux0rb,uy0rb,gx0rb,gy0rb,wz0rb
   real(8),dimension(nx,ny,nz) :: ux0,uy0,gx0,gy0
   real(8),dimension(nx,ny,nz) :: ux0b,uy0b,gx0b,gy0b
   real(8),dimension(nx,ny,nz) :: uxtemp,uytemp,gxtemp,gytemp
   real(8),dimension(nx,ny,nz) :: uxtempE,uytempE,gxtempE,gytempE
   real(8),dimension(nx,ny,nz) :: uxtempS,uytempS,gxtempS,gytempS
   real(8),dimension(nx,ny,nz) :: uxitemp,uyitemp,gxitemp,gyitemp,wzitemp
   real(8),dimension(nx,ny,nz) :: uxftemp,uyftemp,gxftemp,gyftemp,wzftemp
   real(8),dimension(ny,ditime) :: bx1_ini,by1_ini,bx1_fin,by1_fin
   real(8),dimension(nyr,ditime) :: bx1_obs,by1_obs
   real(8),dimension(nyr,ditime) :: bx1r,by1r,bx1r_av,by1r_av
   real(8),dimension(nyr,ditime) :: bx1rb,by1rb
   real(8),dimension(ny,ditime) :: bx1,by1
   real(8),dimension(ny,ditime) :: bx1b,by1b
   real(8),dimension(nx,ditime) :: bx2_fin,by2_fin
   real(8),dimension(ny,ditime) :: bx3_ini,by3_ini,bx3_fin,by3_fin
   real(8),dimension(nx) :: ffx,fsx,fwx,ffxp,fsxp,fwxp
   real(8),dimension(ny) :: ffy,fsy,fwy,ffyp,fsyp,fwyp
   real(8),dimension(nxr) :: ffx_r,fcx_r,fbx_r,fsx_r,fwx_r 
   real(8),dimension(nxr) :: ffxp_r,fsxp_r,fwxp_r
   real(8),dimension(nyr) :: ffy_r,fcy_r,fby_r,fsy_r,fwy_r 
   real(8),dimension(nyr) :: ffyp_r,fsyp_r,fwyp_r
   real(8),dimension(nxrg) :: ffx_rg,fcx_rg,fbx_rg,fsx_rg,fwx_rg 
   real(8),dimension(nxrg) :: ffxp_rg,fsxp_rg,fwxp_rg
   real(8),dimension(nyrg) :: ffy_rg,fcy_rg,fby_rg,fsy_rg,fwy_rg 
   real(8),dimension(nyrg) :: ffyp_rg,fsyp_rg,fwyp_rg
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
!  ********************* Filtrage *********************
   real(8),dimension(nx) :: fifx,ficx,fibx,fiffx,fibbx,fiz1x,fiz2x
   real(8),dimension(nx,2) :: filax,filaxp
   real(8),dimension(nx) :: fifxp,ficxp,fibxp,fiffxp,fibbxp
   real(8),dimension(ny) :: fify,ficy,fiby,fiffy,fibby,fiz1y,fiz2y
   real(8),dimension(ny,2) :: filay,filayp
   real(8),dimension(ny) :: fifyp,ficyp,fibyp,fiffyp,fibbyp
!  ********************* STRET *********************
   real(8),dimension(ny) :: pp4y
   real(8),dimension(nx) :: d5
   real(8),dimension(ny/2,m1+m2+1,nxm) :: a_tot,a2_tot
   real(8),dimension(ny/2,m1,nxm) :: al_tot,al1_tot
   real(8),dimension(ny/2,nxm) :: indx_tot,indx1_tot
!  *************************************************
   real(8) :: adt,bdt,gdt,zero
   integer :: itime,i,j,f,n,ired,jred,irec,jrec,j0
   integer :: i_redg,j_redg,i_rec,j_rec,nodo
!  ********************* Obs. *********************
   real(8),dimension(nx,ny) :: rx,ry
   integer,dimension(ditimev) :: itimev
   real(8),dimension(nyr,ditime+1) :: bx1p,by1p
   integer,dimension(nxrg,nyrg,n_nodos_max) :: iw,jw
   real(8),dimension(nxrg,nyrg,n_nodos_max) :: w_xy
   integer,dimension(nxrg,nyrg) :: n_nodos
   real(8),dimension(nxrg,nyrg) :: sumw_xy
   real(8) :: sigma_r,r,w_r,sumw_r
   integer,dimension(ditime,n_nodos_t_max) :: itw
   real(8),dimension(ditime,n_nodos_t_max) :: w_t
   integer,dimension(ditime) :: n_nodos_t
   real(8),dimension(ditime) :: sumw_t
   real(8),dimension(n_nodos_t_max) :: wbx1r,wby1r
   real(8) :: sigma_t,sumw_d,d,w_d
   real(8) :: ruidox,ruidoy,pendx,pendy
   integer :: longueur,iob,iobi,iobf,iobfilt,rimod,t0,tf,iv
   character(len=4) suffix
   character(len=80) nfichier,ncampo,nvector
   real(8) :: snorm4,snorm5,snorm6,error,norma_wz
!  ********************* LBFGS *********************
   double precision xp(ndim),gr(ndim),diag(ndim),w(nwork)
   double precision fc,eps,xtol
   integer iprint(2),iflag,icall,iter
   logical diagco
!  *************************************************
!
!
!  Chargement des parametres
   call datosinc3d (ffx,fsx,fwx,ffy,fsy,fwy,ffxp,fsxp,fwxp,ffyp,fsyp,fwyp,&
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
!  ******************* Parametros trayectoria Inicial *******************
!
!  Lee estimacion inicial del vector solucion 
   nvector='xp_iter0000'
   open(90,file=nvector,form='unformatted')
   read(90) xp
   close(90)
!
!  Redimensiona xp
   ux0r=reshape(xp(1:nxyr),(/nxr,nyr,nz/))
   uy0r=reshape(xp(nxyr+1:2*nxyr),(/nxr,nyr,nz/))
   gx0r=reshape(xp(2*nxyr+1:3*nxyr),(/nxr,nyr,nz/))
   gy0r=reshape(xp(3*nxyr+1:4*nxyr),(/nxr,nyr,nz/))
   bx1r=reshape(xp(4*nxyr+1:4*nxyr+nytr),(/nyr,ditime/))
   by1r=reshape(xp(4*nxyr+nytr+1:4*nxyr+2*nytr),(/nyr,ditime/))
!
!  Entrada: ux0r, uy0r, gx0r, gy0r, bx1r, by1r
!  Salida: ux0, uy0, gx0, gy0, bx1, by1
   call reconstruye4 (ux0r,uy0r,gx0r,gy0r,bx1r,by1r,&
          ux0,uy0,gx0,gy0,bx1,by1,&
          jyc,jyt,nxr,nyr,nx,ny,nz,ditime)
!
!  Cond. inicial para la trayectoria Inicial
   uxitemp=ux0
   uyitemp=uy0
   gxitemp=gx0 
   gyitemp=gy0
!
!  Cond. de entrada para la trayectoria Inicial
   bx1_ini=bx1
   by1_ini=by1
!
!  Entrada: uxitemp, uyitemp. Salida: wzitemp
   call vortz(uxitemp,uyitemp,wzitemp,ffxp,fsxp,fwxp,ffyp,fsyp,fwyp,nx,ny,nz)
!
!  Guarda cond. inicial para la trayectoria Inicial en dominio de asimilacion
   ncampo='wzini'
   iobi=iobs
   call numcar (iobi,suffix)
   longueur=index(ncampo,' ')-1
   nfichier=ncampo(1:longueur)//suffix
   longueur=index(nfichier,' ')-1
   open(100,file=nfichier(1:longueur),form='formatted')
   do j=j_rec1,j_rec2
   do i=i_rec1,i_rec2
      write(100,*) wzitemp(i,j,1)
   enddo
   enddo 
   close(100)
!
!  ******************* Parametros trayectoria Asimilada *******************
!
!  Lee estimacion del vector solucion en la evaluacion icall de fc y gr (best point found by LBFGS)
   nvector='xp_iter'
   iter=420
   call numcar (iter,suffix)
   longueur=index(nvector,' ')-1
   nfichier=nvector(1:longueur)//suffix  
   longueur=index(nfichier,' ')-1
   open(90,file=nfichier(1:longueur),form='unformatted')
   read(90) xp
   close(90)
!
!  Redimensiona xp
   ux0r=reshape(xp(1:nxyr),(/nxr,nyr,nz/))
   uy0r=reshape(xp(nxyr+1:2*nxyr),(/nxr,nyr,nz/))
   gx0r=reshape(xp(2*nxyr+1:3*nxyr),(/nxr,nyr,nz/))
   gy0r=reshape(xp(3*nxyr+1:4*nxyr),(/nxr,nyr,nz/))
   bx1r=reshape(xp(4*nxyr+1:4*nxyr+nytr),(/nyr,ditime/))
   by1r=reshape(xp(4*nxyr+nytr+1:4*nxyr+2*nytr),(/nyr,ditime/))
!
!  Entrada: ux0r, uy0r, gx0r, gy0r, bx1r, by1r
!  Salida: ux0, uy0, gx0, gy0, bx1, by1
   call reconstruye4 (ux0r,uy0r,gx0r,gy0r,bx1r,by1r,&
          ux0,uy0,gx0,gy0,bx1,by1,&
          jyc,jyt,nxr,nyr,nx,ny,nz,ditime)
!
!  Cond. inicial para la trayectoria Asimilada
   uxftemp=ux0
   uyftemp=uy0
   gxftemp=gx0 
   gyftemp=gy0
!
!  Cond. de entrada para la trayectoria Asimilada
   bx1_fin=bx1
   by1_fin=by1
!
!  Entrada: uxftemp, uyftemp. Salida: wzftemp
   call vortz(uxftemp,uyftemp,wzftemp,ffxp,fsxp,fwxp,ffyp,fsyp,fwyp,nx,ny,nz)
!
!  Guarda cond. inicial para la trayectoria Asimilada en dominio de asimilacion
   ncampo='wzfin'
   iobf=iobs
   call numcar (iobf,suffix)
   longueur=index(ncampo,' ')-1
   nfichier=ncampo(1:longueur)//suffix
   longueur=index(nfichier,' ')-1
   open(100,file=nfichier(1:longueur),form='formatted')
   do j=j_rec1,j_rec2
   do i=i_rec1,i_rec2
      write(100,*) wzftemp(i,j,1)
   enddo
   enddo 
   close(100)
!
!  ******************* Parametros trayectoria de Obs. *******************
!
!  Lee campos PIV que han sido reducidos y normalizados en
!  /home/ale/Programas/fiuba-2012/reduce_normaliza_suaviza_obs_asim82b_video
   ncampo='cmcobsg'
   iob=iobs
   do f=1,fenetre_obs
      call numcar (iob,suffix)
      longueur=index(ncampo,' ')-1
      nfichier=ncampo(1:longueur)//suffix
      longueur=index(nfichier,' ')-1
      open(100,file=nfichier(1:longueur),form='formatted')
      do j=1,nyrg
      do i=1,nxrg
         read(100,*) ux_obsg(f,i,j,1),uy_obsg(f,i,j,1)
      enddo
      enddo 
      close(100)
      iob=iob+1
   enddo
!
!  Entrada: ux_obsg, uy_obsg. Salida: wz_obsg
   do f=1,fenetre_obs
      call vortz_rg(ux_obsg(f,:,:,:),uy_obsg(f,:,:,:),wz_obsg(f,:,:,:),&
             ffxp_rg,fsxp_rg,fwxp_rg,ffyp_rg,fsyp_rg,fwyp_rg,nxrg,nyrg,nz)
   enddo
!
!  Guarda trayectoria de Obs. en dominio de asimilacion
   ncampo='wzobs'
   iob=iobs
   do f=1,fenetre_obs
      call numcar (iob,suffix)
      longueur=index(ncampo,' ')-1
      nfichier=ncampo(1:longueur)//suffix
      longueur=index(nfichier,' ')-1
      open(100,file=nfichier(1:longueur),form='formatted')
      do j=j_redg1,j_redg2
      do i=i_redg1,i_redg2
         write(100,*) wz_obsg(f,i,j,1)
      enddo
      enddo 
      close(100)
      iob=iob+1
   enddo
!
!  ******************* Calculo de trayectorias *******************
!
   do f=1,fenetrem
!
      do n=1,imodulob
!
      itime=n+(f-1)*imodulob
!
      call inc3dorigin (uxitemp,uyitemp,gxitemp,gyitemp,&
             uxitemp,uyitemp,gxitemp,gyitemp,&
             bx1_ini,by1_ini,itime,&
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
      call inc3dorigin (uxftemp,uyftemp,gxftemp,gyftemp,&
             uxftemp,uyftemp,gxftemp,gyftemp,&
             bx1_fin,by1_fin,itime,&
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
      if (mod(itime,imodulo_video).eq.0) then 
!
!     Entrada: uxitemp, uyitemp. Salida: wzitemp
      call vortz(uxitemp,uyitemp,wzitemp,ffxp,fsxp,fwxp,ffyp,fsyp,fwyp,nx,ny,nz)
!
!     Entrada: uxftemp, uyftemp. Salida: wzftemp
      call vortz(uxftemp,uyftemp,wzftemp,ffxp,fsxp,fwxp,ffyp,fsyp,fwyp,nx,ny,nz)
!
!     Guarda campo trayectoria Inicial en dominio de asimilacion
      ncampo='wzini'
      iobi=iobi+1
      call numcar (iobi,suffix)
      longueur=index(ncampo,' ')-1
      nfichier=ncampo(1:longueur)//suffix
      longueur=index(nfichier,' ')-1
      open(100,file=nfichier(1:longueur),form='formatted')
      do j=j_rec1,j_rec2
      do i=i_rec1,i_rec2
         write(100,*) wzitemp(i,j,1)
      enddo
      enddo 
      close(100)
!
!     Guarda campo trayectoria Asimilada en dominio de asimilacion
      ncampo='wzfin'
      iobf=iobf+1
      call numcar (iobf,suffix)
      longueur=index(ncampo,' ')-1
      nfichier=ncampo(1:longueur)//suffix
      longueur=index(nfichier,' ')-1
      open(100,file=nfichier(1:longueur),form='formatted')
      do j=j_rec1,j_rec2
      do i=i_rec1,i_rec2
         write(100,*) wzftemp(i,j,1)
      enddo
      enddo 
      close(100)
!
      endif
!
      enddo
!
   enddo
!
!
end PROGRAM trayectoria1
!
!********************************************************************


