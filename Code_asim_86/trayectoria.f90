!********************************************************************
!
PROGRAM trayectoria
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
   integer,parameter :: fenetre=30 ! taille de la sequence des obs. para asimilacion
   integer,parameter :: fenetrem=fenetre-1
   integer,parameter :: imodulob=128 ! intervalle temporel entre obs.
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
   real(8),dimension(fenetre,nxr,nyr,nz) :: ux_obs,uy_obs
   real(8),dimension(fenetre,nxrg,nyrg,nz) :: ux_obsg,uy_obsg,wz_obsg
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
   real(8),dimension(nx,ny,nz) :: uxmt_fin,uymt_fin,wzmt_fin,uxux_fin,uyuy_fin,uxuy_fin,wzwz_fin
   real(8),dimension(nxrg,nyrg,nz) :: uxmt_obs,uymt_obs,wzmt_obs,uxux_obs,uyuy_obs,uxuy_obs,wzwz_obs
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
   real(8) :: snorm4,snorm5,snorm6,error,norma_obs,norma1_ini,norma1_fin
   real(8),dimension(ditime) :: normav_ini,normav_fin
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
!  ******** Gaussian weights, coordenada espacial, malla grosera ********
   sigma_r=2.
   jw=0. ; iw=0. ; w_xy=0.
   do j_redg = 1,nyrg
   do i_redg = 1,nxrg
      j_rec = jyt+rdxy*j_redg-rdxy
      i_rec = rdxy*i_redg-rdxy+1
!
      sumw_r=0.
      nodo=0
      do j=jyt,jyt+nyr-1
      do i=1,nx
         r=sqrt(real((i-i_rec)**2+(j-j_rec)**2))
         if (r<=r_g) then
            w_r=exp(-(r**2)/(2*sigma_r**2))
            sumw_r = sumw_r + w_r
            nodo=nodo+1
            jw(i_redg,j_redg,nodo)=j
            iw(i_redg,j_redg,nodo)=i
            w_xy(i_redg,j_redg,nodo)=w_r
         endif
      enddo
      enddo
      sumw_xy(i_redg,j_redg)=sumw_r
      n_nodos(i_redg,j_redg)=nodo
   enddo
   enddo 
!  **********************************************************************
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
!  Determina norma vort. para la trayectoria Inicial en dominio de asimilacion  
   norma1_ini=snorm5(wzitemp,nx,ny,nz,nxrg,nyrg,&
             jw,iw,w_xy,sumw_xy,n_nodos,n_nodos_max,&
             i_redg1,i_redg2,j_redg1,j_redg2,nxrg_asim,nyrg_asim)
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
!  Guarda cond. de entrada para la trayectoria Inicial en dominio reducido
   open(80,file='entrada_ini.dat',form='formatted')
   do itime=1,ditime
   do j=1,nyr
      write(80,*) bx1r(j,itime),by1r(j,itime)
   enddo
   enddo  
   close(80)
!
!  ******************* Parametros trayectoria Asimilada *******************
!
!  Lee estimacion del vector solucion en la evaluacion icall de fc y gr (best point found by LBFGS)
   nvector='xp_iter'
   iter=540
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
!  Determina norma vort. para la trayectoria Asimilada en dominio de asimilacion  
   norma1_fin=snorm5(wzftemp,nx,ny,nz,nxrg,nyrg,&
             jw,iw,w_xy,sumw_xy,n_nodos,n_nodos_max,&
             i_redg1,i_redg2,j_redg1,j_redg2,nxrg_asim,nyrg_asim)
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
!  Guarda cond. de entrada para la trayectoria Asimilada en dominio reducido
   open(80,file='entrada_fin.dat',form='formatted')
   do itime=1,ditime
   do j=1,nyr
      write(80,*) bx1r(j,itime),by1r(j,itime)
   enddo
   enddo  
   close(80)
!
!  ******************* Parametros trayectoria de Obs. *******************
!
!  Lee campos PIV que han sido reducidos y normalizados en
!  /home/ale/Programas/fiuba-2014/reduce_normaliza_suaviza_interpola_obs_asim86.m
   ncampo='cmcobsg'
   iob=iobs
   do f=1,fenetre
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
   do f=1,fenetre
      call vortz_rg(ux_obsg(f,:,:,:),uy_obsg(f,:,:,:),wz_obsg(f,:,:,:),&
             ffxp_rg,fsxp_rg,fwxp_rg,ffyp_rg,fsyp_rg,fwyp_rg,nxrg,nyrg,nz)
   enddo
!
!  Determina y guarda norma vort. para la trayectoria de Obs. en dominio de asimilacion  
   open(40,file='Norma_obs',form='formatted')
   do f=1,fenetre
      norma_obs=snorm6(wz_obsg(f,:,:,:),nxrg,nyrg,nz,&
                i_redg1,i_redg2,j_redg1,j_redg2,nxrg_asim,nyrg_asim)
      itime=(f-1)*imodulob
      write(40,*) itime,norma_obs
   enddo
   close(40)
!
!  Guarda trayectoria de Obs. en dominio de asimilacion
   ncampo='wzobs'
   iob=iobs
   do f=1,fenetre
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
   uxmt_obs=0. ; uymt_obs=0. ; wzmt_obs=0. ; uxux_obs=0. ; uyuy_obs=0. ; uxuy_obs=0. ; wzwz_obs=0.
!
!  Campos medios para la trayectoria de Obs. en dominio espacial reducido
   do f=2,fenetre ! dominio temporal de asimilacion
   do j=1,nyrg
   do i=1,nxrg
      uxmt_obs(i,j,1)=uxmt_obs(i,j,1)+ux_obsg(f,i,j,1)
      uymt_obs(i,j,1)=uymt_obs(i,j,1)+uy_obsg(f,i,j,1)
      wzmt_obs(i,j,1)=wzmt_obs(i,j,1)+wz_obsg(f,i,j,1)
      uxux_obs(i,j,1)=uxux_obs(i,j,1)+ux_obsg(f,i,j,1)*ux_obsg(f,i,j,1)
      uyuy_obs(i,j,1)=uyuy_obs(i,j,1)+uy_obsg(f,i,j,1)*uy_obsg(f,i,j,1)
      uxuy_obs(i,j,1)=uxuy_obs(i,j,1)+ux_obsg(f,i,j,1)*uy_obsg(f,i,j,1)
      wzwz_obs(i,j,1)=wzwz_obs(i,j,1)+wz_obsg(f,i,j,1)*wz_obsg(f,i,j,1)
   enddo
   enddo 
   enddo
!
!  Guarda campos medios para la trayectoria de Obs. en dominio de asimilacion
   open(10,file='uxmt_obs',form='formatted')
   open(20,file='uymt_obs',form='formatted')
   open(30,file='wzmt_obs',form='formatted')
   open(40,file='uxux_obs',form='formatted')
   open(50,file='uyuy_obs',form='formatted')
   open(60,file='uxuy_obs',form='formatted')
   open(70,file='wzwz_obs',form='formatted')
   do j=j_redg1,j_redg2
   do i=i_redg1,i_redg2
      write(10,*) uxmt_obs(i,j,1)/fenetrem
      write(20,*) uymt_obs(i,j,1)/fenetrem
      write(30,*) wzmt_obs(i,j,1)/fenetrem
      write(40,*) uxux_obs(i,j,1)/fenetrem - uxmt_obs(i,j,1)/fenetrem*uxmt_obs(i,j,1)/fenetrem
      write(50,*) uyuy_obs(i,j,1)/fenetrem - uymt_obs(i,j,1)/fenetrem*uymt_obs(i,j,1)/fenetrem
      write(60,*) uxuy_obs(i,j,1)/fenetrem - uxmt_obs(i,j,1)/fenetrem*uymt_obs(i,j,1)/fenetrem
      write(70,*) wzwz_obs(i,j,1)/fenetrem - wzmt_obs(i,j,1)/fenetrem*wzmt_obs(i,j,1)/fenetrem
   enddo
   enddo 
   close(10)
   close(20)
   close(30)
   close(40)
   close(50)
   close(60)
   close(70)
!
!  ******************* Calculo de trayectorias *******************
!
   uxmt_fin=0. ; uymt_fin=0. ; wzmt_fin=0. ; uxux_fin=0. ; uyuy_fin=0. ; uxuy_fin=0. ; wzwz_fin=0.
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
!     Entrada: uxitemp, uyitemp. Salida: wzitemp
      call vortz(uxitemp,uyitemp,wzitemp,ffxp,fsxp,fwxp,ffyp,fsyp,fwyp,nx,ny,nz)
!
!     Entrada: uxftemp, uyftemp. Salida: wzftemp
      call vortz(uxftemp,uyftemp,wzftemp,ffxp,fsxp,fwxp,ffyp,fsyp,fwyp,nx,ny,nz)
!
!     Determina norma vort. para la trayectoria Inicial en dominio de asimilacion 
      normav_ini(itime)=snorm5(wzitemp,nx,ny,nz,nxrg,nyrg,&
                        jw,iw,w_xy,sumw_xy,n_nodos,n_nodos_max,&
                        i_redg1,i_redg2,j_redg1,j_redg2,nxrg_asim,nyrg_asim)
!
!     Determina norma vort. para la trayectoria Asimilada en dominio de asimilacion 
      normav_fin(itime)=snorm5(wzftemp,nx,ny,nz,nxrg,nyrg,&
                        jw,iw,w_xy,sumw_xy,n_nodos,n_nodos_max,&
                        i_redg1,i_redg2,j_redg1,j_redg2,nxrg_asim,nyrg_asim)
!
!     Plano y=0 para la trayectoria Asimilada en dominio reconstruido
      j0=floor(ny/2.)
      do i=1,nx
         bx2_fin(i,itime)=uxftemp(i,j0,1)
         by2_fin(i,itime)=uyftemp(i,j0,1)
      enddo
!
!     Plano x=0 para la trayectoria Inicial en dominio reconstruido
      do j=1,ny
         bx3_ini(j,itime)=uxitemp(i_rec1,j,1)
         by3_ini(j,itime)=uyitemp(i_rec1,j,1)
      enddo
!
!     Plano x=0 para la trayectoria Asimilada en dominio reconstruido
      do j=1,ny
         bx3_fin(j,itime)=uxftemp(i_rec1,j,1)
         by3_fin(j,itime)=uyftemp(i_rec1,j,1)
      enddo
!
!     Campos medios para la trayectoria Asimilada en dominio espacial reconstruido
      if (itime.gt.imodulob) then ! dominio temporal de asimilacion
         do j=1,ny
         do i=1,nx
            uxmt_fin(i,j,1)=uxmt_fin(i,j,1)+uxftemp(i,j,1)
            uymt_fin(i,j,1)=uymt_fin(i,j,1)+uyftemp(i,j,1)
            wzmt_fin(i,j,1)=wzmt_fin(i,j,1)+wzftemp(i,j,1)
            uxux_fin(i,j,1)=uxux_fin(i,j,1)+uxftemp(i,j,1)*uxftemp(i,j,1)
            uyuy_fin(i,j,1)=uyuy_fin(i,j,1)+uyftemp(i,j,1)*uyftemp(i,j,1)
            uxuy_fin(i,j,1)=uxuy_fin(i,j,1)+uxftemp(i,j,1)*uyftemp(i,j,1)
            wzwz_fin(i,j,1)=wzwz_fin(i,j,1)+wzftemp(i,j,1)*wzftemp(i,j,1)
         enddo
         enddo 
      endif
!
      enddo
!
!     Entrada: uxitemp, uyitemp. Salida: wzitemp
      call vortz(uxitemp,uyitemp,wzitemp,ffxp,fsxp,fwxp,ffyp,fsyp,fwyp,nx,ny,nz)
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
!     Entrada: uxftemp, uyftemp. Salida: wzftemp
      call vortz(uxftemp,uyftemp,wzftemp,ffxp,fsxp,fwxp,ffyp,fsyp,fwyp,nx,ny,nz)
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
   enddo
!
!  Guarda norma vort. para la trayectoria Inicial en dominio de asimilacion
   open(40,file='Norma_ini',form='formatted')
   itime=0  
   write(40,*) itime,norma1_ini
   do f=1,fenetrem
   do n=1,imodulob
      itime=n+(f-1)*imodulob
      write(40,*) itime,normav_ini(itime)
   enddo
   enddo
   close(40)
!
!  Guarda norma vort. para la trayectoria Asimilada en dominio de asimilacion
   open(40,file='Norma_fin',form='formatted')
   itime=0  
   write(40,*) itime,norma1_fin
   do f=1,fenetrem
   do n=1,imodulob
      itime=n+(f-1)*imodulob
      write(40,*) itime,normav_fin(itime)
   enddo
   enddo
   close(40)
!
!  Guarda plano y=0 para la trayectoria Asimilada en dominio de asimilacion
   open(80,file='plano_y0_fin.dat',form='formatted')
   do itime=imodulob+1,ditime
   do i=i_rec1,i_rec2
      write(80,*) bx2_fin(i,itime),by2_fin(i,itime)
   enddo
   enddo  
   close(80)
!
!  Guarda plano x=0 para la trayectoria Inicial en dominio de asimilacion
   open(80,file='plano_x0_ini.dat',form='formatted')
   do itime=imodulob+1,ditime
   do j=j_rec1,j_rec2
      write(80,*) bx3_ini(j,itime),by3_ini(j,itime)
   enddo
   enddo  
   close(80)
!
!  Guarda plano x=0 para la trayectoria Asimilada en dominio de asimilacion
   open(80,file='plano_x0_fin.dat',form='formatted')
   do itime=imodulob+1,ditime
   do j=j_rec1,j_rec2
      write(80,*) bx3_fin(j,itime),by3_fin(j,itime)
   enddo
   enddo  
   close(80)
!
!  Guarda campos medios para la trayectoria Asimilada en dominio de asimilacion
   open(10,file='uxmt_fin',form='formatted')
   open(20,file='uymt_fin',form='formatted')
   open(30,file='wzmt_fin',form='formatted')
   open(40,file='uxux_fin',form='formatted')
   open(50,file='uyuy_fin',form='formatted')
   open(60,file='uxuy_fin',form='formatted')
   open(70,file='wzwz_fin',form='formatted')
   do j=j_rec1,j_rec2
   do i=i_rec1,i_rec2
      write(10,*) uxmt_fin(i,j,1)/ditime_asim
      write(20,*) uymt_fin(i,j,1)/ditime_asim
      write(30,*) wzmt_fin(i,j,1)/ditime_asim
      write(40,*) uxux_fin(i,j,1)/ditime_asim - uxmt_fin(i,j,1)/ditime_asim*uxmt_fin(i,j,1)/ditime_asim
      write(50,*) uyuy_fin(i,j,1)/ditime_asim - uymt_fin(i,j,1)/ditime_asim*uymt_fin(i,j,1)/ditime_asim
      write(60,*) uxuy_fin(i,j,1)/ditime_asim - uxmt_fin(i,j,1)/ditime_asim*uymt_fin(i,j,1)/ditime_asim
      write(70,*) wzwz_fin(i,j,1)/ditime_asim - wzmt_fin(i,j,1)/ditime_asim*wzmt_fin(i,j,1)/ditime_asim
   enddo
   enddo 
   close(10)
   close(20)
   close(30)
   close(40)
   close(50)
   close(60)
   close(70)
!
!
end PROGRAM trayectoria
!
!********************************************************************


