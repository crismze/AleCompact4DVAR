!********************************************************************
!
PROGRAM asimila2d
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
   integer,parameter :: nxrg_asim=i_redg2-i_redg1+1,nyrg_asim=j_redg2-j_redg1+1 ! malla asimilada grosera
   integer,parameter :: nxm=nx-1,nym=ny-1
   integer,parameter :: mx=nx+1,my=ny+1,mz=nz+2
   integer,parameter :: myr=nyr+1
   integer,parameter :: myrg=nyrg+1
   integer,parameter :: m1=2,m2=2
   integer,parameter :: fenetre=30 ! taille de la sequence des obs.
   integer,parameter :: fenetrem=fenetre-1
   integer,parameter :: imodulob=128 ! intervalle temporel entre obs.
   integer,parameter :: ditime=fenetrem*imodulob
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
   real(8),dimension(fenetre,nxr,nyr,nz) :: ux_obs,uy_obs,wz_obs
   real(8),dimension(fenetre,nxrg,nyrg,nz) :: ux_obsg,uy_obsg,wz_obsg
   real(8),dimension(nxr,nyr,nz) :: ux0r,uy0r,gx0r,gy0r
   real(8),dimension(nxr,nyr,nz) :: ux0r_obs,uy0r_obs,gx0r_obs,gy0r_obs
   real(8),dimension(nxr,nyr,nz) :: ux0rb,uy0rb,gx0rb,gy0rb
   real(8),dimension(nx,ny,nz) :: ux0,uy0,gx0,gy0,wz0
   real(8),dimension(nx,ny,nz) :: ux0b,uy0b,gx0b,gy0b
   real(8),dimension(nx,ny,nz) :: uxtemp,uytemp,gxtemp,gytemp,wztemp
   real(8),dimension(nx,ny,nz) :: uxtempE,uytempE,gxtempE,gytempE
   real(8),dimension(nx,ny,nz) :: uxtempS,uytempS,gxtempS,gytempS
   real(8),dimension(nyr,ditime) :: bx1_obs,by1_obs
   real(8),dimension(nyr,ditime) :: bx1r,by1r
   real(8),dimension(nyr,ditime) :: bx1rb,by1rb
   real(8),dimension(ny,ditime) :: bx1,by1
   real(8),dimension(ny,ditime) :: bx1b,by1b
   real(8),dimension(ny,ditime) :: bx1temp,by1temp
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
   integer :: itime,i,j,f,n,ired,jred,irec,jrec
   integer :: i_redg,j_redg,i_red,j_red,i_rec,j_rec,nodo
!  ********************* Obs. *********************
   real(8),dimension(nx,ny) :: rx,ry
   integer,dimension(ditimev) :: itimev
   real(8),dimension(nyr,ditime+1) :: bx1p,by1p
   integer,dimension(nxrg,nyrg,n_nodos_max) :: iw,jw
   real(8),dimension(nxrg,nyrg,n_nodos_max) :: w_xy
   integer,dimension(nxrg,nyrg) :: n_nodos
   real(8),dimension(nxrg,nyrg) :: sumw_xy
   integer,dimension(nxr,nyr,n_nodos_q_max) :: iwq,jwq
   real(8),dimension(nxr,nyr,n_nodos_q_max) :: w_xyq
   integer,dimension(nxr,nyr) :: n_nodos_q
   real(8),dimension(nxr,nyr) :: sumw_xyq
   integer,dimension(ditime,n_nodos_t_max) :: itw
   real(8),dimension(ditime,n_nodos_t_max) :: w_t
   integer,dimension(ditime) :: n_nodos_t
   real(8),dimension(ditime) :: sumw_t
   real(8) :: ruidox,ruidoy,pendx,pendy
   real(8),dimension(nxrg,nyrg) :: cov_obs
   real(8) :: cov_iniu,cov_inig,cov_suav,cov_tay,cov_gau
   real(8) :: sigma_r,r,w_r,sumw_r
   real(8) :: sigma_q,rq,w_q,sumw_q
   real(8) :: sigma_t,d,w_d,sumw_d
   integer :: longueur,iob,iobi,iobf,rimod,tiv,iv
   character(len=4) suffix
   character(len=80) nfichier,ncampo,nvector
   real(8) :: snorm4,snorm5,snorm6,error,norma_ini,norma_obs
   real(8),dimension(ditime) :: normav
!  ********************* LBFGS *********************
   double precision xp(ndim),gr(ndim),diag(ndim),w(nwork)
   double precision fc,eps,xtol,ddot,gr_norm0,gr_norm,fc_ant
   double precision fc_obs,fc_iniu,fc_inig,fc_suav,fc_tay,fc_gau
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
!  Parametres de LBFGS
   iprint(1)=0 ! specifies the frequency of the output
   iprint(2)=0 ! specifies the type of output
   diagco=.false. ! we do not provide the diagonal matrix Hk0 at each iteration
   eps=1.0e-4 ! determines the accuracy with which the solution is to be found
   xtol=zero ! estimate of the machine precision
   icall=0 ! control the number of evaluations of fc and gr
   iflag=0 ! control the calling to LBFGS
   iter=-1 ! control the number of LBFGS iterations
   fc_ant=1.0e8 ! valor de fc en la evaluacion icall-1
!
!  Lee campos PIV que han sido reducidos, normalizados y suavizados en
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
!  Lee campos PIV que han sido reducidos, normalizados, suavizados 
!  e interpolados a fin de aumentar la resolucion espacial en
!  /home/ale/Programas/fiuba-2014/reduce_normaliza_suaviza_interpola_obs_asim86.m
   ncampo='cmcobs'
   iob=iobs
   do f=1,fenetre
      call numcar (iob,suffix)
      longueur=index(ncampo,' ')-1
      nfichier=ncampo(1:longueur)//suffix
      longueur=index(nfichier,' ')-1
      open(100,file=nfichier(1:longueur),form='formatted')
      do j=1,nyr
      do i=1,nxr
         read(100,*) ux_obs(f,i,j,1),uy_obs(f,i,j,1)
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
!  Determina y guarda norma de la vorticidad en el dominio asimilado
   open(40,file='Norma_Obs',form='formatted')
   do f=1,fenetre
      norma_obs=snorm6(wz_obsg(f,:,:,:),nxrg,nyrg,nz,&
                i_redg1,i_redg2,j_redg1,j_redg2,nxrg_asim,nyrg_asim)
      itime=(f-1)*imodulob
      write(40,*) itime,norma_obs
   enddo
   close(40)
!
!  ************* Gaussian weights, coordenada temporal *************
   sigma_t=42.6
   do itime=1,ditime
!
      sumw_d=0. 
      nodo=0
      do i=1,ditime
         d=sqrt(real((i-itime)**2))
         if (d<=d_g) then
            w_d=exp(-(d**2)/(2*sigma_t**2))
            sumw_d = sumw_d + w_d
            nodo=nodo+1
            itw(itime,nodo)=i
            w_t(itime,nodo)=w_d
         endif
      enddo
      sumw_t(itime)=sumw_d
      n_nodos_t(itime)=nodo
!
!     Guarda weights empleados para calcular bx1r_av, by1r_av en Subrot. snorm
      if (itime==floor(ditime/2.)) then
         open(20,file='weights_t.dat',form='formatted')
         do n=1,nodo
            write(20,*) itw(itime,n),w_t(itime,n)/sumw_d
         enddo
         close(20)
      endif
!
   enddo 
!  *****************************************************************
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
!
!     Guarda weights empleados para calcular wz_av en Subrot. snorm
      if (j_redg==floor(nyrg/2.) .and. i_redg==floor(nxrg/2.)) then
         open(20,file='weights_xy.dat',form='formatted')
         do n=1,nodo
            write(20,*) jw(i_redg,j_redg,n),iw(i_redg,j_redg,n),w_xy(i_redg,j_redg,n)/sumw_r
         enddo
         close(20)
      endif
!
   enddo
   enddo 
!  **********************************************************************
!
!  ******** Gaussian weights, coordenada espacial, malla fina ********
   sigma_q=2.
   jwq=0. ; iwq=0. ; w_xyq=0.
   do j_red = 1,nyr
   do i_red = 1,nxr
!
      sumw_q=0.
      nodo=0
      do j=1,nyr
      do i=1,nxr
         rq=sqrt(real((i-i_red)**2+(j-j_red)**2))
         if (rq<=rq_g) then
            w_q=exp(-(rq**2)/(2*sigma_q**2))
            sumw_q = sumw_q + w_q
            nodo=nodo+1
            jwq(i_red,j_red,nodo)=j
            iwq(i_red,j_red,nodo)=i
            w_xyq(i_red,j_red,nodo)=w_q
         endif
      enddo
      enddo
      sumw_xyq(i_red,j_red)=sumw_q
      n_nodos_q(i_red,j_red)=nodo
!
!     Guarda weights empleados para calcular ux0r_av, uy0r_av en Subrot. snorm
      if (j_red==floor(nyr/2.) .and. i_red==floor(nxr/2.)) then
         open(20,file='weights_xyq.dat',form='formatted')
         do n=1,nodo
            write(20,*) jwq(i_red,j_red,n),iwq(i_red,j_red,n),w_xyq(i_red,j_red,n)/sumw_q
         enddo
         close(20)
      endif
!
   enddo
   enddo 
!  *******************************************************************
!
!  ************************************************************************
!  ******* Construye cond. de entrada a partir de obs. interpoladas *******
!  ************************************************************************
!
!  Cond. de entrada de las obs.
   do f=1,fenetre
      itime=(f-1)*imodulob+1
      i=1
      bx1p(:,itime)=ux_obs(f,i,:,1)
      by1p(:,itime)=uy_obs(f,i,:,1)
   enddo
!
!  Region corresp. al desfasaje espacial entre obs.
   rimod=imodulob/imoduloe
   do f=1,fenetrem
   do n=1,imoduloe-1
      itime=(f-1)*imodulob+n*rimod+1
      i=imoduloe+1-n
      bx1p(:,itime)=ux_obs(f+1,i,:,1)
      by1p(:,itime)=uy_obs(f+1,i,:,1)
   enddo
   enddo
!
!  ************* Guarda cond. de entrada *************
   iv=0
   do f=1,fenetrem
   do n=1,imoduloe
      iv=iv+1
      itimev(iv)=(f-1)*imodulob+(n-1)*rimod+1
   enddo
   enddo
   f=fenetre ; n=1
   iv=iv+1
   itimev(iv)=(f-1)*imodulob+(n-1)*rimod+1
!
   open(80,file='entrada_obs.dat',form='formatted')
   do iv=1,ditimev
   tiv=itimev(iv) 
   do j=1,nyr
      write(80,*) bx1p(j,tiv),by1p(j,tiv)
   enddo
   enddo  
   close(80)
!  ***************************************************
!
!  Lee cond. de entrada que ha sido interpolada a fin de cubrir 
!  el espacio temporal, y suavizada aplicando un filtro en   
!  /home/ale/Programas/fiuba-2014/interpola_suaviza_entrada_asim86.m
   open(80,file='entrada_filt.dat',form='formatted')
   do itime=1,ditime
   do j=1,nyr
      read(80,*) bx1_obs(j,itime),by1_obs(j,itime)
   enddo
   enddo  
   close(80)
!  ************************************************************************
!
!  ************************************************************************
!  ***** Construye cond. inicial a partir de obs. inicial interpolada *****
!  ************************************************************************
!
!  Guarda cond. inicial interpolada
   open(100,file='inicial_obs.dat',form='formatted')
   do j=1,nyr
   do i=1,nxr
      write(100,*) ux_obs(1,i,j,1),uy_obs(1,i,j,1)
   enddo
   enddo 
   close(100)
!
!  Lee cond. inicial interpolada que ha sido suavizada aplicando un filtro en   
!  /home/ale/Programas/fiuba-2014/suaviza_inicial_asim86.m
   open(100,file='inicial_filt.dat',form='formatted')
   do j=1,nyr
   do i=1,nxr
      read(100,*) ux0r_obs(i,j,1),uy0r_obs(i,j,1)
   enddo
   enddo 
   close(100)
!  ************************************************************************
!
!  ************************************************************************
!  ************ Constantes de covarianza en la funcion costo F ************
!  ************************************************************************
!
!  Covariance matrix associated to the observations
   cov_obs=1.
!
!  Covariance matrix associated to the initialization
   cov_suav=100.
!
!  Covariance matrix associated to the dynamical model
   cov_gau=10.
!
   cov_iniu=0.
   cov_inig=0.
   cov_tay=0.
!  ************************************************************************
!
!  Initial estimate of the inflow condition
   bx1r=bx1_obs
   by1r=by1_obs
!
!  Initial estimate of the initial condition 
   ux0r=ux0r_obs
   uy0r=uy0r_obs
   gx0r=0.
   gy0r=0.
!
!  Redimensiona la estimacion inicial del vector solucion para entrar a LBFGS 
   xp(1:nxyr)=reshape(ux0r,(/nxyr/))
   xp(nxyr+1:2*nxyr)=reshape(uy0r,(/nxyr/))
   xp(2*nxyr+1:3*nxyr)=reshape(gx0r,(/nxyr/))
   xp(3*nxyr+1:4*nxyr)=reshape(gy0r,(/nxyr/))
   xp(4*nxyr+1:4*nxyr+nytr)=reshape(bx1r,(/nytr/))
   xp(4*nxyr+nytr+1:4*nxyr+2*nytr)=reshape(by1r,(/nytr/))
!
   do ! minimization while loop
!
!  *****************************************************************************
!  **** Calculate the function value fc and its gradient gr at the point xp ****
!  *****************************************************************************
!
!  ****************************** Run Y=F(X0) ******************************
!
!  Entrada: ux0r, uy0r, gx0r, gy0r, bx1r, by1r
!  Salida: ux0, uy0, gx0, gy0, bx1, by1
   call reconstruye4 (ux0r,uy0r,gx0r,gy0r,bx1r,by1r,&
          ux0,uy0,gx0,gy0,bx1,by1,&
          jyc,jyt,nxr,nyr,nx,ny,nz,ditime)
!
!  Entrada: ux0, uy0, gx0, gy0, bx1, by1. Salida: ux, uy
!
   uxtemp=ux0 ; uytemp=uy0 ; gxtemp=gx0 ; gytemp=gy0
!
!  Determina norma de la vorticidad en el dominio asimilado  
   call vortz(uxtemp,uytemp,wztemp,ffxp,fsxp,fwxp,ffyp,fsyp,fwyp,nx,ny,nz)
   norma_ini=snorm5(wztemp,nx,ny,nz,nxrg,nyrg,&
             jw,iw,w_xy,sumw_xy,n_nodos,n_nodos_max,&
             i_redg1,i_redg2,j_redg1,j_redg2,nxrg_asim,nyrg_asim)
!
   do f=1,fenetrem
!
      do n=1,imodulob
      itime=n+(f-1)*imodulob
!
      call inc3dorigin (uxtemp,uytemp,gxtemp,gytemp,&
             uxtemp,uytemp,gxtemp,gytemp,&
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
!     Determina norma de la vorticidad en el dominio asimilado 
      call vortz(uxtemp,uytemp,wztemp,ffxp,fsxp,fwxp,ffyp,fsyp,fwyp,nx,ny,nz)
      normav(itime)=snorm5(wztemp,nx,ny,nz,nxrg,nyrg,&
                    jw,iw,w_xy,sumw_xy,n_nodos,n_nodos_max,&
                    i_redg1,i_redg2,j_redg1,j_redg2,nxrg_asim,nyrg_asim)
!
      enddo
!
      ux(f,:,:,:)=uxtemp(:,:,:)
      uy(f,:,:,:)=uytemp(:,:,:)
      gx(f,:,:,:)=gxtemp(:,:,:)
      gy(f,:,:,:)=gytemp(:,:,:)
!
   enddo
!
!  Entrada: ux, uy. Salida: wz
   do f=1,fenetrem
      call vortz(ux(f,:,:,:),uy(f,:,:,:),wz(f,:,:,:),&
             ffxp,fsxp,fwxp,ffyp,fsyp,fwyp,nx,ny,nz)
   enddo
!
!  Entrada: wz, ux0r, uy0r, gx0r, gy0r, bx1r, by1r. Salida: fc
   call snorm (fc,wz,wz_obsg,ux0r,uy0r,ux0r_obs,uy0r_obs,gx0r,gy0r,&
          gx0r_obs,gy0r_obs,bx1r,by1r,bx1_obs,by1_obs,&
          fenetre,fenetrem,nx,ny,nz,nxr,nyr,nxrg,nyrg,jyt,rdxy,ditime,&
          cov_obs,cov_iniu,cov_inig,cov_suav,cov_tay,cov_gau,&
          jw,iw,w_xy,sumw_xy,n_nodos,n_nodos_max,&
          jwq,iwq,w_xyq,sumw_xyq,n_nodos_q,n_nodos_q_max,&
          itw,w_t,sumw_t,n_nodos_t,n_nodos_t_max,&
          i_redg1,i_redg2,j_redg1,j_redg2,nxrg_asim,nyrg_asim,&
          fc_obs,fc_iniu,fc_inig,fc_suav,fc_tay,fc_gau)
!
!  ******************* Run the adjoint code X0b=F_B(X0,Y) ******************
!
!  Entrada: wz, ux0r, uy0r, gx0r, gy0r, bx1r, by1r, fc
!  Salida: wzb, ux0rb, uy0rb, gx0rb, gy0rb, bx1rb, by1rb
!
   call SNORM_B(fc,wz,wzb,wz_obsg,ux0r,ux0rb,uy0r,uy0rb,&
          ux0r_obs,uy0r_obs,gx0r,gx0rb,gy0r,gy0rb,&
          gx0r_obs,gy0r_obs,bx1r,bx1rb,by1r,by1rb,bx1_obs,by1_obs,&
          fenetre,fenetrem,nx,ny,nz,nxr,nyr,nxrg,nyrg,jyt,rdxy,ditime,&
          cov_obs,cov_iniu,cov_inig,cov_suav,cov_tay,cov_gau,&
          jw,iw,w_xy,sumw_xy,n_nodos,n_nodos_max,&
          jwq,iwq,w_xyq,sumw_xyq,n_nodos_q,n_nodos_q_max,&
          itw,w_t,sumw_t,n_nodos_t,n_nodos_t_max,&
          i_redg1,i_redg2,j_redg1,j_redg2,nxrg_asim,nyrg_asim)
!
!  Entrada: ux, uy, wzb. Salida: uxb, uyb
   do f=fenetrem,1,-1
      call VORTZ_B(ux(f,:,:,:),uxb(f,:,:,:),uy(f,:,:,:),uyb(f,:,:,:),wz(f,:,:,:),&
             wzb(f,:,:,:),ffxp,fsxp,fwxp,ffyp,fsyp,fwyp,nx,ny,nz)
   enddo
!
!  Entrada: ux0, uy0, gx0, gy0, bx1, by1, uxb, uyb
!  Salida: ux0b, uy0b, gx0b, gy0b, bx1b, by1b 
!
   uxtempE=0. ; uytempE=0. ; gxtempE=0. ; gytempE=0.
!
   do f=fenetrem,1,-1
!
!     ************ recalcul des valeurs intermediaires ************
      if (f.gt.1) then
         ux_rec(1,:,:,:)=ux(f-1,:,:,:)
         uy_rec(1,:,:,:)=uy(f-1,:,:,:)
         gx_rec(1,:,:,:)=gx(f-1,:,:,:)
         gy_rec(1,:,:,:)=gy(f-1,:,:,:)
      else
         ux_rec(1,:,:,:)=ux0 
         uy_rec(1,:,:,:)=uy0 
         gx_rec(1,:,:,:)=gx0 
         gy_rec(1,:,:,:)=gy0
      endif
!
      do n=2,imodulob 

      print *,'Pas recalcul n = ',n-1
      itime=n-1+(f-1)*imodulob
!
      call inc3dorigin (ux_rec(n-1,:,:,:),uy_rec(n-1,:,:,:),gx_rec(n-1,:,:,:),gy_rec(n-1,:,:,:),&
             ux_rec(n,:,:,:),uy_rec(n,:,:,:),gx_rec(n,:,:,:),gy_rec(n,:,:,:),&
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
      enddo 
!     *************************************************************
!
      uxtempE(:,:,:)=uxtempE(:,:,:)+uxb(f,:,:,:)
      uytempE(:,:,:)=uytempE(:,:,:)+uyb(f,:,:,:)
      gxtempE(:,:,:)=gxtempE(:,:,:)+gxb(f,:,:,:)
      gytempE(:,:,:)=gytempE(:,:,:)+gyb(f,:,:,:)
!
      do n=imodulob,1,-1
!
      print *,'Pas adjoint n = ',n
      itime=n+(f-1)*imodulob
!
      uxtempS=0. ; uytempS=0. ; gxtempS=0. ; gytempS=0.
!
      call inc3dadjoint (ux_rec(n,:,:,:),uy_rec(n,:,:,:),gx_rec(n,:,:,:),gy_rec(n,:,:,:),&
             uxtempE,uytempE,gxtempE,gytempE,&
             uxtempS,uytempS,gxtempS,gytempS,&
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
      uxtempE=uxtempS ; uytempE=uytempS ; gxtempE=gxtempS ; gytempE=gytempS
!
      enddo
!
   enddo
!
   ux0b = ux0b + uxtempE
   uy0b = uy0b + uytempE
   gx0b = gx0b + gxtempE
   gy0b = gy0b + gytempE
!
!  Entrada: ux0r, uy0r, gx0r, gy0r, bx1r, by1r, ux0b, uy0b, gx0b, gy0b, bx1b, by1b,
!           ux0rb, uy0rb, gx0rb, gy0rb, bx1rb, by1rb 
!  Salida: ux0rb, uy0rb, gx0rb, gy0rb, bx1rb, by1rb 
!
   call RECONSTRUYE4_B(ux0r,ux0rb,uy0r,uy0rb,gx0r,gx0rb,gy0r,gy0rb,bx1r,bx1rb,&
          by1r,by1rb,ux0,ux0b,uy0,uy0b,gx0,gx0b,gy0,gy0b,bx1,bx1b,by1,by1b,&
          jyc,jyt,nxr,nyr,nx,ny,nz,ditime)
!
!  ******************** Redimensiona para entrar a LBFGS *******************
!
   gr(1:nxyr)=reshape(ux0rb,(/nxyr/))
   gr(nxyr+1:2*nxyr)=reshape(uy0rb,(/nxyr/))
   gr(2*nxyr+1:3*nxyr)=reshape(gx0rb,(/nxyr/))
   gr(3*nxyr+1:4*nxyr)=reshape(gy0rb,(/nxyr/))
   gr(4*nxyr+1:4*nxyr+nytr)=reshape(bx1rb,(/nytr/))
   gr(4*nxyr+nytr+1:4*nxyr+2*nytr)=reshape(by1rb,(/nytr/))
!
!  *****************************************************************************
!  *****************************************************************************
!  *****************************************************************************
!
   icall=icall+1
!
!  We allow at most 2000 evaluations of fc and gr
   if (icall.gt.2000) exit
!
   if (icall==1) then
!     Calcula norma del gradiente inicial
      gr_norm0=dsqrt(ddot(ndim,gr,1,gr,1))
!
!     Abre archivo nuevo para guardar resultados de la optimizacion 
      open(80,file='param_lbfgs.dat',status='replace')
      open(90,file='param_lbfgs1.dat',status='replace') 
   endif
!
   if (fc<fc_ant) then ! Pasa a la proxima iteracion
      iter=iter+1
      fc_ant=fc
      gr_norm=dsqrt(ddot(ndim,gr,1,gr,1))
!
!     Guarda resultados de la optimizacion 
      open(80,file='param_lbfgs.dat',access='append')
      write(80,1000) iter,icall,fc,gr_norm/gr_norm0
      close(80)
1000  format(2i5,2f20.13)
!
!     Guarda estimacion del vector solucion 
      nvector='xp_iter'
      call numcar (iter,suffix)
      longueur=index(nvector,' ')-1
      nfichier=nvector(1:longueur)//suffix  
      longueur=index(nfichier,' ')-1
      open(90,file=nfichier(1:longueur),form='unformatted')
      write(90) xp
      close(90)
!
!     Guarda estimacion del error
      nvector='Error_iter'
      call numcar (iter,suffix)
      longueur=index(nvector,' ')-1
      nfichier=nvector(1:longueur)//suffix  
      longueur=index(nfichier,' ')-1
      open(40,file=nfichier(1:longueur),form='formatted')
      do f=1,fenetrem
         itime=f*imodulob
         error=snorm4(wz_obsg(f+1,:,:,:),wz(f,:,:,:),nx,ny,nz,nxrg,nyrg,&
                 jw,iw,w_xy,sumw_xy,n_nodos,n_nodos_max,&
                 i_redg1,i_redg2,j_redg1,j_redg2,nxrg_asim,nyrg_asim)
         write(40,*) itime,error
      enddo
      close(40)
!
!     Guarda estimacion de la norma de la vorticidad en el dominio asimilado
      nvector='Norma_iter'
      call numcar (iter,suffix)
      longueur=index(nvector,' ')-1
      nfichier=nvector(1:longueur)//suffix  
      longueur=index(nfichier,' ')-1
      open(40,file=nfichier(1:longueur),form='formatted')
      itime=0  
      write(40,*) itime,norma_ini
      do f=1,fenetrem
      do n=1,imodulob
         itime=n+(f-1)*imodulob
         write(40,*) itime,normav(itime)
      enddo
      enddo
      close(40)
!
   endif
!
!  Guarda resultados de la optimizacion 
   open(90,file='param_lbfgs1.dat',access='append')
   write(90,1100) iter,icall,fc,fc_obs,fc_iniu,fc_inig,fc_suav,fc_tay,fc_gau
   close(90)
1100  format(2i5,7f20.13)
!
!  *****************************************************************************
!  *****************************************************************************
!  *****************************************************************************
!
!  Entrada: xp, fc, gr. Salida: xp, iflag 
   call LBFGS (ndim,msave,xp,fc,gr,diagco,diag,iprint,eps,xtol,w,iflag)
!
!  Redimensiona xp
   ux0r=reshape(xp(1:nxyr),(/nxr,nyr,nz/))
   uy0r=reshape(xp(nxyr+1:2*nxyr),(/nxr,nyr,nz/))
   gx0r=reshape(xp(2*nxyr+1:3*nxyr),(/nxr,nyr,nz/))
   gy0r=reshape(xp(3*nxyr+1:4*nxyr),(/nxr,nyr,nz/))
   bx1r=reshape(xp(4*nxyr+1:4*nxyr+nytr),(/nyr,ditime/))
   by1r=reshape(xp(4*nxyr+nytr+1:4*nxyr+2*nytr),(/nyr,ditime/))
!
   if (iflag.le.0) exit
!
   enddo ! minimization while loop
!
!
end PROGRAM asimila2d
!
!********************************************************************


