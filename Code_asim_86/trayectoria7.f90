!********************************************************************
! Guarda campos de veloc. y presion para journal FDR-14
! Calcula drag, lift a partir de informacion disponible solo en inflow
! (cara CD del volumen de control) --> Subrot. aerof_inflow1
!
PROGRAM trayectoria7
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
   integer,parameter :: i_coe1=1,i_coe2=nx ! guardo campos, calculo normas y coef. aerod. en dominio COEx=REC
!    integer,parameter :: i_coe1=i_rec1,i_coe2=i_rec2 ! guardo campos, calculo normas y coef. aerod. en dominio COEx=ASIM
!    integer,parameter :: i_coe1=i_rec1+30,i_coe2=i_rec2 ! guardo campos, calculo normas y coef. aerod. en dominio COEx=ASIM-30i   
!    integer,parameter :: j_coe1=1,j_coe2=ny ! guardo campos, calculo normas y coef. aerod. en dominio COEy=REC
!    integer,parameter :: j_coe1=j_rec1,j_coe2=j_rec2 ! guardo campos, calculo normas y coef. aerod. en dominio COEy=ASIM
   integer,parameter :: j_coe1=j_rec1+30,j_coe2=j_rec2-30 ! guardo campos, calculo normas y coef. aerod. en dominio COEy=ASIM-30j
!    integer,parameter :: j_coe1=j_rec1+20,j_coe2=j_rec2-20 ! guardo campos, calculo normas y coef. aerod. en dominio COEy=ASIM-20j
!    integer,parameter :: j_coe1=j_rec1+10,j_coe2=j_rec2-10 ! guardo campos, calculo normas y coef. aerod. en dominio COEy=ASIM-10j
   integer,parameter :: nxr_coe=i_coe2-i_coe1+1,nyr_coe=j_coe2-j_coe1+1 ! dominio espacial COE, malla fina
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
   real(8),dimension(nx,ny,nz) :: uxitemp,uyitemp,gxitemp,gyitemp,wzitemp,ppitemp
   real(8),dimension(nx,ny,nz) :: uxftemp,uyftemp,gxftemp,gyftemp,wzftemp,ppftemp
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
   integer :: longueur,longueur1,iob,iobi,iobf,iobfilt,rimod,t0,tf,iv
   character(len=4) suffix
   character(len=80) nfichier,ncampo,nvector
   character(len=120) dir
   real(8) :: snorm4,snorm5,snorm6,error,norma_obs,norma1_ini,norma1_fin
   real(8),dimension(ditime) :: normav_ini,normav_fin
!  ****************** aerof_inflow *****************
   real(8),dimension(nx,ny,nz) :: uxm1i,uym1i,uxm2i,uym2i
   real(8),dimension(nx,ny,nz) :: uxm1f,uym1f,uxm2f,uym2f   
   real(8),dimension(ditime) :: uxi_ref,uyi_ref,ppi_ref 
   real(8),dimension(ditime) :: uxf_ref,uyf_ref,ppf_ref   
   real(8) :: norma_uxi,norma_uyi,norma_ppi  
   real(8) :: norma_uxf,norma_uyf,norma_ppf 
   real(8) :: xDrag,yLift  
   integer :: it
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
!  Directorio para guardar campos
   dir = '/home/ale/Asim_86/aerof1/'
   longueur1=index(dir,' ')-1
!
!  Dominio COE
   print *,'i_coe1=',i_coe1
   print *,'i_coe2=',i_coe2
   print *,'j_coe1=',j_coe1
   print *,'j_coe2=',j_coe2   
!
!  Dominio ASIM
   print *,'i_rec1=',i_rec1
   print *,'i_rec2=',i_rec2
   print *,'j_rec1=',j_rec1
   print *,'j_rec2=',j_rec2   
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
   call reconstruye5 (ux0r,uy0r,gx0r,gy0r,bx1r,by1r,&
          ux0,uy0,gx0,gy0,bx1,by1,&
          jyt,nxr,nyr,nx,ny,nz,ditime)
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
!  Guarda cond. inicial para la trayectoria Inicial en dominio COE
   iobi=iobs
   call numcar (iobi,suffix)
   ncampo='uxini'
   longueur=index(ncampo,' ')-1
   nfichier=ncampo(1:longueur)//suffix
   longueur=index(nfichier,' ')-1
   open(100,file=dir(1:longueur1)//nfichier(1:longueur),form='formatted')
   do j=j_coe1,j_coe2
   do i=i_coe1,i_coe2
      write(100,*) uxitemp(i,j,1)
   enddo
   enddo 
   close(100)
   ncampo='uyini'
   longueur=index(ncampo,' ')-1
   nfichier=ncampo(1:longueur)//suffix
   longueur=index(nfichier,' ')-1
   open(100,file=dir(1:longueur1)//nfichier(1:longueur),form='formatted')
   do j=j_coe1,j_coe2
   do i=i_coe1,i_coe2
      write(100,*) uyitemp(i,j,1)
   enddo
   enddo 
   close(100)
!
!  ******************* Parametros trayectoria Asimilada *******************
!
!  Lee estimacion del vector solucion en la evaluacion icall de fc y gr (best point found by LBFGS)
   nvector='xp_iter'
   iter=540
!    iter=280
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
   call reconstruye5 (ux0r,uy0r,gx0r,gy0r,bx1r,by1r,&
          ux0,uy0,gx0,gy0,bx1,by1,&
          jyt,nxr,nyr,nx,ny,nz,ditime)
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
!  Guarda cond. inicial para la trayectoria Asimilada en dominio COE
   iobf=iobs
   call numcar (iobf,suffix)
   ncampo='uxfin'   
   longueur=index(ncampo,' ')-1
   nfichier=ncampo(1:longueur)//suffix
   longueur=index(nfichier,' ')-1
   open(100,file=dir(1:longueur1)//nfichier(1:longueur),form='formatted')
   do j=j_coe1,j_coe2
   do i=i_coe1,i_coe2
      write(100,*) uxftemp(i,j,1)
   enddo
   enddo 
   close(100)
   ncampo='uyfin'   
   longueur=index(ncampo,' ')-1
   nfichier=ncampo(1:longueur)//suffix
   longueur=index(nfichier,' ')-1
   open(100,file=dir(1:longueur1)//nfichier(1:longueur),form='formatted')
   do j=j_coe1,j_coe2
   do i=i_coe1,i_coe2
      write(100,*) uyftemp(i,j,1)
   enddo
   enddo 
   close(100)
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
!  Guarda trayectoria de Obs. en dominio de asimilacion
   ncampo='wzobs'
   iob=iobs
   do f=1,fenetre
      call numcar (iob,suffix)
      longueur=index(ncampo,' ')-1
      nfichier=ncampo(1:longueur)//suffix
      longueur=index(nfichier,' ')-1
      open(100,file=dir(1:longueur1)//nfichier(1:longueur),form='formatted')
      do j=j_redg1,j_redg2
      do i=i_redg1,i_redg2
         write(100,*) wz_obsg(f,i,j,1)
      enddo
      enddo 
      close(100)
      iob=iob+1
   enddo
!  **********************************************************************
!
!  Lee valores de referencia para calculo drag, lift --> generados en una corrida previa de trayectoria7.f90
   open(20,file=dir(1:longueur1)//'Norma_ASIM_ini')
   open(30,file=dir(1:longueur1)//'Norma_ASIM_fin')
   do it=1,ditime
      read(20,1100) itime,uxi_ref(it),uyi_ref(it),ppi_ref(it)
      read(30,1100) itime,uxf_ref(it),uyf_ref(it),ppf_ref(it)
   enddo
   close(20)
   close(30)   
!
!  Abre archivos nuevos
   open(40,file=dir(1:longueur1)//'Norma_ini',status='replace')
   open(50,file=dir(1:longueur1)//'Norma_fin',status='replace') 
   open(60,file=dir(1:longueur1)//'aerof_ini',status='replace')
   open(70,file=dir(1:longueur1)//'aerof_fin',status='replace') 
   open(80,file=dir(1:longueur1)//'rmomentum_ini',status='replace') 
   open(90,file=dir(1:longueur1)//'rmomentum_fin',status='replace')   
   open(100,file=dir(1:longueur1)//'sforce_ini',status='replace') 
   open(110,file=dir(1:longueur1)//'sforce_fin',status='replace')   
!
!  ******************* Calculo de trayectorias *******************
!
   do f=1,fenetrem
!
      do n=1,imodulob
!
      itime=n+(f-1)*imodulob
!
      call inc3dorigin_tray (uxitemp,uyitemp,gxitemp,gyitemp,&
             uxitemp,uyitemp,ppitemp,gxitemp,gyitemp,&
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
      call inc3dorigin_tray (uxftemp,uyftemp,gxftemp,gyftemp,&
             uxftemp,uyftemp,ppftemp,gxftemp,gyftemp,&
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
!     Determina y guarda normas para la trayectoria Inicial en dominio COE 
      norma_uxi=snorm6(uxitemp,nx,ny,nz,i_coe1,i_coe2,j_coe1,j_coe2,nxr_coe,nyr_coe)
      norma_uyi=snorm6(uyitemp,nx,ny,nz,i_coe1,i_coe2,j_coe1,j_coe2,nxr_coe,nyr_coe)
      norma_ppi=snorm6(ppitemp,nx,ny,nz,i_coe1,i_coe2,j_coe1,j_coe2,nxr_coe,nyr_coe)
!
      open(40,file=dir(1:longueur1)//'Norma_ini',access='append')
      write(40,1100) itime,norma_uxi,norma_uyi,norma_ppi
      close(40)
 1100 format(I5,3F16.10) 
!
!     Determina y guarda normas para la trayectoria Final en dominio COE 
      norma_uxf=snorm6(uxftemp,nx,ny,nz,i_coe1,i_coe2,j_coe1,j_coe2,nxr_coe,nyr_coe)
      norma_uyf=snorm6(uyftemp,nx,ny,nz,i_coe1,i_coe2,j_coe1,j_coe2,nxr_coe,nyr_coe)
      norma_ppf=snorm6(ppftemp,nx,ny,nz,i_coe1,i_coe2,j_coe1,j_coe2,nxr_coe,nyr_coe)
!
      open(50,file=dir(1:longueur1)//'Norma_fin',access='append')
      write(50,1100) itime,norma_uxf,norma_uyf,norma_ppf
      close(50)
!
!     Determina y guarda coef. aerod. para la trayectoria Inicial en dominio COE
      call aerof_inflow1(uxitemp,uyitemp,ppitemp,uxm1i,uym1i,uxm2i,uym2i,&
             ffxp,fsxp,fwxp,ffyp,fsyp,fwyp,&
             uxi_ref(itime),uyi_ref(itime),ppi_ref(itime),xDrag,yLift,nx,ny,nz,&
             i_coe1,i_coe2,j_coe1,j_coe2,itime,1)
!
      open(60,file=dir(1:longueur1)//'aerof_ini',access='append')
      write(60,1200) itime,xDrag,yLift
      close(60)
 1200 format(I5,2F16.10) 
!
!     Determina y guarda coef. aerod. para la trayectoria Final en dominio COE
      call aerof_inflow1(uxftemp,uyftemp,ppftemp,uxm1f,uym1f,uxm2f,uym2f,&
             ffxp,fsxp,fwxp,ffyp,fsyp,fwyp,&
             uxf_ref(itime),uyf_ref(itime),ppf_ref(itime),xDrag,yLift,nx,ny,nz,&
             i_coe1,i_coe2,j_coe1,j_coe2,itime,2)
!
      open(70,file=dir(1:longueur1)//'aerof_fin',access='append')
      write(70,1200) itime,xDrag,yLift
      close(70)
!
      enddo
!
!     Entrada: uxitemp, uyitemp. Salida: wzitemp
      call vortz(uxitemp,uyitemp,wzitemp,ffxp,fsxp,fwxp,ffyp,fsyp,fwyp,nx,ny,nz)
!
!     Entrada: uxftemp, uyftemp. Salida: wzftemp
      call vortz(uxftemp,uyftemp,wzftemp,ffxp,fsxp,fwxp,ffyp,fsyp,fwyp,nx,ny,nz)
!
!     ******* Guarda campos trayectoria Inicial en dominio COE *******
      iobi=iobi+1
      call numcar (iobi,suffix)
!
!     Velocidad      
      ncampo='uxini'
      longueur=index(ncampo,' ')-1
      nfichier=ncampo(1:longueur)//suffix
      longueur=index(nfichier,' ')-1
      open(100,file=dir(1:longueur1)//nfichier(1:longueur),form='formatted')
      do j=j_coe1,j_coe2
      do i=i_coe1,i_coe2
         write(100,*) uxitemp(i,j,1)
      enddo
      enddo 
      close(100)
      ncampo='uyini'
      longueur=index(ncampo,' ')-1
      nfichier=ncampo(1:longueur)//suffix
      longueur=index(nfichier,' ')-1
      open(100,file=dir(1:longueur1)//nfichier(1:longueur),form='formatted')
      do j=j_coe1,j_coe2
      do i=i_coe1,i_coe2
         write(100,*) uyitemp(i,j,1)
      enddo
      enddo 
      close(100)
!
!     Presion
      ncampo='ppini'
      longueur=index(ncampo,' ')-1
      nfichier=ncampo(1:longueur)//suffix
      longueur=index(nfichier,' ')-1
      open(100,file=dir(1:longueur1)//nfichier(1:longueur),form='formatted')
      do j=j_coe1,j_coe2
      do i=i_coe1,i_coe2
         write(100,*) ppitemp(i,j,1)
      enddo
      enddo 
      close(100)
!      
!     ******* Guarda campos trayectoria Asimilada en dominio COE *******
      iobf=iobf+1
      call numcar (iobf,suffix)
!
!     Velocidad      
      ncampo='uxfin'
      longueur=index(ncampo,' ')-1
      nfichier=ncampo(1:longueur)//suffix
      longueur=index(nfichier,' ')-1
      open(100,file=dir(1:longueur1)//nfichier(1:longueur),form='formatted')
      do j=j_coe1,j_coe2
      do i=i_coe1,i_coe2
         write(100,*) uxftemp(i,j,1)
      enddo
      enddo 
      close(100)
      ncampo='uyfin'
      longueur=index(ncampo,' ')-1
      nfichier=ncampo(1:longueur)//suffix
      longueur=index(nfichier,' ')-1
      open(100,file=dir(1:longueur1)//nfichier(1:longueur),form='formatted')
      do j=j_coe1,j_coe2
      do i=i_coe1,i_coe2
         write(100,*) uyftemp(i,j,1)
      enddo
      enddo 
      close(100)
!
!     Presion
      ncampo='ppfin'
      longueur=index(ncampo,' ')-1
      nfichier=ncampo(1:longueur)//suffix
      longueur=index(nfichier,' ')-1
      open(100,file=dir(1:longueur1)//nfichier(1:longueur),form='formatted')
      do j=j_coe1,j_coe2
      do i=i_coe1,i_coe2
         write(100,*) ppftemp(i,j,1)
      enddo
      enddo 
      close(100)
!
   enddo
!
!
end PROGRAM trayectoria7
!
!********************************************************************


