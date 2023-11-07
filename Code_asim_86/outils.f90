!******************************************************************
!
subroutine maxmin(ux,uy,uz,nx,ny,nz)
!
!******************************************************************
!
   implicit none
!
   integer :: k, j, i, imax, jmax, kmax, imin, jmin, kmin 
   real(8) :: uxmax, uymax, uzmax, uxmin, uymin, uzmin 
!
   integer , intent(in) :: nx 
   integer , intent(in) :: ny 
   integer , intent(in) :: nz   
   real(8),dimension(nx,ny,nz) :: ux,uy,uz
!
   uxmax=ux(1,1,1)
   uymax=uy(1,1,1)
   uzmax=uz(1,1,1)
   uxmin=ux(1,1,1)
   uymin=uy(1,1,1)
   uzmin=uz(1,1,1)
!
   do k=1,nz
   do j=1,ny
   do i=1,nx
      if (ux(i,j,k).gt.uxmax) then
         uxmax=ux(i,j,k)
         imax=i
         jmax=j
         kmax=k
      endif
      if (uy(i,j,k).gt.uymax) uymax=uy(i,j,k)
      if (uz(i,j,k).gt.uzmax) uzmax=uz(i,j,k)
      if (ux(i,j,k).lt.uxmin) then
         uxmin=ux(i,j,k)
         imin=i
         jmin=j
         kmin=k
      endif
      if (uy(i,j,k).lt.uymin) uymin=uy(i,j,k)
      if (uz(i,j,k).lt.uzmin) uzmin=uz(i,j,k)
   enddo
   enddo
   enddo
!
   print *,'Ux,Uy,Uz max',uxmax,uymax,uzmax
   print *,'i ,j ,k umax',imax,jmax,kmax
   print *,'Ux,Uy,Uz min',uxmin,uymin,uzmin
   print *,'i ,j ,k umin',imin,jmin,kmin
   
   return
end subroutine maxmin
!
!
subroutine annule(a,n) 
!...Translated by PSUITE Trans90                  4.3ZH 15:01:47   1/27/ 4  
!...Switches: -nbejknoc -rl -xl   
  implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  integer , intent(in) :: n 
  real(8) , intent(out) :: a(n) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer :: i 
!-----------------------------------------------
!
!******************************************************************
!
   a=0. 
!
   return  
end subroutine annule
!
!
real(8) function alea (r) 
!...Translated by PSUITE Trans90                  4.3ZH 15:01:47   1/27/ 4  
!...Switches: -nbejknoc -rl -xl   
  implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  real(8) , intent(in) :: r 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer :: ia1, ia0, ia1ma0, ic, ix1, ix0, iy0, iy1 
!-----------------------------------------------
!
!******************************************************************
!
   data ia1,ia0,ia1ma0/1536,1029,507/ 
   data ic/1731/ 
   data ix1,ix0/0,0/ 
!
!      PRINT *,'alea !!!'
!
   if (r >= 0.) then 
      if (r > 0.) go to 20 
!
!     A*X = 2**22*IA1*IX1 + 2**11*(IA1*IX1 + (IA1-IA0)*(IX0-IX1)
!         + IA0*IX0) + IA0*IX0
!
      iy0=ia0*ix0 
      iy1=ia1*ix1 + ia1ma0*(ix0 - ix1) + iy0 
      iy0=iy0 + ic 
      ix0=mod(iy0,2048) 
      iy1=iy1 + (iy0 - ix0)/2048 
      ix1=mod(iy1,2048) 
!
   endif
10 continue 
   alea=ix1*2048 + ix0 
   alea=alea/4194304. 
   return  
!
20 continue 
!   ix1=amod(r,1.)*4194304. + 0.5 
   ix0=mod(ix1,2048) 
   ix1=(ix1 - ix0)/2048 
   go to 10 
!
end function alea
!
!******************************************************************
!
subroutine numcar (num,car)
!
!******************************************************************!
!
   character(len=4) car
!
   if (num.ge.1000) then
      write (car,1) num
1     format (i4)
   else 
      if (num.ge.100) then
         write (car,2) num
2        format ('0',i3)
      else
        if (num.ge.10) then
         write (car,3) num
3        format ('00',i2)
         else
         write (car,4) num
4        format ('000',i1)
         endif 
      endif
   endif
!
   return
 end subroutine numcar
!
!********************************************************************
!
real(8) function epmach (dum)
!
!********************************************************************
!
   USE value_m
!
   implicit none
!
   real(8) :: dum,eps
   eps=1.
   eps=eps/10. 
   call store (eps + 1.)
   do while(v - 1. > 0.) 
      eps=eps/10. 
      call store (eps + 1.) 
   enddo
   epmach=100.*eps
   return
end function epmach
!
!********************************************************************
!
subroutine store(x) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
   USE value_m 
!...Translated by PSUITE Trans90                  4.3ZH 15:01:47   1/27/ 4  
!...Switches: -nbejknoc -rl -xl   
   implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
   real , intent(in) :: x 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!-----------------------------------------------
!
!********************************************************************
!
   v=x 
!
   return  
end subroutine store
!
!******************************************************************

subroutine maxmin1 (ux,nx,ny,nz)

!******************************************************************
!
   implicit none
!
   integer :: k, j, i, imax, jmax, kmax, imin, jmin, kmin 
   real(8) :: uxmax, uymax, uzmax, uxmin, uymin, uzmin 
!
   integer , intent(in) :: nx 
   integer , intent(in) :: ny 
   integer , intent(in) :: nz 
   real(8),dimension(nx,ny,nz) :: ux

   uxmax=ux(1,1,1)
   uxmin=ux(1,1,1)

   do k=1,nz
   do j=1,ny
   do i=1,nx
      if (ux(i,j,k).gt.uxmax) uxmax=ux(i,j,k)
      if (ux(i,j,k).lt.uxmin) uxmin=ux(i,j,k)
   enddo
   enddo
   enddo
   
   print *,'Uxmax,Uxmin',uxmax,uxmin
   
   return
end subroutine maxmin1
!****************************************************
!
real(8) function rand2 (idum) 
!
! ****************************************************
!
!
! =========================================================
!
!  "Long" period ( > 2.E+18 ) random number generator
!  of L'Ecuyer with Bays-Durham shuffle and added safeguards.
!  Returns a uniform random deviate between 0. and 1.
!  (exclusive of the endpoint values).
!
!  Call with IDUM, a negative integer to initialize;
!  thereafter, do not alter IDUM between succesive deviates
!  in a sequence.  RNMAX should approximate the largest
!  floating value that  is less than 1.
!
!  REF: NUMERICAL RECEPIES -  pp 271
!
! =========================================================
!
!...Translated by PSUITE Trans90                  4.3ZH 15:01:47   1/27/ 4  
!...Switches: -nbejknoc -rl -xl   
  implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  integer, intent(inout) :: idum 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
  integer, parameter :: im1 = 2147483563 
  integer, parameter :: im2 = 2147483399 
  integer, parameter :: imm1 = im1 - 1 
  integer, parameter :: ia1 = 40014 
  integer, parameter :: ia2 = 40692 
  integer, parameter :: iq1 = 53668 
  integer, parameter :: iq2 = 52774 
  integer, parameter :: ir1 = 12211 
  integer, parameter :: ir2 = 3791 
  integer, parameter :: ntab = 32 
  integer, parameter :: ndiv = 1 + imm1/ntab 
  real(8), parameter :: am = 1./im1 
  real(8), parameter :: eps = 1.2E-14 
  real(8), parameter :: rnmx = 1. - eps 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer :: idum2, j, k 
  integer , dimension(ntab) :: iv 
  integer :: iy 
!
  save idum2, iv, iy 
!-----------------------------------------------
!
!
!
   data idum2/123456789/ 
   data iv/ntab*0/ 
   data iy/0/ 
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
   if (idum <= 0) then                        !Initialize. 
      idum=max((-idum),1) 
      idum2=idum 
      do j=ntab + 8,1,-1 
         k=idum/iq1 
         idum=ia1*(idum - k*iq1) - k*ir1 
         if (idum < 0) idum=idum + im1 
         if (j > ntab) cycle  
         iv(j)=idum 
      enddo
   endif
!
   k=idum/iq1 
   idum=ia1*(idum - k*iq1) - k*ir1 
   if (idum < 0) idum=idum + im1 
!
   k=idum2/iq2 
   idum2=ia2*(idum2 - k*iq2) - k*ir2 
   if (idum2 < 0) idum2=idum2 + im2 
!
   j=1 + iy/ndiv 
   iy=iv(j) - idum2 
   iv(j)=idum 
   if (iy < 1) iy=iy + imm1 
!===========================
   rand2=min(am*iy,rnmx) 
!===========================
!
   return  
end function rand2
!
!******************************************************************

subroutine gasdev_s(harvest)

!******************************************************************
!
! =========================================================
!
!  Returns in harvest a normally distributed deviate with zero mean and unit variance, 
!  using rand2 as the source of uniform deviates.
!
!  REF: NUMERICAL RECIPES - pp 1152
!
! =========================================================
!
   USE aleatoire_m
!
   implicit none
!
   real(8),intent(out) :: harvest
   real(8) :: rsq,v1,v2
   real(8),save :: g
   logical,save :: gaus_stored=.false.
   real(8),external :: rand2
!                                               
   if (gaus_stored) then
      harvest=g
      gaus_stored=.false.
   else
      do
        v1=rand2(idum)
        v2=rand2(idum)
        v1=2.0*v1-1.0
        v2=2.0*v2-1.0
        rsq=v1**2+v2**2
        if (rsq > 0.0 .and. rsq < 1.0) exit
      enddo
      rsq=sqrt(-2.0*log(rsq)/rsq)  ! Now make the Box-Muller transformation
      harvest=v1*rsq
      g=v2*rsq
      gaus_stored=.true.
   end if
! 
    return  
end subroutine gasdev_s
!
!******************************************************************

subroutine barron_error (error,uxtest,uytest,uxtrue,uytrue,fenetre,nx,ny,nz)

!******************************************************************
!
   implicit none
!
   real(8),dimension(fenetre,nx,ny,nz) :: uxtest,uytest,uxtrue,uytrue
   real(8),dimension(nx,ny,nz) :: cs,ns1,ns2,ns,acs,acs1
   real(8),dimension(fenetre) :: error1
   real(8) :: pi,error,valor
   integer :: f,j,i,fenetre,nx,ny,nz,nxy
!
   pi=acos(-1.)
   nxy=nx*ny
!
   do f=1,fenetre
!
   do j=1,ny
   do i=1,nx
      cs(i,j,1) = 1 + uxtest(f,i,j,1)*uxtrue(f,i,j,1) + uytest(f,i,j,1)*uytrue(f,i,j,1)
      ns1(i,j,1) = sqrt(1 + uxtest(f,i,j,1)*uxtest(f,i,j,1) + uytest(f,i,j,1)*uytest(f,i,j,1))
      ns2(i,j,1) = sqrt(1 + uxtrue(f,i,j,1)*uxtrue(f,i,j,1) + uytrue(f,i,j,1)*uytrue(f,i,j,1))
      ns(i,j,1) = ns1(i,j,1)*ns2(i,j,1)
      acs(i,j,1) = cs(i,j,1)/ns(i,j,1)
   enddo
   enddo 
!
   valor=maxval(dabs(acs))
   do j=1,ny
   do i=1,nx
      acs1(i,j,1) = acos(acs(i,j,1)/valor)*180/pi
   enddo
   enddo 
!
!  promedio espacial del error angular (en grados)
   error1(f)=sum(acs1)/nxy 
!
   enddo 
!
!  promedio temporal de errorp(f)
   error=sum(error1)/fenetre 
! 
   return  
end subroutine barron_error
!
!******************************************************************

subroutine barron (error,uxtest,uytest,uxtrue,uytrue,fenetre,fenetrem,nx,ny,nz)

!******************************************************************
!
   implicit none
!
   real(8),dimension(fenetre,nx,ny,nz) :: uxtrue,uytrue
   real(8),dimension(fenetrem,nx,ny,nz) :: uxtest,uytest
   real(8),dimension(nx,ny,nz) :: acsm
   real(8),dimension(fenetrem) :: errorv
   real(8) :: cs,argu1,argu2,ns1,ns2,ns,acs,pi,error,valor
   integer :: f,j,i,fenetre,fenetrem,nx,ny,nz,nxy
!
   pi=acos(-1.)
   nxy=nx*ny
!
   do f=1,fenetrem
!
   do j=1,ny
   do i=1,nx
      cs = 1 + uxtest(f,i,j,1)*uxtrue(f+1,i,j,1) + uytest(f,i,j,1)*uytrue(f+1,i,j,1)
      argu1 = 1 + uxtest(f,i,j,1)*uxtest(f,i,j,1) + uytest(f,i,j,1)*uytest(f,i,j,1)
      argu2 = 1 + uxtrue(f+1,i,j,1)*uxtrue(f+1,i,j,1) + uytrue(f+1,i,j,1)*uytrue(f+1,i,j,1)
      ns1 = sqrt(argu1)
      ns2 = sqrt(argu2)
      ns = ns1*ns2
      acs = cs/ns
      acsm(i,j,1) = acos(acs)*180/pi
   enddo
   enddo 
!
!  promedio espacial del error angular (en grados)
   errorv(f)=sum(acsm)/nxy 
!
   enddo 
!
!  promedio temporal de errorv(f)
   error=sum(errorv)/fenetrem 
! 
   return  
end subroutine barron
!
!******************************************************************

subroutine snorm (error,wz,wz_obsg,ux0r,uy0r,ux0r_obs,uy0r_obs,gx0r,gy0r,&
             gx0r_obs,gy0r_obs,bx1r,by1r,bx1_obs,by1_obs,&
             fenetre,fenetrem,nx,ny,nz,nxr,nyr,nxrg,nyrg,jyt,rdxy,ditime,&
             cov_obs,cov_iniu,cov_inig,cov_suav,cov_tay,cov_gau,&
             jw,iw,w_xy,sumw_xy,n_nodos,n_nodos_max,&
             jwq,iwq,w_xyq,sumw_xyq,n_nodos_q,n_nodos_q_max,&
             itw,w_t,sumw_t,n_nodos_t,n_nodos_t_max,&
             i_redg1,i_redg2,j_redg1,j_redg2,nxrg_asim,nyrg_asim,&
             error_obs,error_iniu,error_inig,error_suav,error_tay,error_gau)

!******************************************************************
!
   implicit none
!
   real(8),dimension(fenetre,nxrg,nyrg,nz) :: wz_obsg
   real(8),dimension(fenetrem,nxrg,nyrg,nz) :: wz_av
   real(8),dimension(fenetrem,nx,ny,nz) :: wz
   real(8),dimension(nxr,nyr,nz) :: ux0r,uy0r,gx0r,gy0r
   real(8),dimension(nxr,nyr,nz) :: ux0r_av,uy0r_av
   real(8),dimension(nxr,nyr,nz) :: ux0r_obs,uy0r_obs,gx0r_obs,gy0r_obs
   real(8),dimension(nyr,ditime) :: bx1r,by1r,bx1r_av,by1r_av
   real(8),dimension(nyr,ditime) :: bx1_obs,by1_obs
   real(8),dimension(nxrg,nyrg) :: errorm
   real(8),dimension(nxr,nyr) :: errorlu,errorlg,errormq
   real(8),dimension(nyr) :: errorp
   real(8),dimension(ditime) :: errorq
   real(8),dimension(fenetrem) :: errorv
   real(8),dimension(nyr) :: errorn
   real(8),dimension(ditime) :: errorw
   integer,dimension(nxrg,nyrg,n_nodos_max) :: iw,jw
   real(8),dimension(nxrg,nyrg,n_nodos_max) :: w_xy
   integer,dimension(nxrg,nyrg) :: n_nodos
   real(8),dimension(nxrg,nyrg) :: sumw_xy
   real(8),dimension(n_nodos_max) :: wwz
   integer,dimension(nxr,nyr,n_nodos_q_max) :: iwq,jwq
   real(8),dimension(nxr,nyr,n_nodos_q_max) :: w_xyq
   integer,dimension(nxr,nyr) :: n_nodos_q
   real(8),dimension(nxr,nyr) :: sumw_xyq
   real(8),dimension(n_nodos_q_max) :: wux,wuy
   integer,dimension(ditime,n_nodos_t_max) :: itw
   real(8),dimension(ditime,n_nodos_t_max) :: w_t
   integer,dimension(ditime) :: n_nodos_t
   real(8),dimension(ditime) :: sumw_t
   real(8),dimension(n_nodos_t_max) :: wbx1r,wby1r
   real(8) :: deltwz_obs
   real(8) :: deltux_ini,deltuy_ini,deltgx_ini,deltgy_ini
   real(8) :: deltux_suav,deltuy_suav
   real(8) :: deltux_tay,deltuy_tay,deltux_gau,deltuy_gau
   real(8) :: error,error_obs,error_iniu,error_inig,error_suav,error_tay,error_gau
   real(8),dimension(nxrg,nyrg) :: cov_obs
   real(8) :: cov_iniu,cov_inig,cov_suav,cov_tay,cov_gau
   real(8) :: w_r,sumw_r,w_q,sumw_q,w_d,sumw_d
   integer :: f,fenetre,fenetrem,nx,ny,nz
   integer :: itime,n_nodos_max,n_nodos_q_max,n_nodos_t_max,nodo
   integer :: j_redg,i_redg,nxrg,nyrg,nxyrg,jyt,rdxy,ditime
   integer :: j_red,i_red,nxr,nyr,nxyr,j_rec,i_rec,j,i
   integer :: j_redg1,j_redg2,i_redg1,i_redg2,nxrg_asim,nyrg_asim,nxyrg_asim
!
   nxyr = nxr*nyr
   nxyrg = nxrg*nyrg
   nxyrg_asim = nxrg_asim*nyrg_asim
!
!  ********************* error_obs *********************
   do f=1,fenetrem
!
   errorm=0.
   do j_redg = j_redg1,j_redg2
   do i_redg = i_redg1,i_redg2
!
!     /w_r(nodo)*wz_r(nodo)
      wwz=0.
      do nodo=1,n_nodos(i_redg,j_redg)
         i=iw(i_redg,j_redg,nodo)
         j=jw(i_redg,j_redg,nodo)
         w_r=w_xy(i_redg,j_redg,nodo)
         wwz(nodo) = w_r*wz(f,i,j,1)
      enddo
!
!     Weighted average
      sumw_r=sumw_xy(i_redg,j_redg)
      wz_av(f,i_redg,j_redg,1) = sum(wwz)/sumw_r
!
      deltwz_obs = wz_av(f,i_redg,j_redg,1)-wz_obsg(f+1,i_redg,j_redg,1)
      errorm(i_redg,j_redg) = cov_obs(i_redg,j_redg)*deltwz_obs**2
   enddo
   enddo 
!
!  promedio espacial de errorm
   errorv(f)=sum(errorm)/nxyrg_asim 
!
   enddo 
!
!  promedio temporal de errorv
   error_obs=sum(errorv)/fenetrem 
! 
!  ********************* error_iniu *********************
   do j = 1,nyr
   do i = 1,nxr
      deltux_ini = ux0r(i,j,1)-ux0r_obs(i,j,1)
      deltuy_ini = uy0r(i,j,1)-uy0r_obs(i,j,1)
      errorlu(i,j) = cov_iniu*(deltux_ini**2 + deltuy_ini**2)
   enddo
   enddo 
!
!  promedio espacial de errorlu
   error_iniu=sum(errorlu)/nxyr 
!
!  ********************* error_inig *********************
   do j = 1,nyr
   do i = 1,nxr
      deltgx_ini = gx0r(i,j,1)-gx0r_obs(i,j,1)
      deltgy_ini = gy0r(i,j,1)-gy0r_obs(i,j,1)
      errorlg(i,j) = cov_inig*(deltgx_ini**2 + deltgy_ini**2)
   enddo
   enddo 
!
!  promedio espacial de errorlg
   error_inig=sum(errorlg)/nxyr 
!
!  ********************* error_suav *********************
   do j_red = 1,nyr
   do i_red = 1,nxr
!
!     /w_r(nodo)*u_r(nodo)
      wux=0. ; wuy=0.
      do nodo=1,n_nodos_q(i_red,j_red)
         i=iwq(i_red,j_red,nodo)
         j=jwq(i_red,j_red,nodo)
         w_q=w_xyq(i_red,j_red,nodo)
         wux(nodo) = w_q*ux0r(i,j,1)
         wuy(nodo) = w_q*uy0r(i,j,1)
      enddo
!
!     Weighted average
      sumw_q=sumw_xyq(i_red,j_red)
      ux0r_av(i_red,j_red,1) = sum(wux)/sumw_q
      uy0r_av(i_red,j_red,1) = sum(wuy)/sumw_q
!
      deltux_suav = ux0r(i_red,j_red,1)-ux0r_av(i_red,j_red,1)
      deltuy_suav = uy0r(i_red,j_red,1)-uy0r_av(i_red,j_red,1)
      errormq(i_red,j_red) = cov_suav*(deltux_suav**2 + deltuy_suav**2)
   enddo
   enddo 
!
!  promedio espacial de errormq
   error_suav=sum(errormq)/nxyr 
! 
!  ********************* error_tay *********************
   do itime=1,ditime
!
   do j = 1,nyr
      deltux_tay = bx1r(j,itime)-bx1_obs(j,itime)
      deltuy_tay = by1r(j,itime)-by1_obs(j,itime)
      errorn(j) = cov_tay*(deltux_tay**2 + deltuy_tay**2)
   enddo
!
!  promedio espacial de errorn
   errorw(itime)=sum(errorn)/nyr 
!
   enddo 
!
!  promedio temporal de errorw
   error_tay=sum(errorw)/ditime 
!
!  ********************* error_gau *********************
   do itime=1,ditime
!
   do j=1,nyr
!
!     /w_t(nodo)*b1r_t(nodo)
      wbx1r=0.
      wby1r=0.
      do nodo=1,n_nodos_t(itime)
         i=itw(itime,nodo)
         w_d=w_t(itime,nodo)
         wbx1r(nodo) = w_d*bx1r(j,i)
         wby1r(nodo) = w_d*by1r(j,i)
      enddo
!
!     Weighted average
      sumw_d=sumw_t(itime)
      bx1r_av(j,itime) = sum(wbx1r)/sumw_d
      by1r_av(j,itime) = sum(wby1r)/sumw_d
!
      deltux_gau = bx1r(j,itime)-bx1r_av(j,itime)
      deltuy_gau = by1r(j,itime)-by1r_av(j,itime)
      errorp(j) = cov_gau*(deltux_gau**2 + deltuy_gau**2)
   enddo
!
!  promedio espacial de errorp
   errorq(itime)=sum(errorp)/nyr 
!
   enddo 
!
!  promedio temporal de errorq
   error_gau=sum(errorq)/ditime 
!
!  *********************** error ***********************
   error = error_obs + error_iniu + error_inig + error_suav + error_tay + error_gau
!
   return  
end subroutine snorm
!
!******************************************************************

subroutine reconstruye (ux0r,uy0r,bx1r,by1r,ux0,uy0,bx1,by1,&
             jyc,jyt,nxr,nyr,nx,ny,nz,ditime)

!******************************************************************
!
!  reconstruye obs. y cond. de entrada reducidas para satisfacer cond. borde
!  lateral requerida por la dns
!
   implicit none
!
   real(8),dimension(nxr,nyr,nz) :: ux0r,uy0r
   real(8),dimension(nx,ny,nz) :: ux0,uy0
   real(8),dimension(nyr,ditime) :: bx1r,by1r
   real(8),dimension(ny,ditime) :: bx1,by1
   real(8) :: pi,umc,RAMP_base,RAMP_exp,RAMP
   integer :: j0,jred,deltaj,pre
   integer :: jyc,jyt,nxr,nyr,nx,ny,nz,ditime
   integer :: itime,i,j
!
   pi = acos(-1.)
   umc=1. ! veloc. de conveccion
   deltaj=jyt-jyc
   pre=9999 ! parametro de la funcion RAMP
!
!  zona veloc. constante inferior
   do j=1,jyc
      do i=1,nx
         ux0(i,j,1)=umc
         uy0(i,j,1)=0.
      enddo
      do itime=1,ditime
         bx1(j,itime)=umc
         by1(j,itime)=0.
      enddo
   enddo
!  
!  zona veloc. de transicion inferior
   do j=jyc+1,jyt
      j0 = j-jyc
      RAMP_base = 0.5*sin((pi*j0/deltaj)-(pi/2))+0.5
      RAMP_exp = 1/(((j0/deltaj)**pre)+1)
      RAMP = RAMP_base**RAMP_exp 
      do i=1,nx
         ux0(i,j,1)=umc-(umc-ux0r(i,1,1))*RAMP
         uy0(i,j,1)=uy0r(i,1,1)*RAMP
      enddo
      do itime=1,ditime
         bx1(j,itime)=umc-(umc-bx1r(1,itime))*RAMP
         by1(j,itime)=by1r(1,itime)*RAMP
      enddo
   enddo
!    
!  zona obs. reducida
   jred=1
   do j=jyt+1,jyt+nyr-2
      jred=jred+1 ! jred va de 2 a nyr-1
      do i=1,nx
         ux0(i,j,1)=ux0r(i,jred,1)
         uy0(i,j,1)=uy0r(i,jred,1)
      enddo
      do itime=1,ditime
         bx1(j,itime)=bx1r(jred,itime)
         by1(j,itime)=by1r(jred,itime)
      enddo
   enddo
!    
!  zona veloc. de transicion superior
   do j=jyt+nyr-1,2*jyt+nyr-jyc-1
      j0 = j-(jyt+nyr-2)
      RAMP_base = 0.5*sin((pi*j0/deltaj)-(pi/2))+0.5
      RAMP_exp = 1/(((j0/deltaj)**pre)+1)
      RAMP = RAMP_base**RAMP_exp 
      do i=1,nx
         ux0(i,j,1)=ux0r(i,nyr,1)-(ux0r(i,nyr,1)-umc)*RAMP
         uy0(i,j,1)=uy0r(i,nyr,1)*(1-RAMP)
      enddo
      do itime=1,ditime
         bx1(j,itime)=bx1r(nyr,itime)-(bx1r(nyr,itime)-umc)*RAMP
         by1(j,itime)=by1r(nyr,itime)*(1-RAMP)
      enddo
   enddo
!    
!  zona veloc. constante superior
   do j=2*jyt+nyr-jyc,2*jyt+nyr
      do i=1,nx
         ux0(i,j,1)=umc
         uy0(i,j,1)=0.
      enddo
      do itime=1,ditime
         bx1(j,itime)=umc
         by1(j,itime)=0.
      enddo
   enddo
! 
   return  
end subroutine reconstruye
!
!******************************************************************

subroutine reconstruye1 (uy0r,uy0,&
             jyc,jyt,nxr,nyr,nx,ny,nz)

!******************************************************************
!
!  reconstruye campo de obs. reducida para satisfacer cond. borde
!  lateral requerida por la dns
!
   implicit none
!
   real(8),dimension(nxr,nyr,nz) :: uy0r
   real(8),dimension(nx,ny,nz) :: uy0
   real(8) :: pi,RAMP_base,RAMP_exp,RAMP
   integer :: j0,jred,deltaj,pre
   integer :: jyc,jyt,nxr,nyr,nx,ny,nz
   integer :: i,j
!
   pi = acos(-1.)
   deltaj=jyt-jyc
   pre=9999 ! parametro de la funcion RAMP
!
!  zona veloc. constante inferior
   do j=1,jyc
      do i=1,nx
         uy0(i,j,1)=0.
      enddo
   enddo
!  
!  zona veloc. de transicion inferior
   do j=jyc+1,jyt
      j0 = j-jyc
      RAMP_base = 0.5*sin((pi*j0/deltaj)-(pi/2))+0.5
      RAMP_exp = 1/(((j0/deltaj)**pre)+1)
      RAMP = RAMP_base**RAMP_exp 
      do i=1,nx
         uy0(i,j,1)=uy0r(i,1,1)*RAMP
      enddo
   enddo
!    
!  zona obs. reducida
   jred=1
   do j=jyt+1,jyt+nyr-2
      jred=jred+1 ! jred va de 2 a nyr-1
      do i=1,nx
         uy0(i,j,1)=uy0r(i,jred,1)
      enddo
   enddo
!    
!  zona veloc. de transicion superior
   do j=jyt+nyr-1,2*jyt+nyr-jyc-1
      j0 = j-(jyt+nyr-2)
      RAMP_base = 0.5*sin((pi*j0/deltaj)-(pi/2))+0.5
      RAMP_exp = 1/(((j0/deltaj)**pre)+1)
      RAMP = RAMP_base**RAMP_exp 
      do i=1,nx
         uy0(i,j,1)=uy0r(i,nyr,1)*(1-RAMP)
      enddo
   enddo
!    
!  zona veloc. constante superior
   do j=2*jyt+nyr-jyc,2*jyt+nyr
      do i=1,nx
         uy0(i,j,1)=0.
      enddo
   enddo
! 
   return  
end subroutine reconstruye1
!
!******************************************************************

subroutine reconstruye2 (ux0r,uy0r,ux0,uy0,&
             jyc,jyt,nxr,nyr,nx,ny,nz)

!******************************************************************
!
!  reconstruye obs. reducida para satisfacer cond. borde
!  lateral requerida por la dns
!
   implicit none
!
   real(8),dimension(nxr,nyr,nz) :: ux0r,uy0r
   real(8),dimension(nx,ny,nz) :: ux0,uy0
   real(8) :: pi,umc,RAMP_base,RAMP_exp,RAMP
   integer :: j0,jred,deltaj,pre
   integer :: jyc,jyt,nxr,nyr,nx,ny,nz
   integer :: i,j
!
   pi = acos(-1.)
   umc=1. ! veloc. de conveccion
   deltaj=jyt-jyc
   pre=9999 ! parametro de la funcion RAMP
!
!  zona veloc. constante inferior
   do j=1,jyc
      do i=1,nx
         ux0(i,j,1)=umc
         uy0(i,j,1)=0.
      enddo
   enddo
!  
!  zona veloc. de transicion inferior
   do j=jyc+1,jyt
      j0 = j-jyc
      RAMP_base = 0.5*sin((pi*j0/deltaj)-(pi/2))+0.5
      RAMP_exp = 1/(((j0/deltaj)**pre)+1)
      RAMP = RAMP_base**RAMP_exp 
      do i=1,nx
         ux0(i,j,1)=umc-(umc-ux0r(i,1,1))*RAMP
         uy0(i,j,1)=uy0r(i,1,1)*RAMP
      enddo
   enddo
!    
!  zona obs. reducida
   jred=1
   do j=jyt+1,jyt+nyr-2
      jred=jred+1 ! jred va de 2 a nyr-1
      do i=1,nx
         ux0(i,j,1)=ux0r(i,jred,1)
         uy0(i,j,1)=uy0r(i,jred,1)
      enddo
   enddo
!    
!  zona veloc. de transicion superior
   do j=jyt+nyr-1,2*jyt+nyr-jyc-1
      j0 = j-(jyt+nyr-2)
      RAMP_base = 0.5*sin((pi*j0/deltaj)-(pi/2))+0.5
      RAMP_exp = 1/(((j0/deltaj)**pre)+1)
      RAMP = RAMP_base**RAMP_exp 
      do i=1,nx
         ux0(i,j,1)=ux0r(i,nyr,1)-(ux0r(i,nyr,1)-umc)*RAMP
         uy0(i,j,1)=uy0r(i,nyr,1)*(1-RAMP)
      enddo
   enddo
!    
!  zona veloc. constante superior
   do j=2*jyt+nyr-jyc,2*jyt+nyr
      do i=1,nx
         ux0(i,j,1)=umc
         uy0(i,j,1)=0.
      enddo
   enddo
! 
   return  
end subroutine reconstruye2
!
!******************************************************************

subroutine reconstruye3 (bx1r,by1r,bx1,by1,&
             jyc,jyt,nyr,ny,ditime)

!******************************************************************
!
!  reconstruye cond. de entrada reducida para satisfacer cond. borde
!  lateral requerida por la dns
!
   implicit none
!
   real(8),dimension(nyr,ditime) :: bx1r,by1r
   real(8),dimension(ny,ditime) :: bx1,by1
   real(8) :: pi,umc,RAMP_base,RAMP_exp,RAMP
   integer :: j0,jred,deltaj,pre
   integer :: jyc,jyt,nyr,ny,ditime
   integer :: itime,j
!
   pi = acos(-1.)
   umc=0. ! veloc. de conveccion
   deltaj=jyt-jyc
   pre=9999 ! parametro de la funcion RAMP
!
!  zona veloc. constante inferior
   do j=1,jyc
      do itime=1,ditime
         bx1(j,itime)=umc
         by1(j,itime)=0.
      enddo
   enddo
!  
!  zona veloc. de transicion inferior
   do j=jyc+1,jyt
      j0 = j-jyc
      RAMP_base = 0.5*sin((pi*j0/deltaj)-(pi/2))+0.5
      RAMP_exp = 1/(((j0/deltaj)**pre)+1)
      RAMP = RAMP_base**RAMP_exp 
      do itime=1,ditime
         bx1(j,itime)=umc-(umc-bx1r(1,itime))*RAMP
         by1(j,itime)=by1r(1,itime)*RAMP
      enddo
   enddo
!    
!  zona obs. reducida
   jred=1
   do j=jyt+1,jyt+nyr-2
      jred=jred+1 ! jred va de 2 a nyr-1
      do itime=1,ditime
         bx1(j,itime)=bx1r(jred,itime)
         by1(j,itime)=by1r(jred,itime)
      enddo
   enddo
!    
!  zona veloc. de transicion superior
   do j=jyt+nyr-1,2*jyt+nyr-jyc-1
      j0 = j-(jyt+nyr-2)
      RAMP_base = 0.5*sin((pi*j0/deltaj)-(pi/2))+0.5
      RAMP_exp = 1/(((j0/deltaj)**pre)+1)
      RAMP = RAMP_base**RAMP_exp 
      do itime=1,ditime
         bx1(j,itime)=bx1r(nyr,itime)-(bx1r(nyr,itime)-umc)*RAMP
         by1(j,itime)=by1r(nyr,itime)*(1-RAMP)
      enddo
   enddo
!    
!  zona veloc. constante superior
   do j=2*jyt+nyr-jyc,2*jyt+nyr
      do itime=1,ditime
         bx1(j,itime)=umc
         by1(j,itime)=0.
      enddo
   enddo
! 
   return  
end subroutine reconstruye3
!
!******************************************************************

subroutine reconstruye4 (ux0r,uy0r,gx0r,gy0r,bx1r,by1r,&
             ux0,uy0,gx0,gy0,bx1,by1,&
             jyc,jyt,nxr,nyr,nx,ny,nz,ditime)

!******************************************************************
!
!  reconstruye obs. y cond. de entrada reducidas para satisfacer cond. borde
!  lateral requerida por la dns
!
   implicit none
!
   real(8),dimension(nxr,nyr,nz) :: ux0r,uy0r,gx0r,gy0r
   real(8),dimension(nx,ny,nz) :: ux0,uy0,gx0,gy0
   real(8),dimension(nyr,ditime) :: bx1r,by1r
   real(8),dimension(ny,ditime) :: bx1,by1
   real(8) :: pi,umc,RAMP_base,RAMP_exp,RAMP
   integer :: j0,jred,deltaj,pre
   integer :: jyc,jyt,nxr,nyr,nx,ny,nz,ditime
   integer :: itime,i,j
!
   pi = acos(-1.)
   umc=1. ! veloc. de conveccion
   deltaj=jyt-jyc
   pre=9999 ! parametro de la funcion RAMP
!
!  zona veloc. constante inferior
   do j=1,jyc
      do i=1,nx
         ux0(i,j,1)=umc
         uy0(i,j,1)=0.
         gx0(i,j,1)=0.
         gy0(i,j,1)=0.
      enddo
      do itime=1,ditime
         bx1(j,itime)=umc
         by1(j,itime)=0.
      enddo
   enddo
!  
!  zona veloc. de transicion inferior
   do j=jyc+1,jyt
      j0 = j-jyc
      RAMP_base = 0.5*sin((pi*j0/deltaj)-(pi/2))+0.5
      RAMP_exp = 1/(((j0/deltaj)**pre)+1)
      RAMP = RAMP_base**RAMP_exp 
      do i=1,nx
         ux0(i,j,1)=umc-(umc-ux0r(i,1,1))*RAMP
         uy0(i,j,1)=uy0r(i,1,1)*RAMP
         gx0(i,j,1)=gx0r(i,1,1)*RAMP
         gy0(i,j,1)=gy0r(i,1,1)*RAMP
      enddo
      do itime=1,ditime
         bx1(j,itime)=umc-(umc-bx1r(1,itime))*RAMP
         by1(j,itime)=by1r(1,itime)*RAMP
      enddo
   enddo
!    
!  zona obs. reducida
   jred=1
   do j=jyt+1,jyt+nyr-2
      jred=jred+1 ! jred va de 2 a nyr-1
      do i=1,nx
         ux0(i,j,1)=ux0r(i,jred,1)
         uy0(i,j,1)=uy0r(i,jred,1)
         gx0(i,j,1)=gx0r(i,jred,1)
         gy0(i,j,1)=gy0r(i,jred,1)
      enddo
      do itime=1,ditime
         bx1(j,itime)=bx1r(jred,itime)
         by1(j,itime)=by1r(jred,itime)
      enddo
   enddo
!    
!  zona veloc. de transicion superior
   do j=jyt+nyr-1,2*jyt+nyr-jyc-1
      j0 = j-(jyt+nyr-2)
      RAMP_base = 0.5*sin((pi*j0/deltaj)-(pi/2))+0.5
      RAMP_exp = 1/(((j0/deltaj)**pre)+1)
      RAMP = RAMP_base**RAMP_exp 
      do i=1,nx
         ux0(i,j,1)=ux0r(i,nyr,1)-(ux0r(i,nyr,1)-umc)*RAMP
         uy0(i,j,1)=uy0r(i,nyr,1)*(1-RAMP)
         gx0(i,j,1)=gx0r(i,nyr,1)*(1-RAMP)
         gy0(i,j,1)=gy0r(i,nyr,1)*(1-RAMP)
      enddo
      do itime=1,ditime
         bx1(j,itime)=bx1r(nyr,itime)-(bx1r(nyr,itime)-umc)*RAMP
         by1(j,itime)=by1r(nyr,itime)*(1-RAMP)
      enddo
   enddo
!    
!  zona veloc. constante superior
   do j=2*jyt+nyr-jyc,2*jyt+nyr
      do i=1,nx
         ux0(i,j,1)=umc
         uy0(i,j,1)=0.
         gx0(i,j,1)=0.
         gy0(i,j,1)=0.
      enddo
      do itime=1,ditime
         bx1(j,itime)=umc
         by1(j,itime)=0.
      enddo
   enddo
! 
   return  
end subroutine reconstruye4
!
!******************************************************************

subroutine reconstruye5 (ux0r,uy0r,gx0r,gy0r,bx1r,by1r,&
             ux0,uy0,gx0,gy0,bx1,by1,&
             jyt,nxr,nyr,nx,ny,nz,ditime)

!******************************************************************
!
!  reconstruye obs. y cond. de entrada reducidas para satisfacer cond. borde
!  lateral requerida por la dns
!
   implicit none
!
   real(8),dimension(nxr,nyr,nz) :: ux0r,uy0r,gx0r,gy0r
   real(8),dimension(nx,ny,nz) :: ux0,uy0,gx0,gy0
   real(8),dimension(nyr,ditime) :: bx1r,by1r
   real(8),dimension(ny,ditime) :: bx1,by1
   real(8) :: pi,umc
   integer :: j0,jred,deltaj,pre
   integer :: jyt,nxr,nyr,nx,ny,nz,ditime
   integer :: itime,i,j
!
   pi = acos(-1.)
   umc=1. ! veloc. de conveccion
!
!  zona reconstruida inferior
   do j=1,jyt
      do i=1,nx
         ux0(i,j,1)=umc
         uy0(i,j,1)=0.
         gx0(i,j,1)=0.
         gy0(i,j,1)=0.
      enddo
      do itime=1,ditime
         bx1(j,itime)=umc
         by1(j,itime)=0.
      enddo
   enddo
!  
!  zona obs. reducida
   jred=1
   do j=jyt+1,jyt+nyr-2
      jred=jred+1 ! jred va de 2 a nyr-1
      do i=1,nx
         ux0(i,j,1)=ux0r(i,jred,1)
         uy0(i,j,1)=uy0r(i,jred,1)
         gx0(i,j,1)=gx0r(i,jred,1)
         gy0(i,j,1)=gy0r(i,jred,1)
      enddo
      do itime=1,ditime
         bx1(j,itime)=bx1r(jred,itime)
         by1(j,itime)=by1r(jred,itime)
      enddo
   enddo
!    
!  zona reconstruida superior
   do j=jyt+nyr-1,ny
      do i=1,nx
         ux0(i,j,1)=umc
         uy0(i,j,1)=0.
         gx0(i,j,1)=0.
         gy0(i,j,1)=0.
      enddo
      do itime=1,ditime
         bx1(j,itime)=umc
         by1(j,itime)=0.
      enddo
   enddo
! 
   return  
end subroutine reconstruye5
!
!********************************************************************
!
subroutine vortz(ux,uy,wz,ffxp,fsxp,fwxp,ffyp,fsyp,fwyp,nx,ny,nz)
!
!********************************************************************
!
   implicit none
!
   real(8),dimension(nx,ny,nz) :: ux,uy,tx,ty,wz
   real(8),dimension(nx) :: ffxp,fsxp,fwxp
   real(8),dimension(ny) :: ffyp,fsyp,fwyp
!
   integer :: nx,ny,nz 
   integer :: k,j,i 
!
   call derx (tx,uy,ffxp,fsxp,fwxp,nx,ny,nz,1)
   call dery (ty,ux,ffyp,fsyp,fwyp,nx,ny,nz,1)
!
   do j=1,ny
   do i=1,nx
      wz(i,j,1)=tx(i,j,1)-ty(i,j,1)
   enddo
   enddo
!
   return
end subroutine vortz
!
!********************************************************************
!
subroutine vortz_r(ux,uy,wz,ffxp,fsxp,fwxp,ffyp,fsyp,fwyp,nx,ny,nz)
!
!********************************************************************
!
   implicit none
!
   real(8),dimension(nx,ny,nz) :: ux,uy,tx,ty,wz
   real(8),dimension(nx) :: ffxp,fsxp,fwxp
   real(8),dimension(ny) :: ffyp,fsyp,fwyp
!
   integer :: nx,ny,nz 
   integer :: k,j,i 
!
   call derx_r (tx,uy,ffxp,fsxp,fwxp,nx,ny,nz,1)
   call dery_r (ty,ux,ffyp,fsyp,fwyp,nx,ny,nz,1)
!
   do j=1,ny
   do i=1,nx
      wz(i,j,1)=tx(i,j,1)-ty(i,j,1)
   enddo
   enddo
!
   return
end subroutine vortz_r
!
!********************************************************************
!
subroutine vortz_rg(ux,uy,wz,ffxp,fsxp,fwxp,ffyp,fsyp,fwyp,nx,ny,nz)
!
!********************************************************************
!
   implicit none
!
   real(8),dimension(nx,ny,nz) :: ux,uy,tx,ty,wz
   real(8),dimension(nx) :: ffxp,fsxp,fwxp
   real(8),dimension(ny) :: ffyp,fsyp,fwyp
!
   integer :: nx,ny,nz 
   integer :: k,j,i 
!
   call derx_rg (tx,uy,ffxp,fsxp,fwxp,nx,ny,nz,1)
   call dery_rg (ty,ux,ffyp,fsyp,fwyp,nx,ny,nz,1)
!
   do j=1,ny
   do i=1,nx
      wz(i,j,1)=tx(i,j,1)-ty(i,j,1)
   enddo
   enddo
!
   return
end subroutine vortz_rg
!
!******************************************************************

real(8) function snorm1 (uxtest,uytest,uxtrue,uytrue,nx,ny,nz,jyt,nyr)

!******************************************************************
!
   implicit none
!
   real(8),dimension(nx,ny,nz) :: uxtest,uytest,uxtrue,uytrue
   real(8),dimension(nx,ny,nz) :: errorm
   real(8) :: deltux,deltuy
   integer :: j,i,nx,ny,nz,nxyr,jyt,nyr
!
   nxyr=nx*nyr
!
   do j=jyt,jyt+nyr-1
   do i=1,nx
      deltux = uxtest(i,j,1)-uxtrue(i,j,1)
      deltuy = uytest(i,j,1)-uytrue(i,j,1)
      errorm(i,j,1) = deltux**2 + deltuy**2
   enddo
   enddo 
!
!  promedio espacial de errorm
   snorm1=sum(errorm)/nxyr 
!
   return  
end function snorm1
!
!******************************************************************

real(8) function snorm2 (wztest,wztrue,nx,ny,nz)

!******************************************************************
!
   implicit none
!
   real(8),dimension(nx,ny,nz) :: wztest,wztrue
   real(8),dimension(nx,ny,nz) :: errorm
   real(8) :: deltwz
   integer :: j,i,nx,ny,nz,nxy
!
   nxy=nx*ny
!
   do j=1,ny
   do i=1,nx
      deltwz = wztest(i,j,1)-wztrue(i,j,1)
      errorm(i,j,1) = deltwz**2
   enddo
   enddo 
!
!  promedio espacial de errorm
   snorm2=sum(errorm)/nxy 
!
   return  
end function snorm2
!
!******************************************************************

real(8) function snorm3 (uxtest,uytest,uxtrue,uytrue,nx,ny,nz,&
             nxrg,nyrg,jw,iw,w_xy,sumw_xy,n_nodos,n_nodos_max)

!******************************************************************
!
   implicit none
!
   real(8),dimension(nx,ny,nz) :: uxtest,uytest,uxtrue,uytrue
   real(8),dimension(nxrg,nyrg) :: errorm
   integer,dimension(nxrg,nyrg,n_nodos_max) :: iw,jw
   real(8),dimension(nxrg,nyrg,n_nodos_max) :: w_xy
   integer,dimension(nxrg,nyrg) :: n_nodos
   real(8),dimension(nxrg,nyrg) :: sumw_xy
   real(8),dimension(n_nodos_max) :: wuxtest,wuytest,wuxtrue,wuytrue
   real(8) :: uxtest_av,uytest_av,uxtrue_av,uytrue_av
   real(8) :: deltux,deltuy
   real(8) :: w_r,sumw_r
   integer :: j_redg,i_redg,nxrg,nyrg,nxyrg,nodo
   integer :: nx,ny,nz,n_nodos_max,j,i
!
   nxyrg=nxrg*nyrg
!
   do j_redg = 1,nyrg
   do i_redg = 1,nxrg
!
!     /w_r(nodo)*u_r(nodo)
      wuxtest=0. ; wuytest=0. ; wuxtrue=0. ; wuytrue=0.
      do nodo=1,n_nodos(i_redg,j_redg)
         i=iw(i_redg,j_redg,nodo)
         j=jw(i_redg,j_redg,nodo)
         w_r=w_xy(i_redg,j_redg,nodo)
         wuxtest(nodo) = w_r*uxtest(i,j,1)
         wuytest(nodo) = w_r*uytest(i,j,1)
         wuxtrue(nodo) = w_r*uxtrue(i,j,1)
         wuytrue(nodo) = w_r*uytrue(i,j,1)
      enddo
!
!     Weighted average
      sumw_r=sumw_xy(i_redg,j_redg)
      uxtest_av = sum(wuxtest)/sumw_r
      uytest_av = sum(wuytest)/sumw_r
      uxtrue_av = sum(wuxtrue)/sumw_r
      uytrue_av = sum(wuytrue)/sumw_r
!
      deltux = uxtest_av-uxtrue_av
      deltuy = uytest_av-uytrue_av
      errorm(i_redg,j_redg) = deltux**2 + deltuy**2
   enddo
   enddo 
!
!  promedio espacial de errorm
   snorm3=sum(errorm)/nxyrg 
!
   return  
end function snorm3
!
!******************************************************************

real(8) function snorm4 (utrue,utest,nx,ny,nz,nxrg,nyrg,&
             jw,iw,w_xy,sumw_xy,n_nodos,n_nodos_max,&
             i_redg1,i_redg2,j_redg1,j_redg2,nxrg_asim,nyrg_asim)

!******************************************************************
!
   implicit none
!
   real(8),dimension(nxrg,nyrg,nz) :: utrue
   real(8),dimension(nx,ny,nz) :: utest
   real(8),dimension(nxrg,nyrg) :: errorm
   integer,dimension(nxrg,nyrg,n_nodos_max) :: iw,jw
   real(8),dimension(nxrg,nyrg,n_nodos_max) :: w_xy
   integer,dimension(nxrg,nyrg) :: n_nodos
   real(8),dimension(nxrg,nyrg) :: sumw_xy
   real(8),dimension(n_nodos_max) :: wutest
   real(8) :: utest_av
   real(8) :: deltu
   real(8) :: w_r,sumw_r
   integer :: j_redg,i_redg,nxrg,nyrg,nodo
   integer :: nx,ny,nz,n_nodos_max,j,i
   integer :: j_redg1,j_redg2,i_redg1,i_redg2,nxrg_asim,nyrg_asim,nxyrg_asim
!
   nxyrg_asim = nxrg_asim*nyrg_asim
!
   errorm=0.
   do j_redg = j_redg1,j_redg2
   do i_redg = i_redg1,i_redg2
!
!     /w_r(nodo)*u_r(nodo)
      wutest=0.
      do nodo=1,n_nodos(i_redg,j_redg)
         i=iw(i_redg,j_redg,nodo)
         j=jw(i_redg,j_redg,nodo)
         w_r=w_xy(i_redg,j_redg,nodo)
         wutest(nodo) = w_r*utest(i,j,1)
      enddo
!
!     Weighted average
      sumw_r=sumw_xy(i_redg,j_redg)
      utest_av = sum(wutest)/sumw_r
!
      deltu = utrue(i_redg,j_redg,1)-utest_av
      errorm(i_redg,j_redg) = deltu**2
   enddo
   enddo 
!
!  promedio espacial de errorm
   snorm4=sum(errorm)/nxyrg_asim 
!
   return  
end function snorm4
!
!******************************************************************

real(8) function snorm5 (utest,nx,ny,nz,nxrg,nyrg,&
             jw,iw,w_xy,sumw_xy,n_nodos,n_nodos_max,&
             i_redg1,i_redg2,j_redg1,j_redg2,nxrg_asim,nyrg_asim)

!******************************************************************
!
   implicit none
!
   real(8),dimension(nx,ny,nz) :: utest
   real(8),dimension(nxrg,nyrg,nz) :: utest_av
   integer,dimension(nxrg,nyrg,n_nodos_max) :: iw,jw
   real(8),dimension(nxrg,nyrg,n_nodos_max) :: w_xy
   integer,dimension(nxrg,nyrg) :: n_nodos
   real(8),dimension(nxrg,nyrg) :: sumw_xy
   real(8),dimension(n_nodos_max) :: wutest
   real(8) :: w_r,sumw_r
   integer :: j_redg,i_redg,nxrg,nyrg,nodo
   integer :: nx,ny,nz,n_nodos_max,j,i
   real(8),dimension(nxrg,nyrg) :: norma
   integer :: j_redg1,j_redg2,i_redg1,i_redg2,nxrg_asim,nyrg_asim,nxyrg_asim
!
   nxyrg_asim = nxrg_asim*nyrg_asim
!
   norma=0.
   do j_redg = j_redg1,j_redg2
   do i_redg = i_redg1,i_redg2
!
!     /w_r(nodo)*utest_r(nodo)
      wutest=0.
      do nodo=1,n_nodos(i_redg,j_redg)
         i=iw(i_redg,j_redg,nodo)
         j=jw(i_redg,j_redg,nodo)
         w_r=w_xy(i_redg,j_redg,nodo)
         wutest(nodo) = w_r*utest(i,j,1)
      enddo
!
!     Weighted average
      sumw_r=sumw_xy(i_redg,j_redg)
      utest_av(i_redg,j_redg,1) = sum(wutest)/sumw_r
!
      norma(i_redg,j_redg) = utest_av(i_redg,j_redg,1)
   enddo
   enddo 
!
!  promedio espacial de norma
   snorm5=sum(norma)/nxyrg_asim 
!
   return  
end function snorm5
!
!******************************************************************

real(8) function snorm6 (utrue,nxrg,nyrg,nz,&
             i_redg1,i_redg2,j_redg1,j_redg2,nxrg_asim,nyrg_asim)

!******************************************************************
!
   implicit none
!
   real(8),dimension(nxrg,nyrg,nz) :: utrue
   real(8),dimension(nxrg,nyrg) :: norma
   integer :: j_redg,i_redg,nxrg,nyrg,nz
   integer :: j_redg1,j_redg2,i_redg1,i_redg2,nxrg_asim,nyrg_asim,nxyrg_asim
!
   nxyrg_asim = nxrg_asim*nyrg_asim
!
   norma=0.
   do j_redg = j_redg1,j_redg2
   do i_redg = i_redg1,i_redg2
      norma(i_redg,j_redg) = utrue(i_redg,j_redg,1)
   enddo
   enddo 
!
!  promedio espacial de norma
   snorm6=sum(norma)/nxyrg_asim 
!
   return  
end function snorm6
!
!******************************************************************

subroutine snorm4_campo (utrue,utest,errorm,&
             nx,ny,nz,nxrg,nyrg,jw,iw,w_xy,sumw_xy,n_nodos,n_nodos_max,&
             i_redg1,i_redg2,j_redg1,j_redg2,nxrg_asim,nyrg_asim)

!******************************************************************
!
   implicit none
!
   real(8),dimension(nxrg,nyrg,nz) :: utrue
   real(8),dimension(nx,ny,nz) :: utest
   real(8),dimension(nxrg,nyrg) :: errorm
   integer,dimension(nxrg,nyrg,n_nodos_max) :: iw,jw
   real(8),dimension(nxrg,nyrg,n_nodos_max) :: w_xy
   integer,dimension(nxrg,nyrg) :: n_nodos
   real(8),dimension(nxrg,nyrg) :: sumw_xy
   real(8),dimension(n_nodos_max) :: wutest
   real(8) :: utest_av
   real(8) :: deltu
   real(8) :: w_r,sumw_r
   integer :: j_redg,i_redg,nxrg,nyrg,nodo
   integer :: nx,ny,nz,n_nodos_max,j,i
   integer :: j_redg1,j_redg2,i_redg1,i_redg2,nxrg_asim,nyrg_asim,nxyrg_asim
!
   nxyrg_asim = nxrg_asim*nyrg_asim
!
   errorm=0.
   do j_redg = j_redg1,j_redg2
   do i_redg = i_redg1,i_redg2
!
!     /w_r(nodo)*u_r(nodo)
      wutest=0.
      do nodo=1,n_nodos(i_redg,j_redg)
         i=iw(i_redg,j_redg,nodo)
         j=jw(i_redg,j_redg,nodo)
         w_r=w_xy(i_redg,j_redg,nodo)
         wutest(nodo) = w_r*utest(i,j,1)
      enddo
!
!     Weighted average
      sumw_r=sumw_xy(i_redg,j_redg)
      utest_av = sum(wutest)/sumw_r
!
      deltu = utrue(i_redg,j_redg,1)-utest_av
      errorm(i_redg,j_redg) = deltu**2
   enddo
   enddo 
!
   return  
end subroutine snorm4_campo
!
!******************************************************************

subroutine prueba1 (u0,u,fenetrem,imodulob)

!******************************************************************
!
!  reemplaza a Run Y=F(X0) 
!
   implicit none
!
   integer :: n,f,fenetrem,imodulob
   real(8) :: u0,utemp 
!
   real(8),dimension(fenetrem) :: u

   utemp=u0

   do f=1,fenetrem
      do n=1,imodulob
         utemp=3*utemp  ! reemplaza a call inc3dorigin
      enddo
      u(f)=utemp
   enddo
   
   return
end subroutine prueba1
!
!******************************************************************

subroutine prueba2 (u0,b1,u,fenetrem,imodulob,ditime)

!******************************************************************
!
!  reemplaza a Run Y=F(X0) 
!
   implicit none
!
   integer :: n,f,fenetrem,imodulob,itime,ditime
   real(8) :: u0,utemp 
!
   real(8),dimension(fenetrem) :: u
   real(8),dimension(ditime) :: b1

   utemp=u0

   do f=1,fenetrem
      do n=1,imodulob
         itime=n+(f-1)*imodulob
         utemp=3*utemp*b1(itime)  ! reemplaza a call inc3dorigin
      enddo
      u(f)=utemp
   enddo
   
   return
end subroutine prueba2
!
!********************************************************************
!
subroutine ludcmp(a,n,indx,d)
!
!********************************************************************
!
! Given a matrix a(1:n,1:n), this routine replaces it by the LU decomposition
! of a rowwise permutation of itself. a and n are input; a is output; indx(1:n)
! is an output vector that records the row permutation effected by the partial
! pivoting; d is output as +-1 depending on whether the number of row
! interchanges was even or odd, respectively. This routine is used in
! combination with lubksb to solve linear equations or invert a matrix
!
      implicit none
!
      integer,parameter :: NMAX=1000
      real(8),parameter :: TINY=1.0e-20
!
      real(8),dimension(n,n) :: a
      integer,dimension(n) :: indx
      real(8),dimension(NMAX) :: vv
      real(8) :: d,aamax,dum,sum
      integer :: n,i,imax,j,k
!
      d=1.
      do 12 i=1,n
        aamax=0.
        do 11 j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11      continue
        if (aamax.eq.0.) then
           print *,'singular matrix in ludcmp'
           stop
        endif
        vv(i)=1./aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.
        do 16 i=j,n
          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax) then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if (a(j,j).eq.0.) a(j,j)=TINY
        if (j.ne.n) then
          dum=1./a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
!
      return
end subroutine ludcmp
!
!********************************************************************
!
subroutine lubksb(a,n,indx,b)
!
!********************************************************************
!
! Solves the set of n linear equations AX = B. Here a is input, not as the matrix A
! but rather as its LU decomposition, determined by the routine ludcmp. indx is
! input as the permutation vector returned by ludcmp. b(1:n) is input as the
! right-hand side vector B, and returns with the solution vector X. a, n, and indx
! are not modified by this routine and can be left in place for successive
! calls with different right-hand sides b. This routine takes into account the
! possibility that b will begin with many zero elements, so it is efficient for
! use in matrix inversion
!
      implicit none
!
      real(8),dimension(n,n) :: a
      real(8),dimension(n) :: b
      integer,dimension(n) :: indx
      real(8) :: sum
      integer :: n,i,ii,j,ll
!
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0) then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue
!
      return
end subroutine lubksb
!
!********************************************************************
!
subroutine luinv(a,y,n)
!
!********************************************************************
!
! Using the above LU decomposition and backsubstitution routines, it is completely
! straightforward to find the inverse of a matrix a(1:n,1:n) column by column.
! The matrix y will contain the inverse of the original matrix a, which will have
! been destroyed
!
   implicit none
!
   real(8),dimension(n,n) :: a,y
   integer,dimension(n) :: indx
   real(8) :: d
   integer :: n,i,j
!
   do i=1,n 
      do j=1,n
         y(i,j)=0.
      enddo 
      y(i,i)=1.
   enddo 
   call ludcmp(a,n,indx,d)
   do j=1,n
     call lubksb(a,n,indx,y(1,j))
   enddo
!
   return
end subroutine luinv
!
!********************************************************************
!
subroutine interpx_gap(i,j,xm,ym,i1,i2,j1,j2,step_i,step_j,ux_obsg,f,&
             fenetre,nxrg,nyrg,nz,nvt_max,nip_max,nff_max)
!
!********************************************************************
!
      USE paramx_m 
      USE paramy_m 
!
      implicit none
!
      real(8),dimension(fenetre,nxrg,nyrg,nz) :: ux_obsg
      real(8),dimension(nvt_max) :: uxv,uyv,gamma
      real(8),dimension(nff_max,nff_max) :: vaf,vafinv,betax,xmvf
      real(8),dimension(nff_max) :: vafip,alfay,ymipv,nvf
      real(8) :: cexr,ceyr,xm,ym,uxip,uyip
      integer :: i1,i2,j1,j2,nip,nv,nvt,iv,jv
      integer :: inv,jnv,nff,nvv,f,i,j,k,p,step_i,step_j
      integer :: fenetre,nxrg,nyrg,nz,nvt_max,nip_max,nff_max
!
!     ********** Posicion y velocidad de vecinos de IP_n **********
!     ********** Nomenclatura vecinos: V_nv (general); ************
!     ********** V_inv,jnv (relativa a la formacion) **************
      nv=0
      jnv=0
!
!     **** fila inferior de la formacion ****
      jv=j1-1
      jnv=jnv+1
      ymipv(jnv)=(jv-1)*dyg ! posicion y/D del vecino virtual jnv de la formacion 
      inv=0
      do iv=i1-1,i2+1,step_i
         if (iv.eq.i) cycle
         inv=inv+1
         xmvf(inv,jnv)=(iv-1)*dxg ! posicion x/D del vecino inv en la fila jnv de la formacion 
         nv=nv+1
         uxv(nv)=ux_obsg(f,iv,jv,1) ! velocidad del vecino nv corresp. a IP_n, comp. x 
      enddo
      nvf(jnv)=inv ! numero de vecinos en la fila jnv de la formacion
!
!     **** filas intermedias de la formacion ****
      do jv=j1,j2,step_j
         if (jv.eq.j) cycle 
         jnv=jnv+1
         ymipv(jnv)=(jv-1)*dyg 
         inv=0
         do iv=i1-1,i2+1,i2-i1+2
            inv=inv+1
            xmvf(inv,jnv)=(iv-1)*dxg  
            nv=nv+1
            uxv(nv)=ux_obsg(f,iv,jv,1) 
         enddo
         nvf(jnv)=inv
      enddo
!
!     **** fila superior de la formacion ****
      jv=j2+1
      jnv=jnv+1
      ymipv(jnv)=(jv-1)*dyg 
      inv=0
      do iv=i1-1,i2+1,step_i
         if (iv.eq.i) cycle
         inv=inv+1
         xmvf(inv,jnv)=(iv-1)*dxg 
         nv=nv+1
         uxv(nv)=ux_obsg(f,iv,jv,1) 
      enddo
      nvf(jnv)=inv
!
      nvt=nv ! numero de vecinos totales en la formacion corresp. a IP_n
      if (nvt.gt.nvt_max) then
         print *,'nvt > nv_max ----> STOP' 
         stop
      endif
!
      nff=jnv ! numero de filas de la formacion
      if (nff.gt.nff_max) then
         print *,'nff > nff_max ----> STOP' 
         stop
      endif
!
!     ******************************************************
!     ******** Coef. esquema interp. segun x, betax ********
!     ******************************************************
      betax(:,:)=0.
!
      do jnv=1,nff
!
         nvv=nvf(jnv)
!
!        ****** Vandermonde matrix corresponding to the interpolating polynomial for IPV_(jnv) ******
         do inv=1,nvv
            k=0
            do p=1,nvv
               k=k+1
               vaf(inv,k)=xmvf(inv,jnv)**(p-1) ! basis function k evaluated at xmvf(inv,jnv)
            enddo
         enddo
!
!        *********** Inverse of the Vandermonde matrix ***********
         call luinv(vaf(1:nvv,1:nvv),vafinv(1:nvv,1:nvv),nvv)
!
!        ******** Coef. que ligan a IPV_(jnv) con sus vecinos (inv,jnv) ********
         k=0
         do p=1,nvv
            k=k+1
            vafip(k)=xm**(p-1) ! basis function k evaluated at xmip(n)
         enddo
!
         do inv=1,nvv
!           coef. que liga IPV_(jnv) con el vecino (inv,jnv)
            betax(inv,jnv)=dot_product(vafip(1:nvv),vafinv(1:nvv,inv))
         enddo
!
      enddo ! do jnv=1,nff
!
!     *****************************************************
!     ******** Coef. esquema interp. segun y, alfay ********
!     *****************************************************
      alfay(:)=0.
!
!     ****** Vandermonde matrix corresponding to the interpolating polynomial for IP ******
      do jnv=1,nff
         k=0
         do p=1,nff
            k=k+1
            vaf(jnv,k)=ymipv(jnv)**(p-1) ! basis function k evaluated at ymipv(jnv)
         enddo
      enddo
!
!     *********** Inverse of the Vandermonde matrix ***********
      call luinv(vaf(1:nff,1:nff),vafinv(1:nff,1:nff),nff)
!
!     ******** Coef. que ligan a IP con sus vecinos IPV_(jnv) ********
      k=0
      do p=1,nff
         k=k+1
         vafip(k)=ym**(p-1) ! basis function k evaluated at ymip(n)
      enddo
!
      do jnv=1,nff
!        coef. que liga IP con el vecino IPV_(jnv)
         alfay(jnv)=dot_product(vafip(1:nff),vafinv(1:nff,jnv))
      enddo
!
!     *****************************************************
!     ********* Coef. esquema interp. total, gamma ********
!     *****************************************************
!
!     ******** Coef. que ligan a IP_n con sus vecinos nv ********
      nv=0
      jnv=0
!
!     **** fila inferior de la formacion ****
      jnv=jnv+1
      inv=0
      do iv=i1-1,i2+1,step_i
         if (iv.eq.i) cycle
         inv=inv+1
         nv=nv+1
!        el vecino inv en la fila jnv de la formacion se corresponde con el vecino nv de IP
         gamma(nv)=betax(inv,jnv)*alfay(jnv)
      enddo
!
!     **** filas intermedias de la formacion ****
      do jv=j1,j2,step_j
         if (jv.eq.j) cycle 
         jnv=jnv+1
         inv=0
         do iv=i1-1,i2+1,i2-i1+2
            inv=inv+1
            nv=nv+1
            gamma(nv)=betax(inv,jnv)*alfay(jnv)
         enddo
      enddo
!
!     **** fila superior de la formacion ****
      jnv=jnv+1
      inv=0
      do iv=i1-1,i2+1,step_i
         if (iv.eq.i) cycle
         inv=inv+1
         nv=nv+1
         gamma(nv)=betax(inv,jnv)*alfay(jnv)
      enddo
!
!     ************************************************************
!     ******* Velocidad de nodos IP (Interpolated Points) ********
!     ************************************************************
      uxip=0.
      do nv=1,nvt
         uxip=uxip+gamma(nv)*uxv(nv) ! velocidad del IP-node n, comp. x
      enddo
!
!     ************************************************************
!     ************ Reconstruye campo veloc. en el gap ************
!     ************************************************************
      ux_obsg(f,i,j,1)=uxip
!
      return
end subroutine interpx_gap
!
!********************************************************************
!
subroutine interpy_gap(i,j,xm,ym,i1,i2,j1,j2,step_i,step_j,uy_obsg,f,&
             fenetre,nxrg,nyrg,nz,nvt_max,nip_max,ncf_max)
!
!********************************************************************
!
      USE paramx_m 
      USE paramy_m 
!
      implicit none
!
      real(8),dimension(fenetre,nxrg,nyrg,nz) :: uy_obsg
      real(8),dimension(nvt_max) :: uxv,uyv,gamma
      real(8),dimension(ncf_max,ncf_max) :: vac,vacinv,betay,ymvf
      real(8),dimension(ncf_max) :: vacip,alfax,xmipv,nvc
      real(8) :: cexr,ceyr,xm,ym,uxip,uyip
      integer :: i1,i2,j1,j2,nip,nv,nvt,iv,jv
      integer :: inv,jnv,ncf,nvv,f,i,j,k,p,step_i,step_j
      integer :: fenetre,nxrg,nyrg,nz,nvt_max,nip_max,ncf_max
!
!     ********** Posicion y velocidad de vecinos de IP_n **********
!     ********** nomenclatura vecinos: V_nv (general); ************
!     ********** V_inv,jnv (relativa a la formacion) **************
      nv=0
      inv=0
!
!     **** columna izquierda de la formacion ****
      iv=i1-1
      inv=inv+1
      xmipv(inv)=(iv-1)*dxg ! posicion x/D del vecino virtual inv de la formacion 
      jnv=0
      do jv=j1-1,j2+1,step_j
         if (jv.eq.j) cycle
         jnv=jnv+1
         ymvf(inv,jnv)=(jv-1)*dyg ! posicion y/D del vecino jnv en la columna inv de la formacion 
         nv=nv+1
         uyv(nv)=uy_obsg(f,iv,jv,1) ! velocidad del vecino nv corresp. a IP_n, comp. y
      enddo
      nvc(inv)=jnv ! numero de vecinos en la columna inv de la formacion
!
!     **** columnas intermedias de la formacion ****
      do iv=i1,i2,step_i
         if (iv.eq.i) cycle 
         inv=inv+1
         xmipv(inv)=(iv-1)*dxg 
         jnv=0
         do jv=j1-1,j2+1,j2-j1+2
            jnv=jnv+1
            ymvf(inv,jnv)=(jv-1)*dyg  
            nv=nv+1
            uyv(nv)=uy_obsg(f,iv,jv,1)
         enddo
         nvc(inv)=jnv
      enddo
!
!     **** columna derecha de la formacion ****
      iv=i2+1
      inv=inv+1
      xmipv(inv)=(iv-1)*dxg 
      jnv=0
      do jv=j1-1,j2+1,step_j
         if (jv.eq.j) cycle
         jnv=jnv+1
         ymvf(inv,jnv)=(jv-1)*dyg 
         nv=nv+1
         uyv(nv)=uy_obsg(f,iv,jv,1)
      enddo
      nvc(inv)=jnv
!
      nvt=nv ! numero de vecinos totales en la formacion corresp. a IP_n
      if (nvt.gt.nvt_max) then
         print *,'nvt > nv_max ----> STOP' 
         stop
      endif
!
      ncf=inv ! numero de columnas de la formacion
      if (ncf.gt.ncf_max) then
         print *,'ncf > ncf_max ----> STOP' 
         stop
      endif
!
!     ******************************************************
!     ******** Coef. esquema interp. segun y, betay ********
!     ******************************************************
      betay(:,:)=0.
!
      do inv=1,ncf
!
         nvv=nvc(inv)
!
!        ****** Vandermonde matrix corresponding to the interpolating polynomial for IPVc_(inv) ******
         do jnv=1,nvv
            k=0
            do p=1,nvv
               k=k+1
               vac(jnv,k)=ymvf(inv,jnv)**(p-1) ! basis function k evaluated at ymvf(inv,jnv)
            enddo
         enddo
!
!        *********** Inverse of the Vandermonde matrix ***********
         call luinv(vac(1:nvv,1:nvv),vacinv(1:nvv,1:nvv),nvv)
!
!        ******** Coef. que ligan a IPVc_(inv) con sus vecinos (inv,jnv) ********
         k=0
         do p=1,nvv
            k=k+1
            vacip(k)=ym**(p-1) ! basis function k evaluated at ymip(n)
         enddo
!
         do jnv=1,nvv
!           coef. que liga IPVc_(inv) con el vecino (inv,jnv)
            betay(inv,jnv)=dot_product(vacip(1:nvv),vacinv(1:nvv,jnv))
         enddo
!
      enddo ! do inv=1,ncf
!
!     ******************************************************
!     ******** Coef. esquema interp. segun x, alfax ********
!     ******************************************************
      alfax(:)=0.
!
!     ****** Vandermonde matrix corresponding to the interpolating polynomial for IP ******
      do inv=1,ncf
         k=0
         do p=1,ncf
            k=k+1
            vac(inv,k)=xmipv(inv)**(p-1) ! basis function k evaluated at xmipv(inv)
         enddo
      enddo
!
!     *********** Inverse of the Vandermonde matrix ***********
      call luinv(vac(1:ncf,1:ncf),vacinv(1:ncf,1:ncf),ncf)
!
!     ******** Coef. que ligan a IP con sus vecinos IPVc_(inv) ********
      k=0
      do p=1,ncf
         k=k+1
         vacip(k)=xm**(p-1) ! basis function k evaluated at xmip(n)
      enddo
!
      do inv=1,ncf
!        coef. que liga IP con el vecino IPVc_(inv)
         alfax(inv)=dot_product(vacip(1:ncf),vacinv(1:ncf,inv))
      enddo
!
!     *****************************************************
!     ********* Coef. esquema interp. total, gamma ********
!     *****************************************************
!
!     ******** Coef. que ligan a IP_n con sus vecinos nv ********
      nv=0
      inv=0
!
!     **** columna izquierda de la formacion ****
      inv=inv+1
      jnv=0
      do jv=j1-1,j2+1,step_j
         if (jv.eq.j) cycle
         jnv=jnv+1
         nv=nv+1
!        el vecino jnv en la columna inv de la formacion se corresponde con el vecino nv de IP
         gamma(nv)=betay(inv,jnv)*alfax(inv)
      enddo
!
!     **** columnas intermedias de la formacion ****
      do iv=i1,i2,step_i
         if (iv.eq.i) cycle 
         inv=inv+1
         jnv=0
         do jv=j1-1,j2+1,j2-j1+2
            jnv=jnv+1
            nv=nv+1
            gamma(nv)=betay(inv,jnv)*alfax(inv)
         enddo
      enddo
!
!     **** columna derecha de la formacion ****
      inv=inv+1
      jnv=0
      do jv=j1-1,j2+1,step_j
         if (jv.eq.j) cycle
         jnv=jnv+1
         nv=nv+1
         gamma(nv)=betay(inv,jnv)*alfax(inv)
      enddo
!
!     ************************************************************
!     ******* Velocidad de nodos IP (Interpolated Points) ********
!     ************************************************************
      uyip=0.
      do nv=1,nvt
         uyip=uyip+gamma(nv)*uyv(nv) ! velocidad del IP-node n, comp. y
      enddo
!
!     ************************************************************
!     ************ Reconstruye campo veloc. en el gap ************
!     ************************************************************
      uy_obsg(f,i,j,1)=uyip
!
      return
end subroutine interpy_gap
!


