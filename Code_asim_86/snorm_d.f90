!        Generated by TAPENADE     (INRIA, Tropics team)
!  Tapenade 3.4 (r3342) - 27 Jan 2010 15:25
!
!  Differentiation of snorm in forward (tangent) mode:
!   variations   of useful results: error
!   with respect to varying inputs: uy0r by1r ux uy ux0r bx1r
!   RW status of diff variables: error:out uy0r:in by1r:in ux:in
!                uy:in ux0r:in bx1r:in
!
!  Entrada: wz, ux0r, uy0r, gx0r, gy0r, bx1r, by1r, 
!           wzd, ux0rd, uy0rd, gx0rd, gy0rd, bx1rd, by1rd
!  Salida: errord 
!  
!********************************************************************
!
SUBROUTINE SNORM_D(errord, wz, wzd, wz_obsg, ux0r, ux0rd, uy0r, &
&  uy0rd, ux0r_obs, uy0r_obs, gx0r, gx0rd, gy0r, gy0rd, gx0r_obs, &
&  gy0r_obs, bx1r, bx1rd, by1r, by1rd, bx1_obs, by1_obs, fenetre, &
&  fenetrem, nx, ny, nz, nxr, nyr, nxrg, nyrg, jyt, rdxy, ditime, &
&  cov_obs, cov_iniu, cov_inig, cov_suav, cov_tay, cov_gau, jw, iw, w_xy, &
&  sumw_xy, n_nodos, n_nodos_max, jwq, iwq, w_xyq, sumw_xyq, n_nodos_q, &
&  n_nodos_q_max, itw, w_t, sumw_t, n_nodos_t, n_nodos_t_max, &
&  i_redg1, i_redg2, j_redg1, j_redg2, nxrg_asim, nyrg_asim)
!
  implicit none
!
  INTEGER :: f, fenetre, fenetrem, nx, ny, nz
  INTEGER :: j_redg, i_redg, nxrg, nyrg, nxyrg, jyt, rdxy, ditime
  INTEGER :: j_redg1, j_redg2, i_redg1, i_redg2, nxrg_asim, nyrg_asim, nxyrg_asim
!
  REAL*8, DIMENSION(fenetre, nxrg, nyrg, nz) :: wz_obsg
  REAL*8, DIMENSION(fenetrem, nxrg, nyrg, nz) :: wz_av
  REAL*8, DIMENSION(fenetrem, nxrg, nyrg, nz) :: wz_avd
  REAL*8, DIMENSION(fenetrem, nx, ny, nz) :: wz
  REAL*8, DIMENSION(fenetrem, nx, ny, nz) :: wzd
  INTEGER :: j_red, i_red, nxr, nyr, nxyr, j_rec, i_rec, j, i
  REAL*8, DIMENSION(nxr, nyr, nz) :: ux0r, uy0r, gx0r, gy0r
  REAL*8, DIMENSION(nxr, nyr, nz) :: ux0rd, uy0rd, gx0rd, gy0rd
  REAL*8, DIMENSION(nxr, nyr, nz) :: ux0r_av, uy0r_av
  REAL*8, DIMENSION(nxr, nyr, nz) :: ux0r_avd, uy0r_avd
  REAL*8, DIMENSION(nxr, nyr, nz) :: ux0r_obs, uy0r_obs, gx0r_obs, gy0r_obs
  REAL*8, DIMENSION(nyr, ditime) :: bx1r, by1r, bx1r_av, by1r_av
  REAL*8, DIMENSION(nyr, ditime) :: bx1rd, by1rd, bx1r_avd, by1r_avd
  REAL*8, DIMENSION(nyr, ditime) :: bx1_obs, by1_obs
  REAL*8, DIMENSION(nxrg, nyrg) :: errorm
  REAL*8, DIMENSION(nxrg, nyrg) :: errormd
  REAL*8, DIMENSION(nxr, nyr) :: errorlu, errorlg, errormq
  REAL*8, DIMENSION(nxr, nyr) :: errorlud, errorlgd, errormqd
  REAL*8, DIMENSION(nyr) :: errorp
  REAL*8, DIMENSION(nyr) :: errorpd
  REAL*8, DIMENSION(ditime) :: errorq
  REAL*8, DIMENSION(ditime) :: errorqd
  REAL*8, DIMENSION(fenetrem) :: errorv
  REAL*8, DIMENSION(fenetrem) :: errorvd
  REAL*8, DIMENSION(nyr) :: errorn
  REAL*8, DIMENSION(nyr) :: errornd
  REAL*8, DIMENSION(ditime) :: errorw
  REAL*8, DIMENSION(ditime) :: errorwd
  INTEGER :: itime, n_nodos_max, n_nodos_q_max, n_nodos_t_max, nodo
  INTEGER, DIMENSION(nxrg, nyrg, n_nodos_max) :: iw, jw
  REAL*8, DIMENSION(nxrg, nyrg, n_nodos_max) :: w_xy
  INTEGER, DIMENSION(nxrg, nyrg) :: n_nodos
  REAL*8, DIMENSION(nxrg, nyrg) :: sumw_xy
  REAL*8, DIMENSION(n_nodos_max) :: wwz
  REAL*8, DIMENSION(n_nodos_max) :: wwzd
  INTEGER, DIMENSION(nxr, nyr, n_nodos_q_max) :: iwq, jwq
  REAL*8, DIMENSION(nxr, nyr, n_nodos_q_max) :: w_xyq
  INTEGER, DIMENSION(nxr, nyr) :: n_nodos_q
  REAL*8, DIMENSION(nxr, nyr) :: sumw_xyq
  REAL*8, DIMENSION(n_nodos_q_max) :: wux, wuy
  REAL*8, DIMENSION(n_nodos_q_max) :: wuxd, wuyd
  INTEGER, DIMENSION(ditime, n_nodos_t_max) :: itw
  REAL*8, DIMENSION(ditime, n_nodos_t_max) :: w_t
  INTEGER, DIMENSION(ditime) :: n_nodos_t
  REAL*8, DIMENSION(ditime) :: sumw_t
  REAL*8, DIMENSION(n_nodos_t_max) :: wbx1r, wby1r
  REAL*8, DIMENSION(n_nodos_t_max) :: wbx1rd, wby1rd
  REAL*8 :: deltwz_obs
  REAL*8 :: deltwz_obsd
  REAL*8 :: deltux_ini, deltuy_ini, deltgx_ini, deltgy_ini
  REAL*8 :: deltux_inid, deltuy_inid, deltgx_inid, deltgy_inid
  REAL*8 :: deltux_suav, deltuy_suav
  REAL*8 :: deltux_suavd, deltuy_suavd
  REAL*8 :: deltux_tay, deltuy_tay, deltux_gau, deltuy_gau
  REAL*8 :: deltux_tayd, deltuy_tayd, deltux_gaud, deltuy_gaud
  REAL*8 :: error, error_obs, error_iniu, error_inig, error_suav, error_tay, error_gau
  REAL*8 :: errord, error_obsd, error_iniud, error_inigd, error_suavd, error_tayd, error_gaud
  REAL*8 :: cov_obs, cov_iniu, cov_inig, cov_suav, cov_tay, cov_gau
  REAL*8 :: w_r, sumw_r, w_q, sumw_q, w_d, sumw_d
!
  nxyr = nxr*nyr
  nxyrg = nxrg*nyrg
  nxyrg_asim = nxrg_asim*nyrg_asim
!
! ********************* error_obs *********************
  errorvd = 0.0_8
  wz_avd = 0.0_8
!
  DO f=1,fenetrem
!
    errorm = 0.
    errormd = 0.0_8
    DO j_redg=j_redg1,j_redg2
      DO i_redg=i_redg1,i_redg2
!
!       /w_r(nodo)*wz_r(nodo)
        wwz = 0.
        wwzd = 0.0_8
        DO nodo=1,n_nodos(i_redg, j_redg)
          i = iw(i_redg, j_redg, nodo)
          j = jw(i_redg, j_redg, nodo)
          w_r = w_xy(i_redg, j_redg, nodo)
          wwzd(nodo) = w_r*wzd(f, i, j, 1)
          wwz(nodo) = w_r*wz(f, i, j, 1)
        END DO
!
!       Weighted average
        sumw_r = sumw_xy(i_redg, j_redg)
        wz_avd(f, i_redg, j_redg, 1) = SUM(wwzd)/sumw_r
        wz_av(f, i_redg, j_redg, 1) = SUM(wwz)/sumw_r
!
        deltwz_obsd = wz_avd(f, i_redg, j_redg, 1)
        deltwz_obs = wz_av(f, i_redg, j_redg, 1) - wz_obsg(f+1, i_redg, j_redg, 1)
        errormd(i_redg, j_redg) = cov_obs*2*deltwz_obs*deltwz_obsd
        errorm(i_redg, j_redg) = cov_obs*deltwz_obs**2
      END DO
    END DO
!
!   promedio espacial de errorm
    errorvd(f) = SUM(errormd)/nxyrg_asim
    errorv(f) = SUM(errorm)/nxyrg_asim
  END DO
!
! promedio temporal de errorv
  error_obsd = SUM(errorvd)/fenetrem
  error_obs = SUM(errorv)/fenetrem
!
! ********************* error_iniu *********************
  errorlud = 0.0_8
! 
  DO j=1,nyr
    DO i=1,nxr
      deltux_inid = ux0rd(i, j, 1)
      deltux_ini = ux0r(i, j, 1) - ux0r_obs(i, j, 1)
      deltuy_inid = uy0rd(i, j, 1)
      deltuy_ini = uy0r(i, j, 1) - uy0r_obs(i, j, 1)
      errorlud(i, j) = cov_iniu*(2*deltux_ini*deltux_inid+2*deltuy_ini*deltuy_inid)
      errorlu(i, j) = cov_iniu*(deltux_ini**2+deltuy_ini**2)
    END DO
  END DO
!
! promedio espacial de errorlu
  error_iniud = SUM(errorlud)/nxyr
  error_iniu = SUM(errorlu)/nxyr
!
! ********************* error_inig *********************
  errorlgd = 0.0_8
! 
  DO j=1,nyr
    DO i=1,nxr
      deltgx_inid = gx0rd(i, j, 1)
      deltgx_ini = gx0r(i, j, 1) - gx0r_obs(i, j, 1)
      deltgy_inid = gy0rd(i, j, 1)
      deltgy_ini = gy0r(i, j, 1) - gy0r_obs(i, j, 1)
      errorlgd(i, j) = cov_inig*(2*deltgx_ini*deltgx_inid+2*deltgy_ini*deltgy_inid)
      errorlg(i, j) = cov_inig*(deltgx_ini**2+deltgy_ini**2)
    END DO
  END DO
!
! promedio espacial de errorlg
  error_inigd = SUM(errorlgd)/nxyr
  error_inig = SUM(errorlg)/nxyr
!
! ********************* error_suav *********************
  errormqd = 0.0_8
  uy0r_avd = 0.0_8
  ux0r_avd = 0.0_8
!
  DO j_red=1,nyr
    DO i_red=1,nxr
!
!     /w_r(nodo)*u_r(nodo)
      wux = 0.
      wuy = 0.
      wuxd = 0.0_8
      wuyd = 0.0_8
      DO nodo=1,n_nodos_q(i_red, j_red)
        i = iwq(i_red, j_red, nodo)
        j = jwq(i_red, j_red, nodo)
        w_q = w_xyq(i_red, j_red, nodo)
        wuxd(nodo) = w_q*ux0rd(i, j, 1)
        wux(nodo) = w_q*ux0r(i, j, 1)
        wuyd(nodo) = w_q*uy0rd(i, j, 1)
        wuy(nodo) = w_q*uy0r(i, j, 1)
      END DO
!
!     Weighted average
      sumw_q = sumw_xyq(i_red, j_red)
      ux0r_avd(i_red, j_red, 1) = SUM(wuxd)/sumw_q
      ux0r_av(i_red, j_red, 1) = SUM(wux)/sumw_q
      uy0r_avd(i_red, j_red, 1) = SUM(wuyd)/sumw_q
      uy0r_av(i_red, j_red, 1) = SUM(wuy)/sumw_q
!
      deltux_suavd = ux0rd(i_red, j_red, 1) - ux0r_avd(i_red, j_red, 1)
      deltux_suav = ux0r(i_red, j_red, 1) - ux0r_av(i_red, j_red, 1)
      deltuy_suavd = uy0rd(i_red, j_red, 1) - uy0r_avd(i_red, j_red, 1)
      deltuy_suav = uy0r(i_red, j_red, 1) - uy0r_av(i_red, j_red, 1)
      errormqd(i_red, j_red) = cov_suav*(2*deltux_suav*deltux_suavd+2*deltuy_suav*deltuy_suavd)
      errormq(i_red, j_red) = cov_suav*(deltux_suav**2+deltuy_suav**2)
    END DO
  END DO
!
! promedio espacial de errormq
  error_suavd = SUM(errormqd)/nxyr
  error_suav = SUM(errormq)/nxyr
! 
! ********************* error_tay *********************
  errornd = 0.0_8
  errorwd = 0.0_8
!
  DO itime=1,ditime
!
    DO j=1,nyr
      deltux_tayd = bx1rd(j, itime)
      deltux_tay = bx1r(j, itime) - bx1_obs(j, itime)
      deltuy_tayd = by1rd(j, itime)
      deltuy_tay = by1r(j, itime) - by1_obs(j, itime)
      errornd(j) = cov_tay*(2*deltux_tay*deltux_tayd + 2*deltuy_tay*deltuy_tayd)
      errorn(j) = cov_tay*(deltux_tay**2 + deltuy_tay**2)
    END DO
!
!   promedio espacial de errorn
    errorwd(itime) = SUM(errornd)/nyr
    errorw(itime) = SUM(errorn)/nyr
  END DO
!
! promedio temporal de errorw
  error_tayd = SUM(errorwd)/ditime
  error_tay = SUM(errorw)/ditime
! 
! ********************* error_gau *********************
  errorpd = 0.0_8
  errorqd = 0.0_8
  bx1r_avd = 0.0_8
  by1r_avd = 0.0_8
!
  DO itime=1,ditime
!
    DO j=1,nyr
!
!     /w_t(nodo)*b1r_t(nodo)
      wbx1r = 0.
      wby1r = 0.
      wby1rd = 0.0_8
      wbx1rd = 0.0_8
      DO nodo=1,n_nodos_t(itime)
        i = itw(itime, nodo)
        w_d = w_t(itime, nodo)
        wbx1rd(nodo) = w_d*bx1rd(j, i)
        wbx1r(nodo) = w_d*bx1r(j, i)
        wby1rd(nodo) = w_d*by1rd(j, i)
        wby1r(nodo) = w_d*by1r(j, i)
      END DO
!
!     Weighted average
      sumw_d = sumw_t(itime)
      bx1r_avd(j, itime) = SUM(wbx1rd)/sumw_d
      bx1r_av(j, itime) = SUM(wbx1r)/sumw_d
      by1r_avd(j, itime) = SUM(wby1rd)/sumw_d
      by1r_av(j, itime) = SUM(wby1r)/sumw_d
!
      deltux_gaud = bx1rd(j, itime) - bx1r_avd(j, itime)
      deltux_gau = bx1r(j, itime) - bx1r_av(j, itime)
      deltuy_gaud = by1rd(j, itime) - by1r_avd(j, itime)
      deltuy_gau = by1r(j, itime) - by1r_av(j, itime)
      errorpd(j) = cov_gau*(2*deltux_gau*deltux_gaud+2*deltuy_gau*deltuy_gaud)
      errorp(j) = cov_gau*(deltux_gau**2+deltuy_gau**2)
    END DO
!
!   promedio espacial de errorp
    errorqd(itime) = SUM(errorpd)/nyr
    errorq(itime) = SUM(errorp)/nyr
  END DO
!
! promedio temporal de errorq
  error_gaud = SUM(errorqd)/ditime
  error_gau = SUM(errorq)/ditime
!
! *********************** error ***********************
  errord = error_obsd + error_iniud + error_inigd + error_suavd + error_tayd + error_gaud
  error = error_obs + error_iniu + error_inig + error_suav + error_tay + error_gau
!
  RETURN
END SUBROUTINE SNORM_D