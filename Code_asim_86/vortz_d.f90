!        Generated by TAPENADE     (INRIA, Tropics team)
!  Tapenade 3.4 (r3342) - 27 Jan 2010 15:25
!
!  Differentiation of vortz in forward (tangent) mode:
!   variations   of useful results: wz ux uy
!   with respect to varying inputs: ux uy
!   RW status of diff variables: wz:out ux:in-out uy:in-out
!
!  Entrada: ux, uy, uxd, uyd. Salida: wzd
!  
!********************************************************************
!
SUBROUTINE VORTZ_D(ux, uxd, uy, uyd, wz, wzd, ffxp, fsxp, fwxp, ffyp, &
&  fsyp, fwyp, nx, ny, nz)
!
  implicit none
!
  INTEGER :: nx, ny, nz
!
  REAL*8, DIMENSION(nx, ny, nz) :: ux, uy, tx, ty, wz
  REAL*8, DIMENSION(nx, ny, nz) :: uxd, uyd, txd, tyd, wzd
  REAL*8, DIMENSION(nx) :: ffxp, fsxp, fwxp
  REAL*8, DIMENSION(ny) :: ffyp, fsyp, fwyp
  INTEGER :: k, j, i
!
  txd = 0.0_8
  tyd = 0.0_8
  wzd = 0.0_8
!
  CALL DERX_D(tx, txd, uy, uyd, ffxp, fsxp, fwxp, nx, ny, nz, 1)
  CALL DERY_D(ty, tyd, ux, uxd, ffyp, fsyp, fwyp, nx, ny, nz, 1)
!
  DO j=1,ny
    DO i=1,nx
      wzd(i, j, 1) = txd(i, j, 1) - tyd(i, j, 1)
      wz(i, j, 1) = tx(i, j, 1) - ty(i, j, 1)
    END DO
  END DO
!
  RETURN
END SUBROUTINE VORTZ_D