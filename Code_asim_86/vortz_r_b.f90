!        Generated by TAPENADE     (INRIA, Tropics team)
!  Tapenade 3.4 (r3342) - 27 Jan 2010 15:25
!
!  Differentiation of vortz in reverse (adjoint) mode:
!   gradient     of useful results: wz
!   with respect to varying inputs: wz ux uy
!   RW status of diff variables: wz:in-out ux:out uy:out
!
!  Entrada: ux, uy, wzb. Salida: uxb, uyb
!
!********************************************************************
!
SUBROUTINE VORTZ_R_B(ux, uxb, uy, uyb, wz, wzb, ffxp, fsxp, fwxp, ffyp, &
&  fsyp, fwyp, nx, ny, nz)
!
  implicit none
!
  INTEGER :: nx, ny, nz
!
  REAL*8, DIMENSION(nx, ny, nz) :: ux, uy, tx, ty, wz
  REAL*8, DIMENSION(nx, ny, nz) :: uxb, uyb, txb, tyb, wzb
  REAL*8, DIMENSION(nx) :: ffxp, fsxp, fwxp
  REAL*8, DIMENSION(ny) :: ffyp, fsyp, fwyp
  INTEGER :: k, j, i
!
! ********************** Forward sweep **********************
!
  CALL DERX_R(tx, uy, ffxp, fsxp, fwxp, nx, ny, nz, 1)
  CALL DERY_R(ty, ux, ffyp, fsyp, fwyp, nx, ny, nz, 1)

! ********************** Backward sweep **********************
!
  txb = 0.0_8
  tyb = 0.0_8
!
  DO j=ny,1,-1
    DO i=nx,1,-1
      txb(i, j, 1) = txb(i, j, 1) + wzb(i, j, 1)
      tyb(i, j, 1) = tyb(i, j, 1) - wzb(i, j, 1)
      wzb(i, j, 1) = 0.0_8
    END DO
  END DO
!
  uxb = 0.0_8
  uyb = 0.0_8
!
  CALL DERY_R_B(ty, tyb, ux, uxb, ffyp, fsyp, fwyp, nx, ny, nz, 1)
  CALL DERX_R_B(tx, txb, uy, uyb, ffxp, fsxp, fwxp, nx, ny, nz, 1)
!
END SUBROUTINE VORTZ_R_B
