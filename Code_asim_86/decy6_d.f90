!        Generated by TAPENADE     (INRIA, Tropics team)
!  Tapenade 3.4 (r3342) - 27 Jan 2010 15:25
!
!  Differentiation of decy6 in forward (tangent) mode:
!   variations   of useful results: ty
!   with respect to varying inputs: uy
!
!********************************************************************
!
SUBROUTINE DECY6_D(ty, tyd, uy, uyd, cfy6, csy6, cwy6, nx, nxm, ny, nym&
&  , nz, npaire)
!
  USE paramy_m
  implicit none
!
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  INTEGER :: nx, nxm, ny, nym, nz, npaire
  REAL*8 :: ty(nxm, nym, nz)
  REAL*8 :: tyd(nxm, nym, nz)
  REAL*8, INTENT(IN) :: uy(nxm, ny, nz)
  REAL*8, INTENT(IN) :: uyd(nxm, ny, nz)
  REAL*8 :: ry(nxm, nz, nym)
  REAL*8 :: ryd(nxm, nz, nym)
  REAL*8 :: di(nxm, nz, ny)
  REAL*8 :: did(nxm, nz, ny)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  INTEGER :: k, j, i
  REAL*8, DIMENSION(nym) :: cfy6, csy6, cwy6
  did = 0.0_8
!-----------------------------------------------
!
!********************************************************************
!
  DO k=1,nz
    DO j=1,ny
      DO i=1,nxm
        did(i, k, j) = uyd(i, j, k)
        di(i, k, j) = uy(i, j, k)
      END DO
    END DO
  END DO
!
  CALL DECY16_D(ry, ryd, di, did, cfy6, csy6, cwy6, nx, nxm, ny, nym, nz&
&          , npaire)
  tyd = 0.0_8
!
  DO k=1,nz
    DO j=1,nym
      DO i=1,nxm
        tyd(i, j, k) = yly*ryd(i, k, j)
        ty(i, j, k) = ry(i, k, j)*yly
      END DO
    END DO
  END DO
!
  RETURN
END SUBROUTINE DECY6_D