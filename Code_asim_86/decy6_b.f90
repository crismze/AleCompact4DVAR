!        Generated by TAPENADE     (INRIA, Tropics team)
!  Tapenade 3.4 (r3342) - 27 Jan 2010 15:25
!
!  Differentiation of decy6 in reverse (adjoint) mode:
!   gradient     of useful results: ty
!   with respect to varying inputs: uy
!
!********************************************************************
!
SUBROUTINE DECY6_B(ty, tyb, uy, uyb, cfy6, csy6, cwy6, nx, nxm, ny, nym&
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
  REAL*8 :: tyb(nxm, nym, nz)
  REAL*8, INTENT(IN) :: uy(nxm, ny, nz)
  REAL*8 :: uyb(nxm, ny, nz)
  REAL*8 :: ry(nxm, nz, nym)
  REAL*8 :: ryb(nxm, nz, nym)
  REAL*8 :: di(nxm, nz, ny)
  REAL*8 :: dib(nxm, nz, ny)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  INTEGER :: k, j, i
  REAL*8, DIMENSION(nym) :: cfy6, csy6, cwy6
  ryb = 0.0_8
!
  DO k=nz,1,-1
    DO j=nym,1,-1
      DO i=nxm,1,-1
        ryb(i, k, j) = ryb(i, k, j) + yly*tyb(i, j, k)
        tyb(i, j, k) = 0.0_8
      END DO
    END DO
  END DO
  CALL DECY16_B(ry, ryb, di, dib, cfy6, csy6, cwy6, nx, nxm, ny, nym, nz&
&          , npaire)
  uyb = 0.0_8
  DO k=nz,1,-1
    DO j=ny,1,-1
      DO i=nxm,1,-1
        uyb(i, j, k) = uyb(i, j, k) + dib(i, k, j)
        dib(i, k, j) = 0.0_8
      END DO
    END DO
  END DO
!
END SUBROUTINE DECY6_B
