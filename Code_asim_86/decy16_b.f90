!        Generated by TAPENADE     (INRIA, Tropics team)
!  Tapenade 3.4 (r3342) - 27 Jan 2010 15:25
!
!  Differentiation of decy16 in reverse (adjoint) mode:
!   gradient     of useful results: ty
!   with respect to varying inputs: uy
!
!********************************************************************
!
SUBROUTINE DECY16_B(ty, tyb, uy, uyb, cfy6, csy6, cwy6, nx, nxm, ny, nym&
&  , nz, npaire)
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
  USE dericey6_m 
  USE paramy_m 
  implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  INTEGER :: nx, nxm, ny, nym, nz, npaire
  REAL*8, INTENT(INOUT) :: ty(nxm, nz, nym)
  REAL*8 :: tyb(nxm, nz, nym)
  REAL*8, INTENT(IN) :: uy(nxm, nz, ny)
  REAL*8 :: uyb(nxm, nz, ny)
  REAL*8, INTENT(IN) :: cfy6(nym)
  REAL*8, INTENT(IN) :: csy6(nym)
  REAL*8, INTENT(IN) :: cwy6(nym)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  INTEGER :: k, i, j
  REAL*8 :: tempb0
  REAL*8 :: tempb
!-----------------------------------------------
!
!********************************************************************
!
  IF (npaire .EQ. 0) THEN
    DO j=1,nym-1,1
      DO k=nz,1,-1
        DO i=nxm,1,-1
          tyb(i, k, j+1) = tyb(i, k, j+1) - cwy6(j)*cfy6(j)*tyb(i, k, j)
          tyb(i, k, j) = cwy6(j)*tyb(i, k, j)
        END DO
      END DO
    END DO
    DO k=nz,1,-1
      DO i=nxm,1,-1
        tyb(i, k, nym) = cwy6(nym)*tyb(i, k, nym)
      END DO
    END DO
    DO j=nym,2,-1
      DO k=nz,1,-1
        DO i=nxm,1,-1
          tyb(i, k, j-1) = tyb(i, k, j-1) - csy6(j)*tyb(i, k, j)
        END DO
      END DO
    END DO
    uyb = 0.0_8
    DO k=nz,1,-1
      DO i=nxm,1,-1
        tempb0 = bciy6*tyb(i, k, nym)
        uyb(i, k, ny) = uyb(i, k, ny) + 2.*tempb0 + aciy6*tyb(i, k, nym)
        uyb(i, k, nym) = uyb(i, k, nym) - tempb0 - aciy6*tyb(i, k, nym)
        uyb(i, k, nym-1) = uyb(i, k, nym-1) - tempb0
        tyb(i, k, nym) = 0.0_8
        uyb(i, k, nym) = uyb(i, k, nym) + aciy6*tyb(i, k, nym-1)
        uyb(i, k, nym-1) = uyb(i, k, nym-1) - aciy6*tyb(i, k, nym-1)
        uyb(i, k, ny) = uyb(i, k, ny) + bciy6*tyb(i, k, nym-1)
        uyb(i, k, nym-2) = uyb(i, k, nym-2) - bciy6*tyb(i, k, nym-1)
        tyb(i, k, nym-1) = 0.0_8
      END DO
    END DO
    DO j=nym-2,3,-1
      DO k=nz,1,-1
        DO i=nxm,1,-1
          uyb(i, k, j+1) = uyb(i, k, j+1) + aciy6*tyb(i, j, k)
          uyb(i, k, j) = uyb(i, k, j) - aciy6*tyb(i, j, k)
          uyb(i, k, j+2) = uyb(i, k, j+2) + bciy6*tyb(i, j, k)
          uyb(i, k, j-1) = uyb(i, k, j-1) - bciy6*tyb(i, j, k)
          tyb(i, j, k) = 0.0_8
        END DO
      END DO
    END DO
    DO k=nz,1,-1
      DO i=nxm,1,-1
        uyb(i, k, 3) = uyb(i, k, 3) + aciy6*tyb(i, k, 2)
        uyb(i, k, 2) = uyb(i, k, 2) - aciy6*tyb(i, k, 2)
        uyb(i, k, 4) = uyb(i, k, 4) + bciy6*tyb(i, k, 2)
        uyb(i, k, 1) = uyb(i, k, 1) - bciy6*tyb(i, k, 2)
        tyb(i, k, 2) = 0.0_8
        tempb = bciy6*tyb(i, k, 1)
        uyb(i, k, 2) = uyb(i, k, 2) + tempb + aciy6*tyb(i, k, 1)
        uyb(i, k, 1) = uyb(i, k, 1) - 2.*tempb - aciy6*tyb(i, k, 1)
        uyb(i, k, 3) = uyb(i, k, 3) + tempb
        tyb(i, k, 1) = 0.0_8
      END DO
    END DO
  ELSE
    uyb = 0.0_8
  END IF
END SUBROUTINE DECY16_B
