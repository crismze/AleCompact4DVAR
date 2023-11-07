!        Generated by TAPENADE     (INRIA, Tropics team)
!  Tapenade 3.4 (r3342) - 27 Jan 2010 15:25
!
!  Differentiation of deryy in reverse (adjoint) mode:
!   gradient     of useful results: ty uy
!   with respect to varying inputs: uy
!
!********************************************************************
!
SUBROUTINE DERYY_B(ty, tyb, uy, uyb, sfy, ssy, swy, nx, ny, nz, npaire)
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
  USE paramy_m 
  USE derivey_m 
  implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  INTEGER, INTENT(IN) :: nx
  INTEGER, INTENT(IN) :: ny
  INTEGER, INTENT(IN) :: nz
  INTEGER, INTENT(IN) :: npaire
  REAL*8, INTENT(INOUT) :: ty(nx, nz, ny)
  REAL*8 :: tyb(nx, nz, ny)
  REAL*8, INTENT(IN) :: uy(nx, nz, ny)
  REAL*8 :: uyb(nx, nz, ny)
  REAL*8, INTENT(IN) :: sfy(ny)
  REAL*8, INTENT(IN) :: ssy(ny)
  REAL*8, INTENT(IN) :: swy(ny)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  INTEGER :: k, j, i
  INTEGER :: branch
  REAL*8 :: tempb6
  REAL*8 :: tempb5
  REAL*8 :: tempb4
  REAL*8 :: tempb3
  REAL*8 :: tempb2
  REAL*8 :: tempb1
  REAL*8 :: tempb0
  REAL*8 :: tempb
!-----------------------------------------------
!
!********************************************************************
!
  IF (ncly .EQ. 1) THEN
    IF (npaire .EQ. 1) THEN
      CALL PUSHINTEGER4(1)
    ELSE
      CALL PUSHINTEGER4(0)
    END IF
    IF (npaire .EQ. 0) THEN
      DO j=1,ny-1,1
        DO k=nz,1,-1
          DO i=nx,1,-1
            tyb(i, k, j+1) = tyb(i, k, j+1) - swy(j)*sfy(j)*tyb(i, k, j)
            tyb(i, k, j) = swy(j)*tyb(i, k, j)
          END DO
        END DO
      END DO
      DO k=nz,1,-1
        DO i=nx,1,-1
          tyb(i, k, ny) = swy(ny)*tyb(i, k, ny)
        END DO
      END DO
      DO j=ny,2,-1
        DO k=nz,1,-1
          DO i=nx,1,-1
            tyb(i, k, j-1) = tyb(i, k, j-1) - ssy(j)*tyb(i, k, j)
          END DO
        END DO
      END DO
      DO k=nz,1,-1
        DO i=nx,1,-1
          tyb(i, k, ny) = 0.0_8
          tempb6 = asjy*tyb(i, k, ny-1)
          uyb(i, k, ny) = uyb(i, k, ny) + tempb6
          uyb(i, k, ny-1) = uyb(i, k, ny-1) - bsjy*3*tyb(i, k, ny-1) - 2&
&            *tempb6
          uyb(i, k, ny-2) = uyb(i, k, ny-2) + tempb6
          uyb(i, k, ny-3) = uyb(i, k, ny-3) + bsjy*tyb(i, k, ny-1)
          tyb(i, k, ny-1) = 0.0_8
        END DO
      END DO
      DO k=nz,1,-1
        DO i=nx,1,-1
          DO j=ny-2,3,-1
            tempb4 = asjy*tyb(i, k, j)
            tempb5 = bsjy*tyb(i, k, j)
            uyb(i, k, j+1) = uyb(i, k, j+1) + tempb4
            uyb(i, k, j) = uyb(i, k, j) - 2*tempb5 - 2*tempb4
            uyb(i, k, j-1) = uyb(i, k, j-1) + tempb4
            uyb(i, k, j+2) = uyb(i, k, j+2) + tempb5
            uyb(i, k, j-2) = uyb(i, k, j-2) + tempb5
            tyb(i, k, j) = 0.0_8
          END DO
        END DO
      END DO
      DO k=nz,1,-1
        DO i=nx,1,-1
          tempb3 = asjy*tyb(i, k, 2)
          uyb(i, k, 3) = uyb(i, k, 3) + tempb3
          uyb(i, k, 2) = uyb(i, k, 2) - bsjy*3*tyb(i, k, 2) - 2*tempb3
          uyb(i, k, 1) = uyb(i, k, 1) + tempb3
          uyb(i, k, 4) = uyb(i, k, 4) + bsjy*tyb(i, k, 2)
          tyb(i, k, 2) = 0.0_8
          tyb(i, k, 1) = 0.0_8
        END DO
      END DO
    END IF
    CALL POPINTEGER4(branch)
    IF (.NOT.branch .LT. 1) THEN
      DO j=1,ny-1,1
        DO k=nz,1,-1
          DO i=nx,1,-1
            tyb(i, k, j+1) = tyb(i, k, j+1) - swy(j)*sfy(j)*tyb(i, k, j)
            tyb(i, k, j) = swy(j)*tyb(i, k, j)
          END DO
        END DO
      END DO
      DO k=nz,1,-1
        DO i=nx,1,-1
          tyb(i, k, ny) = swy(ny)*tyb(i, k, ny)
        END DO
      END DO
      DO j=ny,2,-1
        DO k=nz,1,-1
          DO i=nx,1,-1
            tyb(i, k, j-1) = tyb(i, k, j-1) - ssy(j)*tyb(i, k, j)
          END DO
        END DO
      END DO
      DO k=nz,1,-1
        DO i=nx,1,-1
          uyb(i, k, ny-1) = uyb(i, k, ny-1) + asjy*2*tyb(i, k, ny)
          uyb(i, k, ny) = uyb(i, k, ny) + ((-2)*bsjy-2*asjy)*tyb(i, k, &
&            ny)
          uyb(i, k, ny-2) = uyb(i, k, ny-2) + bsjy*2*tyb(i, k, ny)
          tyb(i, k, ny) = 0.0_8
          tempb2 = asjy*tyb(i, k, ny-1)
          uyb(i, k, ny) = uyb(i, k, ny) + tempb2
          uyb(i, k, ny-1) = uyb(i, k, ny-1) - bsjy*tyb(i, k, ny-1) - 2*&
&            tempb2
          uyb(i, k, ny-2) = uyb(i, k, ny-2) + tempb2
          uyb(i, k, ny-3) = uyb(i, k, ny-3) + bsjy*tyb(i, k, ny-1)
          tyb(i, k, ny-1) = 0.0_8
        END DO
      END DO
      DO k=nz,1,-1
        DO i=nx,1,-1
          DO j=ny-2,3,-1
            tempb0 = asjy*tyb(i, k, j)
            tempb1 = bsjy*tyb(i, k, j)
            uyb(i, k, j+1) = uyb(i, k, j+1) + tempb0
            uyb(i, k, j) = uyb(i, k, j) - 2*tempb1 - 2*tempb0
            uyb(i, k, j-1) = uyb(i, k, j-1) + tempb0
            uyb(i, k, j+2) = uyb(i, k, j+2) + tempb1
            uyb(i, k, j-2) = uyb(i, k, j-2) + tempb1
            tyb(i, k, j) = 0.0_8
          END DO
        END DO
      END DO
      DO k=nz,1,-1
        DO i=nx,1,-1
          tempb = asjy*tyb(i, k, 2)
          uyb(i, k, 3) = uyb(i, k, 3) + tempb
          uyb(i, k, 2) = uyb(i, k, 2) - bsjy*tyb(i, k, 2) - 2*tempb
          uyb(i, k, 1) = uyb(i, k, 1) + tempb
          uyb(i, k, 4) = uyb(i, k, 4) + bsjy*tyb(i, k, 2)
          tyb(i, k, 2) = 0.0_8
          uyb(i, k, 2) = uyb(i, k, 2) + asjy*2*tyb(i, k, 1)
          uyb(i, k, 1) = uyb(i, k, 1) + ((-2)*bsjy-2*asjy)*tyb(i, k, 1)
          uyb(i, k, 3) = uyb(i, k, 3) + bsjy*2*tyb(i, k, 1)
          tyb(i, k, 1) = 0.0_8
        END DO
      END DO
    END IF
  END IF
END SUBROUTINE DERYY_B
