!        Generated by TAPENADE     (INRIA, Tropics team)
!  Tapenade 3.4 (r3342) - 27 Jan 2010 15:25
!
!  Differentiation of decx6 in reverse (adjoint) mode:
!   gradient     of useful results: tx
!   with respect to varying inputs: ux
!
!********************************************************************
!
SUBROUTINE DECX6_B(tx, txb, ux, uxb, cfx6, csx6, cwx6, nx, nxm, ny, nym&
&  , nz, npaire)
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
  USE paramx_m 
  USE dericex6_m 
  implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  INTEGER, INTENT(IN) :: nx
  INTEGER, INTENT(IN) :: nxm
  INTEGER :: ny
  INTEGER, INTENT(IN) :: nym
  INTEGER, INTENT(IN) :: nz
  INTEGER, INTENT(IN) :: npaire
  REAL, INTENT(INOUT) :: tx(nxm, nym, nz)
  REAL :: txb(nxm, nym, nz)
  REAL, INTENT(IN) :: ux(nx, nym, nz)
  REAL :: uxb(nx, nym, nz)
  REAL, INTENT(IN) :: cfx6(nxm)
  REAL, INTENT(IN) :: csx6(nxm)
  REAL, INTENT(IN) :: cwx6(nxm)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  INTEGER :: k, j, i
  INTEGER :: branch
  REAL :: tempb0
  REAL :: tempb
!-----------------------------------------------
!
  IF (npaire .EQ. 1) THEN
    CALL PUSHINTEGER4(1)
  ELSE
    CALL PUSHINTEGER4(0)
  END IF
  IF (npaire .EQ. 0) THEN
    DO i=1,nxm-1,1
      DO k=nz,1,-1
        DO j=nym,1,-1
          txb(i+1, j, k) = txb(i+1, j, k) - cwx6(i)*cfx6(i)*txb(i, j, k)
          txb(i, j, k) = cwx6(i)*txb(i, j, k)
        END DO
      END DO
    END DO
    DO k=nz,1,-1
      DO j=nym,1,-1
        txb(nxm, j, k) = cwx6(nxm)*txb(nxm, j, k)
      END DO
    END DO
    DO i=nxm,2,-1
      DO k=nz,1,-1
        DO j=nym,1,-1
          txb(i-1, j, k) = txb(i-1, j, k) - csx6(i)*txb(i, j, k)
        END DO
      END DO
    END DO
    uxb = 0.0
    DO k=nz,1,-1
      DO j=nym,1,-1
        tempb0 = bcix6*txb(nxm, j, k)
        uxb(nx, j, k) = uxb(nx, j, k) + 2.*tempb0 + acix6*txb(nxm, j, k)
        uxb(nxm, j, k) = uxb(nxm, j, k) - tempb0 - acix6*txb(nxm, j, k)
        uxb(nxm-1, j, k) = uxb(nxm-1, j, k) - tempb0
        txb(nxm, j, k) = 0.0
        uxb(nxm, j, k) = uxb(nxm, j, k) + acix6*txb(nxm-1, j, k)
        uxb(nxm-1, j, k) = uxb(nxm-1, j, k) - acix6*txb(nxm-1, j, k)
        uxb(nx, j, k) = uxb(nx, j, k) + bcix6*txb(nxm-1, j, k)
        uxb(nxm-2, j, k) = uxb(nxm-2, j, k) - bcix6*txb(nxm-1, j, k)
        txb(nxm-1, j, k) = 0.0
      END DO
    END DO
    DO i=nxm-2,3,-1
      DO k=nz,1,-1
        DO j=nym,1,-1
          uxb(i+1, j, k) = uxb(i+1, j, k) + acix6*txb(i, j, k)
          uxb(i, j, k) = uxb(i, j, k) - acix6*txb(i, j, k)
          uxb(i+2, j, k) = uxb(i+2, j, k) + bcix6*txb(i, j, k)
          uxb(i-1, j, k) = uxb(i-1, j, k) - bcix6*txb(i, j, k)
          txb(i, j, k) = 0.0
        END DO
      END DO
    END DO
    DO k=nz,1,-1
      DO j=nym,1,-1
        uxb(3, j, k) = uxb(3, j, k) + acix6*txb(2, j, k)
        uxb(2, j, k) = uxb(2, j, k) - acix6*txb(2, j, k)
        uxb(4, j, k) = uxb(4, j, k) + bcix6*txb(2, j, k)
        uxb(1, j, k) = uxb(1, j, k) - bcix6*txb(2, j, k)
        txb(2, j, k) = 0.0
        tempb = bcix6*txb(1, j, k)
        uxb(2, j, k) = uxb(2, j, k) + tempb + acix6*txb(1, j, k)
        uxb(1, j, k) = uxb(1, j, k) - 2.*tempb - acix6*txb(1, j, k)
        uxb(3, j, k) = uxb(3, j, k) + tempb
        txb(1, j, k) = 0.0
      END DO
    END DO
  ELSE
    uxb = 0.0
  END IF
  CALL POPINTEGER4(branch)
  IF (.NOT.branch .LT. 1) THEN
    DO i=1,nxm-1,1
      DO k=nz,1,-1
        DO j=nym,1,-1
          txb(i+1, j, k) = txb(i+1, j, k) - cwx6(i)*cfx6(i)*txb(i, j, k)
          txb(i, j, k) = cwx6(i)*txb(i, j, k)
        END DO
      END DO
    END DO
    DO k=nz,1,-1
      DO j=nym,1,-1
        txb(nxm, j, k) = cwx6(nxm)*txb(nxm, j, k)
      END DO
    END DO
    DO i=nxm,2,-1
      DO k=nz,1,-1
        DO j=nym,1,-1
          txb(i-1, j, k) = txb(i-1, j, k) - csx6(i)*txb(i, j, k)
        END DO
      END DO
    END DO
    DO k=nz,1,-1
      DO j=nym,1,-1
        uxb(nx, j, k) = uxb(nx, j, k) + acix6*txb(nxm, j, k)
        uxb(nxm, j, k) = uxb(nxm, j, k) + (bcix6-acix6)*txb(nxm, j, k)
        uxb(nxm-1, j, k) = uxb(nxm-1, j, k) - bcix6*txb(nxm, j, k)
        txb(nxm, j, k) = 0.0
        uxb(nxm, j, k) = uxb(nxm, j, k) + acix6*txb(nxm-1, j, k)
        uxb(nxm-1, j, k) = uxb(nxm-1, j, k) - acix6*txb(nxm-1, j, k)
        uxb(nx, j, k) = uxb(nx, j, k) + bcix6*txb(nxm-1, j, k)
        uxb(nxm-2, j, k) = uxb(nxm-2, j, k) - bcix6*txb(nxm-1, j, k)
        txb(nxm-1, j, k) = 0.0
      END DO
    END DO
    DO i=nxm-2,3,-1
      DO k=nz,1,-1
        DO j=nym,1,-1
          uxb(i+1, j, k) = uxb(i+1, j, k) + acix6*txb(i, j, k)
          uxb(i, j, k) = uxb(i, j, k) - acix6*txb(i, j, k)
          uxb(i+2, j, k) = uxb(i+2, j, k) + bcix6*txb(i, j, k)
          uxb(i-1, j, k) = uxb(i-1, j, k) - bcix6*txb(i, j, k)
          txb(i, j, k) = 0.0
        END DO
      END DO
    END DO
    DO k=nz,1,-1
      DO j=nym,1,-1
        uxb(3, j, k) = uxb(3, j, k) + acix6*txb(2, j, k)
        uxb(2, j, k) = uxb(2, j, k) - acix6*txb(2, j, k)
        uxb(4, j, k) = uxb(4, j, k) + bcix6*txb(2, j, k)
        uxb(1, j, k) = uxb(1, j, k) - bcix6*txb(2, j, k)
        txb(2, j, k) = 0.0
        uxb(2, j, k) = uxb(2, j, k) + (acix6-bcix6)*txb(1, j, k)
        uxb(1, j, k) = uxb(1, j, k) - acix6*txb(1, j, k)
        uxb(3, j, k) = uxb(3, j, k) + bcix6*txb(1, j, k)
        txb(1, j, k) = 0.0
      END DO
    END DO
  END IF
END SUBROUTINE DECX6_B
