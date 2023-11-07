!        Generated by TAPENADE     (INRIA, Tropics team)
!  Tapenade 3.4 (r3342) - 27 Jan 2010 15:25
!
!  Differentiation of inter6 in forward (tangent) mode:
!   variations   of useful results: tx
!   with respect to varying inputs: ux
!
!********************************************************************
!
SUBROUTINE INTER6_D(tx, txd, ux, uxd, cifx6, cisx6, ciwx6, nx, nxm, ny, &
&  nz, npaire)
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
  USE paramx_m 
  USE interpol6_m 
  implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  INTEGER, INTENT(IN) :: nx, nxm, ny, nz, npaire
  REAL*8, INTENT(INOUT) :: tx(nxm, ny, nz)
  REAL*8, INTENT(INOUT) :: txd(nxm, ny, nz)
  REAL*8, INTENT(IN) :: ux(nx, ny, nz)
  REAL*8, INTENT(IN) :: uxd(nx, ny, nz)
  REAL*8, INTENT(IN) :: cifx6(nxm)
  REAL*8, INTENT(IN) :: cisx6(nxm)
  REAL*8, INTENT(IN) :: ciwx6(nxm)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  INTEGER :: k, j, i
!-----------------------------------------------
!
!********************************************************************
!   
  IF (npaire .EQ. 1) THEN
    txd = 0.0_8
    DO k=1,nz
      DO j=1,ny
        txd(1, j, k) = aicix6*(uxd(2, j, k)+uxd(1, j, k)) + bicix6*(uxd(&
&          3, j, k)+uxd(2, j, k))
        tx(1, j, k) = aicix6*(ux(2, j, k)+ux(1, j, k)) + bicix6*(ux(3, j&
&          , k)+ux(2, j, k))
        txd(2, j, k) = aicix6*(uxd(3, j, k)+uxd(2, j, k)) + bicix6*(uxd(&
&          4, j, k)+uxd(1, j, k))
        tx(2, j, k) = aicix6*(ux(3, j, k)+ux(2, j, k)) + bicix6*(ux(4, j&
&          , k)+ux(1, j, k))
      END DO
    END DO
    DO i=3,nxm-2
      DO k=1,nz
        DO j=1,ny
          txd(i, j, k) = aicix6*(uxd(i+1, j, k)+uxd(i, j, k)) + bicix6*(&
&            uxd(i+2, j, k)+uxd(i-1, j, k))
          tx(i, j, k) = aicix6*(ux(i+1, j, k)+ux(i, j, k)) + bicix6*(ux(&
&            i+2, j, k)+ux(i-1, j, k))
        END DO
      END DO
    END DO
    DO k=1,nz
      DO j=1,ny
        txd(nxm-1, j, k) = aicix6*(uxd(nxm, j, k)+uxd(nxm-1, j, k)) + &
&          bicix6*(uxd(nx, j, k)+uxd(nxm-2, j, k))
        tx(nxm-1, j, k) = aicix6*(ux(nxm, j, k)+ux(nxm-1, j, k)) + &
&          bicix6*(ux(nx, j, k)+ux(nxm-2, j, k))
        txd(nxm, j, k) = aicix6*(uxd(nx, j, k)+uxd(nxm, j, k)) + bicix6*&
&          (uxd(nxm, j, k)+uxd(nxm-1, j, k))
        tx(nxm, j, k) = aicix6*(ux(nx, j, k)+ux(nxm, j, k)) + bicix6*(ux&
&          (nxm, j, k)+ux(nxm-1, j, k))
      END DO
    END DO
    DO i=2,nxm
      DO k=1,nz
        DO j=1,ny
          txd(i, j, k) = txd(i, j, k) - cisx6(i)*txd(i-1, j, k)
          tx(i, j, k) = tx(i, j, k) - tx(i-1, j, k)*cisx6(i)
        END DO
      END DO
    END DO
    DO k=1,nz
      DO j=1,ny
        txd(nxm, j, k) = ciwx6(nxm)*txd(nxm, j, k)
        tx(nxm, j, k) = tx(nxm, j, k)*ciwx6(nxm)
      END DO
    END DO
    DO i=nxm-1,1,-1
      DO k=1,nz
        DO j=1,ny
          txd(i, j, k) = ciwx6(i)*(txd(i, j, k)-cifx6(i)*txd(i+1, j, k))
          tx(i, j, k) = (tx(i, j, k)-cifx6(i)*tx(i+1, j, k))*ciwx6(i)
        END DO
      END DO
    END DO
  ELSE
    txd = 0.0_8
  END IF
END SUBROUTINE INTER6_D