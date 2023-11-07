!        Generated by TAPENADE     (INRIA, Tropics team)
!  Tapenade 3.4 (r3342) - 27 Jan 2010 15:25
!
!  Differentiation of derxx in forward (tangent) mode:
!   variations   of useful results: tx
!   with respect to varying inputs: ux
!
!********************************************************************
!
SUBROUTINE DERXX_D(tx, txd, ux, uxd, sfx, ssx, swx, nx, ny, nz, npaire)
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
  USE paramx_m 
  USE derivex_m 
  implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  INTEGER, INTENT(IN) :: nx
  INTEGER, INTENT(IN) :: ny
  INTEGER, INTENT(IN) :: nz
  INTEGER, INTENT(IN) :: npaire
  REAL*8, INTENT(INOUT) :: tx(nx, ny, nz)
  REAL*8, INTENT(INOUT) :: txd(nx, ny, nz)
  REAL*8, INTENT(IN) :: ux(nx, ny, nz)
  REAL*8, INTENT(IN) :: uxd(nx, ny, nz)
  REAL*8, INTENT(IN) :: sfx(nx)
  REAL*8, INTENT(IN) :: ssx(nx)
  REAL*8, INTENT(IN) :: swx(nx)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  INTEGER :: k, j, i
!-----------------------------------------------
!
!********************************************************************
!
  IF (nclx .EQ. 2) THEN
    txd = 0.0_8
    DO k=1,nz
      DO j=1,ny
        txd(1, j, k) = as1x*uxd(1, j, k) + bs1x*uxd(2, j, k) + cs1x*uxd(&
&          3, j, k) + ds1x*uxd(4, j, k)
        tx(1, j, k) = as1x*ux(1, j, k) + bs1x*ux(2, j, k) + cs1x*ux(3, j&
&          , k) + ds1x*ux(4, j, k)
        txd(2, j, k) = as2x*(uxd(3, j, k)-2*uxd(2, j, k)+uxd(1, j, k))
        tx(2, j, k) = as2x*(ux(3, j, k)-ux(2, j, k)-ux(2, j, k)+ux(1, j&
&          , k))
      END DO
    END DO
!
    DO i=3,nx-2
      DO k=1,nz
        DO j=1,ny
          txd(i, j, k) = asix*(uxd(i+1, j, k)+uxd(i-1, j, k)) + bsix*(&
&            uxd(i+2, j, k)+uxd(i-2, j, k)) - (asix*2.+bsix*2.)*uxd(i, j&
&            , k)
          tx(i, j, k) = asix*(ux(i+1, j, k)+ux(i-1, j, k)) + bsix*(ux(i+&
&            2, j, k)+ux(i-2, j, k)) - (asix*2.+bsix*2.)*ux(i, j, k)
        END DO
      END DO
    END DO
!
    DO k=1,nz
      DO j=1,ny
        txd(nx-1, j, k) = asmx*(uxd(nx, j, k)-2*uxd(nx-1, j, k)+uxd(nx-2&
&          , j, k))
        tx(nx-1, j, k) = asmx*(ux(nx, j, k)-ux(nx-1, j, k)-ux(nx-1, j, k&
&          )+ux(nx-2, j, k))
        txd(nx, j, k) = asnx*uxd(nx, j, k) + bsnx*uxd(nx-1, j, k) + csnx&
&          *uxd(nx-2, j, k) + dsnx*uxd(nx-3, j, k)
        tx(nx, j, k) = asnx*ux(nx, j, k) + bsnx*ux(nx-1, j, k) + csnx*ux&
&          (nx-2, j, k) + dsnx*ux(nx-3, j, k)
      END DO
    END DO
    DO i=2,nx
      DO k=1,nz
        DO j=1,ny
          txd(i, j, k) = txd(i, j, k) - ssx(i)*txd(i-1, j, k)
          tx(i, j, k) = tx(i, j, k) - tx(i-1, j, k)*ssx(i)
        END DO
      END DO
    END DO
    DO k=1,nz
      DO j=1,ny
        txd(nx, j, k) = swx(nx)*txd(nx, j, k)
        tx(nx, j, k) = tx(nx, j, k)*swx(nx)
      END DO
    END DO
    DO i=nx-1,1,-1
      DO k=1,nz
        DO j=1,ny
          txd(i, j, k) = swx(i)*(txd(i, j, k)-sfx(i)*txd(i+1, j, k))
          tx(i, j, k) = (tx(i, j, k)-sfx(i)*tx(i+1, j, k))*swx(i)
        END DO
      END DO
    END DO
  ELSE
    txd = 0.0_8
  END IF
!
  RETURN
END SUBROUTINE DERXX_D