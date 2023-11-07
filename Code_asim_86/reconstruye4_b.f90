!        Generated by TAPENADE     (INRIA, Tropics team)
!  Tapenade 3.4 (r3342) - 27 Jan 2010 15:25
!
!  Differentiation of reconstruye in reverse (adjoint) mode:
!   gradient     of useful results: ux0 uy0 bx1 by1
!   with respect to varying inputs: ux0 uy0r by1r uy0 bx1 ux0r
!                bx1r by1
!   RW status of diff variables: ux0:in-out uy0r:out by1r:out uy0:in-out
!                bx1:in-out ux0r:out bx1r:out by1:in-out
!
!  Entrada: ux0r, uy0r, gx0r, gy0r, bx1r, by1r, ux0b, uy0b, gx0b, gy0b, bx1b, by1b 
!  Salida: ux0rb, uy0rb, gx0rb, gy0rb, bx1rb, by1rb 
!  
!********************************************************************
!
SUBROUTINE RECONSTRUYE4_B(ux0r, ux0rb, uy0r, uy0rb, gx0r, gx0rb, gy0r, gy0rb, bx1r, bx1rb, &
&  by1r, by1rb, ux0, ux0b, uy0, uy0b, gx0, gx0b, gy0, gy0b, bx1, bx1b, by1, by1b, &
&  jyc, jyt, nxr, nyr, nx, ny, nz, ditime)
!
  implicit none
!
  INTEGER :: jyc, jyt, nxr, nyr, nx, ny, nz, ditime
!
  REAL*8, DIMENSION(nxr, nyr, nz) :: ux0r, uy0r, gx0r, gy0r
  REAL*8, DIMENSION(nxr, nyr, nz) :: ux0rb, uy0rb, gx0rb, gy0rb
  REAL*8, DIMENSION(nx, ny, nz) :: ux0, uy0, gx0, gy0
  REAL*8, DIMENSION(nx, ny, nz) :: ux0b, uy0b, gx0b, gy0b
  REAL*8, DIMENSION(nyr, ditime) :: bx1r, by1r
  REAL*8, DIMENSION(nyr, ditime) :: bx1rb, by1rb
  REAL*8, DIMENSION(ny, ditime) :: bx1, by1
  REAL*8, DIMENSION(ny, ditime) :: bx1b, by1b
  REAL*8 :: pi, umc, ramp_base, ramp_exp, ramp
  INTEGER :: j0, jred, deltaj, pre
  INTEGER :: itime, i, j
!
  pi = ACOS(-1.)
  deltaj = jyt - jyc
  pre = 9999
!  
  DO j=jyc+1,jyt
    j0 = j - jyc
    ramp_base = 0.5*SIN(pi*j0/deltaj-pi/2) + 0.5
    ramp_exp = 1/((j0/deltaj)**pre+1)
    CALL PUSHREAL8(ramp)
    ramp = ramp_base**ramp_exp
  END DO
!    
  jred = 1
  DO j=jyt+1,jyt+nyr-2
    CALL PUSHINTEGER4(jred)
    jred = jred + 1
  END DO
!    
  DO j=jyt+nyr-1,2*jyt+nyr-jyc-1
    j0 = j - (jyt+nyr-2)
    ramp_base = 0.5*SIN(pi*j0/deltaj-pi/2) + 0.5
    ramp_exp = 1/((j0/deltaj)**pre+1)
    CALL PUSHREAL8(ramp)
    ramp = ramp_base**ramp_exp
  END DO
!
  DO j=2*jyt+nyr,2*jyt+nyr-jyc,-1
    DO itime=ditime,1,-1
      by1b(j, itime) = 0.0_8
      bx1b(j, itime) = 0.0_8
    END DO
    DO i=nx,1,-1
      uy0b(i, j, 1) = 0.0_8
      ux0b(i, j, 1) = 0.0_8
      gy0b(i, j, 1) = 0.0_8
      gx0b(i, j, 1) = 0.0_8
    END DO
  END DO
!
!a  uy0rb = 0.0_8
!a  by1rb = 0.0_8
!a  ux0rb = 0.0_8
!a  bx1rb = 0.0_8
!
  DO j=2*jyt+nyr-jyc-1,jyt+nyr-1,-1
    DO itime=ditime,1,-1
      by1rb(nyr, itime) = by1rb(nyr, itime) + (1-ramp)*by1b(j, itime)
      by1b(j, itime) = 0.0_8
      bx1rb(nyr, itime) = bx1rb(nyr, itime) + (1.0-ramp)*bx1b(j, itime)
      bx1b(j, itime) = 0.0_8
    END DO
    DO i=nx,1,-1
      uy0rb(i, nyr, 1) = uy0rb(i, nyr, 1) + (1-ramp)*uy0b(i, j, 1)
      uy0b(i, j, 1) = 0.0_8
      ux0rb(i, nyr, 1) = ux0rb(i, nyr, 1) + (1.0-ramp)*ux0b(i, j, 1)
      ux0b(i, j, 1) = 0.0_8
      gy0rb(i, nyr, 1) = gy0rb(i, nyr, 1) + (1-ramp)*gy0b(i, j, 1)
      gy0b(i, j, 1) = 0.0_8
      gx0rb(i, nyr, 1) = gx0rb(i, nyr, 1) + (1.0-ramp)*gx0b(i, j, 1)
      gx0b(i, j, 1) = 0.0_8
    END DO
    CALL POPREAL8(ramp)
  END DO
!
  DO j=jyt+nyr-2,jyt+1,-1
    DO itime=ditime,1,-1
      by1rb(jred, itime) = by1rb(jred, itime) + by1b(j, itime)
      by1b(j, itime) = 0.0_8
      bx1rb(jred, itime) = bx1rb(jred, itime) + bx1b(j, itime)
      bx1b(j, itime) = 0.0_8
    END DO
    DO i=nx,1,-1
      uy0rb(i, jred, 1) = uy0rb(i, jred, 1) + uy0b(i, j, 1)
      uy0b(i, j, 1) = 0.0_8
      ux0rb(i, jred, 1) = ux0rb(i, jred, 1) + ux0b(i, j, 1)
      ux0b(i, j, 1) = 0.0_8
      gy0rb(i, jred, 1) = gy0rb(i, jred, 1) + gy0b(i, j, 1)
      gy0b(i, j, 1) = 0.0_8
      gx0rb(i, jred, 1) = gx0rb(i, jred, 1) + gx0b(i, j, 1)
      gx0b(i, j, 1) = 0.0_8
    END DO
    CALL POPINTEGER4(jred)
  END DO
!
  DO j=jyt,jyc+1,-1
    DO itime=ditime,1,-1
      by1rb(1, itime) = by1rb(1, itime) + ramp*by1b(j, itime)
      by1b(j, itime) = 0.0_8
      bx1rb(1, itime) = bx1rb(1, itime) + ramp*bx1b(j, itime)
      bx1b(j, itime) = 0.0_8
    END DO
    DO i=nx,1,-1
      uy0rb(i, 1, 1) = uy0rb(i, 1, 1) + ramp*uy0b(i, j, 1)
      uy0b(i, j, 1) = 0.0_8
      ux0rb(i, 1, 1) = ux0rb(i, 1, 1) + ramp*ux0b(i, j, 1)
      ux0b(i, j, 1) = 0.0_8
      gy0rb(i, 1, 1) = gy0rb(i, 1, 1) + ramp*gy0b(i, j, 1)
      gy0b(i, j, 1) = 0.0_8
      gx0rb(i, 1, 1) = gx0rb(i, 1, 1) + ramp*gx0b(i, j, 1)
      gx0b(i, j, 1) = 0.0_8
    END DO
    CALL POPREAL8(ramp)
  END DO
!
  DO j=jyc,1,-1
    DO itime=ditime,1,-1
      by1b(j, itime) = 0.0_8
      bx1b(j, itime) = 0.0_8
    END DO
    DO i=nx,1,-1
      uy0b(i, j, 1) = 0.0_8
      ux0b(i, j, 1) = 0.0_8
      gy0b(i, j, 1) = 0.0_8
      gx0b(i, j, 1) = 0.0_8
    END DO
  END DO
!
END SUBROUTINE RECONSTRUYE4_B
