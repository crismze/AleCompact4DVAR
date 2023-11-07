!        Generated by TAPENADE     (INRIA, Tropics team)
!  Tapenade 3.4 (r3342) - 27 Jan 2010 15:25
!
!  Differentiation of funcion5 in forward (tangent) mode:
!   variations   of useful results: us2d2
!   with respect to varying inputs: us2d1
!   RW status of diff variables: us2d1:in us2d2:out
!
!  Entrada: us2d1d
!  Salida: us2d2d
! 
!********************************************************************
!
SUBROUTINE FUNCION5_D(us2d1d, us2d2d, mx, my, nx, ny, nxm&
&  , nym)
!
  implicit none
!
  INTEGER :: i, j, mx, my, nx, ny, nxm, nym
!
  REAL*8, DIMENSION(mx, my) :: us2d1, us2d2
  REAL*8, DIMENSION(mx, my) :: us2d1d, us2d2d
  REAL*8, DIMENSION(mx, my) :: ur1, ue1, w1, w11
  REAL*8, DIMENSION(mx, my) :: ur1d, ue1d, w1d, w11d
  REAL*8 :: pi, xx1, xx2, xx3, xx4
  REAL*8 :: xx1d, xx2d, xx3d, xx4d
  REAL*8 :: arg1
  REAL*8 :: arg2
  REAL :: result1
  INTRINSIC COS
  INTRINSIC SIN
  INTRINSIC ACOS
  INTRINSIC SQRT
!
! ******************* temporellebis16STR ********************
! *************** p^d(mx,my) ---> p^(mx,my) *****************
!
  pi = ACOS(-1.)
!
  DO j=1,my
    DO i=1,mx
      ur1d(i, j) = 0.0_8
      ur1(i, j) = 0.
      ue1d(i, j) = 0.0_8
      ue1(i, j) = 0.
    END DO
  END DO
  w1d = 0.0_8
!
  DO j=1,ny
    DO i=1,nxm/2+1
      w1d(i, j) = us2d1d(i, j)
      w1(i, j) = us2d1(i, j)
    END DO
  END DO
  w11d = 0.0_8
!
  DO j=1,ny
    DO i=1,nxm/2+1
      w11d(i, j) = us2d1d(nxm-i+2, j)
      w11(i, j) = us2d1(nxm-i+2, j)
    END DO
  END DO
  ue1d = 0.0_8
  ur1d = 0.0_8
!
  DO j=1,ny
    DO i=1,nxm/2+1
      arg1 = (i-1)/2.*pi/nxm
      arg2 = (nxm-i+1)/2.*pi/nxm
      ur1d(i, j) = COS(arg1)*w1d(i, j) + COS(arg2)*w11d(i, j)
      ur1(i, j) = w1(i, j)*COS(arg1) + w11(i, j)*COS(arg2)
      arg1 = (i-1)/2.*pi/nxm
      arg2 = (nxm-i+1)/2.*pi/nxm
      ue1d(i, j) = SIN(arg1)*w1d(i, j) - SIN(arg2)*w11d(i, j)
      ue1(i, j) = w1(i, j)*SIN(arg1) - w11(i, j)*SIN(arg2)
    END DO
  END DO
!
  DO j=1,ny
    DO i=2,nxm/2
      ur1d(i+nxm/2, j) = ur1d(nxm/2-i+2, j)
      ur1(i+nxm/2, j) = ur1(nxm/2-i+2, j)
    END DO
  END DO
!
  DO j=1,ny
    DO i=2,nxm/2
      ue1d(i+nxm/2, j) = -ue1d(nxm/2-i+2, j)
      ue1(i+nxm/2, j) = -ue1(nxm/2-i+2, j)
    END DO
  END DO
  us2d2d = 0.0_8
!
  DO j=1,ny
    DO i=1,nxm/2+1
      arg1 = (j-1)/2.*pi/nym
      arg2 = (nym-j+1)/2.*pi/nym
      xx1d = COS(arg1)*ur1d(i, j) - COS(arg2)*ue1d(i, j)
      xx1 = ur1(i, j)*COS(arg1) - ue1(i, j)*COS(arg2)
      arg1 = (j-1)/2.*pi/nym
      arg2 = (j-1)/2.*pi/nym
      xx2d = SIN(arg1)*ur1d(i, nym-j+2) + COS(arg2)*ue1d(i, nym-j+2)
      xx2 = ur1(i, nym-j+2)*SIN(arg1) + ue1(i, nym-j+2)*COS(arg2)
      arg1 = (nym-j+1)/2.*pi/nym
      arg2 = (j-1)/2.*pi/nym
      xx3d = SIN(arg2)*ue1d(i, nym-j+2) - SIN(arg1)*ur1d(i, nym-j+2)
      xx3 = -(ur1(i, nym-j+2)*SIN(arg1)) + ue1(i, nym-j+2)*SIN(arg2)
      arg1 = (j-1)/2.*pi/nym
      arg2 = (j-1)/2.*pi/nym
      xx4d = SIN(arg1)*ur1d(i, j) + COS(arg2)*ue1d(i, j)
      xx4 = ur1(i, j)*SIN(arg1) + ue1(i, j)*COS(arg2)
      us2d2d(2*i-1, j) = xx1d + xx2d
      us2d2(2*i-1, j) = xx1 + xx2
      us2d2d(2*i, j) = xx3d + xx4d
      us2d2(2*i, j) = xx3 + xx4
    END DO
  END DO
!
  DO i=1,mx
    us2d2d(i, ny) = 0.0_8
    us2d2(i, ny) = 0.
  END DO
!
  DO i=1,mx
    us2d2d(i, ny) = us2d2d(i, (ny+1)/2)
    us2d2(i, ny) = us2d2(i, (ny+1)/2)
  END DO
!
! (CN). Asi recupero la entrada de csfft2d corresp. a funcion6
  DO i=3,nxm
    DO j=1,my
      result1 = SQRT(2.)
      us2d2d(i, j) = result1*us2d2d(i, j)
      us2d2(i, j) = result1*us2d2(i, j)
    END DO
  END DO
!
END SUBROUTINE FUNCION5_D