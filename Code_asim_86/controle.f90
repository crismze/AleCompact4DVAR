!**********************************************************************
!
subroutine control (ut1,ut2,t)
!
!**********************************************************************
!
   USE controle_m
   USE ecoulement_m
!
  implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  real(8)  :: ut1 
  real(8)  :: ut2 
  real(8)  :: t 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer :: i
  real(8) :: pi
!-----------------------------------------------
!
   pi=acos(-1.)
   uampl=0.5
   if (t.le.5.) then
      ut1=u2+uampl*(1.+cos(pi*t/5.+pi))
   else
      ut1=u2+uampl*2.
   endif
!
   return
end subroutine control
