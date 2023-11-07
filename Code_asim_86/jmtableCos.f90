! $Header: /home/teuler/cvsroot/lib/jmtable.f90,v 6.2 2000/03/10 09:44:45 teuler Exp $
! JMFFTLIB : A library of portable fourier transform subroutines
!            emulating Cray SciLib
! Author   : Jean-Marie Teuler, CNRS-IDRIS (teuler@idris.fr)
!
! Permission is granted to copy and distribute this file or modified
! versions of this file for no fee, provided the copyright notice and
! this permission notice are preserved on all copies.

subroutine jmtableCos(table,ntable,itable,n)

  implicit none

  ! Arguments
  integer, intent(in) :: ntable, itable, n
  real(kind=8), intent(out), dimension(0:ntable-1) :: table

  ! Variables locales
  real(kind=8) :: pi

  integer :: i

  pi = acos(-1.0_8)

  ! Calcul des cosinus (necessaire pour les transformee en sinus et en cosinus)
  do i = 0,n-1
    table(itable+i) = cos(pi*real(i,kind=8)/real(n,kind=8))
  end do

end subroutine jmtableCos
