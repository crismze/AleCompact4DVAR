! $Header: /home/teuler/cvsroot/lib/jmsetstop.f90,v 6.2 2000/03/10 11:57:59 teuler Exp $
! JMFFTLIB : A library of portable fourier transform subroutines
!            emulating Cray SciLib
! Author   : Jean-Marie Teuler, CNRS-IDRIS (teuler@idris.fr)
!
! Permission is granted to copy and distribute this file or modified
! versions of this file for no fee, provided the copyright notice and
! this permission notice are preserved on all copies.

subroutine jmsetstop(arret)

  implicit none

  ! Arguments
  integer, intent(in) :: arret

  ! Variables locales
  integer :: arret2

  arret2 = arret
  call jmgetsetstop(arret2,'s')

end subroutine jmsetstop
