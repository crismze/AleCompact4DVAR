! $Header: /home/teuler/cvsroot/lib/jmgetstop.f90,v 6.2 2000/03/10 11:57:58 teuler Exp $
! JMFFTLIB : A library of portable fourier transform subroutines
!            emulating Cray SciLib
! Author   : Jean-Marie Teuler, CNRS-IDRIS (teuler@idris.fr)
!
! Permission is granted to copy and distribute this file or modified
! versions of this file for no fee, provided the copyright notice and
! this permission notice are preserved on all copies.

subroutine jmgetstop(arret)

  implicit none

  ! Arguments
  integer, intent(out) :: arret

  ! Variables locales

  call jmgetsetstop(arret,'g')

end subroutine jmgetstop
