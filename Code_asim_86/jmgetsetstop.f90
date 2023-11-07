! $Header: /home/teuler/cvsroot/lib/jmgetsetstop.f90,v 6.2 2000/03/10 11:57:58 teuler Exp $
! JMFFTLIB : A library of portable fourier transform subroutines
!            emulating Cray SciLib
! Author   : Jean-Marie Teuler, CNRS-IDRIS (teuler@idris.fr)
!
! Permission is granted to copy and distribute this file or modified
! versions of this file for no fee, provided the copyright notice and
! this permission notice are preserved on all copies.

subroutine jmgetsetstop(arret,type)

  ! Subroutine qui permet de stocker une valeur statique
  ! Ceci evite de recourir a un common ...

  implicit none

  ! Arguments
  integer, intent(inout) :: arret
  character(len=1), intent(in) :: type

  ! Variables locales

  ! Variable statique
  integer, save :: arret_last = 1

  if (type == 's') then
    arret_last = arret
  else if (type == 'g') then 
    arret = arret_last
  end if

end subroutine jmgetsetstop
