! $Header: /home/teuler/cvsroot/lib/jmscfft2d.f90,v 6.5 2000/02/22 17:25:26 teuler Exp $
! JMFFTLIB : A library of portable fourier transform subroutines
!            emulating Cray SciLib
! Author   : Jean-Marie Teuler, CNRS-IDRIS (teuler@idris.fr)
!
! Permission is granted to copy and distribute this file or modified
! versions of this file for no fee, provided the copyright notice and
! this permission notice are preserved on all copies.

subroutine cosperio(isign,a,d,nx,mx,ny,my,vlen,table,work,isys)

  implicit none

  ! Arguments
  integer, intent(in) :: isign, vlen
  integer, intent(in) :: nx, ny, mx, my
  real(kind=8), intent(in), dimension(nx,my) :: a
  real(kind=8), intent(out), dimension(nx,my) :: d
  real(kind=8), intent(inout), dimension(100+2*(nx+ny)) :: table
  real(kind=8), intent(inout), dimension(4*vlen*max(nx,ny)) :: work
  integer, intent(in) :: isys
!---------------------------------------------------------
  !Variables locales
  real(kind=8), dimension(my,nx) :: b,c
  integer :: ifail, i, j,jdeb,jfin,nwork_temp
  integer :: nwork
  integer :: debx, deby
  character(len=1) :: init
  logical :: debut, fin,nxpaire,nypaire
!---------------------------------------------------------
  nxpaire= (mod(nx,2)==0)
  nypaire= (mod(ny,2)==0)

  debx=nx
  deby=0
  init = 'i'
  call mycos(1, nx-1, a, init, debx, deby, table, work, ifail)
!
  do j=1,ny
  do i=1,nx
    a(i,j)=a(i,j)
!*sqrt(2.)*2./(nx-1)
    print *,a(i,j)
 enddo        
 pause
  enddo
!---------------------------------------------------------  
  do j=1,my
  do i=1,nx
    b(j,i)=a(i,j)
  enddo
  enddo
!---------------------------------------------------------
 call scfftm(0, ny, 1, 1.   , b, ny, c, ny, table, work, 0)
 call scfftm(1, ny, 1, 1./ny   , b, ny, c, ny, table, work, 0)
!---------------------------------------------------------
 do i=1,nx
  do j=1,ny
     print *,c(j,i)
  enddo
  pause
  enddo
 do j=1,my
 do i=1,nx
    d(i,j)=c(j,i)
    print '(i4,i4,f15.8)',i,j,d(i,j)
 enddo
 pause
 enddo
!---------------------------------------------------------

end subroutine cosperio
