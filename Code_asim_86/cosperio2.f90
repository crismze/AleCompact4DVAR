! $Header: /home/teuler/cvsroot/lib/jmscfft2d.f90,v 6.5 2000/02/22 17:25:26 teuler Exp $
! JMFFTLIB : A library of portable fourier transform subroutines
!            emulating Cray SciLib
! Author   : Jean-Marie Teuler, CNRS-IDRIS (teuler@idris.fr)
!
! Permission is granted to copy and distribute this file or modified
! versions of this file for no fee, provided the copyright notice and
! this permission notice are preserved on all copies.

subroutine cosperio(isign,a,d,mx,my,vlen,table,work,isys)

  implicit none

  ! Arguments
  integer, intent(in) :: isign, vlen
  integer, intent(in) :: mx, my
  real(kind=8), intent(in), dimension(mx-1,my) :: a
  real(kind=8), intent(out), dimension(mx-1,my) :: d
  real(kind=8), intent(inout), dimension(100+2*(mx+my)) :: table
  real(kind=8), intent(inout), dimension(4*vlen*max(mx,my)) :: work
  integer, intent(in) :: isys
!---------------------------------------------------------
  !Variables locales
  real(kind=8), dimension(my,mx) :: b,c
  integer :: ifail, i, j,jdeb,jfin,nwork_temp
  integer :: nwork, nx,ny
  integer :: debx, deby, jumpx, jumpy, incx, incy
  character(len=1) :: init
  logical :: debut, fin,nxpaire,nypaire
!---------------------------------------------------------

  debx=0
  deby=0
  jumpx=nx
  jumpy=2*(nx/2+1)
  incx=my
  incy=1
  d(:,:)=a(:,:)
  init = 'i'
  call mycos(1,mx-2,a,init,debx,deby,jumpx,jumpy,incx,incy,table,work,ifail)
!
 do j=1,my
 do i=1,mx
   print '(i4,i4,f15.8)',i,j,a(j,i)
 enddo        
 pause
 enddo
 
!
  debx=nx
  deby=0
  jumpx=ny
  jumpy=2*(ny/2+1)
  incx=1
  incy=1
 call scfftm(0,ny,1,1.,b,ny,c,ny,table,work,0,debx,deby,jumpx,jumpy,incx,incy)
 call scfftm(1,ny,1,1./ny,b,ny,c,ny,table,work,0,debx,deby,jumpx,jumpy,incx,incy) 
  do j=1,ny
  do i=1,nx
    a(i,j)=a(i,j)
!*sqrt(2.)*2./(nx-1)
    print *,a(i,j),d(i,j)
 enddo        
 pause
  enddo
!---------------------------------------------------------  

!---------------------------------------------------------

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
