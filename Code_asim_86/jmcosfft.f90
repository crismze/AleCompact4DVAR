! $Header: /home/teuler/cvsroot/lib/jmscfft.f90,v 6.4 2000/02/22 17:25:26 teuler Exp $
! JMFFTLIB : A library of portable fourier transform subroutines
!            emulating Cray SciLib
! Author   : Jean-Marie Teuler, CNRS-IDRIS (teuler@idris.fr)
!
! Permission is granted to copy and distribute this file or modified
! versions of this file for no fee, provided the copyright notice and
! this permission notice are preserved on all copies.

subroutine cosfft(isign,n,scale,x,y,table,work,isys)

  implicit none

  ! Arguments
  integer, intent(in) :: isign
  integer, intent(in) :: n
  real(kind=8), intent(in) :: scale
  real(kind=8), intent(in), dimension(0:n-1) :: x
  real(kind=8), intent(out), dimension(0:2*(n/2)+1) :: y
  real(kind=8), intent(inout), dimension(0:100+2*n-1) :: table
  real(kind=8), intent(inout), dimension(0:2*n-1) :: work
  integer, intent(in) :: isys

  ! Variables locales
  real :: fac1, fac2
  integer :: i
  integer :: ntable, nwork
  integer :: nfact
  integer, dimension(0:99) :: fact
  integer :: dimx, debx, incx, jumpx, jump1x
  integer :: dimy, deby, incy, jumpy
  character(len=*), parameter :: nomsp = 'SCFFT'
  real(kind=8), dimension (0:n+1) :: a2

  ! Positionnement a 0 du code de retour
  call jmsetcode(0)

  ! Verification des conditions
  if (isign /= 0 .and. isign /=-1 .and. isign /= 1) &
  & call jmerreur1(nomsp,2,isign)
  if (n < 1) call jmerreur1(nomsp,23,n)
  if (mod(n,2) /= 0) call jmerreur1(nomsp,24,n)

  ! Gestion de table
  ntable = 100+2*n

  ! Gestion de work
  nwork = 2*n
  !=======================================
  dimx = n         ; debx = 0 ; incx = 1 ; jumpx = 0 ; jump1x = n+1
  dimy = 2*(n/2)+2 ; deby = 0 ; incy = 1 ; jumpy = 0
!
  ! Test sur isign
  if (isign == 0) then
    ! Pour la factorisation
    call jmfact(n,fact,100,    0,nfact)
    table(0:nfact-1) = fact(0:nfact-1)
    ! Pour les sinus et cosinus
    call jmtable(table,ntable,100+0,n)
    return
  else
    nfact = nint(table(0))
    fact(0:nfact-1) = nint(table(0:nfact-1))
  end if
  !
  !=====================================
  if (isign == -1) then
     fac1=1./n
     fac2=0.25
  else
     fac1=1.
     fac2=.5
  endif
  !
  do j=0,n+1
     jj=(j)*jumpx+1
     work(0)=x(jj)*fac1
     work(n)=x(jj+n*incx)*fac1
     work(1)=0.
     work(n+1)=0.
     !
     sum = 0.
     do i=1,n-1,2
        sum=sum+x((i)*incx+jj)
     enddo
     a2(j)=2.*fac1*sum
     do i=2,n-1,2
        work(i)=fac1*x(jj+i*incx)
        work(i+1)=fac1*x(jj+(i-1)*inc-x(jj+(i+1)*incx))
     enddo
     do i=0,jump1-1
        x(jj+i*incx)=work(i+1)
     enddo
  !=====================================

  ! On appelle le sous-programme principal

  call jmscm1dxy(1,n,fact,100,0,table,ntable,100+0,work,nwork, &
  & x,dimx,debx,incx,jumpx,y,dimy,deby,incy,jumpy,isign,scale)

  !======================================
  do j=0,n-1
     jj=j*jump+1
     do i=0,jump1-1
        work(i)=x(jj+i*incx)
     enddo
     do i=1,n-1
        ii=i*incx+jj
        x(ii)=fac2*((work(i+1)+work(n+3-i) - &
          0.5*(work(i+1)-work(n+3-i))/trig   
  enddo
  !======================================

end subroutine cosfft
