subroutine SLFFT3D(x,nx,ny,nz,mx,my,mz,nxm,nym,mzm,nwork,ntrigsX,ntrigsY,ntrigsZ,clx,ncly,nclz,isign)

  implicit none

  ! Arguments
  integer, intent(in) :: nx,ny,nz,mx,my,mz,nxm,nym,nzm,nwork,ntrigsX,ntrigsY,ntrigsZ
  integer, intent(in) :: nclx, ncly, nclz, isign
  real(8),  intent(inout), dimension(mx,my,mz) :: x

  ! Variables locales
  ! Le tableau de travail
  real(8), dimension(nwork) :: work

  ! Les tableaux pour les factorisations et les sinus et cosinus
  real(8), dimension(0:ntrigsX-1) :: trigsX
  real(8), dimension(0:ntrigsY-1) :: trigsY
  real(8), dimension(0:ntrigsZ-1) :: trigsZ
  integer, parameter :: nfax = 19
  integer, dimension(0:nfax-1) :: ifaxX, ifaxY, ifaxZ

  ! Autres variables locales
  integer :: inc, jump
  integer :: i, j, k
  real(8) :: pi
 
  ! La constante pi
  pi=acos(-1.0_8)
!************************FFT3D***(0-0-0)***********************************
  if ((nclx==0).and.(ncly==0).and.(nclz==0)) then
!
  ! Initialisation des factorisations et tableaux trigonométrique
  ! Note : L'appel a jmcftfax calcule les cos et sin(2*pi/n), tandis
  ! que les appels suivants calculent les cos(pi/n)
  ! Il y a un peu de gaspillage, mais c'est plus simple comme ca
  call fftfax(nxm,ifaxX,trigsX)
  call fftfax(nym,ifaxY,trigsY)
  call fftfax(nzm,ifaxZ,trigsZ)
!
  ! On applique une TF reelle->complexe 
  inc = 1
  jump = mx
  call rfftmlt(x,work,trigsX,ifaxX,inc,jump,nx,ny,isign)

  ! On applique une TF reelle->complexe
  inc = mx
  jump = 1
  call rfftmlt(x,work,trigsY,ifaxY,inc,jump,ny,nx,isign)

  ! On applique une TF reelle->complexe
  inc = mx
  jump = 1
  call rfftmlt(x,work,trigsY,ifaxY,inc,jump,ny,nx,isign)
!
  endif
!*************************************************************************

!
end subroutine SLFFT3D
