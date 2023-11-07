subroutine sinfftmlt(x,work,trigs,ifax,inc,jump,n,m)

  implicit none

  ! Arguments
  integer, intent(in) ::  m, n
  real(8), dimension(0:(m-1)*jump+(n-1)*inc):: x
  real(8), dimension(0:2*m*n+m*(n+2)-1) :: work
  real(8), dimension(0:3*n-1) :: trigs
  integer, dimension(0:18) :: ifax
  integer :: inc, jump

  ! Variables locales
  integer :: ntrigs, nwork
  real(8) :: pi
  character(len=*), parameter :: nomsp = 'SINFFTMLT'
  integer :: incx, jumpx, isign
  integer :: i, j
  real(8) :: s, t, u

  ! Gestion de pi
  pi = acos(-1.0_8)

  ! Gestion de table
  ntrigs = 3*n

  ! Gestion de work (dimension pour jmrfftmlt)
  nwork = 2*m*n

  ! On prepare le tableau d'entree
  work(:) = 0
  do j = 0,m-1
    work(nwork+j*(n+2)) = 0
  end do
  do i = 1,n-1
    do j = 0,m-1
      s = trigs(2*n+i)
      t = x(j*jump+(i-1)*inc)
      u = x(j*jump+(n-i-1)*inc)
      work(nwork+i+j*(n+2)) = s*(t+u)+0.5_8*(t-u)
    end do
  end do

  ! On appelle le sous-programme de transformee de Fourier
  isign = -1
  incx = 1
  jumpx = n+2
  call rfftmlt(work(nwork),work,trigs,ifax,incx,jumpx,n,m,isign)

  ! On reconstitue x
  ! Note : Il faut tenir compte des particularites de rfftmlt, qui met un
  !        facteur 1/n par defaut et qui prend une exponentielle negative
  !        Ceci ne s'applique pas bien sur au terme sauvegarde

  ! D'abord les indices impairs
  do i = 0,n/2-1
    do j = 0,m-1
      x(j*jump+(2*i+1)*inc) = -n*work(nwork+2*i+3+2*j*(n/2+1))
    end do
  end do

  ! Ensuite les indices pairs
  ! Cas particulier indice 0
  do j = 0,m-1
    x(j*jump) = n*work(nwork+2*j*(n/2+1))/2.0_8
  end do
 ! Cas general : recurrence
  do i = 1,n/2-1
    do j = 0,m-1
      x(j*jump+2*i*inc) = x(j*jump+(2*i-2)*inc) + n*work(nwork+j*2*(n/2+1)+2*i)
    end do
  end do

end subroutine sinfftmlt
