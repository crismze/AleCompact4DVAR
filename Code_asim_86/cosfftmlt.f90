subroutine cosfftmlt(x,work,trigs,ifax,inc,jump,n,m)

  implicit none

  ! Arguments
  integer, intent(in) ::  m, n
  real(8), dimension(0:(m-1)*jump+n*inc):: x
  real(8), dimension(0:2*m*n+m+m*(n+2)-1) :: work
  real(8), dimension(0:4*n-1) :: trigs
  integer, dimension(0:18) :: ifax
  integer :: inc, jump

  ! Variables locales
  integer :: ntrigs, nwork
  real(8) :: pi
  character(len=*), parameter :: nomsp = 'COSFFTMLT'
  integer :: incx, jumpx
  integer :: i, j, isign

  ! Gestion de pi
  pi = acos(-1.0_8)

  ! Gestion de table
  ntrigs = 4*n

  ! Gestion de work (dimension pour jmrfftmlt)
  nwork = 2*m*n

  ! On calcule par anticipation le premier terme de rang impair
  do j = 0, m-1
    work(nwork+j) = 0.5_8 * ( x(j*jump) - x(j*jump+n*inc) )
  end do
  do i = 1, n-1
    do j = 0, m-1
      work(nwork+j) = work(nwork+j) + trigs(3*n+i) * x(j*jump+i*inc)
    end do
  end do

  ! On prepare le tableau d'entree
  do i = 0, n-1
    do j = 0, m-1
      work(nwork+m+i+j*(n+2)) = &
         0.5_8 * (x(j*jump+i*inc) + x(j*jump+(n-i)*inc) ) &
        - trigs(2*n+i) * ( x(j*jump+i*inc) - x(j*jump+(n-i)*inc) )
    end do
  end do

  ! On appelle le sous-programme de transformee de Fourier
  isign = -1
  incx = 1
  jumpx = n+2
  call rfftmlt(work(nwork+m),work,trigs,ifax,incx,jumpx,n,m,isign)

  ! On reconstitue x
  ! Note : Il faut tenir compte des particularites de rfftmlt, qui met un
  !        facteur 1/n par defaut et qui prend une exponentielle negative
  !        Ceci ne s'applique pas bien sur au terme sauvegarde

  ! Traitement des indices pairs
  do i = 0, n, 2
    do j = 0, m-1
!      x(j*jump+i*inc) =work(nwork+m+i+2*j*(n/2+1))
       x(j*jump+i*inc) =n*work(nwork+m+i+2*j*(n/2+1))
    end do
  end do

  ! Traitement des indices impairs
  do j = 0, m-1
    x(j*jump+inc) = work(nwork+j)
  end do
  do i = 3, n, 2
    do j = 0, m-1
!      x(j*jump+i*inc) =x(j*jump+(i-2)*inc) - work(nwork+m+i+2*j*(n/2+1))
     x(j*jump+i*inc) = x(j*jump+(i-2)*inc) - n*work(nwork+m+i+2*j*(n/2+1))
    end do
  end do

end subroutine cosfftmlt
