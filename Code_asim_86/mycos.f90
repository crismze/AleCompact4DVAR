!***********************************
  subroutine mycos(m,n,x,init,debx,deby,jumpx,jumpy,incx,incy,table,work,ifail)

  implicit none

  ! Arguments
  integer, intent(in) ::  m, n
  integer, intent(in) :: debx, deby, jumpx, jumpy, incx, incy
  real(8), intent(inout), dimension( 0:m-1, 0:n ) :: x
  character( len=1 ), intent( in ) :: init
  real(8), dimension( 0:100+4*n-1 ), intent( inout ) :: table
  real(8), dimension( 0:m*(3*n+7)-1 ), intent( inout ) :: work
  integer, intent(inout) :: ifail

  ! Variables locales
  integer, save :: lastn = -1
  logical, save :: first = .true.
  integer :: ntable, nwork
  real(8), save :: pi
  integer :: isign
  real(8) :: scale
  integer :: i, j, nfact, ltoto
  integer ,dimension(0:99) :: fact
  integer :: dimx
  integer :: dimy
  character(len=*), parameter :: nomsp = 'SCFFTM'

  ! Gestion de pi
  if ( first ) pi = acos( -1.0_8 )

  ! Gestion de table
  ntable = 100 + 4*n

  ! Gestion de work
  nwork = m*(2*n+4)

  ! Gestion de ifail
  ifail = 0

  ! Test sur M et N
  if (m < 1) then
    ifail = 1
    return
  else if (n < 1) then
    ifail = 2
    return
  end if

  ! Test sur isign
  if (init == 'i' .or. init == 'I') then

    first = .false.
    lastn = n

    ! On initialise les coefficients trigonometriques
    isign = 0

  ! Positionnement a 0 du code de retour
  call jmsetcode(0)

  ! Gestion de table
  ntable = 100+2*n
  ! Gestion de work
  nwork = 2*n*m

  ! Test sur isign

    ! Pour la factorisation
    call jmfact(n,fact,100,    0,nfact)
    table(0:nfact-1) = fact(0:nfact-1)
    ! Pour les sinus et cosinus
     call jmtable(table,ntable,100+0,n)
!------------------------------------------------------------------------

    ! Il en manque...
!dir$ ivdep
!ocl novrec
!CDIR NODEP
    do i = 0, n-1
      table( 100 + 2*n + i ) = cos( ( pi * i ) / n )
      table( 100 + 3*n + i ) = sin( ( pi * i ) / n )
    end do

  else if (init == 's' .or. init == 'S') then

    ! On verifie que le precedent appel etait bon
    if (first) then
      ifail = 4
      return
    else if (lastn /= n) then
      ifail = 5
      return
    end if

  else if (init == 'r' .or. init == 'R') then

   ! On admet que l'on recupere le bon "table" (avec fact et trig)
   ! On ne fait donc aucune verification

  else

    ! Valeurs erronnees de ifail
    ifail = 3
    return

  end if

  ! On calcule par anticipation le premier terme de rang impair
!dir$ ivdep
!ocl novrec
!CDIR NODEP
  do j = 0, m-1
    work( nwork + j ) = 0.5 * ( x( j, 0 ) - x( j, n ) )
  end do
  do i = 1, n-1
    do j = 0, m-1
      work( nwork + j ) = work( nwork + j ) + &
                          table( 100 + 2*n + i ) * x( j, i )
    end do
  end do

  ! On prepare le tableau d'entree
  ! Note : On utilise work( nwork+m:nwork+m+m*n-1)
  ! Attention : en plus, il faut transposer !
  do i = 0, n-1
!dir$ ivdep
!ocl novrec
!CDIR NODEP
    do j = 0, m-1
      work( nwork + m + i+j*n ) = 0.5_8 * ( x( j, i ) + x( j, n-i ) )    &
                                  - table( 100 + 3*n + i ) * ( x( j, i ) &
                                  - x( j, n-i ) )
    end do
  end do

  ! On appelle le sous-programme de transformee de Fourier
  isign = 1
  scale = sqrt( 2._8 / real(n,kind=8) )
!
  nfact = nint(table(0))
  fact(0:nfact-1) = nint(table(0:nfact-1))
!
  ! On appelle le sous-programme principal
  dimx = n*m   
  dimy = 2*(n/2+1)*m
  call jmscm1dxy(m,n,fact,100,0,table,ntable,100+0,work,nwork, &
  & work(nwork+m),dimx,debx,incx,jumpx,work(nwork+m),dimy,deby,incy,jumpy,isign,scale)
!----------------------------------------------------------------------
  ! On reconstitue x
  ! Traitement des indices pairs
  do j = 0, m-1
!dir$ ivdep
!ocl novrec
!CDIR NODEP
    do i = 0, n, 2
      x( j, i ) = work( nwork + m + i + 2*j*(n/2+1) )
    end do
  end do

  ! Traitement des indices impairs
!dir$ ivdep
!ocl novrec
!CDIR NODEP
  do j = 0, m-1
    x( j, 1 ) = scale*work( nwork + j )
  end do
  do j = 0, m-1
    do i = 3, n, 2
      x( j, i ) = x( j, i-2 ) + work( nwork + m + i + 2*j*(n/2+1) )
    end do
  end do

end subroutine mycos
