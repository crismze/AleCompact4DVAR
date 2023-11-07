
subroutine SLFFT2D(x,y,mx,my,nxm,nym,nwork,ntrigsX,ntrigsY,isign)
!
  implicit none
!
  integer, intent(in) :: mx,my,nxm,nym,nwork,ntrigsX,ntrigsY
  integer, intent(in) :: isign
  real(8), dimension(mx,my) :: x,y
  real(8), dimension(nwork) :: work
  real(8), dimension(mx,my) :: x1,x2
  real(8), dimension(0:ntrigsX-1) :: trigsX
  real(8), dimension(0:ntrigsY-1) :: trigsY
  integer, parameter :: nfax = 19
  integer, dimension(0:nfax-1) :: ifaxX, ifaxY
  integer :: inc, jump
  integer :: i, j
!
  call fftfax(nxm,ifaxX,trigsX)
  call cftfax(nym,ifaxY,trigsY)
!
  if (isign==-1) then
     inc = 1
     jump = mx
     call rfftmlt(x,work,trigsX,ifaxX,inc,jump,nxm,nym,isign)
     inc = mx
     jump = 1
     do j=1,my
     do i=1,mx/2
        x1(i,j)=x(2*i-1,j)
        x2(i,j)=x(2*i,j)
     enddo
     enddo
     call cfftmlt(x1,x2,work,trigsY,ifaxY,inc,jump,nym,nxm,isign)
     do j=1,my
     do i=1,mx/2
        x(2*i-1,j)=x1(i,j)
        x(2*i,j)=x2(i,j)
     enddo
     enddo
     do j=1,my
     do i=1,mx
        y(i,j)=x(i,j)/nym
     enddo
     enddo
  endif
!
  if (isign==1) then
     inc = mx
     jump = 1
     do j=1,my
     do i=1,mx/2
        x1(i,j)=x(2*i-1,j)
        x2(i,j)=x(2*i,j)
     enddo
     enddo
     call cfftmlt(x1,x2,work,trigsY,ifaxY,inc,jump,nym,nxm,isign)
     do j=1,my
     do i=1,mx/2+2
        x(2*i-1,j)=x1(i,j)
        x(2*i,j)=x2(i,j)
     enddo
     enddo
     inc = 1
     jump = mx
     call rfftmlt(x,work,trigsX,ifaxX,inc,jump,nxm,nym,isign)
     do j=1,my
     do i=1,mx
        y(i,j)=x(i,j)
     enddo
     enddo
  endif
!
end subroutine SLFFT2D
