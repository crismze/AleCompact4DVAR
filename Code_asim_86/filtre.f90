!*********************************************************************
!
subroutine filx(tx,ux,rx,sx,vx,fiffx,fifx,ficx,fibx,fibbx,filax,fiz1x,&
     fiz2x,nx,ny,nz,npaire) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
  USE paramx_m 
  USE parfix_m 
!...Translated by PSUITE Trans90                  4.3ZH  9:53:12   1/27/ 4  
!...Switches: -yv INDDO=0 -nbejkno
  implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  integer , intent(in) :: nx 
  integer , intent(in) :: ny 
  integer , intent(in) :: nz 
  integer , intent(in) :: npaire 
  real(8) , intent(inout) :: tx(nx,ny,nz) 
  real(8) , intent(in) :: ux(nx,ny,nz) 
  real(8) , intent(inout) :: rx(nx,ny,nz) 
  real(8) , intent(inout) :: sx(ny,nz) 
  real(8) , intent(inout) :: vx(ny,nz) 
  real(8) , intent(in) :: fiffx(nx) 
  real(8) , intent(in) :: fifx(nx) 
  real(8) , intent(in) :: ficx(nx) 
  real(8) , intent(in) :: fibx(nx) 
  real(8) , intent(in) :: fibbx(nx) 
  real(8) , intent(in) :: filax(nx,2) 
  real(8) , intent(in) :: fiz1x(nx) 
  real(8) , intent(in) :: fiz2x(nx) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer :: k, j, i 
  real(8) :: xcoef 
!-----------------------------------------------
!
!*********************************************************************
!
   if (nclx==0) then
      do k=1,nz
      do j=1,ny
         rx(1,j,k)=fiaix*ux(1,j,k)+&
                   fibix*(ux(2,j,k)+ux(nx,j,k))+&
                   ficix*(ux(3,j,k)+ux(nx-1,j,k))+&
                   fidix*(ux(4,j,k)+ux(nx-2,j,k)) 
         rx(2,j,k)=fiaix*ux(2,j,k)+&
                   fibix*(ux(3,j,k)+ux(1,j,k))+&
                   ficix*(ux(4,j,k)+ux(nx,j,k))+&
                   fidix*(ux(5,j,k)+ux(nx-1,j,k)) 
         rx(3,j,k)=fiaix*ux(3,j,k)+&
                   fibix*(ux(4,j,k)+ux(2,j,k))+&
                   ficix*(ux(5,j,k)+ux(1,j,k))+&
                   fidix*(ux(6,j,k)+ux(nx,j,k)) 
      enddo
      enddo
!
      do i=4,nx-3
      do k=1,nz
      do j=1,ny
         rx(i,j,k)=fiaix*ux(i,j,k)+&
                   fibix*(ux(i+1,j,k)+ux(i-1,j,k))+& 
                   ficix*(ux(i+2,j,k)+ux(i-2,j,k))+& 
                   fidix*(ux(i+3,j,k)+ux(i-3,j,k)) 
      enddo
      enddo
      enddo
!
      do k=1,nz
      do j=1,ny
         rx(nx,j,k)=fiaix*ux(nx,j,k)+&
                    fibix*(ux(1,j,k)+ux(nx-1,j,k))+& 
                    ficix*(ux(2,j,k)+ux(nx-2,j,k))+&
                    fidix*(ux(3,j,k)+ux(nx-3,j,k)) 
         rx(nx-1,j,k)=fiaix*ux(nx-1,j,k)+&
                      fibix*(ux(nx,j,k)+ux(nx-2,j,k))+&
                      ficix*(ux(1,j,k)+ux(nx-3,j,k))+&
                      fidix*(ux(2,j,k)+ux(nx-4,j,k)) 
         rx(nx-2,j,k)=fiaix*ux(nx-2,j,k)+&
                      fibix*(ux(nx-1,j,k)+ux(nx-3,j,k))+& 
                      ficix*(ux(nx,j,k)+ux(nx-4,j,k))+& 
                      fidix*(ux(1,j,k)+ux(nx-5,j,k)) 
      enddo
      enddo
!
      do k=1,nz
      do j=1,ny
      do i=1,nx-2
         rx(i+1,j,k)=rx(i+1,j,k)-filax(i,1)*rx(i,j,k)
         rx(i+2,j,k)=rx(i+2,j,k)-filax(i,2)*rx(i,j,k)
      enddo
      enddo
      enddo
!
      do k=1,nz
      do j=1,ny
         rx(nx,j,k)=rx(nx,j,k)-filax(nx-1,1)*rx(nx-1,j,k)
      enddo
      enddo
!
      do k=1,nz
      do j=1,ny
         rx(nx,j,k)=rx(nx,j,k)*fiffx(nx)
         rx(nx-1,j,k)=(rx(nx-1,j,k)-fifx(nx-1)*rx(nx,j,k))*&
                      fiffx(nx-1)
         rx(nx-2,j,k)=(rx(nx-2,j,k)-fifx(nx-2)*rx(nx-1,j,k)-&
                      ficx(nx-2)*rx(nx,j,k))*fiffx(nx-2)
         rx(nx-3,j,k)=(rx(nx-3,j,k)-fifx(nx-3)*rx(nx-2,j,k)-&
                       ficx(nx-3)*rx(nx-1,j,k)-&
                       fibx(nx-3)*rx(nx,j,k))*fiffx(nx-3)
      enddo
      enddo
      do k=1,nz
      do j=1,ny
      do i=nx-4,1,-1
         rx(i,j,k)=(rx(i,j,k)-fifx(i)*rx(i+1,j,k)-&
                    ficx(i)*rx(i+2,j,k)-&
                    fibx(i)*rx(i+3,j,k)-&
                    fibbx(i)*rx(i+4,j,k))*fiffx(i)
      enddo
      enddo
      enddo
!
      xcoef=1./2.
!
      do k=1,nz
      do j=1,ny
         sx(j,k)=fih1x*(-fibex*rx(1,j,k)+fibex*rx(nx-1,j,k)*xcoef+&
                    fialx*rx(nx,j,k)*xcoef)+&
             fih2x*(fialx*rx(1,j,k)*xcoef+fibex*rx(2,j,k)*xcoef-&
                    fibex*rx(nx,j,k))
         vx(j,k)=fih3x*(-fibex*rx(1,j,k)+fibex*rx(nx-1,j,k)*xcoef+&
                    fialx*rx(nx,j,k)*xcoef)+&
             fih4x*(fialx*rx(1,j,k)*xcoef+fibex*rx(2,j,k)*xcoef-&
                    fibex*rx(nx,j,k))
      enddo
      enddo
      do k=1,nz
      do j=1,ny
      do i=1,nx
         tx(i,j,k)=rx(i,j,k)-fiz1x(i)*sx(j,k)-fiz2x(i)*vx(j,k)
      enddo
      enddo
      enddo
   endif
!
   if (nclx==1) then
      if (npaire==1) then
         do k=1,nz
         do j=1,ny
            tx(1,j,k)=fiaix*ux(1,j,k)+&
                      fibix*(ux(2,j,k)+ux(2,j,k))+&
                      ficix*(ux(3,j,k)+ux(3,j,k))+&
                      fidix*(ux(4,j,k)+ux(4,j,k))
            tx(2,j,k)=fiaix*ux(2,j,k)+&
                      fibix*(ux(3,j,k)+ux(1,j,k))+&
                      ficix*(ux(4,j,k)+ux(2,j,k))+&
                      fidix*(ux(5,j,k)+ux(3,j,k))
            tx(3,j,k)=fiaix*ux(3,j,k)+&
                      fibix*(ux(4,j,k)+ux(2,j,k))+&
                      ficix*(ux(5,j,k)+ux(1,j,k))+&
                      fidix*(ux(6,j,k)+ux(2,j,k))
         enddo
         enddo
         do i=4,nx-3
         do k=1,nz
         do j=1,ny
            tx(i,j,k)=fiaix*ux(i,j,k)+&
                      fibix*(ux(i+1,j,k)+ux(i-1,j,k))+&
                      ficix*(ux(i+2,j,k)+ux(i-2,j,k))+&
                      fidix*(ux(i+3,j,k)+ux(i-3,j,k))
         enddo
         enddo
         enddo
         do k=1,nz
         do j=1,ny
            tx(nx,j,k)=fiaix*ux(nx,j,k)+&
                       fibix*(ux(nx-1,j,k)+ux(nx-1,j,k))+&
                       ficix*(ux(nx-2,j,k)+ux(nx-2,j,k))+&
                       fidix*(ux(nx-3,j,k)+ux(nx-3,j,k))
            tx(nx-1,j,k)=fiaix*ux(nx-1,j,k)+&
                          fibix*(ux(nx,j,k)+ux(nx-2,j,k))+&
                          ficix*(ux(nx-1,j,k)+ux(nx-3,j,k))+&
                          fidix*(ux(nx-2,j,k)+ux(nx-4,j,k))
            tx(nx-2,j,k)=fiaix*ux(nx-2,j,k)+&
                         fibix*(ux(nx-1,j,k)+ux(nx-3,j,k))+&
                         ficix*(ux(nx,j,k)+ux(nx-4,j,k))+&
                         fidix*(ux(nx-1,j,k)+ux(nx-5,j,k))
         enddo
         enddo
!
         do k=1,nz
         do j=1,ny
            do i=1,nx-2
               tx(i+1,j,k)=tx(i+1,j,k)-filax(i,1)*tx(i,j,k)
               tx(i+2,j,k)=tx(i+2,j,k)-filax(i,2)*tx(i,j,k)
            enddo
            tx(nx,j,k)=tx(nx,j,k)-filax(nx-1,1)*tx(nx-1,j,k)
         enddo
         enddo
!
         do k=1,nz
         do j=1,ny
            tx(nx,j,k)=tx(nx,j,k)*fiffx(nx)
            tx(nx-1,j,k)=(tx(nx-1,j,k)-fifx(nx-1)*tx(nx,j,k))*&
                         fiffx(nx-1)
            tx(nx-2,j,k)=(tx(nx-2,j,k)-fifx(nx-2)*tx(nx-1,j,k)-&
                         ficx(nx-2)*tx(nx,j,k))*fiffx(nx-2)
            tx(nx-3,j,k)=(tx(nx-3,j,k)-fifx(nx-3)*tx(nx-2,j,k)-&
                          ficx(nx-3)*tx(nx-1,j,k)-&
                          fibx(nx-3)*tx(nx,j,k))*fiffx(nx-3)
            do i=nx-4,1,-1
               tx(i,j,k)=(tx(i,j,k)-fifx(i)*tx(i+1,j,k)-&
                          ficx(i)*tx(i+2,j,k)-&
                          fibx(i)*tx(i+3,j,k)-&
                          fibbx(i)*tx(i+4,j,k))*fiffx(i)
            enddo
         enddo
         enddo
      endif
      if (npaire==0) then
         do k=1,nz
         do j=1,ny
            tx(1,j,k)=fiaix*ux(1,j,k)+&
                      fibix*(ux(2,j,k)-ux(2,j,k))+&
                      ficix*(ux(3,j,k)-ux(3,j,k))+&
                      fidix*(ux(4,j,k)-ux(4,j,k))
            tx(2,j,k)=fiaix*ux(2,j,k)+&
                      fibix*(ux(3,j,k)+ux(1,j,k))+&
                      ficix*(ux(4,j,k)-ux(2,j,k))+&
                      fidix*(ux(5,j,k)-ux(3,j,k))
            tx(3,j,k)=fiaix*ux(3,j,k)+&
                      fibix*(ux(4,j,k)+ux(2,j,k))+&
                      ficix*(ux(5,j,k)+ux(1,j,k))+&
                      fidix*(ux(6,j,k)-ux(2,j,k))
         enddo
         enddo
         do i=4,nx-3
         do k=1,nz
         do j=1,ny
            tx(i,j,k)=fiaix*ux(i,j,k)+&
                      fibix*(ux(i+1,j,k)+ux(i-1,j,k))+&
                      ficix*(ux(i+2,j,k)+ux(i-2,j,k))+&
                      fidix*(ux(i+3,j,k)+ux(i-3,j,k))
         enddo
         enddo
         enddo
         do k=1,nz
         do j=1,ny
            tx(nx,j,k)=fiaix*ux(nx,j,k)+&
                       fibix*(ux(nx-1,j,k)-ux(nx-1,j,k))+&
                       ficix*(ux(nx-2,j,k)-ux(nx-2,j,k))+&
                       fidix*(ux(nx-3,j,k)-ux(nx-3,j,k))
            tx(nx-1,j,k)=fiaix*ux(nx-1,j,k)+&
                         fibix*(ux(nx,j,k)+ux(nx-2,j,k))+&
                         ficix*(-ux(nx-1,j,k)+ux(nx-3,j,k))+&
                         fidix*(-ux(nx-2,j,k)+ux(nx-4,j,k))
            tx(nx-2,j,k)=fiaix*ux(nx-2,j,k)+&
                         fibix*(ux(nx-1,j,k)+ux(nx-3,j,k))+&
                         ficix*(ux(nx,j,k)+ux(nx-4,j,k))+&
                         fidix*(-ux(nx-1,j,k)+ux(nx-5,j,k))
         enddo
         enddo
!
         do k=1,nz
         do j=1,ny
            do i=1,nx-2
               tx(i+1,j,k)=tx(i+1,j,k)-filax(i,1)*tx(i,j,k)
               tx(i+2,j,k)=tx(i+2,j,k)-filax(i,2)*tx(i,j,k)
            enddo
            tx(nx,j,k)=tx(nx,j,k)-filax(nx-1,1)*tx(nx-1,j,k)
         enddo
         enddo
!
         do k=1,nz
         do j=1,ny
            tx(nx,j,k)=tx(nx,j,k)*fiffx(nx)
            tx(nx-1,j,k)=(tx(nx-1,j,k)-fifx(nx-1)*tx(nx,j,k))*&
                         fiffx(nx-1)
            tx(nx-2,j,k)=(tx(nx-2,j,k)-fifx(nx-2)*tx(nx-1,j,k)-&
                         ficx(nx-2)*tx(nx,j,k))*fiffx(nx-2)
            tx(nx-3,j,k)=(tx(nx-3,j,k)-fifx(nx-3)*tx(nx-2,j,k)-&
                          ficx(nx-3)*tx(nx-1,j,k)-&
                          fibx(nx-3)*tx(nx,j,k))*fiffx(nx-3)
            do i=nx-4,1,-1
               tx(i,j,k)=(tx(i,j,k)-fifx(i)*tx(i+1,j,k)-&
                          ficx(i)*tx(i+2,j,k)-&
                          fibx(i)*tx(i+3,j,k)-&
                          fibbx(i)*tx(i+4,j,k))*fiffx(i)
            enddo
         enddo
         enddo
      endif
   endif
!
   if (nclx==2) then
      do k=1,nz
      do j=1,ny
         tx(1,j,k)=fia1x*ux(1,j,k)+fib1x*ux(2,j,k)+&
                   fic1x*ux(3,j,k)+fid1x*ux(4,j,k)+&
                   fie1x*ux(5,j,k)
         tx(2,j,k)=fia2x*ux(2,j,k)+fib2x*ux(1,j,k)+&
                   fic2x*ux(3,j,k)+fid2x*ux(4,j,k)+&
                   fie2x*ux(5,j,k)
         tx(3,j,k)=fia3x*ux(3,j,k)+fib3x*ux(1,j,k)+&
                   fic3x*ux(2,j,k)+fid3x*ux(4,j,k)+&
                   fie3x*ux(5,j,k)
      enddo
      enddo
      do i=4,nx-3
      do k=1,nz
      do j=1,ny
         tx(i,j,k)=fiaix*ux(i,j,k)+&
                   fibix*(ux(i+1,j,k)+ux(i-1,j,k))+& 
                   ficix*(ux(i+2,j,k)+ux(i-2,j,k))+& 
                   fidix*(ux(i+3,j,k)+ux(i-3,j,k)) 
      enddo
      enddo
      enddo
      do k=1,nz
      do j=1,ny
         tx(nx,j,k)=fianx*ux(nx,j,k)+fibnx*ux(nx-1,j,k)+&
                    ficnx*ux(nx-2,j,k)+fidnx*ux(nx-3,j,k)+&
                    fienx*ux(nx-4,j,k)
         tx(nx-1,j,k)=fiamx*ux(nx-1,j,k)+fibmx*ux(nx,j,k)+&
                      ficmx*ux(nx-2,j,k)+fidmx*ux(nx-3,j,k)+&
                      fiemx*ux(nx-4,j,k)
         tx(nx-2,j,k)=fiapx*ux(nx-2,j,k)+fibpx*ux(nx,j,k)+&
                      ficpx*ux(nx-1,j,k)+fidpx*ux(nx-3,j,k)+&
                      fiepx*ux(nx-4,j,k)
      enddo
      enddo
!
      do i=1,nx-2
      do k=1,nz
      do j=1,ny
         tx(i+1,j,k)=tx(i+1,j,k)-filax(i,1)*tx(i,j,k)
         tx(i+2,j,k)=tx(i+2,j,k)-filax(i,2)*tx(i,j,k)
      enddo
      enddo
      enddo
!
      do k=1,nz
      do j=1,ny
         tx(nx,j,k)=tx(nx,j,k)-filax(nx-1,1)*tx(nx-1,j,k)
      enddo
      enddo
!
      do k=1,nz
      do j=1,ny
         tx(nx,j,k)=tx(nx,j,k)*fiffx(nx)
         tx(nx-1,j,k)=(tx(nx-1,j,k)-fifx(nx-1)*tx(nx,j,k))*&
                      fiffx(nx-1)
         tx(nx-2,j,k)=(tx(nx-2,j,k)-fifx(nx-2)*tx(nx-1,j,k)-&
                      ficx(nx-2)*tx(nx,j,k))*fiffx(nx-2)
         tx(nx-3,j,k)=(tx(nx-3,j,k)-fifx(nx-3)*tx(nx-2,j,k)-&
                       ficx(nx-3)*tx(nx-1,j,k)-&
                       fibx(nx-3)*tx(nx,j,k))*fiffx(nx-3)
      enddo
      enddo
      do i=nx-4,1,-1
      do k=1,nz
      do j=1,ny
         tx(i,j,k)=(tx(i,j,k)-fifx(i)*tx(i+1,j,k)-&
                    ficx(i)*tx(i+2,j,k)-&
                    fibx(i)*tx(i+3,j,k)-&
                    fibbx(i)*tx(i+4,j,k))*fiffx(i)
      enddo
      enddo
      enddo
!
   endif
!
   return  
end subroutine filx
!
!*********************************************************************
!
subroutine fily(ty,uy,ry,sy,vy,fiffy,fify,ficy,fiby,fibby,filay,fiz1y,&
     fiz2y,nx,ny,nz,npaire) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
  USE paramy_m 
  USE parfiy_m 
!...Translated by PSUITE Trans90                  4.3ZH  9:53:12   1/27/ 4  
!...Switches: -yv INDDO=0 -nbejkno
  implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  integer , intent(in) :: nx 
  integer , intent(in) :: ny 
  integer , intent(in) :: nz 
  integer , intent(in) :: npaire 
  real(8) , intent(inout) :: ty(nx,ny,nz) 
  real(8) , intent(in) :: uy(nx,ny,nz) 
  real(8) , intent(inout) :: ry(nx,ny,nz) 
  real(8) , intent(inout) :: sy(nx,nz) 
  real(8) , intent(inout) :: vy(nx,nz) 
  real(8) , intent(in) :: fiffy(ny) 
  real(8) , intent(in) :: fify(ny) 
  real(8) , intent(in) :: ficy(ny) 
  real(8) , intent(in) :: fiby(ny) 
  real(8) , intent(in) :: fibby(ny) 
  real(8) , intent(in) :: filay(ny,2) 
  real(8) , intent(in) :: fiz1y(ny) 
  real(8) , intent(in) :: fiz2y(ny) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer :: k, i, j 
  real(8) :: xcoef 
!-----------------------------------------------
!
!*********************************************************************
!
   if (ncly==0) then
      do k=1,nz
      do i=1,nx
         ry(i,1,k)=fiaiy*uy(i,1,k)+&
                   fibiy*(uy(i,2,k)+uy(i,ny,k))+& 
                   ficiy*(uy(i,3,k)+uy(i,ny-1,k))+& 
                   fidiy*(uy(i,4,k)+uy(i,ny-2,k)) 
         ry(i,2,k)=fiaiy*uy(i,2,k)+&
                   fibiy*(uy(i,3,k)+uy(i,1,k))+& 
                   ficiy*(uy(i,4,k)+uy(i,ny,k))+ &
                   fidiy*(uy(i,5,k)+uy(i,ny-1,k)) 
         ry(i,3,k)=fiaiy*uy(i,3,k)+&
                   fibiy*(uy(i,4,k)+uy(i,2,k))+&
                   ficiy*(uy(i,5,k)+uy(i,1,k))+& 
                   fidiy*(uy(i,6,k)+uy(i,ny,k)) 
      enddo
      enddo
      do j=4,ny-3
      do k=1,nz
      do i=1,nx
         ry(i,j,k)=fiaiy*uy(i,j,k)+&
                   fibiy*(uy(i,j+1,k)+uy(i,j-1,k))+& 
                   ficiy*(uy(i,j+2,k)+uy(i,j-2,k))+& 
                   fidiy*(uy(i,j+3,k)+uy(i,j-3,k)) 
      enddo
      enddo
      enddo
      do k=1,nz
      do i=1,nx
         ry(i,ny,k)=fiaiy*uy(i,ny,k)+&
                    fibiy*(uy(i,1,k)+uy(i,ny-1,k))+& 
                    ficiy*(uy(i,2,k)+uy(i,ny-2,k))+& 
                    fidiy*(uy(i,3,k)+uy(i,ny-3,k)) 
         ry(i,ny-1,k)=fiaiy*uy(i,ny-1,k)+&
                      fibiy*(uy(i,ny,k)+uy(i,ny-2,k))+& 
                      ficiy*(uy(i,1,k)+uy(i,ny-3,k))+& 
                      fidiy*(uy(i,2,k)+uy(i,ny-4,k)) 
         ry(i,ny-2,k)=fiaiy*uy(i,ny-2,k)+&
                      fibiy*(uy(i,ny-1,k)+uy(i,ny-3,k))+& 
                      ficiy*(uy(i,ny,k)+uy(i,ny-4,k))+& 
                      fidiy*(uy(i,1,k)+uy(i,ny-5,k)) 
      enddo
      enddo
!
      do k=1,nz
      do i=1,nx
      do j=1,ny-2
         ry(i,j+1,k)=ry(i,j+1,k)-filay(j,1)*ry(i,j,k)
         ry(i,j+2,k)=ry(i,j+2,k)-filay(j,2)*ry(i,j,k)
      enddo
      enddo
      enddo
      do k=1,nz
      do i=1,nx
         ry(i,ny,k)=ry(i,ny,k)-filay(ny-1,1)*ry(i,ny-1,k)
      enddo
      enddo
!
      do k=1,nz
      do i=1,nx
         ry(i,ny,k)=ry(i,ny,k)*fiffy(ny)
         ry(i,ny-1,k)=(ry(i,ny-1,k)-fify(ny-1)*ry(i,ny,k))*&
                      fiffy(ny-1)
         ry(i,ny-2,k)=(ry(i,ny-2,k)-fify(ny-2)*ry(i,ny-1,k)-&
                      ficy(ny-2)*ry(i,ny,k))*fiffy(ny-2)
         ry(i,ny-3,k)=(ry(i,ny-3,k)-fify(ny-3)*ry(i,ny-2,k)-&
                       ficy(ny-3)*ry(i,ny-1,k)-&
                       fiby(ny-3)*ry(i,ny,k))*fiffy(ny-3)
      enddo
      enddo
      do k=1,nz
      do i=1,nx
      do j=ny-4,1,-1
         ry(i,j,k)=(ry(i,j,k)-fify(j)*ry(i,j+1,k)-&
                    ficy(j)*ry(i,j+2,k)-&
                    fiby(j)*ry(i,j+3,k)-&
                    fibby(j)*ry(i,j+4,k))*fiffy(j)
      enddo
      enddo
      enddo
!
      xcoef=1./2.
!
      do k=1,nz
      do i=1,nx
         sy(i,k)=fih1y*(-fibey*ry(i,1,k)+fibey*ry(i,ny-1,k)*xcoef+&
                    fialy*ry(i,ny,k)*xcoef)+&
             fih2y*(fialy*ry(i,1,k)*xcoef+fibey*ry(i,2,k)*xcoef-&
                    fibey*ry(i,ny,k))
         vy(i,k)=fih3y*(-fibey*ry(i,1,k)+fibey*ry(i,ny-1,k)*xcoef+&
                    fialy*ry(i,ny,k)*xcoef)+&
             fih4y*(fialy*ry(i,1,k)*xcoef+fibey*ry(i,2,k)*xcoef-&
                    fibey*ry(i,ny,k))
      enddo
      enddo
      do k=1,nz
      do j=1,ny
      do i=1,nx
         ty(i,j,k)=ry(i,j,k)-fiz1y(j)*sy(i,k)-fiz2y(j)*vy(i,k)
      enddo
      enddo
      enddo
   endif
!
   if (ncly==1) then
      if (npaire==1) then
         do k=1,nz
         do i=1,nx
            ty(i,1,k)=fiaiy*uy(i,1,k)+&
                      fibiy*(uy(i,2,k)+uy(i,2,k))+& 
                      ficiy*(uy(i,3,k)+uy(i,3,k))+& 
                      fidiy*(uy(i,4,k)+uy(i,4,k)) 
            ty(i,2,k)=fiaiy*uy(i,2,k)+&
                      fibiy*(uy(i,3,k)+uy(i,1,k))+& 
                      ficiy*(uy(i,4,k)+uy(i,2,k))+& 
                      fidiy*(uy(i,5,k)+uy(i,3,k)) 
            ty(i,3,k)=fiaiy*uy(i,3,k)+&
                      fibiy*(uy(i,4,k)+uy(i,2,k))+& 
                      ficiy*(uy(i,5,k)+uy(i,1,k))+& 
                      fidiy*(uy(i,6,k)+uy(i,2,k)) 
         enddo
         enddo
         do j=4,ny-3
         do k=1,nz
         do i=1,nx
            ty(i,j,k)=fiaiy*uy(i,j,k)+&
                      fibiy*(uy(i,j+1,k)+uy(i,j-1,k))+& 
                      ficiy*(uy(i,j+2,k)+uy(i,j-2,k))+& 
                      fidiy*(uy(i,j+3,k)+uy(i,j-3,k)) 
         enddo
         enddo
         enddo
         do k=1,nz
         do i=1,nx
            ty(i,ny,k)=fiaiy*uy(i,ny,k)+&
                       fibiy*(uy(i,ny-1,k)+uy(i,ny-1,k))+&
                       ficiy*(uy(i,ny-2,k)+uy(i,ny-2,k))+&
                       fidiy*(uy(i,ny-3,k)+uy(i,ny-3,k))
            ty(i,ny-1,k)=fiaiy*uy(i,ny-1,k)+&
                         fibiy*(uy(i,ny,k)+uy(i,ny-2,k))+&
                         ficiy*(uy(i,ny-1,k)+uy(i,ny-3,k))+&
                         fidiy*(uy(i,ny-2,k)+uy(i,ny-4,k))
            ty(i,ny-2,k)=fiaiy*uy(i,ny-2,k)+&
                         fibiy*(uy(i,ny-1,k)+uy(i,ny-3,k))+&
                         ficiy*(uy(i,ny,k)+uy(i,ny-4,k))+&
                         fidiy*(uy(i,ny-1,k)+uy(i,ny-5,k))
         enddo
         enddo
         do k=1,nz
         do i=1,nx
            do j=1,ny-2
               ty(i,j+1,k)=ty(i,j+1,k)-filay(j,1)*ty(i,j,k)
               ty(i,j+2,k)=ty(i,j+2,k)-filay(j,2)*ty(i,j,k)
            enddo
            ty(i,ny,k)=ty(i,ny,k)-filay(ny-1,1)*ty(i,ny-1,k)
         enddo
         enddo
         do k=1,nz
         do i=1,nx
            ty(i,ny,k)=ty(i,ny,k)*fiffy(ny)
            ty(i,ny-1,k)=(ty(i,ny-1,k)-fify(ny-1)*ty(i,ny,k))*&
                         fiffy(ny-1)
            ty(i,ny-2,k)=(ty(i,ny-2,k)-fify(ny-2)*ty(i,ny-1,k)-&
                         ficy(ny-2)*ty(i,ny,k))*fiffy(ny-2)
            ty(i,ny-3,k)=(ty(i,ny-3,k)-fify(ny-3)*ty(i,ny-2,k)-&
                          ficy(ny-3)*ty(i,ny-1,k)-&
                          fiby(ny-3)*ty(i,ny,k))*fiffy(ny-3)
            do j=ny-4,1,-1
               ty(i,j,k)=(ty(i,j,k)-fify(j)*ty(i,j+1,k)-&
                          ficy(j)*ty(i,j+2,k)-&
                          fiby(j)*ty(i,j+3,k)-&
                          fibby(j)*ty(i,j+4,k))*fiffy(j)
            enddo
         enddo
         enddo
      endif
      if (npaire==0) then
         do k=1,nz
         do i=1,nx
            ty(i,1,k)=fiaiy*uy(i,1,k)+&
                      fibiy*(uy(i,2,k)-uy(i,2,k))+& 
                      ficiy*(uy(i,3,k)-uy(i,3,k))+& 
                      fidiy*(uy(i,4,k)-uy(i,4,k)) 
            ty(i,2,k)=fiaiy*uy(i,2,k)+&
                      fibiy*(uy(i,3,k)+uy(i,1,k))+& 
                      ficiy*(uy(i,4,k)-uy(i,2,k))+& 
                      fidiy*(uy(i,5,k)-uy(i,3,k)) 
            ty(i,3,k)=fiaiy*uy(i,3,k)+&
                      fibiy*(uy(i,4,k)+uy(i,2,k))+&
                      ficiy*(uy(i,5,k)+uy(i,1,k))+& 
                      fidiy*(uy(i,6,k)-uy(i,2,k)) 
         enddo
         enddo
         do j=4,ny-3
         do k=1,nz
         do i=1,nx
            ty(i,j,k)=fiaiy*uy(i,j,k)+&
                      fibiy*(uy(i,j+1,k)+uy(i,j-1,k))+& 
                      ficiy*(uy(i,j+2,k)+uy(i,j-2,k))+&
                      fidiy*(uy(i,j+3,k)+uy(i,j-3,k)) 
         enddo
         enddo
         enddo
         do k=1,nz
         do i=1,nx
            ty(i,ny,k)=fiaiy*uy(i,ny,k)+&
                       fibiy*(uy(i,ny-1,k)-uy(i,ny-1,k))+&
                       ficiy*(uy(i,ny-2,k)-uy(i,ny-2,k))+&
                       fidiy*(uy(i,ny-3,k)-uy(i,ny-3,k))
            ty(i,ny-1,k)=fiaiy*uy(i,ny-1,k)+&
                         fibiy*(uy(i,ny,k)+uy(i,ny-2,k))+&
                         ficiy*(-uy(i,ny-1,k)+uy(i,ny-3,k))+&
                         fidiy*(-uy(i,ny-2,k)+uy(i,ny-4,k))
            ty(i,ny-2,k)=fiaiy*uy(i,ny-2,k)+&
                         fibiy*(uy(i,ny-1,k)+uy(i,ny-3,k))+&
                         ficiy*(uy(i,ny,k)+uy(i,ny-4,k))+&
                         fidiy*(-uy(i,ny-1,k)+uy(i,ny-5,k))
         enddo
         enddo
         do k=1,nz
         do i=1,nx
            do j=1,ny-2
               ty(i,j+1,k)=ty(i,j+1,k)-filay(j,1)*ty(i,j,k)
               ty(i,j+2,k)=ty(i,j+2,k)-filay(j,2)*ty(i,j,k)
            enddo
            ty(i,ny,k)=ty(i,ny,k)-filay(ny-1,1)*ty(i,ny-1,k)
         enddo
         enddo
         do k=1,nz
         do i=1,nx
            ty(i,ny,k)=ty(i,ny,k)*fiffy(ny)
            ty(i,ny-1,k)=(ty(i,ny-1,k)-fify(ny-1)*ty(i,ny,k))*&
                         fiffy(ny-1)
            ty(i,ny-2,k)=(ty(i,ny-2,k)-fify(ny-2)*ty(i,ny-1,k)-&
                         ficy(ny-2)*ty(i,ny,k))*fiffy(ny-2)
            ty(i,ny-3,k)=(ty(i,ny-3,k)-fify(ny-3)*ty(i,ny-2,k)-&
                          ficy(ny-3)*ty(i,ny-1,k)-&
                          fiby(ny-3)*ty(i,ny,k))*fiffy(ny-3)
            do j=ny-4,1,-1
               ty(i,j,k)=(ty(i,j,k)-fify(j)*ty(i,j+1,k)-&
                          ficy(j)*ty(i,j+2,k)-&
                          fiby(j)*ty(i,j+3,k)-&
                          fibby(j)*ty(i,j+4,k))*fiffy(j)
            enddo
         enddo
         enddo
      endif
   endif
!
   if (ncly.eq.2) then
      do k=1,nz
      do i=1,nx
         ty(i,1,k)=fia1y*uy(i,1,k)+fib1y*uy(i,2,k)+&
                   fic1y*uy(i,3,k)+fid1y*uy(i,4,k)+&
                   fie1y*uy(i,5,k)
         ty(i,2,k)=fia2y*uy(i,2,k)+fib2y*uy(i,1,k)+&
                   fic2y*uy(i,3,k)+fid2y*uy(i,4,k)+&
                   fie2y*uy(i,5,k)
         ty(i,3,k)=fia3y*uy(i,3,k)+fib3y*uy(i,1,k)+&
                   fic3y*uy(i,2,k)+fid3y*uy(i,4,k)+&
                   fie3y*uy(i,5,k)
      enddo
      enddo
      do j=4,ny-3
      do k=1,nz
      do i=1,nx
         ty(i,j,k)=fiaiy*uy(i,j,k)+&
                   fibiy*(uy(i,j+1,k)+uy(i,j-1,k))+&
                   ficiy*(uy(i,j+2,k)+uy(i,j-2,k))+&
                   fidiy*(uy(i,j+3,k)+uy(i,j-3,k)) 
      enddo
      enddo
      enddo
      do k=1,nz
      do i=1,nx
         ty(i,ny,k)=fiany*uy(i,ny,k)+fibny*uy(i,ny-1,k)+&
                    ficny*uy(i,ny-2,k)+fidny*uy(i,ny-3,k)+&
                    fieny*uy(i,ny-4,k)
         ty(i,ny-1,k)=fiamy*uy(i,ny-1,k)+fibmy*uy(i,ny,k)+&
                      ficmy*uy(i,ny-2,k)+fidmy*uy(i,ny-3,k)+&
                      fiemy*uy(i,ny-4,k)
         ty(i,ny-2,k)=fiapy*uy(i,ny-2,k)+fibpy*uy(i,ny,k)+&
                      ficpy*uy(i,ny-1,k)+fidpy*uy(i,ny-3,k)+&
                      fiepy*uy(i,ny-4,k)
      enddo
      enddo
!
      do k=1,nz
      do i=1,nx
         do j=1,ny-2
            ty(i,j+1,k)=ty(i,j+1,k)-filay(j,1)*ty(i,j,k)
            ty(i,j+2,k)=ty(i,j+2,k)-filay(j,2)*ty(i,j,k)
         enddo
         ty(i,ny,k)=ty(i,ny,k)-filay(ny-1,1)*ty(i,ny-1,k)
      enddo
      enddo
!
      do k=1,nz
      do i=1,nx
         ty(i,ny,k)=ty(i,ny,k)*fiffy(ny)
         ty(i,ny-1,k)=(ty(i,ny-1,k)-fify(ny-1)*ty(i,ny,k))*&
                      fiffy(ny-1)
         ty(i,ny-2,k)=(ty(i,ny-2,k)-fify(ny-2)*ty(i,ny-1,k)-&
                      ficy(ny-2)*ty(i,ny,k))*fiffy(ny-2)
         ty(i,ny-3,k)=(ty(i,ny-3,k)-fify(ny-3)*ty(i,ny-2,k)-&
                       ficy(ny-3)*ty(i,ny-1,k)-&
                       fiby(ny-3)*ty(i,ny,k))*fiffy(ny-3)
         do j=ny-4,1,-1
            ty(i,j,k)=(ty(i,j,k)-fify(j)*ty(i,j+1,k)-&
                       ficy(j)*ty(i,j+2,k)-&
                       fiby(j)*ty(i,j+3,k)-&
                       fibby(j)*ty(i,j+4,k))*fiffy(j)
         enddo
      enddo
      enddo
!
   endif
!
   return  
end subroutine fily


!
!*********************************************************************
!
subroutine filz(tz,uz,rz,sz,vz,fiffz,fifz,ficz,fibz,fibbz,filaz,fiz1z,&
     fiz2z,nx,ny,nz,npaire) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
  USE paramz_m 
  USE parfiz_m 
!...Translated by PSUITE Trans90                  4.3ZH  9:53:12   1/27/ 4  
!...Switches: -yv INDDO=0 -nbejkno
  implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  integer , intent(in) :: nx 
  integer , intent(in) :: ny 
  integer , intent(in) :: nz 
  integer , intent(in) :: npaire 
  real(8) , intent(inout) :: tz(nx,ny,nz) 
  real(8) , intent(in) :: uz(nx,ny,nz) 
  real(8) , intent(inout) :: rz(nx,ny,nz) 
  real(8) , intent(inout) :: sz(nx,ny) 
  real(8) , intent(inout) :: vz(nx,ny) 
  real(8) , intent(in) :: fiffz(nz) 
  real(8) , intent(in) :: fifz(nz) 
  real(8) , intent(in) :: ficz(nz) 
  real(8) , intent(in) :: fibz(nz) 
  real(8) , intent(in) :: fibbz(nz) 
  real(8) , intent(in) :: filaz(nz,2) 
  real(8) , intent(in) :: fiz1z(nz) 
  real(8) , intent(in) :: fiz2z(nz) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer :: j, i, k 
  real(8) :: xcoef 
!-----------------------------------------------
!
!*********************************************************************
!
   if (nclz==0) then
      do j=1,ny
      do i=1,nx
         rz(i,j,1)=fiaiz*uz(i,j,1)+&
                   fibiz*(uz(i,j,2)+uz(i,j,nz))+&
                   ficiz*(uz(i,j,3)+uz(i,j,nz-1))+&
                   fidiz*(uz(i,j,4)+uz(i,j,nz-2)) 
         rz(i,j,2)=fiaiz*uz(i,j,2)+&
                   fibiz*(uz(i,j,3)+uz(i,j,1))+ &
                   ficiz*(uz(i,j,4)+uz(i,j,nz))+&
                   fidiz*(uz(i,j,5)+uz(i,j,nz-1)) 
         rz(i,j,3)=fiaiz*uz(i,j,3)+&
                   fibiz*(uz(i,j,4)+uz(i,j,2))+&
                   ficiz*(uz(i,j,5)+uz(i,j,1))+&
                   fidiz*(uz(i,j,6)+uz(i,j,nz)) 
      enddo
      enddo
      do k=4,nz-3
      do j=1,ny
      do i=1,nx
         rz(i,j,k)=fiaiz*uz(i,j,k)+&
                   fibiz*(uz(i,j,k+1)+uz(i,j,k-1))+&
                   ficiz*(uz(i,j,k+2)+uz(i,j,k-2))+&
                   fidiz*(uz(i,j,k+3)+uz(i,j,k-3)) 
      enddo
      enddo
      enddo
      do j=1,ny
      do i=1,nx
         rz(i,j,nz)=fiaiz*uz(i,j,nz)+&
                    fibiz*(uz(i,j,1)+uz(i,j,nz-1))+&
                    ficiz*(uz(i,j,2)+uz(i,j,nz-2))+&
                    fidiz*(uz(i,j,3)+uz(i,j,nz-3)) 
         rz(i,j,nz-1)=fiaiz*uz(i,j,nz-1)+&
                      fibiz*(uz(i,j,nz)+uz(i,j,nz-2))+&
                      ficiz*(uz(i,j,1)+uz(i,j,nz-3))+&
                      fidiz*(uz(i,j,2)+uz(i,j,nz-4)) 
         rz(i,j,nz-2)=fiaiz*uz(i,j,nz-2)+&
                      fibiz*(uz(i,j,nz-1)+uz(i,j,nz-3))+&
                      ficiz*(uz(i,j,nz)+uz(i,j,nz-4))+&
                      fidiz*(uz(i,j,1)+uz(i,j,nz-5)) 
      enddo
      enddo
!
      do j=1,ny
      do i=1,nx
         do k=1,nz-2
            rz(i,j,k+1)=rz(i,j,k+1)-filaz(k,1)*rz(i,j,k)
            rz(i,j,k+2)=rz(i,j,k+2)-filaz(k,2)*rz(i,j,k)
         enddo
         rz(i,j,nz)=rz(i,j,nz)-filaz(nz-1,1)*rz(i,j,nz-1)
      enddo
      enddo
!
      do j=1,ny
      do i=1,nx
         rz(i,j,nz)=rz(i,j,nz)*fiffz(nz)
         rz(i,j,nz-1)=(rz(i,j,nz-1)-fifz(nz-1)*rz(i,j,nz))*&
                      fiffz(nz-1)
         rz(i,j,nz-2)=(rz(i,j,nz-2)-fifz(nz-2)*rz(i,j,nz-1)-&
                      ficz(nz-2)*rz(i,j,nz))*fiffz(nz-2)
         rz(i,j,nz-3)=(rz(i,j,nz-3)-fifz(nz-3)*rz(i,j,nz-2)-&
                       ficz(nz-3)*rz(i,j,nz-1)-&
                       fibz(nz-3)*rz(i,j,nz))*fiffz(nz-3)
      enddo
      enddo
      do j=1,ny
      do i=1,nx
      do k=nz-4,1,-1
         rz(i,j,k)=(rz(i,j,k)-fifz(k)*rz(i,j,k+1)-&
                    ficz(k)*rz(i,j,k+2)-&
                    fibz(k)*rz(i,j,k+3)-&
                    fibbz(k)*rz(i,j,k+4))*fiffz(k)
      enddo
      enddo
      enddo
!
      xcoef=1./2.
!
      do j=1,ny
      do i=1,nx
         sz(i,j)=fih1z*(-fibez*rz(i,j,1)+fibez*rz(i,j,nz-1)*xcoef+&
                    fialz*rz(i,j,nz)*xcoef)+&
             fih2z*(fialz*rz(i,j,1)*xcoef+fibez*rz(i,j,2)*xcoef-&
                    fibez*rz(i,j,nz))
         vz(i,j)=fih3z*(-fibez*rz(i,j,1)+fibez*rz(i,j,nz-1)*xcoef+&
                    fialz*rz(i,j,nz)*xcoef)+&
             fih4z*(fialz*rz(i,j,1)*xcoef+fibez*rz(i,j,2)*xcoef-&
                    fibez*rz(i,j,nz))
      enddo
      enddo
      do k=1,nz
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=rz(i,j,k)-fiz1z(k)*sz(i,j)-fiz2z(k)*vz(i,j)
      enddo
      enddo
      enddo
   endif
!
   if (nclz==1) then
      if (npaire==1) then
         do j=1,ny
         do i=1,nx
            tz(i,j,1)=fiaiz*uz(i,j,1)+&
                      fibiz*(uz(i,j,2)+uz(i,j,2))+&
                      ficiz*(uz(i,j,3)+uz(i,j,3))+&
                      fidiz*(uz(i,j,4)+uz(i,j,4)) 
            tz(i,j,2)=fiaiz*uz(i,j,2)+&
                      fibiz*(uz(i,j,3)+uz(i,j,1))+&
                      ficiz*(uz(i,j,4)+uz(i,j,2))+&
                      fidiz*(uz(i,j,5)+uz(i,j,3)) 
            tz(i,j,3)=fiaiz*uz(i,j,3)+&
                      fibiz*(uz(i,j,4)+uz(i,j,2))+&
                      ficiz*(uz(i,j,5)+uz(i,j,1))+&
                      fidiz*(uz(i,j,6)+uz(i,j,2)) 
         enddo
         enddo
         do k=4,nz-3
         do j=1,ny
         do i=1,nx
            tz(i,j,k)=fiaiz*uz(i,j,k)+&
                      fibiz*(uz(i,j,k+1)+uz(i,j,k-1))+&
                      ficiz*(uz(i,j,k+2)+uz(i,j,k-2))+&
                      fidiz*(uz(i,j,k+3)+uz(i,j,k-3)) 
         enddo
         enddo
         enddo
         do j=1,ny
         do i=1,nx
            tz(i,j,nz)=fiaiz*uz(i,j,nz)+&
                      fibiz*(uz(i,j,nz-1)+uz(i,j,nz-1))+&
                      ficiz*(uz(i,j,nz-2)+uz(i,j,nz-2))+&
                      fidiz*(uz(i,j,nz-3)+uz(i,j,nz-3))
            tz(i,j,nz-1)=fiaiz*uz(i,j,nz-1)+&
                      fibiz*(uz(i,j,nz)+uz(i,j,nz-2))+&
                      ficiz*(uz(i,j,nz-1)+uz(i,j,nz-3))+&
                      fidiz*(uz(i,j,nz-2)+uz(i,j,nz-4))
            tz(i,j,nz-2)=fiaiz*uz(i,j,nz-2)+&
                      fibiz*(uz(i,j,nz-1)+uz(i,j,nz-3))+&
                      ficiz*(uz(i,j,nz)+uz(i,j,nz-4))+&
                      fidiz*(uz(i,j,nz-1)+uz(i,j,nz-5))
         enddo
         enddo
         do j=1,ny
         do i=1,nx
            do k=1,nz-2
               tz(i,j,k+1)=tz(i,j,k+1)-filaz(k,1)*tz(i,j,k)
               tz(i,j,k+2)=tz(i,j,k+2)-filaz(k,2)*tz(i,j,k)
            enddo
            tz(i,j,nz)=tz(i,j,nz)-filaz(nz-1,1)*tz(i,j,nz-1)
         enddo
         enddo
         do j=1,ny
         do i=1,nx
            tz(i,j,nz)=tz(i,j,nz)*fiffz(nz)
            tz(i,j,nz-1)=(tz(i,j,nz-1)-fifz(nz-1)*tz(i,j,nz))*&
                         fiffz(nz-1)
            tz(i,j,nz-2)=(tz(i,j,nz-2)-fifz(nz-2)*tz(i,j,nz-1)-&
                         ficz(nz-2)*tz(i,j,nz))*fiffz(nz-2)
            tz(i,j,nz-3)=(tz(i,j,nz-3)-fifz(nz-3)*tz(i,j,nz-2)-&
                          ficz(nz-3)*tz(i,j,nz-1)-&
                          fibz(nz-3)*tz(i,j,nz))*fiffz(nz-3)
            do k=nz-4,1,-1
               tz(i,j,k)=(tz(i,j,k)-fifz(k)*tz(i,j,k+1)-&
                          ficz(k)*tz(i,j,k+2)-&
                          fibz(k)*tz(i,j,k+3)-&
                          fibbz(k)*tz(i,j,k+4))*fiffz(k)
            enddo
         enddo
         enddo
      endif
      if (npaire==0) then
         do j=1,ny
         do i=1,nx
            tz(i,j,1)=fiaiz*uz(i,j,1)+&
                      fibiz*(uz(i,j,2)-uz(i,j,2))+&
                      ficiz*(uz(i,j,3)-uz(i,j,3))+&
                      fidiz*(uz(i,j,4)-uz(i,j,4)) 
            tz(i,j,2)=fiaiz*uz(i,j,2)+&
                      fibiz*(uz(i,j,3)+uz(i,j,1))+&
                      ficiz*(uz(i,j,4)-uz(i,j,2))+&
                      fidiz*(uz(i,j,5)-uz(i,j,3)) 
            tz(i,j,3)=fiaiz*uz(i,j,3)+&
                      fibiz*(uz(i,j,4)+uz(i,j,2))+&
                      ficiz*(uz(i,j,5)+uz(i,j,1))+&
                      fidiz*(uz(i,j,6)-uz(i,j,2)) 
         enddo
         enddo
         do k=4,nz-3
         do j=1,ny
         do i=1,nx
            tz(i,j,k)=fiaiz*uz(i,j,k)+&
                      fibiz*(uz(i,j,k+1)+uz(i,j,k-1))+&
                      ficiz*(uz(i,j,k+2)+uz(i,j,k-2))+&
                      fidiz*(uz(i,j,k+3)+uz(i,j,k-3)) 
         enddo
         enddo
         enddo
         do j=1,ny
         do i=1,nx
            tz(i,j,nz)=fiaiz*uz(i,j,nz)+&
                      fibiz*(uz(i,j,nz-1)-uz(i,j,nz-1))+&
                      ficiz*(uz(i,j,nz-2)-uz(i,j,nz-2))+&
                      fidiz*(uz(i,j,nz-3)-uz(i,j,nz-3))
            tz(i,j,nz-1)=fiaiz*uz(i,j,nz-1)+&
                      fibiz*(uz(i,j,nz)+uz(i,j,nz-2))+&
                      ficiz*(-uz(i,j,nz-1)+uz(i,j,nz-3))+&
                      fidiz*(-uz(i,j,nz-2)+uz(i,j,nz-4))
            tz(i,j,nz-2)=fiaiz*uz(i,j,nz-2)+&
                      fibiz*(uz(i,j,nz-1)+uz(i,j,nz-3))+&
                      ficiz*(uz(i,j,nz)+uz(i,j,nz-4))+&
                      fidiz*(-uz(i,j,nz-1)+uz(i,j,nz-5))
         enddo
         enddo
         do j=1,ny
         do i=1,nx
            do k=1,nz-2
               tz(i,j,k+1)=tz(i,j,k+1)-filaz(k,1)*tz(i,j,k)
               tz(i,j,k+2)=tz(i,j,k+2)-filaz(k,2)*tz(i,j,k)
            enddo
            tz(i,j,nz)=tz(i,j,nz)-filaz(nz-1,1)*tz(i,j,nz-1)
         enddo
         enddo
         do j=1,ny
         do i=1,nx
            tz(i,j,nz)=tz(i,j,nz)*fiffz(nz)
            tz(i,j,nz-1)=(tz(i,j,nz-1)-fifz(nz-1)*tz(i,j,nz))*&
                         fiffz(nz-1)
            tz(i,j,nz-2)=(tz(i,j,nz-2)-fifz(nz-2)*tz(i,j,nz-1)-&
                         ficz(nz-2)*tz(i,j,nz))*fiffz(nz-2)
            tz(i,j,nz-3)=(tz(i,j,nz-3)-fifz(nz-3)*tz(i,j,nz-2)-&
                          ficz(nz-3)*tz(i,j,nz-1)-&
                          fibz(nz-3)*tz(i,j,nz))*fiffz(nz-3)
            do k=nz-4,1,-1
               tz(i,j,k)=(tz(i,j,k)-fifz(k)*tz(i,j,k+1)-&
                          ficz(k)*tz(i,j,k+2)-&
                          fibz(k)*tz(i,j,k+3)-&
                          fibbz(k)*tz(i,j,k+4))*fiffz(k)
            enddo
         enddo
         enddo
      endif
   endif
!
   if (nclz==2) then
      do j=1,ny
      do i=1,nx
         tz(i,j,1)=fia1z*uz(i,j,1)+fib1z*uz(i,j,2)+&
                   fic1z*uz(i,j,3)+fid1z*uz(i,j,4)+&
                   fie1z*uz(i,j,5)
         tz(i,j,2)=fia2z*uz(i,j,2)+fib2z*uz(i,j,1)+&
                   fic2z*uz(i,j,3)+fid2z*uz(i,j,4)+&
                   fie2z*uz(i,j,5)
         tz(i,j,3)=fia3z*uz(i,j,3)+fib3z*uz(i,j,1)+&
                   fic3z*uz(i,j,2)+fid3z*uz(i,j,4)+&
                   fie3z*uz(i,j,5)
      enddo
      enddo
      do k=4,nz-3
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=fiaiz*uz(i,j,k)+&
                   fibiz*(uz(i,j,k+1)+uz(i,j,k-1))+&
                   ficiz*(uz(i,j,k+2)+uz(i,j,k-2))+&
                   fidiz*(uz(i,j,k+3)+uz(i,j,k-3)) 
      enddo
      enddo
      enddo
      do j=1,ny
      do i=1,nx
         tz(i,j,nz)=fianz*uz(i,j,nz)+fibnz*uz(i,j,nz-1)+&
                    ficnz*uz(i,j,nz-2)+fidnz*uz(i,j,nz-3)+&
                    fienz*uz(i,j,nz-4)
         tz(i,j,nz-1)=fiamz*uz(i,j,nz-1)+fibmz*uz(i,j,nz)+&
                      ficmz*uz(i,j,nz-2)+fidmz*uz(i,j,nz-3)+&
                      fiemz*uz(i,j,nz-4)
         tz(i,j,nz-2)=fiapz*uz(i,j,nz-2)+fibpz*uz(i,j,nz)+&
                      ficpz*uz(i,j,nz-1)+fidpz*uz(i,j,nz-3)+&
                      fiepz*uz(i,j,nz-4)
      enddo
      enddo
!
      do j=1,ny
      do i=1,nx
         do k=1,nz-2
            tz(i,j,k+1)=tz(i,j,k+1)-filaz(k,1)*tz(i,j,k)
            tz(i,j,k+2)=tz(i,j,k+2)-filaz(k,2)*tz(i,j,k)
         enddo
!
         tz(i,j,nz)=tz(i,j,nz)-filaz(nz-1,1)*tz(i,j,nz-1)
      enddo
      enddo
!
      do j=1,ny
      do i=1,nx
         tz(i,j,nz)=tz(i,j,nz)*fiffz(nz)
         tz(i,j,nz-1)=(tz(i,j,nz-1)-fifz(nz-1)*tz(i,j,nz))*&
                      fiffz(nz-1)
         tz(i,j,nz-2)=(tz(i,j,nz-2)-fifz(nz-2)*tz(i,j,nz-1)-&
                      ficz(nz-2)*tz(i,j,nz))*fiffz(nz-2)
         tz(i,j,nz-3)=(tz(i,j,nz-3)-fifz(nz-3)*tz(i,j,nz-2)-&
                       ficz(nz-3)*tz(i,j,nz-1)-&
                       fibz(nz-3)*tz(i,j,nz))*fiffz(nz-3)
         do k=nz-4,1,-1
            tz(i,j,k)=(tz(i,j,k)-fifz(k)*tz(i,j,k+1)-&
                       ficz(k)*tz(i,j,k+2)-&
                       fibz(k)*tz(i,j,k+3)-&
                       fibbz(k)*tz(i,j,k+4))*fiffz(k)
         enddo
      enddo
      enddo
!
   endif
!
   return  
end subroutine filz
!
!*********************************************************************
!
subroutine filtres(fiffx,fifx,ficx,fibx,fibbx,filax,fiffxp,fifxp,ficxp,&
     fibxp,fibbxp,filaxp,fiffy,fify,ficy,fiby,fibby,filay,fiffyp,fifyp,&
     ficyp,fibyp,fibbyp,filayp,fiffz,fifz,ficz,fibz,fibbz,filaz,fiffzp,&
     fifzp,ficzp,fibzp,fibbzp,filazp,fiz1x,fiz2x,fiz1y,fiz2y,fiz1z,fiz2z,&
     nx,ny,nz) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
  USE paramx_m 
  USE paramy_m 
  USE paramz_m 
  USE parfix_m 
  USE parfiy_m 
  USE parfiz_m 
!...Translated by PSUITE Trans90                  4.3ZH  9:53:12   1/27/ 4  
!...Switches: -yv INDDO=0 -nbejkno
  implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  integer  :: nx 
  integer  :: ny 
  integer  :: nz 
  real(8)  :: fiffx(nx) 
  real(8)  :: fifx(nx) 
  real(8)  :: ficx(nx) 
  real(8)  :: fibx(nx) 
  real(8)  :: fibbx(nx) 
  real(8)  :: filax(nx,2) 
  real(8)  :: fiffxp(nx) 
  real(8)  :: fifxp(nx) 
  real(8)  :: ficxp(nx) 
  real(8)  :: fibxp(nx) 
  real(8)  :: fibbxp(nx) 
  real(8)  :: filaxp(nx,2) 
  real(8)  :: fiffy(ny) 
  real(8)  :: fify(ny) 
  real(8)  :: ficy(ny) 
  real(8)  :: fiby(ny) 
  real(8)  :: fibby(ny) 
  real(8)  :: filay(ny,2) 
  real(8)  :: fiffyp(ny) 
  real(8)  :: fifyp(ny) 
  real(8)  :: ficyp(ny) 
  real(8)  :: fibyp(ny) 
  real(8)  :: fibbyp(ny) 
  real(8)  :: filayp(ny,2) 
  real(8)  :: fiffz(nz) 
  real(8)  :: fifz(nz) 
  real(8)  :: ficz(nz) 
  real(8)  :: fibz(nz) 
  real(8)  :: fibbz(nz) 
  real(8)  :: filaz(nz,2) 
  real(8)  :: fiffzp(nz) 
  real(8)  :: fifzp(nz) 
  real(8)  :: ficzp(nz) 
  real(8)  :: fibzp(nz) 
  real(8)  :: fibbzp(nz) 
  real(8)  :: filazp(nz,2) 
  real(8)  :: fiz1x(nx) 
  real(8)  :: fiz2x(nx) 
  real(8)  :: fiz1y(ny) 
  real(8)  :: fiz2y(ny) 
  real(8)  :: fiz1z(nz) 
  real(8)  :: fiz2z(nz) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer :: i, j, k 
  real(8) :: xcoef 
!-----------------------------------------------
!
!*********************************************************************
!
   call coefficients()
!
   if (nx>1) then
      xcoef=1./2.
      if (nclx==0) then
         fiffx(1)=0.
         fifx(1)=0.
         ficx(1)=1.+fibex/xcoef
         fibx(1)=fialx
         fibbx(1)=fibex
         fiffx(2)=0.
         fifx(2)=fialx+fialx*xcoef
         ficx(2)=1.+fibex*xcoef
         fibx(2)=fialx
         fibbx(2)=fibex
!
         do i=3,nx-2
            fiffx(i)=fibex
            fifx(i)=fialx
            ficx(i)=1.
            fibx(i)=fialx
            fibbx(i)=fibex
         enddo
!
         fiffx(nx-1)=fibex
         fifx(nx-1)=fialx
         ficx(nx-1)=1.+fibex*xcoef
         fibx(nx-1)=fialx+fialx*xcoef
         fibbx(nx-1)=0.
         fiffx(nx)=fibex
         fifx(nx)=fialx
         ficx(nx)=1.+fibex/xcoef
         fibx(nx)=0.
         fibbx(nx)=0.
      endif
!
      if (nclx==1) then
         fiffx(1)=0.
         fifx(1)=0.
         ficx(1)=1.
         fibx(1)=fialx+fialx
         fibbx(1)=fibex+fibex
         fiffx(2)=0.
         fifx(2)=fialx
         ficx(2)=1.+fibex
         fibx(2)=fialx
         fibbx(2)=fibex
         do i=3,nx-2
            fiffx(i)=fibex
            fifx(i)=fialx
            ficx(i)=1.
            fibx(i)=fialx
            fibbx(i)=fibex
         enddo
         fiffx(nx-1)=fibex
         fifx(nx-1)=fialx
         ficx(nx-1)=1.+fibex
         fibx(nx-1)=fialx
         fibbx(nx-1)=0.
         fiffx(nx)=fibex+fibex
         fifx(nx)=fialx+fialx
         ficx(nx)=1.
         fibx(nx)=0.
         fibbx(nx)=0.
         do i=1,nx
            fiffxp(i)=fiffx(i)
            fifxp(i)=fifx(i)
            ficxp(i)=ficx(i)
            fibxp(i)=fibx(i)
            fibbxp(i)=fibbx(i)
         enddo
         fibx(1)=0.
         fibbx(1)=0.
         ficx(2)=1.-fibex
         ficx(nx-1)=1.-fibex
         fifx(nx)=0.
         fiffx(nx)=0. 
      endif
!
      if (nclx==2) then
         fiffx(1)=0.
         fifx(1)=0.
         ficx(1)=1.
         fibx(1)=0.
         fibbx(1)=0.
         fiffx(2)=0.
         fifx(2)=0.
         ficx(2)=1.
         fibx(2)=0.
         fibbx(2)=0.
         fiffx(3)=0.
         fifx(3)=0.
         ficx(3)=1.
         fibx(3)=0.
         fibbx(3)=0.
!
         do i=4,nx-3
            fiffx(i)=fibex
            fifx(i)=fialx
            ficx(i)=1.
            fibx(i)=fialx
            fibbx(i)=fibex
         enddo
!
         fiffx(nx-2)=0.
         fifx(nx-2)=0.
         ficx(nx-2)=1.
         fibx(nx-2)=0.
         fibbx(nx-2)=0.
         fiffx(nx-1)=0.
         fifx(nx-1)=0.
         ficx(nx-1)=1.
         fibx(nx-1)=0.
         fibbx(nx-1)=0.
         fiffx(nx)=0.
         fifx(nx)=0.
         ficx(nx)=1.
         fibx(nx)=0.
         fibbx(nx)=0.
      endif
!
      call prepare_filtre(fiffx,fifx,ficx,fibx,fibbx,filax,nx)
      if (nclx==1) then
         call prepare_filtre(fiffxp,fifxp,ficxp,fibxp,&
                             fibbxp,filaxp,nx)
         do i=1,nx
            fiffxp(i)=1./fiffxp(i)
         enddo
      endif
!
      if (nclx==0) then
         call cyclix(fiffx,fifx,ficx,fibx,fibbx,&
                     filax,fiz1x,fiz2x,nx)
      endif
      do i=1,nx
         fiffx(i)=1./fiffx(i)
      enddo
      endif
!
      if (ny>1) then
      xcoef=1./2.
      if (ncly==0) then
         fiffy(1)=0.
         fify(1)=0.
         ficy(1)=1.+fibey/xcoef
         fiby(1)=fialy
         fibby(1)=fibey
         fiffy(2)=0.
         fify(2)=fialy+xcoef*fialy
         ficy(2)=1.+xcoef*fibey
         fiby(2)=fialy
         fibby(2)=fibey
!
         do j=3,ny-2
            fiffy(j)=fibey
            fify(j)=fialy
            ficy(j)=1.
            fiby(j)=fialy
            fibby(j)=fibey
         enddo
!
         fiffy(ny-1)=fibey
         fify(ny-1)=fialy
         ficy(ny-1)=1.+xcoef*fibey
         fiby(ny-1)=fialy+xcoef*fialy
         fibby(ny-1)=0.
         fiffy(ny)=fibey
         fify(ny)=fialy
         ficy(ny)=1.+fibey/xcoef
         fiby(ny)=0.
         fibby(ny)=0.
         do j=1,ny
            fiffyp(j)=fiffy(j)
            fifyp(j)=fify(j)
            ficyp(j)=ficy(j)
            fibyp(j)=fiby(j)
            fibbyp(j)=fibby(j)
         enddo
      endif
!
      if (ncly==1) then
         fiffy(1)=0.
         fify(1)=0.
         ficy(1)=1.
         fiby(1)=fialy+fialy
         fibby(1)=fibey+fibey
         fiffy(2)=0.
         fify(2)=fialy
         ficy(2)=1.+fibey
         fiby(2)=fialy
         fibby(2)=fibey
         do j=3,ny-2
            fiffy(j)=fibey
            fify(j)=fialy
            ficy(j)=1.
            fiby(j)=fialy
            fibby(j)=fibey
         enddo
         fiffy(ny-1)=fibey
         fify(ny-1)=fialy
         ficy(ny-1)=1.+fibey
         fiby(ny-1)=fialy
         fibby(ny-1)=0.
         fiffy(ny)=fibey+fibey
         fify(ny)=fialy+fialy
         ficy(ny)=1.
         fiby(ny)=0.
         fibby(ny)=0.
         do j=1,ny
            fiffyp(j)=fiffy(j)
            fifyp(j)=fify(j)
            ficyp(j)=ficy(j)
            fibyp(j)=fiby(j)
            fibbyp(j)=fibby(j)
         enddo
         fiby(1)=0.
         fibby(1)=0.
         ficy(2)=1.-fibey
         ficy(ny-1)=1.-fibey
         fify(ny)=0.
         fiffy(ny)=0. 
      endif
!
      if (ncly==2) then
         fiffy(1)=0.
         fify(1)=0.
         ficy(1)=1.
         fiby(1)=0.
         fibby(1)=0.
         fiffy(2)=0.
         fify(2)=0.
         ficy(2)=1.
         fiby(2)=0.
         fibby(2)=0.
         fiffy(3)=0.
         fify(3)=0.
         ficy(3)=1.
         fiby(3)=0.
         fibby(3)=0.
!
         do j=4,ny-3
            fiffy(j)=fibey
            fify(j)=fialy
            ficy(j)=1.
            fiby(j)=fialy
            fibby(j)=fibey
         enddo
!
         fiffy(ny-2)=0.
         fify(ny-2)=0.
         ficy(ny-2)=1.
         fiby(ny-2)=0.
         fibby(ny-2)=0.
         fiffy(ny-1)=0.
         fify(ny-1)=0.
         ficy(ny-1)=1.
         fiby(ny-1)=0.
         fibby(ny-1)=0.
         fiffy(ny)=0.
         fify(ny)=0.
         ficy(ny)=1.
         fiby(ny)=0.
         fibby(ny)=0.
      endif
!
      call prepare_filtre(fiffy,fify,ficy,fiby,fibby,filay,ny)
      if (ncly==1) then
         call prepare_filtre(fiffyp,fifyp,ficyp,fibyp,&
                             fibbyp,filayp,ny)
         do j=1,ny
            fiffyp(j)=1./fiffyp(j)
         enddo
      endif
!
      if (ncly==0) then
         call cycliy(fiffy,fify,ficy,fiby,fibby,&
                     filay,fiz1y,fiz2y,ny)
      endif
      do j=1,ny
         fiffy(j)=1./fiffy(j)
      enddo
      endif
!
      if (nz>1) then
         xcoef=1./2.
         if (nclz==0) then
         fiffz(1)=0.
         fifz(1)=0.
         ficz(1)=1.+fibez/xcoef
         fibz(1)=fialz
         fibbz(1)=fibez
         fiffz(2)=0.
         fifz(2)=fialz+fialz*xcoef
         ficz(2)=1.+fibez*xcoef
         fibz(2)=fialz
         fibbz(2)=fibez
!
         do k=3,nz-2
            fiffz(k)=fibez
            fifz(k)=fialz
            ficz(k)=1.
            fibz(k)=fialz
            fibbz(k)=fibez
         enddo
!
         fiffz(nz-1)=fibez
         fifz(nz-1)=fialz
         ficz(nz-1)=1.+fibez*xcoef
         fibz(nz-1)=fialz+fialz*xcoef
         fibbz(nz-1)=0.
         fiffz(nz)=fibez
         fifz(nz)=fialz
         ficz(nz)=1.+fibez/xcoef
         fibz(nz)=0.
         fibbz(nz)=0.
         do k=1,nz
            fiffzp(k)=fiffz(k)
            fifzp(k)=fifz(k)
            ficzp(k)=ficz(k)
            fibzp(k)=fibz(k)
            fibbzp(k)=fibbz(k)
         enddo
      endif
!
      if (nclz==1) then
         fiffz(1)=0.
         fifz(1)=0.
         ficz(1)=1.
         fibz(1)=fialz+fialz
         fibbz(1)=fibez+fibez
         fiffz(2)=0.
         fifz(2)=fialz
         ficz(2)=1.+fibez
         fibz(2)=fialz
         fibbz(2)=fibez
         do k=3,nz-2
            fiffz(k)=fibez
            fifz(k)=fialz
            ficz(k)=1.
            fibz(k)=fialz
            fibbz(k)=fibez
         enddo
         fiffz(nz-1)=fibez
         fifz(nz-1)=fialz
         ficz(nz-1)=1.+fibez
         fibz(nz-1)=fialz
         fibbz(nz-1)=0.
         fiffz(nz)=fibez+fibez
         fifz(nz)=fialz+fialz
         ficz(nz)=1.
         fibz(nz)=0.
         fibbz(nz)=0.
         do k=1,nz
            fiffzp(k)=fiffz(k)
            fifzp(k)=fifz(k)
            ficzp(k)=ficz(k)
            fibzp(k)=fibz(k)
            fibbzp(k)=fibbz(k)
         enddo
         fibz(1)=0.
         fibbz(1)=0.
         ficz(2)=1.-fibez
         ficz(nz-1)=1.-fibez
         fifz(nz)=0.
         fiffz(nz)=0. 
      endif
!
      if (nclz==2) then
         fiffz(1)=0.
         fifz(1)=0.
         ficz(1)=1.
         fibz(1)=0.
         fibbz(1)=0.
         fiffz(2)=0.
         fifz(2)=0.
         ficz(2)=1.
         fibz(2)=0.
         fibbz(2)=0.
         fiffz(3)=0.
         fifz(3)=0.
         ficz(3)=1.
         fibz(3)=0.
         fibbz(3)=0.
!
         do k=4,nz-3
            fiffz(k)=fibez
            fifz(k)=fialz
            ficz(k)=1.
            fibz(k)=fialz
            fibbz(k)=fibez
         enddo
!
         fiffz(nz-2)=0.
         fifz(nz-2)=0.
         ficz(nz-2)=1.
         fibz(nz-2)=0.
         fibbz(nz-2)=0.
         fiffz(nz-1)=0.
         fifz(nz-1)=0.
         ficz(nz-1)=1.
         fibz(nz-1)=0.
         fibbz(nz-1)=0.
         fiffz(nz)=0.
         fifz(nz)=0.
         ficz(nz)=1.
         fibz(nz)=0.
         fibbz(nz)=0.
      endif
!
      call prepare_filtre(fiffz,fifz,ficz,fibz,fibbz,&
                          filaz,nz)
      if (nclz==1) then
         call prepare_filtre(fiffzp,fifzp,ficzp,fibzp,&
                             fibbzp,filazp,nz)
         do k=1,nz
            fiffzp(k)=1./fiffzp(k)
         enddo
      endif
      if (nclz==0) then
         call cycliz(fiffz,fifz,ficz,fibz,fibbz,&
                     filaz,fiz1z,fiz2z,nz)
      endif
      do k=1,nz
         fiffz(k)=1./fiffz(k)
      enddo
   endif
!
   return  
end subroutine filtres
!
!*********************************************************************
!
subroutine prepare_filtre(aff,af,a,ab,abb,al,n) 
!...Translated by PSUITE Trans90                  4.3ZH  9:53:12   1/27/ 4  
!...Switches: -yv INDDO=0 -nbejkno
  implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  integer , intent(in) :: n 
  real(8) , intent(inout) :: aff(n) 
  real(8) , intent(inout) :: af(n) 
  real(8) , intent(inout) :: a(n) 
  real(8) , intent(inout) :: ab(n) 
  real(8) , intent(inout) :: abb(n) 
  real(8) , intent(out) :: al(n,2) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer :: l, k, i, j 
  real(8) :: tiny, dum 
!-----------------------------------------------
!
!*********************************************************************
!
   tiny=1.E-10 
!
   print *,'aff(2)=',aff(2) 
!
   aff(1)=a(1) 
   af(1)=ab(1) 
   a(1)=abb(1) 
   ab(1)=0. 
   abb(1)=0. 
!
   aff(2)=af(2) 
   af(2)=a(2) 
   a(2)=ab(2) 
   ab(2)=abb(2) 
   abb(2)=0. 
!
   print *,'aff(2)=',aff(2) 
!
   l=2 
!
   do k=1,n 
      dum=aff(k) 
      i=k 
      if (l<n) l=l+1 
      do j=k+1,l 
         if (abs(aff(j)) <= abs(dum)) cycle  
         print *,'aff(',j,')=',aff(j),'.gt.abs(aff(',k,')=',aff(k) 
         dum=aff(j) 
         i=j 
      enddo
      if (dum==0.) aff(k)=0. 
!
      if (i/=k) then 
         dum=aff(k) 
         aff(k)=aff(i) 
         aff(i)=dum 
!
         dum=af(k) 
         af(k)=af(i) 
         af(i)=dum 
!
         dum=a(k) 
         a(k)=a(i) 
         a(i)=dum 
!
         dum=ab(k) 
         ab(k)=ab(i) 
         ab(i)=dum 
!
         dum=abb(k) 
         abb(k)=abb(i) 
         abb(i)=dum 
!
      endif
      do i=k+1,l 
         dum=aff(i)/aff(k) 
         if (dum<=0.) print *,'Dum<0 pour i=',aff(i),aff(k) 
         al(k,i-k)=dum 
         aff(i)=af(i)-dum*af(k) 
         if (i<=2) print *,'aff(i)=',i,k,aff(i),af(i),af(k),dum 
         af(i)=a(i)-dum*a(k) 
         a(i)=ab(i)-dum*ab(k) 
         ab(i)=abb(i)-dum*abb(k) 
         abb(i)=0. 
      enddo
   enddo
!
   return  
end subroutine prepare_filtre
!
!*********************************************************************
!
subroutine coefficients 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
  USE parfix_m 
  USE parfiy_m 
  USE parfiz_m 
!...Translated by PSUITE Trans90                  4.3ZH  9:53:12   1/27/ 4  
!...Switches: -yv INDDO=0 -nbejkno
  implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!-----------------------------------------------
!
!*********************************************************************
!
!      omega1=1.5
!      v1=0.95
!
!      omega2=2.
!      v2=0.5
!
!      c11=0.5+9.*cos(omega1)/16.
!     1        -cos(3.*omega1)/16.-v1
!      c12=0.75+cos(omega1)+0.25*cos(2.*omega1)
!     1        -2.*v1*cos(omega1)
!      c13=10.*cos(omega1)/16.+cos(2.*omega1)
!     1        +6.*cos(3.*omega1)/16.
!     2        -2.*v1*cos(2.*omega1)
!
!      c21=0.5+9.*cos(omega2)/16.
!     1        -cos(3.*omega2)/16.-v2
!      c22=0.75+cos(omega2)+0.25*cos(2.*omega2)
!     1        -2.*v2*cos(omega2)
!      c23=10.*cos(omega2)/16.+cos(2.*omega2)
!     1        +6.*cos(3.*omega2)/16.
!     2        -2.*v2*cos(2.*omega2)
!
!      fialx=(c11*c23-c21*c13)/(-c12*c23+c22*c13)
!      fibex=-(c11+fialx*c12)/c13
!
!      print *,'Alpha,Beta=',fialx,fibex
!      fialy=fialx
!      fibey=fibex
!      print *,'Alpha,Beta=',fialy,fibey
!
!      fialz=fialx
!      fibez=fibex
!      print *,'Alpha,Beta=',fialz,fibez
!
!      fiaix=0.25*(2.+3.*fialx)
!      fibix=1./16.*(9.+16.*fialx+10.*fibex)/2.
!      ficix=0.25*(fialx+4.*fibex)/2.
!      fidix=1./16.*(6.*fibex-1.)/2.
!
!      print *,'a,b,c,d=',fiaix,fibix,ficix,fidix
!
!      fiaiy=0.25*(2.+3.*fialy)
!      fibiy=1./16.*(9.+16.*fialy+10.*fibey)/2.
!      ficiy=0.25*(fialy+4.*fibey)/2.
!      fidiy=1./16.*(6.*fibey-1.)/2.
!
!      print *,'a,b,c,d=',fiaiy,fibiy,ficiy,fidiy
!
!      fiaiz=0.25*(2.+3.*fialz)
!      fibiz=1./16.*(9.+16.*fialz+10.*fibez)/2.
!      ficiz=0.25*(fialz+4.*fibez)/2.
!      fidiz=1./16.*(6.*fibez-1.)/2.
!
!      fialx=0.15
!      fialx=0.3
      fialx=0.45 
!      fialx=0.6
!      fialx=0.64
      fibex=(3. - 2.*fialx)/10. 
      fiaix=(2. + 3.*fialx)/4. 
      fibix=(6. + 7.*fialx)/8. 
      ficix=(6. + fialx)/20. 
      fidix=(2. - 3.*fialx)/40. 
      fibix=fibix*0.5 
      ficix=ficix*0.5 
      fidix=fidix*0.5 
      print *,'Alpha,Beta=',fialx,fibex 
!      print *,'a,b,c,d=',fiaix,fibix,ficix,fidix
!
!      fialy=0.15
!      fialy=0.3
      fialy=0.45 
!      fialy=0.6
!      fialy=0.64
      fibey=(3. - 2.*fialy)/10. 
      fiaiy=(2. + 3.*fialy)/4. 
      fibiy=(6. + 7.*fialy)/8. 
      ficiy=(6. + fialy)/20. 
      fidiy=(2. - 3.*fialy)/40. 
      fibiy=fibiy*0.5 
      ficiy=ficiy*0.5 
      fidiy=fidiy*0.5 
      print *,'Alpha,Beta=',fialy,fibey 
!      print *,'a,b,c,d=',fiaiy,fibiy,ficiy,fidiy
!
!      fialz=0.15
!      fialz=0.3
      fialz=0.45 
!      fialz=0.6
!      fialz=0.64
      fibez=(3. - 2.*fialz)/10. 
      fiaiz=(2. + 3.*fialz)/4. 
      fibiz=(6. + 7.*fialz)/8. 
      ficiz=(6. + fialz)/20. 
      fidiz=(2. - 3.*fialz)/40. 
      fibiz=fibiz*0.5 
      ficiz=ficiz*0.5 
      fidiz=fidiz*0.5 
      print *,'Alpha,Beta=',fialz,fibez 
!      print *,'a,b,c,d=',fiaiz,fibiz,ficiz,fidiz
!
      fia1x=15./16. 
      fib1x=1./4. 
      fic1x=-3./8. 
      fid1x=1./4. 
      fie1x=-1./16. 
!
      fia1y=15./16. 
      fib1y=1./4. 
      fic1y=-3./8. 
      fid1y=1./4. 
      fie1y=-1./16. 
!
      fia1z=15./16. 
      fib1z=1./4. 
      fic1z=-3./8. 
      fid1z=1./4. 
      fie1z=-1./16. 
!
      fia2x=3./4. 
      fib2x=1./16. 
      fic2x=3./8. 
      fid2x=-1./4. 
      fie2x=1./16. 
!
      fia2y=3./4. 
      fib2y=1./16. 
      fic2y=3./8. 
      fid2y=-1./4. 
      fie2y=1./16. 
!
      fia2z=3./4. 
      fib2z=1./16. 
      fic2z=3./8. 
      fid2z=-1./4. 
      fie2z=1./16. 
!
      fia3x=5./8. 
      fib3x=-1./16. 
      fic3x=1./4. 
      fid3x=4./16. 
      fie3x=-1./16. 
!
      fia3y=5./8. 
      fib3y=-1./16. 
      fic3y=1./4. 
      fid3y=4./16. 
      fie3y=-1./16. 
!
      fia3z=5./8. 
      fib3z=-1./16. 
      fic3z=1./4. 
      fid3z=4./16. 
      fie3z=-1./16. 
!
      fianx=fia1x 
      fibnx=fib1x 
      ficnx=fic1x 
      fidnx=fid1x 
      fienx=fie1x 
!
      fiany=fia1y 
      fibny=fib1y 
      ficny=fic1y 
      fidny=fid1y 
      fieny=fie1y 
!
      fianz=fia1z 
      fibnz=fib1z 
      ficnz=fic1z 
      fidnz=fid1z 
      fienz=fie1z 
!
      fiamx=fia2x 
      fibmx=fib2x 
      ficmx=fic2x 
      fidmx=fid2x 
      fiemx=fie2x 
!
      fiamy=fia2y 
      fibmy=fib2y 
      ficmy=fic2y 
      fidmy=fid2y 
      fiemy=fie2y 
!
      fiamz=fia2z 
      fibmz=fib2z 
      ficmz=fic2z 
      fidmz=fid2z 
      fiemz=fie2z 
!
      fiapx=fia3x 
      fibpx=fib3x 
      ficpx=fic3x 
      fidpx=fid3x 
      fiepx=fie3x 
!
      fiapy=fia3y 
      fibpy=fib3y 
      ficpy=fic3y 
      fidpy=fid3y 
      fiepy=fie3y 
!
      fiapz=fia3z 
      fibpz=fib3z 
      ficpz=fic3z 
      fidpz=fid3z 
      fiepz=fie3z 
!
   return  
end subroutine coefficients
!
!*********************************************************************
!
subroutine cyclix(fiffx,fifx,ficx,fibx,fibbx,filax,fiz1x,fiz2x,nx) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
  USE parfix_m 
!...Translated by PSUITE Trans90                  4.3ZH  9:53:12   1/27/ 4  
!...Switches: -yv INDDO=0 -nbejkno
  implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  integer , intent(in) :: nx 
  real(8) , intent(in) :: fiffx(nx) 
  real(8) , intent(in) :: fifx(nx) 
  real(8) , intent(in) :: ficx(nx) 
  real(8) , intent(in) :: fibx(nx) 
  real(8) , intent(in) :: fibbx(nx) 
  real(8) , intent(in) :: filax(nx,2) 
  real(8) , intent(inout) :: fiz1x(nx) 
  real(8) , intent(inout) :: fiz2x(nx) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer :: i 
  real(8) :: xcoef, xdet 
!-----------------------------------------------
!
!*********************************************************************
!
   do i=1,nx
      fiz1x(i)=0.
      fiz2x(i)=0.
   enddo
!
   xcoef=2.
!
   fiz1x(1)=xcoef
   fiz1x(nx-1)=-1
   fiz2x(2)=-1.
   fiz2x(nx)=xcoef
!
   do i=1,nx-2
      fiz1x(i+1)=fiz1x(i+1)-filax(i,1)*fiz1x(i)
      fiz1x(i+2)=fiz1x(i+2)-filax(i,2)*fiz1x(i)
      fiz2x(i+1)=fiz2x(i+1)-filax(i,1)*fiz2x(i)
      fiz2x(i+2)=fiz2x(i+2)-filax(i,2)*fiz2x(i)
   enddo
!
   fiz1x(nx)=fiz1x(nx)-filax(nx-1,1)*fiz1x(nx-1)
   fiz2x(nx)=fiz2x(nx)-filax(nx-1,1)*fiz2x(nx-1)
!
   fiz1x(nx)=fiz1x(nx)/fiffx(nx)
   fiz2x(nx)=fiz2x(nx)/fiffx(nx)
   fiz1x(nx-1)=(fiz1x(nx-1)-fifx(nx-1)*fiz1x(nx))/fiffx(nx-1)
   fiz2x(nx-1)=(fiz2x(nx-1)-fifx(nx-1)*fiz2x(nx))/fiffx(nx-1)
   fiz1x(nx-2)=(fiz1x(nx-2)-fifx(nx-2)*fiz1x(nx-1)-&
                ficx(nx-2)*fiz1x(nx))/fiffx(nx-2)
   fiz2x(nx-2)=(fiz2x(nx-2)-fifx(nx-2)*fiz2x(nx-1)-&
                ficx(nx-2)*fiz2x(nx))/fiffx(nx-2)
   fiz1x(nx-3)=(fiz1x(nx-3)-fifx(nx-3)*fiz1x(nx-2)-&
                ficx(nx-3)*fiz1x(nx-1)-&
                fibx(nx-3)*fiz1x(nx))/fiffx(nx-3)
   fiz2x(nx-3)=(fiz2x(nx-3)-fifx(nx-3)*fiz2x(nx-2)-&
                ficx(nx-3)*fiz2x(nx-1)-&
                fibx(nx-3)*fiz2x(nx))/fiffx(nx-3)
   do i=nx-4,1,-1
      fiz1x(i)=(fiz1x(i)-fifx(i)*fiz1x(i+1)-&
                ficx(i)*fiz1x(i+2)-&
                fibx(i)*fiz1x(i+3)-&
                fibbx(i)*fiz1x(i+4))/fiffx(i)
      fiz2x(i)=(fiz2x(i)-fifx(i)*fiz2x(i+1)-&
                ficx(i)*fiz2x(i+2)-&
                fibx(i)*fiz2x(i+3)-&
                fibbx(i)*fiz2x(i+4))/fiffx(i)
   enddo
!
   xdet=(1.-fibex*fiz1x(1)+fiz1x(nx-1)*fibex/xcoef+&
                            fialx*fiz1x(nx)/xcoef)*&
        (1.+fialx*fiz2x(1)/xcoef+fibex*fiz2x(2)/xcoef-&
                                     fibex*fiz2x(nx))-&
        (-fibex*fiz2x(1)+fiz2x(nx-1)*fibex/xcoef+&
                          fialx*fiz2x(nx)/xcoef)*&
        (fialx*fiz1x(1)/xcoef+fibex*fiz1x(2)/xcoef-&
                                  fibex*fiz1x(nx))
!
   print *,'Apres cxclix',xdet
   print *,'Fiz1x',fiz1x(1),fiz1x(2),fiz1x(nx-1),fiz1x(nx)
   print *,'Fiz2x',fiz2x(1),fiz2x(2),fiz2x(nx-1),fiz2x(nx)
!
   fih1x=(1.+fialx*fiz2x(1)/xcoef+fibex*fiz2x(2)/xcoef-&
                                  fibex*fiz2x(nx))/xdet
!
   fih2x=-(fialx*fiz1x(1)/xcoef+fibex*fiz1x(2)/xcoef-&
                               fibex*fiz1x(nx))/xdet
!
   fih3x=-(-fibex*fiz2x(1)+fiz2x(nx-1)*fibex/xcoef+&
                           fialx*fiz2x(nx)/xcoef)/xdet
!
   fih4x=(1.-fibex*fiz1x(1)+fiz1x(nx-1)*fibex/xcoef+&
                            fialx*fiz1x(nx)/xcoef)/xdet
!
   print *,'Attention !!!,fihix=',fih1x,fih2x,fih3x,fih4x
!
   return  
end subroutine cyclix
!
!*********************************************************************
!
subroutine cycliy(fiffy,fify,ficy,fiby,fibby,filay,fiz1y,fiz2y,ny) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
  USE parfiy_m 
!...Translated by PSUITE Trans90                  4.3ZH  9:53:12   1/27/ 4  
!...Switches: -yv INDDO=0 -nbejkno
  implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  integer , intent(in) :: ny 
  real(8) , intent(in) :: fiffy(ny) 
  real(8) , intent(in) :: fify(ny) 
  real(8) , intent(in) :: ficy(ny) 
  real(8) , intent(in) :: fiby(ny) 
  real(8) , intent(in) :: fibby(ny) 
  real(8) , intent(in) :: filay(ny,2) 
  real(8) , intent(inout) :: fiz1y(ny) 
  real(8) , intent(inout) :: fiz2y(ny) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer :: j 
  real(8) :: xcoef, ydet 
!-----------------------------------------------
!
!*********************************************************************
!
   do j=1,ny
      fiz1y(j)=0.
      fiz2y(j)=0.
   enddo
!
   xcoef=2.
!
   fiz1y(1)=xcoef
   fiz1y(ny-1)=-1.
   fiz2y(2)=-1.
   fiz2y(ny)=xcoef
!
   do j=1,2
      print *,'Filay',j,'=',filay(j,1),filay(j,2)
   enddo
!
   do j=1,ny-2
      fiz1y(j+1)=fiz1y(j+1)-filay(j,1)*fiz1y(j)
      fiz1y(j+2)=fiz1y(j+2)-filay(j,2)*fiz1y(j)
      fiz2y(j+1)=fiz2y(j+1)-filay(j,1)*fiz2y(j)
      fiz2y(j+2)=fiz2y(j+2)-filay(j,2)*fiz2y(j)
   enddo
!
   fiz1y(ny)=fiz1y(ny)-filay(ny-1,1)*fiz1y(ny-1)
   fiz2y(ny)=fiz2y(ny)-filay(ny-1,1)*fiz2y(ny-1)
!
   fiz1y(ny)=fiz1y(ny)/fiffy(ny)
   fiz2y(ny)=fiz2y(ny)/fiffy(ny)
   fiz1y(ny-1)=(fiz1y(ny-1)-fify(ny-1)*fiz1y(ny))/fiffy(ny-1)
   fiz2y(ny-1)=(fiz2y(ny-1)-fify(ny-1)*fiz2y(ny))/fiffy(ny-1)
   fiz1y(ny-2)=(fiz1y(ny-2)-fify(ny-2)*fiz1y(ny-1)-&
                ficy(ny-2)*fiz1y(ny))/fiffy(ny-2)
   fiz2y(ny-2)=(fiz2y(ny-2)-fify(ny-2)*fiz2y(ny-1)-&
                ficy(ny-2)*fiz2y(ny))/fiffy(ny-2)
   fiz1y(ny-3)=(fiz1y(ny-3)-fify(ny-3)*fiz1y(ny-2)-&
                ficy(ny-3)*fiz1y(ny-1)-&
                fiby(ny-3)*fiz1y(ny))/fiffy(ny-3)
   fiz2y(ny-3)=(fiz2y(ny-3)-fify(ny-3)*fiz2y(ny-2)-&
                ficy(ny-3)*fiz2y(ny-1)-&
                fiby(ny-3)*fiz2y(ny))/fiffy(ny-3)
   do j=ny-4,1,-1
      fiz1y(j)=(fiz1y(j)-fify(j)*fiz1y(j+1)-&
                ficy(j)*fiz1y(j+2)-&
                fiby(j)*fiz1y(j+3)-&
                fibby(j)*fiz1y(j+4))/fiffy(j)
      fiz2y(j)=(fiz2y(j)-fify(j)*fiz2y(j+1)-&
                ficy(j)*fiz2y(j+2)-&
                fiby(j)*fiz2y(j+3)-&
                fibby(j)*fiz2y(j+4))/fiffy(j)
   enddo
!
   ydet=(1.-fibey*fiz1y(1)+fibey*fiz1y(ny-1)/xcoef+&
                              fialy*fiz1y(ny)/xcoef)*&
        (1.+fialy*fiz2y(1)/xcoef+fibey*fiz2y(2)/xcoef-&
                                  fibey*fiz2y(ny))-&
        (-fibey*fiz2y(1)+fibey*fiz2y(ny-1)/xcoef+&
                              fialy*fiz2y(ny)/xcoef)*&
        (fialy*fiz1y(1)/xcoef+fibey*fiz1y(2)/xcoef-&
                                  fibey*fiz1y(ny))
!
   print *,'Apres cycliy',ydet
   print *,'Fiz1y',fiz1y(1),fiz1y(2),fiz1y(ny-1),fiz1y(ny)
   print *,'Fiz2y',fiz2y(1),fiz2y(2),fiz2y(ny-1),fiz2y(ny)
!
   fih1y=(1.+fialy*fiz2y(1)/xcoef+fibey*fiz2y(2)/xcoef-&
          fibey*fiz2y(ny))/ydet
!
   fih2y=-(fialy*fiz1y(1)/xcoef+fibey*fiz1y(2)/xcoef-&
           fibey*fiz1y(ny))/ydet
!
   fih3y=-(-fibey*fiz2y(1)+fibey*fiz2y(ny-1)/xcoef+&
           fialy*fiz2y(ny)/xcoef)/ydet
!   
   fih4y=(1.-fibey*fiz1y(1)+fibey*fiz1y(ny-1)/xcoef+&
          fialy*fiz1y(ny)/xcoef)/ydet
!
   print *,'Attention !!!,fihiy=',fih1y,fih2y,fih3y,fih4y
!
   return  
end subroutine cycliy
!
!*********************************************************************
!
subroutine cycliz(fiffz,fifz,ficz,fibz,fibbz,filaz,fiz1z,fiz2z,nz) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
  USE parfiz_m 
!...Translated by PSUITE Trans90                  4.3ZH  9:53:12   1/27/ 4  
!...Switches: -yv INDDO=0 -nbejkno
  implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  integer , intent(in) :: nz 
  real(8) , intent(in) :: fiffz(nz) 
  real(8) , intent(in) :: fifz(nz) 
  real(8) , intent(in) :: ficz(nz) 
  real(8) , intent(in) :: fibz(nz) 
  real(8) , intent(in) :: fibbz(nz) 
  real(8) , intent(in) :: filaz(nz,2) 
  real(8) , intent(inout) :: fiz1z(nz) 
  real(8) , intent(inout) :: fiz2z(nz) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer :: k 
  real(8) :: xcoef, zdet 
!-----------------------------------------------
!
!*********************************************************************
!
   do k=1,nz
      fiz1z(k)=0.
      fiz2z(k)=0.
   enddo
!
   xcoef=2.
!
   fiz1z(1)=xcoef
   fiz1z(nz-1)=-1.
   fiz2z(2)=-1.
   fiz2z(nz)=xcoef
!
   do k=1,nz-2
      fiz1z(k+1)=fiz1z(k+1)-filaz(k,1)*fiz1z(k)
      fiz1z(k+2)=fiz1z(k+2)-filaz(k,2)*fiz1z(k)
      fiz2z(k+1)=fiz2z(k+1)-filaz(k,1)*fiz2z(k)
      fiz2z(k+2)=fiz2z(k+2)-filaz(k,2)*fiz2z(k)
   enddo
!
   fiz1z(nz)=fiz1z(nz)-filaz(nz-1,1)*fiz1z(nz-1)
   fiz2z(nz)=fiz2z(nz)-filaz(nz-1,1)*fiz2z(nz-1)
!
   fiz1z(nz)=fiz1z(nz)/fiffz(nz)
   fiz2z(nz)=fiz2z(nz)/fiffz(nz)
   fiz1z(nz-1)=(fiz1z(nz-1)-fifz(nz-1)*fiz1z(nz))/fiffz(nz-1)
   fiz2z(nz-1)=(fiz2z(nz-1)-fifz(nz-1)*fiz2z(nz))/fiffz(nz-1)
   fiz1z(nz-2)=(fiz1z(nz-2)-fifz(nz-2)*fiz1z(nz-1)-&
                ficz(nz-2)*fiz1z(nz))/fiffz(nz-2)
   fiz2z(nz-2)=(fiz2z(nz-2)-fifz(nz-2)*fiz2z(nz-1)-&
                ficz(nz-2)*fiz2z(nz))/fiffz(nz-2)
   fiz1z(nz-3)=(fiz1z(nz-3)-fifz(nz-3)*fiz1z(nz-2)-&
                ficz(nz-3)*fiz1z(nz-1)-&
                fibz(nz-3)*fiz1z(nz))/fiffz(nz-3)
   fiz2z(nz-3)=(fiz2z(nz-3)-fifz(nz-3)*fiz2z(nz-2)-&
                ficz(nz-3)*fiz2z(nz-1)-&
                fibz(nz-3)*fiz2z(nz))/fiffz(nz-3)
   do k=nz-4,1,-1
      fiz1z(k)=(fiz1z(k)-fifz(k)*fiz1z(k+1)-&
                ficz(k)*fiz1z(k+2)-&
                fibz(k)*fiz1z(k+3)-&
                fibbz(k)*fiz1z(k+4))/fiffz(k)
      fiz2z(k)=(fiz2z(k)-fifz(k)*fiz2z(k+1)-&
                ficz(k)*fiz2z(k+2)-&
                fibz(k)*fiz2z(k+3)-&
                fibbz(k)*fiz2z(k+4))/fiffz(k)
   enddo
!
   zdet=(1.-fibez*fiz1z(1)+fiz1z(nz-1)*fibez/xcoef+&
                            fialz*fiz1z(nz)/xcoef)*&
        (1.+fialz*fiz2z(1)/xcoef+fibez*fiz2z(2)/xcoef-&
                                     fibez*fiz2z(nz))-&
        (-fibez*fiz2z(1)+fiz2z(nz-1)*fibez/xcoef+&
                          fialz*fiz2z(nz)/xcoef)*&
        (fialz*fiz1z(1)/xcoef+fibez*fiz1z(2)/xcoef-&
                                   fibez*fiz1z(nz))
!
   print *,'Apres cycliz',zdet
!
   fih1z=(1.+fialz*fiz2z(1)/xcoef+fibez*fiz2z(2)/xcoef-&
                                  fibez*fiz2z(nz))/zdet
!
   fih2z=-(fialz*fiz1z(1)/xcoef+fibez*fiz1z(2)/xcoef-&
                                  fibez*fiz1z(nz))/zdet
!
   fih3z=-(-fibez*fiz2z(1)+fiz2z(nz-1)*fibez/xcoef+&
                            fialz*fiz2z(nz)/xcoef)/zdet
!
   fih4z=(1.-fibez*fiz1z(1)+fiz1z(nz-1)*fibez/xcoef+&
                            fialz*fiz1z(nz)/xcoef)/zdet
!
   print *,'Attention !!!,fihiz=',fih1z,fih2z,fih3z,fih4z
!
   return  
end subroutine cycliz
