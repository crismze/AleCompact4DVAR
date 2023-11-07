!********************************************************************
!
subroutine derx(tx,ux,ffx,fsx,fwx,nx,ny,nz,npaire) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
  USE paramx_m
  USE derivex_m
  implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  integer , intent(in) :: nx 
  integer , intent(in) :: ny 
  integer , intent(in) :: nz 
  integer , intent(in) :: npaire 
  real , intent(inout) :: tx(nx,ny,nz) 
  real(8) , intent(in) :: ux(nx,ny,nz) 
  real(8) , intent(in) :: ffx(nx) 
  real(8) , intent(in) :: fsx(nx) 
  real(8) , intent(in) :: fwx(nx) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer :: k, j, i 
!-----------------------------------------------
!
!********************************************************************
!
   if (nclx==2) then 
      do k=1,nz 
      do j=1,ny 
         tx(1,j,k)=af1x*ux(1,j,k)+bf1x*ux(2,j,k)+cf1x*ux(3,j,k) 
         tx(2,j,k)=af2x*(ux(3,j,k)-ux(1,j,k)) 
      enddo 
      enddo 
      do i=3,nx-2 
      do k=1,nz 
      do j=1,ny 
         tx(i,j,k)=afix*(ux(i+1,j,k)-ux(i-1,j,k))&
                  +bfix*(ux(i+2,j,k)-ux(i-2,j,k)) 
      enddo 
      enddo 
      enddo 
      do k=1,nz 
      do j=1,ny 
         tx(nx-1,j,k)=afmx*(ux(nx,j,k)-ux(nx-2,j,k)) 
         tx(nx,j,k)=(-afnx*ux(nx,j,k))-bfnx*ux(nx-1,j,k)-cfnx*ux(nx-2,j,k) 
      enddo
      enddo
      do i=2,nx 
      do k=1,nz 
      do j=1,ny 
         tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*fsx(i) 
      enddo
      enddo
      enddo
      do k=1,nz 
      do j=1,ny 
         tx(nx,j,k)=tx(nx,j,k)*fwx(nx) 
      enddo 
      enddo 
      do i=nx-1,1,-1 
      do k=1,nz 
      do j=1,ny 
         tx(i,j,k)=(tx(i,j,k)-ffx(i)*tx(i+1,j,k))*fwx(i) 
      enddo 
      enddo 
      enddo 
   endif 
!
end subroutine derx 
!
!********************************************************************
!
subroutine dery(ty,uy,ffy,fsy,fwy,nx,ny,nz,npaire) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
  USE paramy_m 
  USE derivey_m 
  implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  integer , intent(in) :: nx 
  integer , intent(in) :: ny 
  integer , intent(in) :: nz 
  integer , intent(in) :: npaire 
  real(8) , intent(inout) :: ty(nx,nz,ny) 
  real(8) , intent(in) :: uy(nx,nz,ny) 
  real(8) , intent(in) :: ffy(ny) 
  real(8) , intent(in) :: fsy(ny) 
  real(8) , intent(in) :: fwy(ny) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer :: k, j, i 
!-----------------------------------------------
!
!********************************************************************
!
   if (ncly==1) then 
      if (npaire==1) then 
         do k=1,nz 
         do i=1,nx 
            ty(i,k,1)=0. 
            ty(i,k,2)=afjy*(uy(i,k,3)-uy(i,k,1))&
                     +bfjy*(uy(i,k,4)-uy(i,k,2)) 
         enddo 
         enddo 
         do k=1,nz 
         do i=1,nx 
         do j=3,ny-2 
            ty(i,k,j)=afjy*(uy(i,k,j+1)-uy(i,k,j-1))&
                     +bfjy*(uy(i,k,j+2)-uy(i,k,j-2)) 
         enddo 
         enddo 
         enddo 
         do k=1,nz 
         do i=1,nx 
            ty(i,k,ny-1)=afjy*(uy(i,k,ny)-uy(i,k,ny-2))&
                        +bfjy*(uy(i,k,ny-1)-uy(i,k,ny-3)) 
            ty(i,k,ny)=0. 
         enddo 
         enddo 
         do j=2,ny 
         do k=1,nz 
         do i=1,nx 
            ty(i,k,j)=ty(i,k,j)-ty(i,k,j-1)*fsy(j) 
         enddo 
         enddo 
         enddo 
         do k=1,nz 
         do i=1,nx 
            ty(i,k,ny)=ty(i,k,ny)*fwy(ny) 
         enddo 
         enddo 
         do j=ny-1,1,-1 
         do k=1,nz 
         do i=1,nx 
            ty(i,k,j)=(ty(i,k,j)-ffy(j)*ty(i,k,j+1))*fwy(j) 
         enddo 
         enddo 
         enddo 
      endif
      if (npaire==0) then 
         do k=1,nz 
         do i=1,nx 
            ty(i,k,1)=afjy*(uy(i,k,2)+uy(i,k,2))&
                     +bfjy*(uy(i,k,3)+uy(i,k,3)) 
            ty(i,k,2)=afjy*(uy(i,k,3)-uy(i,k,1))&
                     +bfjy*(uy(i,k,4)+uy(i,k,2)) 
         enddo 
         enddo 
         do k=1,nz 
         do i=1,nx 
         do j=3,ny-2 
            ty(i,k,j)=afjy*(uy(i,k,j+1)-uy(i,k,j-1))&
                     +bfjy*(uy(i,k,j+2)-uy(i,k,j-2)) 
         enddo 
         enddo 
         enddo 
         do k=1,nz 
         do i=1,nx 
            ty(i,k,ny-1)=afjy*(uy(i,k,ny)-uy(i,k,ny-2))&
                        +bfjy*((-uy(i,k,ny-1))-uy(i,k,ny-3)) 
            ty(i,k,ny)=afjy*((-uy(i,k,ny-1))-uy(i,k,ny-1))&
                      +bfjy*((-uy(i,k,ny-2))-uy(i,k,ny-2)) 
         enddo 
         enddo 
         do j=2,ny 
         do k=1,nz 
         do i=1,nx 
            ty(i,k,j)=ty(i,k,j)-ty(i,k,j-1)*fsy(j) 
         enddo 
         enddo 
         enddo 
         do k=1,nz 
         do i=1,nx 
            ty(i,k,ny)=ty(i,k,ny)*fwy(ny) 
         enddo
         enddo 
         do j=ny-1,1,-1 
         do k=1,nz 
         do i=1,nx 
            ty(i,k,j)=(ty(i,k,j)-ffy(j)*ty(i,k,j+1))*fwy(j) 
         enddo 
         enddo 
         enddo 
      endif
   endif
!
!********************STRET STRET************************
   do k=1,nz 
   do j=1,ny 
   do i=1,nx 
       ty(i,j,k)=ty(i,k,j)*yly
   enddo 
   enddo
   enddo 
!*******************************************************
!
   return  
end subroutine dery
!
!********************************************************************
!
subroutine derz(tz,uz,rz,sz,ffz,fsz,fwz,nx,ny,nz,npaire) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
  USE paramz_m 
  USE derivez_m 
!...Translated by PSUITE Trans90                  4.3ZH  9:45:39   1/26/ 4  
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
  real(8) , intent(in) :: ffz(nz) 
  real(8) , intent(in) :: fsz(nz) 
  real(8) , intent(in) :: fwz(nz) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer :: j, i, k 
!-----------------------------------------------
!
!********************************************************************
!
!
   if (nclz==0) then
      do j=1,ny
      do i=1,nx
         tz(i,j,1)=afkz*(uz(i,j,2)-uz(i,j,nz  ))&
                  +bfkz*(uz(i,j,3)-uz(i,j,nz-1))
         rz(i,j,1)=-1.
         tz(i,j,2)=afkz*(uz(i,j,3)-uz(i,j,1 ))&
                  +bfkz*(uz(i,j,4)-uz(i,j,nz))
         rz(i,j,2)=0.
      enddo
      enddo
      do k=3,nz-2
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=afkz*(uz(i,j,k+1)-uz(i,j,k-1))&
                  +bfkz*(uz(i,j,k+2)-uz(i,j,k-2))
         rz(i,j,k)=0.
      enddo
      enddo
      enddo
      do j=1,ny
      do i=1,nx
         tz(i,j,nz-1)=afkz*(uz(i,j,nz)-uz(i,j,nz-2))&
                     +bfkz*(uz(i,j,1 )-uz(i,j,nz-3))
         rz(i,j,nz-1)=0.
         tz(i,j,nz  )=afkz*(uz(i,j,1)-uz(i,j,nz-1))&
                     +bfkz*(uz(i,j,2)-uz(i,j,nz-2))
         rz(i,j,nz  )=alfakz
      enddo
      enddo
      do k=2,nz
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=tz(i,j,k)-tz(i,j,k-1)*fsz(k)
         rz(i,j,k)=rz(i,j,k)-rz(i,j,k-1)*fsz(k)
      enddo
      enddo
      enddo
      do j=1,ny
      do i=1,nx
         tz(i,j,nz)=tz(i,j,nz)*fwz(nz)
         rz(i,j,nz)=rz(i,j,nz)*fwz(nz)
      enddo
      enddo
      do k=nz-1,1,-1
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=(tz(i,j,k)-ffz(k)*tz(i,j,k+1))*fwz(k)
         rz(i,j,k)=(rz(i,j,k)-ffz(k)*rz(i,j,k+1))*fwz(k)
      enddo
      enddo
      enddo
      do j=1,ny
      do i=1,nx
         sz(i,j)=(   tz(i,j,1)-alfakz*tz(i,j,nz))/&
                 (1.+rz(i,j,1)-alfakz*rz(i,j,nz))
      enddo
      enddo
      do k=1,nz
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=tz(i,j,k)-sz(i,j)*rz(i,j,k)
      enddo
      enddo
      enddo
   endif
!
   if (nclz==1) then
      if (npaire==1) then
         do j=1,ny
         do i=1,nx
            tz(i,j,1)=0.
            tz(i,j,2)=afkz*(uz(i,j,3)-uz(i,j,1))&
                     +bfkz*(uz(i,j,4)-uz(i,j,2))
         enddo
         enddo
         do k=3,nz-2
         do j=1,ny
         do i=1,nx
               tz(i,j,k)=afkz*(uz(i,j,k+1)-uz(i,j,k-1))&
                        +bfkz*(uz(i,j,k+2)-uz(i,j,k-2))
         enddo
         enddo
         enddo
         do j=1,ny
         do i=1,nx
            tz(i,j,nz-1)=afkz*(uz(i,j,nz  )-uz(i,j,nz-2))&
                        +bfkz*(uz(i,j,nz-1)-uz(i,j,nz-3))
            tz(i,j,nz  )=0.
         enddo
         enddo
         do k=2,nz
         do j=1,ny
         do i=1,nx
            tz(i,j,k)=tz(i,j,k)-tz(i,j,k-1)*fsz(k)
         enddo
         enddo
         enddo
         do j=1,ny
         do i=1,nx
            tz(i,j,nz)=tz(i,j,nz)*fwz(nz)
         enddo
         enddo
         do k=nz-1,1,-1
         do j=1,ny
         do i=1,nx
            tz(i,j,k)=(tz(i,j,k)-ffz(k)*tz(i,j,k+1))*fwz(k)
         enddo
         enddo
         enddo
      endif
      if (npaire==0) then
         do j=1,ny
         do i=1,nx
            tz(i,j,1)=afkz*(uz(i,j,2)+uz(i,j,2))&
                     +bfkz*(uz(i,j,3)+uz(i,j,3))
            tz(i,j,2)=afkz*(uz(i,j,3)-uz(i,j,1))&
                     +bfkz*(uz(i,j,4)+uz(i,j,2))
         enddo
         enddo
         do k=3,nz-2
         do j=1,ny
         do i=1,nx
            tz(i,j,k)=afkz*(uz(i,j,k+1)-uz(i,j,k-1))&
                     +bfkz*(uz(i,j,k+2)-uz(i,j,k-2))
         enddo
         enddo
         enddo
         do j=1,ny
         do i=1,nx
            tz(i,j,nz-1)=afkz*( uz(i,j,nz  )-uz(i,j,nz-2))&
                        +bfkz*(-uz(i,j,nz-1)-uz(i,j,nz-3))
            tz(i,j,nz  )=afkz*(-uz(i,j,nz-1)-uz(i,j,nz-1))&
                        +bfkz*(-uz(i,j,nz-2)-uz(i,j,nz-2))
         enddo
         enddo
         do k=2,nz
         do j=1,ny
         do i=1,nx
            tz(i,j,k)=tz(i,j,k)-tz(i,j,k-1)*fsz(k)
         enddo
         enddo
         enddo
         do j=1,ny
         do i=1,nx
            tz(i,j,nz)=tz(i,j,nz)*fwz(nz)
         enddo
         enddo
         do k=nz-1,1,-1
         do j=1,ny
         do i=1,nx
            tz(i,j,k)=(tz(i,j,k)-ffz(k)*tz(i,j,k+1))*fwz(k)
         enddo
         enddo
         enddo
      endif
   endif
!
   if (nclz==2) then
      do j=1,ny
      do i=1,nx
         tz(i,j,1)=af1z*uz(i,j,1)+bf1z*uz(i,j,2)&
                  +cf1z*uz(i,j,3)
         tz(i,j,2)=af2z*(uz(i,j,3)-uz(i,j,1))
      enddo
      enddo
      do k=3,nz-2
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=afkz*(uz(i,j,k+1)-uz(i,j,k-1))&
                  +bfkz*(uz(i,j,k+2)-uz(i,j,k-2))
      enddo
      enddo
      enddo
      do j=1,ny
      do i=1,nx
         tz(i,j,nz-1)= afmz*(uz(i,j,nz)-uz(i,j,nz-2))
         tz(i,j,nz  )=-afnz*uz(i,j,nz)-bfnz*uz(i,j,nz-1)&
                      -cfnz*uz(i,j,nz-2)
      enddo
      enddo
      do k=2,nz
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=tz(i,j,k)-tz(i,j,k-1)*fsz(k)
      enddo
      enddo
      enddo
      do j=1,ny
      do i=1,nx
         tz(i,j,nz)=tz(i,j,nz)*fwz(nz)
      enddo
      enddo
      do k=nz-1,1,-1
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=(tz(i,j,k)-ffz(k)*tz(i,j,k+1))*fwz(k)
      enddo
      enddo
      enddo
   endif
!
   return
end subroutine derz
!
!********************************************************************
!
subroutine derxx(tx,ux,sfx,ssx,swx,nx,ny,nz,npaire) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
  USE paramx_m 
  USE derivex_m 
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
  real(8) , intent(in) :: sfx(nx) 
  real(8) , intent(in) :: ssx(nx) 
  real(8) , intent(in) :: swx(nx) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer :: k, j, i 
!-----------------------------------------------
!
!********************************************************************
!
   if (nclx==2) then
      do k=1,nz
      do j=1,ny
         tx(1,j,k)=as1x*ux(1,j,k)+bs1x*ux(2,j,k)&
                  +cs1x*ux(3,j,k)+ds1x*ux(4,j,k)
         tx(2,j,k)=as2x*(ux(3,j,k)-ux(2,j,k)&
                        -ux(2,j,k)+ux(1,j,k))
      enddo
      enddo
!
      do i=3,nx-2
      do k=1,nz
      do j=1,ny
         tx(i,j,k)=asix*(ux(i+1,j,k)+ux(i-1,j,k))&
                  +bsix*(ux(i+2,j,k)+ux(i-2,j,k))&
                  -(asix*2.+bsix*2.)*ux(i,j,k)
      enddo
      enddo
      enddo
!
      do k=1,nz
      do j=1,ny
         tx(nx-1,j,k)=asmx*(ux(nx  ,j,k)-ux(nx-1,j,k)&
                           -ux(nx-1,j,k)+ux(nx-2,j,k))
         tx(nx  ,j,k)=asnx*ux(nx  ,j,k)+bsnx*ux(nx-1,j,k)&
                     +csnx*ux(nx-2,j,k)+dsnx*ux(nx-3,j,k)
      enddo
      enddo
      do i=2,nx
      do k=1,nz
      do j=1,ny
         tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*ssx(i)
      enddo
      enddo
      enddo
      do k=1,nz
      do j=1,ny
         tx(nx,j,k)=tx(nx,j,k)*swx(nx)
      enddo
      enddo
      do i=nx-1,1,-1
      do k=1,nz
      do j=1,ny
         tx(i,j,k)=(tx(i,j,k)-sfx(i)*tx(i+1,j,k))*swx(i)
      enddo
      enddo
      enddo
   endif
!
   return  
end subroutine derxx
!
!********************************************************************
!
subroutine deryy(ty,uy,sfy,ssy,swy,nx,ny,nz,npaire) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
  USE paramy_m 
  USE derivey_m 
  implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  integer , intent(in) :: nx 
  integer , intent(in) :: ny 
  integer , intent(in) :: nz 
  integer , intent(in) :: npaire 
  real(8) , intent(inout) :: ty(nx,nz,ny) 
  real(8) , intent(in) :: uy(nx,nz,ny) 
  real(8) , intent(in) :: sfy(ny) 
  real(8) , intent(in) :: ssy(ny) 
  real(8) , intent(in) :: swy(ny) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer :: k, j, i 
!-----------------------------------------------
!
!********************************************************************
!
   if (ncly==1) then
      if (npaire==1) then
         do k=1,nz
         do i=1,nx
            ty(i,k,1)=asjy*(uy(i,k,2)-uy(i,k,1)&
                           -uy(i,k,1)+uy(i,k,2))&
                     +bsjy*(uy(i,k,3)-uy(i,k,1)&
                           -uy(i,k,1)+uy(i,k,3))
            ty(i,k,2)=asjy*(uy(i,k,3)-uy(i,k,2)&
                           -uy(i,k,2)+uy(i,k,1))&
                     +bsjy*(uy(i,k,4)-uy(i,k,2)&
                           -uy(i,k,2)+uy(i,k,2))
         enddo
         enddo
         do k=1,nz
         do i=1,nx
         do j=3,ny-2
            ty(i,k,j)=asjy*(uy(i,k,j+1)-uy(i,k,j  )&
                           -uy(i,k,j  )+uy(i,k,j-1))&
                     +bsjy*(uy(i,k,j+2)-uy(i,k,j  )&
                           -uy(i,k,j  )+uy(i,k,j-2))
         enddo
         enddo
         enddo
         do k=1,nz
         do i=1,nx
            ty(i,k,ny-1)=asjy*(uy(i,k,ny  )-uy(i,k,ny-1)&
                              -uy(i,k,ny-1)+uy(i,k,ny-2))&
                        +bsjy*(uy(i,k,ny-1)-uy(i,k,ny-1)&
                              -uy(i,k,ny-1)+uy(i,k,ny-3))
            ty(i,k,ny  )=asjy*(uy(i,k,ny-1)-uy(i,k,ny  )&
                              -uy(i,k,ny  )+uy(i,k,ny-1))&
                        +bsjy*(uy(i,k,ny-2)-uy(i,k,ny  )&
                              -uy(i,k,ny  )+uy(i,k,ny-2))
         enddo
         enddo
         do j=2,ny
         do k=1,nz
         do i=1,nx
            ty(i,k,j)=ty(i,k,j)-ty(i,k,j-1)*ssy(j)
         enddo
         enddo
         enddo
         do k=1,nz
         do i=1,nx
            ty(i,k,ny)=ty(i,k,ny)*swy(ny)
         enddo
         enddo
         do j=ny-1,1,-1
         do k=1,nz
         do i=1,nx
            ty(i,k,j)=(ty(i,k,j)-sfy(j)*ty(i,k,j+1))*swy(j)
         enddo
         enddo
         enddo
      endif
      if (npaire==0) then
         do k=1,nz
         do i=1,nx
            ty(i,k,1)=0.
            ty(i,k,2)=asjy*(uy(i,k,3)-uy(i,k,2)&
                           -uy(i,k,2)+uy(i,k,1))&
                     +bsjy*(uy(i,k,4)-uy(i,k,2)&
                           -uy(i,k,2)-uy(i,k,2))
         enddo
         enddo
         do k=1,nz
         do i=1,nx
         do j=3,ny-2
            ty(i,k,j)=asjy*(uy(i,k,j+1)-uy(i,k,j  )&
                           -uy(i,k,j  )+uy(i,k,j-1))&
                     +bsjy*(uy(i,k,j+2)-uy(i,k,j  )&
                           -uy(i,k,j  )+uy(i,k,j-2))
         enddo
         enddo
         enddo
         do k=1,nz
         do i=1,nx
            ty(i,k,ny-1)=asjy*(uy(i,k,ny  )-uy(i,k,ny-1)&
                              -uy(i,k,ny-1)+uy(i,k,ny-2))&
                        +bsjy*(-uy(i,k,ny-1)-uy(i,k,ny-1)&
                               -uy(i,k,ny-1)+uy(i,k,ny-3))
            ty(i,k,ny  )=0.
         enddo
         enddo
         do j=2,ny
         do k=1,nz
         do i=1,nx
            ty(i,k,j)=ty(i,k,j)-ty(i,k,j-1)*ssy(j)
         enddo
         enddo
         enddo
         do k=1,nz
         do i=1,nx
            ty(i,k,ny)=ty(i,k,ny)*swy(ny)
         enddo
         enddo
         do j=ny-1,1,-1
         do k=1,nz
         do i=1,nx
            ty(i,k,j)=(ty(i,k,j)-sfy(j)*ty(i,k,j+1))*swy(j)
         enddo
         enddo
         enddo
      endif
   endif
!
   return  
end subroutine deryy
!
!********************************************************************
!
subroutine derzz(tz,uz,rz,sz,sfz,ssz,swz,nx,ny,nz,npaire) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
  USE paramz_m 
  USE derivez_m 
!...Translated by PSUITE Trans90                  4.3ZH  9:45:39   1/26/ 4  
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
  real(8) , intent(in) :: sfz(nz) 
  real(8) , intent(in) :: ssz(nz) 
  real(8) , intent(in) :: swz(nz) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer :: j, i, k 
!-----------------------------------------------
!
!********************************************************************
!
   if (nclz==0) then
      do j=1,ny
      do i=1,nx
         tz(i,j,1)=askz*(uz(i,j,2)-uz(i,j,1   )&
                        -uz(i,j,1)+uz(i,j,nz  ))&
                  +bskz*(uz(i,j,3)-uz(i,j,1   )&
                        -uz(i,j,1)+uz(i,j,nz-1))
         rz(i,j,1)=-1.
         tz(i,j,2)=askz*(uz(i,j,3)-uz(i,j,2 )&
                        -uz(i,j,2)+uz(i,j,1 ))&
                  +bskz*(uz(i,j,4)-uz(i,j,2 )&
                        -uz(i,j,2)+uz(i,j,nz))
         rz(i,j,2)=0.
      enddo
      enddo
      do k=3,nz-2
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=askz*(uz(i,j,k+1)-uz(i,j,k  )&
                        -uz(i,j,k  )+uz(i,j,k-1))&
                  +bskz*(uz(i,j,k+2)-uz(i,j,k  )&
                        -uz(i,j,k  )+uz(i,j,k-2))
          rz(i,j,k)=0.
      enddo
      enddo
      enddo
      do j=1,ny
      do i=1,nx
         tz(i,j,nz-1)=askz*(uz(i,j,nz  )-uz(i,j,nz-1)&
                           -uz(i,j,nz-1)+uz(i,j,nz-2))&
                     +bskz*(uz(i,j,1   )-uz(i,j,nz-1)&
                           -uz(i,j,nz-1)+uz(i,j,nz-3))
         rz(i,j,nz-1)=0.
         tz(i,j,nz  )=askz*(uz(i,j,1 )-uz(i,j,nz  )&
                           -uz(i,j,nz)+uz(i,j,nz-1))&
                     +bskz*(uz(i,j,2 )-uz(i,j,nz  )&
                           -uz(i,j,nz)+uz(i,j,nz-2))
         rz(i,j,nz  )=alsakz
      enddo
      enddo
      do k=2,nz
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=tz(i,j,k)-tz(i,j,k-1)*ssz(k)
         rz(i,j,k)=rz(i,j,k)-rz(i,j,k-1)*ssz(k)
      enddo
      enddo
      enddo
      do j=1,ny
      do i=1,nx
         tz(i,j,nz)=tz(i,j,nz)*swz(nz)
         rz(i,j,nz)=rz(i,j,nz)*swz(nz)
      enddo
      enddo
      do k=nz-1,1,-1
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=(tz(i,j,k)-sfz(k)*tz(i,j,k+1))*swz(k)
         rz(i,j,k)=(rz(i,j,k)-sfz(k)*rz(i,j,k+1))*swz(k)
      enddo
      enddo
      enddo
      do j=1,ny
      do i=1,nx
         sz(i,j)=(   tz(i,j,1)-alsakz*tz(i,j,nz))/&
                 (1.+rz(i,j,1)-alsakz*rz(i,j,nz))
      enddo
      enddo
      do k=1,nz
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=tz(i,j,k)-sz(i,j)*rz(i,j,k)
      enddo
      enddo
      enddo
   endif
!
   if (nclz==1) then
      if (npaire==1) then
         do j=1,ny
         do i=1,nx
            tz(i,j,1)=askz*(uz(i,j,2)-uz(i,j,1)&
                           -uz(i,j,1)+uz(i,j,2))&
                     +bskz*(uz(i,j,3)-uz(i,j,1)&
                           -uz(i,j,1)+uz(i,j,3))
            tz(i,j,2)=askz*(uz(i,j,3)-uz(i,j,2)&
                           -uz(i,j,2)+uz(i,j,1))&
                     +bskz*(uz(i,j,4)-uz(i,j,2)&
                           -uz(i,j,2)+uz(i,j,2))
         enddo
         enddo
         do k=3,nz-2
         do j=1,ny
         do i=1,nx
            tz(i,j,k)=askz*(uz(i,j,k+1)-uz(i,j,k  )&
                           -uz(i,j,k  )+uz(i,j,k-1))&
                     +bskz*(uz(i,j,k+2)-uz(i,j,k  )&
                           -uz(i,j,k  )+uz(i,j,k-2))
         enddo
         enddo
         enddo
         do j=1,ny
         do i=1,nx
            tz(i,j,nz-1)=askz*(uz(i,j,nz  )-uz(i,j,nz-1)&
                              -uz(i,j,nz-1)+uz(i,j,nz-2))&
                        +bskz*(uz(i,j,nz-1)-uz(i,j,nz-1)&
                              -uz(i,j,nz-1)+uz(i,j,nz-3))
            tz(i,j,nz  )=askz*(uz(i,j,nz-1)-uz(i,j,nz  )&
                              -uz(i,j,nz  )+uz(i,j,nz-1))&
                        +bskz*(uz(i,j,nz-2)-uz(i,j,nz  )&
                              -uz(i,j,nz  )+uz(i,j,nz-2))
         enddo
         enddo
         do k=2,nz
         do j=1,ny
         do i=1,nx
            tz(i,j,k)=tz(i,j,k)-tz(i,j,k-1)*ssz(k)
         enddo
         enddo
         enddo
         do j=1,ny
         do i=1,nx
            tz(i,j,nz)=tz(i,j,nz)*swz(nz)
         enddo
         enddo
         do k=nz-1,1,-1
         do j=1,ny
         do i=1,nx
            tz(i,j,k)=(tz(i,j,k)-sfz(k)*tz(i,j,k+1))*swz(k)
         enddo
         enddo
         enddo
      endif
      if (npaire==0) then
         do j=1,ny
         do i=1,nx
            tz(i,j,1)=0.
            tz(i,j,2)=askz*(uz(i,j,3)-uz(i,j,2)&
                           -uz(i,j,2)+uz(i,j,1))&
                     +bskz*(uz(i,j,4)-uz(i,j,2)&
                           -uz(i,j,2)-uz(i,j,2))
         enddo
         enddo
         do k=3,nz-2
         do j=1,ny
         do i=1,nx
            tz(i,j,k)=askz*(uz(i,j,k+1)-uz(i,j,k  )&
                           -uz(i,j,k  )+uz(i,j,k-1))&
                     +bskz*(uz(i,j,k+2)-uz(i,j,k  )&
                           -uz(i,j,k  )+uz(i,j,k-2))
         enddo
         enddo
         enddo
         do j=1,ny
         do i=1,nx
            tz(i,j,nz-1)=askz*(uz(i,j,nz  )-uz(i,j,nz-1)&
                              -uz(i,j,nz-1)+uz(i,j,nz-2))&
                        +bskz*(-uz(i,j,nz-1)-uz(i,j,nz-1)&
                              -uz(i,j,nz-1)+uz(i,j,nz-3))
            tz(i,j,nz  )=0.
         enddo
         enddo
         do k=2,nz
         do j=1,ny
         do i=1,nx
            tz(i,j,k)=tz(i,j,k)-tz(i,j,k-1)*ssz(k)
         enddo
         enddo
         enddo
         do j=1,ny
         do i=1,nx
            tz(i,j,nz)=tz(i,j,nz)*swz(nz)
         enddo
         enddo
         do k=nz-1,1,-1
         do j=1,ny
         do i=1,nx
            tz(i,j,k)=(tz(i,j,k)-sfz(k)*tz(i,j,k+1))*swz(k)
         enddo
         enddo
         enddo
      endif
   endif
!
   if (nclz==2) then
      do j=1,ny
      do i=1,nx
         tz(i,j,1)=as1z*uz(i,j,1)+bs1z*uz(i,j,2)&
                  +cs1z*uz(i,j,3)+ds1z*uz(i,j,4)
         tz(i,j,2)=as2z*(uz(i,j,3)-uz(i,j,2)&
                        -uz(i,j,2)+uz(i,j,1))
      enddo
      enddo
      do k=3,nz-2
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=askz*(uz(i,j,k+1)-uz(i,j,k  )&
                        -uz(i,j,k  )+uz(i,j,k-1))&
                  +bskz*(uz(i,j,k+2)-uz(i,j,k  )&
                        -uz(i,j,k  )+uz(i,j,k-2))
      enddo
      enddo
      enddo
      do j=1,ny
      do i=1,nx
         tz(i,j,nz-1)=asmz*(uz(i,j,nz  )-uz(i,j,nz-1)&
                           -uz(i,j,nz-1)+uz(i,j,nz-2))
         tz(i,j,nz  )=asnz*uz(i,j,nz  )+bsnz*uz(i,j,nz-1)&
                     +csnz*uz(i,j,nz-2)+dsnz*uz(i,j,nz-3)
      enddo
      enddo
      do k=2,nz
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=tz(i,j,k)-tz(i,j,k-1)*ssz(k)
      enddo
      enddo
      enddo
      do j=1,ny
      do i=1,nx
         tz(i,j,nz)=tz(i,j,nz)*swz(nz)
      enddo
      enddo
      do k=nz-1,1,-1
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=(tz(i,j,k)-sfz(k)*tz(i,j,k+1))*swz(k)
      enddo
      enddo
      enddo
   endif
!
   return  
end subroutine derzz
!
!********************************************************************
!
subroutine derxb(tx,ux,rx,sx,ffx,fsx,fwx,nx,ny,nz,npaire)
!
!********************************************************************
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
   USE paramx_m
   USE derivex_m
!...Translated by PSUITE Trans90                  4.3ZH  9:45:39   1/26/ 4  
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
   real(8) , intent(in) :: ffx(nx) 
   real(8) , intent(in) :: fsx(nx) 
   real(8) , intent(in) :: fwx(nx) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
   integer :: k, j, i 
!-----------------------------------------------
!
!   if (nclx.eq.1) then
      if (npaire.eq.1) then
         do k=1,nz
         do j=1,ny
            tx(1,j,k)=0.
            tx(2,j,k)=afix*(ux(3,j,k)-ux(1,j,k))&
                 +bfix*(ux(4,j,k)-ux(2,j,k))
         enddo
         enddo
         do i=3,nx-2
         do k=1,nz
         do j=1,ny
            tx(i,j,k)=afix*(ux(i+1,j,k)-ux(i-1,j,k))&
                    +bfix*(ux(i+2,j,k)-ux(i-2,j,k))
         enddo
         enddo
         enddo
         do k=1,nz
         do j=1,ny
            tx(nx-1,j,k)=afix*(ux(nx  ,j,k)-ux(nx-2,j,k))&
                 +bfix*(ux(nx-1,j,k)-ux(nx-3,j,k))
            tx(nx  ,j,k)=0.
         enddo
         enddo
         do i=2,nx
         do k=1,nz
         do j=1,ny
            tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*fsx(i)
         enddo
         enddo
         enddo
         do k=1,nz
         do j=1,ny
            tx(nx,j,k)=tx(nx,j,k)*fwx(nx)
         enddo
         enddo
         do i=nx-1,1,-1
         do k=1,nz
         do j=1,ny
            tx(i,j,k)=(tx(i,j,k)-ffx(i)*tx(i+1,j,k))*fwx(i)
         enddo
         enddo
         enddo
      endif
      if (npaire.eq.0) then
         do k=1,nz
         do j=1,ny
            tx(1,j,k)=afix*(2.*ux(2,j,k)-2.*ux(1,j,k))&
                 +bfix*(2.*ux(3,j,k)-2.*ux(1,j,k))
            tx(2,j,k)=afix*(ux(3,j,k)-ux(1,j,k))&
                 +bfix*(ux(4,j,k)-2.*ux(1,j,k)+ux(2,j,k))
         enddo
         enddo
         do i=3,nx-2
         do k=1,nz
         do j=1,ny
            tx(i,j,k)=afix*(ux(i+1,j,k)-ux(i-1,j,k))&
                 +bfix*(ux(i+2,j,k)-ux(i-2,j,k))
         enddo
         enddo
         enddo
         do k=1,nz
         do j=1,ny
            tx(nx-1,j,k)=afix*(ux(nx  ,j,k)-ux(nx-2,j,k))&
                 +bfix*(2.*ux(nx,j,k)-ux(nx-1,j,k)-ux(nx-3,j,k))
            tx(nx  ,j,k)=afix*(2.*ux(nx,j,k)-2.*ux(nx-1,j,k))&
                 +bfix*(2.*ux(nx,j,k)-2.*ux(nx-2,j,k))
         enddo
         enddo
         do i=2,nx
         do k=1,nz
         do j=1,ny
            tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*fsx(i)
         enddo
         enddo
         enddo
         do k=1,nz
         do j=1,ny
            tx(nx,j,k)=tx(nx,j,k)*fwx(nx)
         enddo
         enddo
         do i=nx-1,1,-1
         do k=1,nz
         do j=1,ny
            tx(i,j,k)=(tx(i,j,k)-ffx(i)*tx(i+1,j,k))*fwx(i)
         enddo
         enddo
         enddo
      endif
!   endif
!
end subroutine derxb
!********************************************************************
!
subroutine derxxb(tx,ux,rx,sx,sfx,ssx,swx,nx,ny,nz,npaire)
!
!********************************************************************
!
!
!   M o d u l e s 
!-----------------------------------------------
  USE paramx_m 
  USE derivex_m 
!...Translated by PSUITE Trans90                  4.3ZH  9:45:39   1/26/ 4  
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
  real(8) , intent(in) :: sfx(nx) 
  real(8) , intent(in) :: ssx(nx) 
  real(8) , intent(in) :: swx(nx) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer :: k, j, i 
!-----------------------------------------------

!  if (nclx.eq.2) then
     if (npaire.eq.1) then
        do k=1,nz
        do j=1,ny
           tx(1,j,k)=asix*(ux(2,j,k)-ux(1,j,k)&
                -ux(1,j,k)+ux(2,j,k))&
                +bsix*(ux(3,j,k)-ux(1,j,k)&
                -ux(1,j,k)+ux(3,j,k))
           tx(2,j,k)=asix*(ux(3,j,k)-ux(2,j,k)&
                -ux(2,j,k)+ux(1,j,k))&
                +bsix*(ux(4,j,k)-ux(2,j,k)&
                -ux(2,j,k)+ux(2,j,k))
        enddo
        enddo
        do i=3,nx-2
        do k=1,nz
        do j=1,ny
           tx(i,j,k)=asix*(ux(i+1,j,k)-ux(i  ,j,k)&
                -ux(i  ,j,k)+ux(i-1,j,k))&
                +bsix*(ux(i+2,j,k)-ux(i  ,j,k)&
                -ux(i  ,j,k)+ux(i-2,j,k))
        enddo
        enddo
        enddo
        do k=1,nz
        do j=1,ny
           tx(nx-1,j,k)=asix*(ux(nx  ,j,k)-ux(nx-1,j,k)&
                -ux(nx-1,j,k)+ux(nx-2,j,k))&
                +bsix*(ux(nx-1,j,k)-ux(nx-1,j,k)&
                -ux(nx-1,j,k)+ux(nx-3,j,k))
           tx(nx  ,j,k)=asix*(ux(nx-1,j,k)-ux(nx  ,j,k)&
                -ux(nx  ,j,k)+ux(nx-1,j,k))&
                +bsix*(ux(nx-2,j,k)-ux(nx  ,j,k)&
                -ux(nx  ,j,k)+ux(nx-2,j,k))
        enddo
        enddo
        do i=2,nx
        do k=1,nz
        do j=1,ny
           tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*ssx(i)
        enddo
        enddo
        enddo
        do k=1,nz
        do j=1,ny
           tx(nx,j,k)=tx(nx,j,k)*swx(nx)
        enddo
        enddo
        do i=nx-1,1,-1
        do k=1,nz
        do j=1,ny
           tx(i,j,k)=(tx(i,j,k)-sfx(i)*tx(i+1,j,k))*swx(i)
        enddo
        enddo
        enddo
     endif
     if (npaire.eq.0) then
        do k=1,nz
        do j=1,ny
           tx(1,j,k)=0.
           tx(2,j,k)=asix*(ux(3,j,k)-ux(2,j,k)&
                -ux(2,j,k)+ux(1,j,k))&
                +bsix*(ux(4,j,k)-3.*ux(2,j,k)&
                -2.*ux(1,j,k))
        enddo
        enddo
        do i=3,nx-2
        do k=1,nz
        do j=1,ny
           tx(i,j,k)=asix*(ux(i+1,j,k)-ux(i  ,j,k)&
                -ux(i  ,j,k)+ux(i-1,j,k))&
                +bsix*(ux(i+2,j,k)-ux(i  ,j,k)&
                -ux(i  ,j,k)+ux(i-2,j,k))
        enddo
        enddo
        enddo
        do k=1,nz
        do j=1,ny
           tx(nx-1,j,k)=asix*(ux(nx  ,j,k)-ux(nx-1,j,k)&
                -ux(nx-1,j,k)+ux(nx-2,j,k))&
                +bsix*(2.*ux(nx,j,k)-ux(nx-1,j,k)&
                -2.*ux(nx-1,j,k)+ux(nx-3,j,k))
           tx(nx  ,j,k)=0.
        enddo
        enddo
        do i=2,nx
        do k=1,nz
        do j=1,ny
           tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*ssx(i)
        enddo
        enddo
        enddo
        do k=1,nz
        do j=1,ny
           tx(nx,j,k)=tx(nx,j,k)*swx(nx)
        enddo
        enddo
        do i=nx-1,1,-1
        do k=1,nz
        do j=1,ny
           tx(i,j,k)=(tx(i,j,k)-sfx(i)*tx(i+1,j,k))*swx(i)
        enddo
        enddo
        enddo
     endif
!
     return
end subroutine derxxb
!
!********************************************************************
!
subroutine decx6(tx,ux,cfx6,csx6,cwx6,nx,nxm,ny,nym,nz,npaire) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
  USE paramx_m 
  USE dericex6_m 
  implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  integer , intent(in) :: nx 
  integer , intent(in) :: nxm 
  integer  :: ny 
  integer , intent(in) :: nym 
  integer , intent(in) :: nz 
  integer , intent(in) :: npaire 
  real , intent(inout) :: tx(nxm,nym,nz) 
  real , intent(in) :: ux(nx,nym,nz) 
  real , intent(in) :: cfx6(nxm) 
  real , intent(in) :: csx6(nxm) 
  real , intent(in) :: cwx6(nxm) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer :: k, j, i 
!-----------------------------------------------
!
      if (npaire==1) then
         do k=1,nz
         do j=1,nym
            tx(1,j,k)=acix6*(ux(2,j,k)-ux(1,j,k))&
                     +bcix6*(ux(3,j,k)-ux(2,j,k))             
            tx(2,j,k)=acix6*(ux(3,j,k)-ux(2,j,k))&
                     +bcix6*(ux(4,j,k)-ux(1,j,k))             
         enddo
         enddo
         do i=3,nxm-2
         do k=1,nz
         do j=1,nym
            tx(i,j,k)=acix6*(ux(i+1,j,k)-ux(i,j,k))&
                     +bcix6*(ux(i+2,j,k)-ux(i-1,j,k))              
         enddo
         enddo
         enddo
         do k=1,nz
         do j=1,nym
            tx(nxm-1,j,k)=acix6*(ux(nxm,j,k)-ux(nxm-1,j,k))&
                         +bcix6*(ux(nx,j,k)-ux(nxm-2,j,k))         
            tx(nxm,j,k)=acix6*(ux(nx,j,k)-ux(nxm,j,k))&
                       +bcix6*(ux(nxm,j,k)-ux(nxm-1,j,k))            
         enddo
         enddo
         do i=2,nxm
         do k=1,nz
         do j=1,nym
            tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*csx6(i)
         enddo
         enddo
         enddo
         do k=1,nz
         do j=1,nym
            tx(nxm,j,k)=tx(nxm,j,k)*cwx6(nxm)
         enddo
         enddo
         do i=nxm-1,1,-1
         do k=1,nz
         do j=1,nym
            tx(i,j,k)=(tx(i,j,k)-cfx6(i)*tx(i+1,j,k))*cwx6(i)
         enddo
         enddo
         enddo
      endif
      if (npaire==0) then
         do k=1,nz
         do j=1,nym
            tx(1,j,k)=acix6*(ux(2,j,k)-ux(1,j,k))&
                     +bcix6*(ux(3,j,k)-2.*ux(1,j,k)+ux(2,j,k)) 
            tx(2,j,k)=acix6*(ux(3,j,k)-ux(2,j,k))&
                     +bcix6*(ux(4,j,k)-ux(1,j,k))
         enddo
         enddo
         do i=3,nxm-2
         do k=1,nz
         do j=1,nym
            tx(i,j,k)=acix6*(ux(i+1,j,k)-ux(i,j,k))&
                     +bcix6*(ux(i+2,j,k)-ux(i-1,j,k))
         enddo
         enddo
         enddo
         do k=1,nz
         do j=1,nym
            tx(nxm-1,j,k)=acix6*(ux(nxm,j,k)-ux(nxm-1,j,k))&
                         +bcix6*(ux(nx,j,k)-ux(nxm-2,j,k)) 
            tx(nxm,j,k)=acix6*(ux(nx,j,k)-ux(nxm,j,k))&
                       +bcix6*(2.*ux(nx,j,k)-ux(nxm,j,k)-ux(nxm-1,j,k)) 
         enddo
         enddo
         do i=2,nxm
         do k=1,nz
         do j=1,nym
            tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*csx6(i)
         enddo
         enddo
         enddo
         do k=1,nz
         do j=1,nym
            tx(nxm,j,k)=tx(nxm,j,k)*cwx6(nxm)
         enddo
         enddo
         do i=nxm-1,1,-1
         do k=1,nz
         do j=1,nym
            tx(i,j,k)=(tx(i,j,k)-cfx6(i)*tx(i+1,j,k))*cwx6(i)
         enddo
         enddo
         enddo
      endif
!
end subroutine decx6
!
!********************************************************************
!
subroutine inter6(tx,ux,cifx6,cisx6,ciwx6,nx,nxm,ny,nz,npaire) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
  USE paramx_m 
  USE interpol6_m 
  implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  integer , intent(in) :: nx,nxm,ny,nz,npaire
  real(8) , intent(inout) :: tx(nxm,ny,nz) 
  real(8) , intent(in) :: ux(nx,ny,nz) 
  real(8) , intent(in) :: cifx6(nxm) 
  real(8) , intent(in) :: cisx6(nxm) 
  real(8) , intent(in) :: ciwx6(nxm) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer :: k, j, i 
!-----------------------------------------------
!
!********************************************************************
!   
      if (npaire==1) then
         do k=1,nz
         do j=1,ny
            tx(1,j,k)=aicix6*(ux(2,j,k)+ux(1,j,k))&
                     +bicix6*(ux(3,j,k)+ux(2,j,k))             
            tx(2,j,k)=aicix6*(ux(3,j,k)+ux(2,j,k))&
                     +bicix6*(ux(4,j,k)+ux(1,j,k))             
         enddo
         enddo
         do i=3,nxm-2
         do k=1,nz
         do j=1,ny
            tx(i,j,k)=aicix6*(ux(i+1,j,k)+ux(i,j,k))&
                     +bicix6*(ux(i+2,j,k)+ux(i-1,j,k))              
         enddo
         enddo
         enddo
         do k=1,nz
         do j=1,ny
            tx(nxm-1,j,k)=aicix6*(ux(nxm,j,k)+ux(nxm-1,j,k))&
                         +bicix6*(ux(nx,j,k)+ux(nxm-2,j,k))         
            tx(nxm,j,k)=aicix6*(ux(nx,j,k)+ux(nxm,j,k))&
                       +bicix6*(ux(nxm,j,k)+ux(nxm-1,j,k))            
         enddo
         enddo
         do i=2,nxm
         do k=1,nz
         do j=1,ny
               tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*cisx6(i)
         enddo
         enddo
         enddo
         do k=1,nz
         do j=1,ny
            tx(nxm,j,k)=tx(nxm,j,k)*ciwx6(nxm)
         enddo
         enddo
         do i=nxm-1,1,-1
         do k=1,nz
         do j=1,ny
            tx(i,j,k)=(tx(i,j,k)-cifx6(i)*tx(i+1,j,k))*ciwx6(i)
         enddo
         enddo
         enddo 
      endif
!  
end subroutine inter6 
!
!********************************************************************
!
subroutine deci6(tx,ux,cfi6,csi6,cwi6,nxm,nx,ny,nz,npaire) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
  USE paramx_m 
  USE dericex6_m 
  implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
   integer :: nx,nxm,ny,nz,npaire
   real(8),dimension(nx,ny,nz) :: tx
   real(8),dimension(nxm,ny,nz) :: ux
   real(8),dimension(nx) :: cfi6,csi6,cwi6
!
   integer :: i,j,k
!
      if (npaire==1) then
         do k=1,nz
         do j=1,ny
            tx(1,j,k)=0.  
            tx(2,j,k)=acix6*(ux(2,j,k)-ux(1,j,k))&
                 +bcix6*(ux(3,j,k)-ux(1,j,k))  
         enddo
         enddo
         do i=3,nx-2
         do k=1,nz
         do j=1,ny
            tx(i,j,k)=acix6*(ux(i,j,k)-ux(i-1,j,k))&
                 +bcix6*(ux(i+1,j,k)-ux(i-2,j,k))
         enddo
         enddo
         enddo
         do k=1,nz
         do j=1,ny
            tx(nx-1,j,k)=acix6*(ux(nx-1,j,k)-ux(nx-2,j,k))&
                        +bcix6*(ux(nx-1,j,k)-ux(nx-3,j,k))
            tx(nx,j,k)=0.
         enddo
         enddo
         do i=2,nx
         do k=1,nz
         do j=1,ny
            tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*csi6(i)
         enddo
         enddo
         enddo
         do k=1,nz
         do j=1,ny
            tx(nx,j,k)=tx(nx,j,k)*cwi6(nx)
         enddo
         enddo
         do i=nx-1,1,-1
         do k=1,nz
         do j=1,ny
            tx(i,j,k)=(tx(i,j,k)-cfi6(i)*tx(i+1,j,k))*cwi6(i)
         enddo
         enddo
         enddo
      endif
!
end subroutine deci6 
!
!********************************************************************
!
subroutine interi6(tx,ux,cifi6,cisi6,ciwi6,nxm,nx,ny,nz,npaire) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
  USE paramx_m 
  USE interpol6_m 
  implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
   integer :: nx,nxm,ny,nym,nz,npaire
   real(8),dimension(nx,ny,nz) :: tx
   real(8),dimension(nxm,ny,nz) :: ux
   real(8),dimension(nx) :: cifi6,cisi6,ciwi6
!
   integer :: i,j,k
!
         if (npaire==1) then
            do k=1,nz
            do j=1,ny
               tx(1,j,k)=aicix6*(ux(1,j,k)+ux(1,j,k))&
                       +bicix6*(ux(2,j,k)+ux(2,j,k))     
               tx(2,j,k)=aicix6*(ux(2,j,k)+ux(1,j,k))&
                       +bicix6*(ux(3,j,k)+ux(1,j,k))  
            enddo
            enddo
           
            do i=3,nx-2
            do k=1,nz
            do j=1,ny
               tx(i,j,k)=aicix6*(ux(i,j,k)+ux(i-1,j,k))&
                       +bicix6*(ux(i+1,j,k)+ux(i-2,j,k))
            enddo
            enddo
            enddo
            do k=1,nz
            do j=1,ny
               tx(nx-1,j,k)=aicix6*(ux(nx-1,j,k)+ux(nx-2,j,k))&
                           +bicix6*(ux(nx-1,j,k)+ux(nx-3,j,k))
               tx(nx,j,k)=aicix6*(ux(nx-1,j,k)+ux(nx-1,j,k))&
                         +bicix6*(ux(nx-2,j,k)+ux(nx-2,j,k)) 
            enddo
            enddo
            do i=2,nx
            do k=1,nz
            do j=1,ny
               tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*cisi6(i)
            enddo
            enddo
            enddo
          
            do k=1,nz
            do j=1,ny
               tx(nx,j,k)=tx(nx,j,k)*ciwi6(nx)
            enddo
            enddo
            do i=nx-1,1,-1
            do k=1,nz
            do j=1,ny
               tx(i,j,k)=(tx(i,j,k)-cifi6(i)*tx(i+1,j,k))*ciwi6(i)
            enddo
            enddo
            enddo 
         endif
!
end subroutine interi6
!
!********************************************************************
!
subroutine intery6(ty,uy,cify6,cisy6,ciwy6,nx,ny,nym,nz,npaire) 
!
  USE paramy_m
  implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  integer  :: nx,ny,nym,nz,npaire
  real(8)  :: ty(nx,nym,nz) 
  real(8) , intent(in) :: uy(nx,ny,nz) 
  real(8)  :: ry(nx,nz,nym) 
  real(8)  :: di(nx,nz,ny) 
  real(8)  :: cify6(nym) 
  real(8)  :: cisy6(nym) 
  real(8)  :: ciwy6(nym) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer :: k, j, i 
!-----------------------------------------------
!
!********************************************************************
!
   do k=1,nz 
   do j=1,ny 
   do i=1,nx 
      di(i,k,j)=uy(i,j,k) 
   enddo 
   enddo 
   enddo 
!
   call intery16(ry,di,cify6,cisy6,ciwy6,nx,ny,nym,nz,npaire) 
!
   do k=1,nz 
   do j=1,nym 
   do i=1,nx 
      ty(i,j,k)=ry(i,k,j) 
   enddo 
   enddo 
   enddo 
!
   return  
end subroutine intery6
!
!********************************************************************
!
subroutine intery16(ty,uy,cify6,cisy6,ciwy6,nx,ny,nym,nz,npaire) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
  USE interpoly6_m 
  USE paramy_m 
  USE derivey_m 
  implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  integer , intent(in) :: nx,ny,nym,nz,npaire 
  real(8) , intent(inout) :: ty(nx,nz,nym) 
  real(8) , intent(in) :: uy(nx,nz,ny) 
  real(8) , intent(in) :: cify6(nym) 
  real(8) , intent(in) :: cisy6(nym) 
  real(8) , intent(in) :: ciwy6(nym) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer :: k, i, j 
!-----------------------------------------------
!
!********************************************************************
!
      if (npaire==1) then
         do k=1,nz
         do i=1,nx
            ty(i,k,1)=aiciy6*(uy(i,k,2)+uy(i,k,1))&
                     +biciy6*(uy(i,k,3)+uy(i,k,2))             
            ty(i,k,2)=aiciy6*(uy(i,k,3)+uy(i,k,2))&
                     +biciy6*(uy(i,k,4)+uy(i,k,1))             
         enddo
         enddo
         do j=3,nym-2
         do k=1,nz
         do i=1,nx
            ty(i,k,j)=aiciy6*(uy(i,k,j+1)+uy(i,k,j))&
                     +biciy6*(uy(i,k,j+2)+uy(i,k,j-1))              
         enddo
         enddo
         enddo
         do k=1,nz
         do i=1,nx
            ty(i,k,nym-1)=aiciy6*(uy(i,k,nym)+uy(i,k,nym-1))&
                         +biciy6*(uy(i,k,ny)+uy(i,k,nym-2))         
            ty(i,k,nym)=aiciy6*(uy(i,k,ny)+uy(i,k,nym))&
                       +biciy6*(uy(i,k,nym)+uy(i,k,nym-1))            
         enddo
         enddo
         do j=2,nym
         do k=1,nz
         do i=1,nx
            ty(i,k,j)=ty(i,k,j)-ty(i,k,j-1)*cisy6(j)
         enddo
         enddo
         enddo
         do k=1,nz
         do i=1,nx
            ty(i,k,nym)=ty(i,k,nym)*ciwy6(nym)
         enddo
         enddo
         do j=nym-1,1,-1
         do k=1,nz
         do i=1,nx
            ty(i,k,j)=(ty(i,k,j)-cify6(j)*ty(i,k,j+1))*ciwy6(j)
         enddo
         enddo
         enddo 
      endif
! 
   return  
end subroutine intery16
!
!********************************************************************
!
subroutine decy6(ty,uy,cfy6,csy6,cwy6,nx,nxm,ny,nym,nz,npaire) 
!
  USE paramy_m
  implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  integer  :: nx,nxm,ny,nym,nz,npaire
  real(8)  :: ty(nxm,nym,nz) 
  real(8) , intent(in) :: uy(nxm,ny,nz) 
  real(8)  :: ry(nxm,nz,nym) 
  real(8)  :: di(nxm,nz,ny) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer :: k, j, i 
  real(8), dimension(nym) :: cfy6, csy6, cwy6 
!-----------------------------------------------
!
!********************************************************************
!
   do k=1,nz 
   do j=1,ny 
   do i=1,nxm 
      di(i,k,j)=uy(i,j,k) 
   enddo 
   enddo 
   enddo 
!
   call decy16 (ry,di,cfy6,csy6,cwy6,nx,nxm,ny,nym,nz,npaire) 
!
   do k=1,nz 
   do j=1,nym 
   do i=1,nxm 
      ty(i,j,k)=ry(i,k,j)*yly 
   enddo 
   enddo 
   enddo 
!
   return  
end subroutine decy6 
!
!********************************************************************
!
subroutine decy16(ty,uy,cfy6,csy6,cwy6,nx,nxm,ny,nym,nz,npaire) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
  USE dericey6_m 
  USE paramy_m 
  implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  integer  :: nx,nxm,ny,nym,nz,npaire
  real(8) , intent(inout) :: ty(nxm,nz,nym) 
  real(8) , intent(in) :: uy(nxm,nz,ny) 
  real(8) , intent(in) :: cfy6(nym) 
  real(8) , intent(in) :: csy6(nym) 
  real(8) , intent(in) :: cwy6(nym) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer :: k, i, j 
!-----------------------------------------------
!
!********************************************************************
!
      if (npaire==0) then
         do k=1,nz
         do i=1,nxm
            ty(i,k,1)=aciy6*(uy(i,k,2)-uy(i,k,1))&
                     +bciy6*(uy(i,k,3)-2.*uy(i,k,1)+uy(i,k,2)) 
            ty(i,k,2)=aciy6*(uy(i,k,3)-uy(i,k,2))&
                     +bciy6*(uy(i,k,4)-uy(i,k,1))
         enddo
         enddo
         do j=3,nym-2
         do k=1,nz
         do i=1,nxm
            ty(i,j,k)=aciy6*(uy(i,k,j+1)-uy(i,k,j))&
                     +bciy6*(uy(i,k,j+2)-uy(i,k,j-1))
         enddo
         enddo
         enddo
         do k=1,nz
         do i=1,nxm
            ty(i,k,nym-1)=aciy6*(uy(i,k,nym)-uy(i,k,nym-1))&
                         +bciy6*(uy(i,k,ny)-uy(i,k,nym-2)) 
            ty(i,k,nym)=aciy6*(uy(i,k,ny)-uy(i,k,nym))&
                +bciy6*(2.*uy(i,k,ny)-uy(i,k,nym)-uy(i,k,nym-1)) 
         enddo
         enddo
         do j=2,nym
         do k=1,nz
         do i=1,nxm
            ty(i,k,j)=ty(i,k,j)-ty(i,k,j-1)*csy6(j)
         enddo
         enddo
         enddo
         do k=1,nz
         do i=1,nxm
            ty(i,k,nym)=ty(i,k,nym)*cwy6(nym)
         enddo
         enddo
         do j=nym-1,1,-1
         do k=1,nz
         do i=1,nxm
            ty(i,k,j)=(ty(i,k,j)-cfy6(j)*ty(i,k,j+1))*cwy6(j)
         enddo
         enddo
         enddo
      endif
!
   return  
end subroutine decy16
!
!********************************************************************
!
subroutine interiy6(ty,uy,cifi6y,cisi6y,ciwi6y,nx,nym,ny,nz,npaire) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
  Use paramy_m
  implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
   integer :: nx,nxm,ny,nym,nz,npaire
   real(8),dimension(nx,ny,nz) :: ty
   real(8),dimension(nx,nym,nz) :: uy
   real(8),dimension(nx,nz,ny) :: ry
   real(8),dimension(nx,nz,nym) :: di
   real(8),dimension(ny) :: cifi6y,cisi6y,ciwi6y
!
   integer :: i,j,k
!
      do k=1,nz
      do j=1,nym
      do i=1,nx
         di(i,k,j)=uy(i,j,k)
      enddo
      enddo
      enddo
! 
      call interiy16(ry,di,cifi6y,cisi6y,ciwi6y,nx,nym,ny,nz,npaire)
!
      do k=1,nz
      do j=1,ny
      do i=1,nx
         ty(i,j,k)=ry(i,k,j)
      enddo
      enddo
      enddo
!
      return 
end subroutine interiy6
!
!********************************************************************
!
subroutine interiy16(ty,uy,cifi6y,cisi6y,ciwi6y,nx,nym,ny,nz,npaire) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
  USE interpoly6_m 
  USE paramy_m 
  implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
   integer :: nx,nxm,ny,nym,nz,npaire
   real(8),dimension(nx,nz,ny) :: ty
   real(8),dimension(nx,nz,nym) :: uy
   real(8),dimension(ny) :: cifi6y,cisi6y,ciwi6y
!
   integer :: i,j,k
!
         if (npaire==1) then
           do k=1,nz
            do i=1,nx
               ty(i,k,1)=aiciy6*(uy(i,k,1)+uy(i,k,1))&
                       +biciy6*(uy(i,k,2)+uy(i,k,2))             
               ty(i,k,2)=aiciy6*(uy(i,k,2)+uy(i,k,1))&
                       +biciy6*(uy(i,k,3)+uy(i,k,1))             
            enddo
            enddo
            do j=3,ny-2
            do k=1,nz
            do i=1,nx
               ty(i,k,j)=aiciy6*(uy(i,k,j)+uy(i,k,j-1))&
                       +biciy6*(uy(i,k,j+1)+uy(i,k,j-2))              
            enddo
            enddo
            enddo
            do k=1,nz
            do i=1,nx
               ty(i,k,ny-1)=aiciy6*(uy(i,k,ny-1)+uy(i,k,ny-2))&
                           +biciy6*(uy(i,k,ny-1)+uy(i,k,ny-3))         
               ty(i,k,ny)=aiciy6*(uy(i,k,ny-1)+uy(i,k,ny-1))&
                         +biciy6*(uy(i,k,ny-2)+uy(i,k,ny-2))            
            enddo
            enddo
            do j=2,ny
            do k=1,nz
            do i=1,nx
               ty(i,k,j)=ty(i,k,j)-ty(i,k,j-1)*cisi6y(j)
            enddo
            enddo
            enddo
            do k=1,nz
            do i=1,nx
               ty(i,k,ny)=ty(i,k,ny)*ciwi6y(ny)
            enddo
            enddo
            do j=ny-1,1,-1
            do k=1,nz
            do i=1,nx
               ty(i,k,j)=(ty(i,k,j)-cifi6y(j)*ty(i,k,j+1))*ciwi6y(j)
            enddo
            enddo
            enddo 
         endif
!
   return  
end subroutine interiy16 
!
!********************************************************************
!
subroutine deciy6(ty,uy,cfi6y,csi6y,cwi6y,nx,nym,ny,nz,npaire) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
  USE paramy_m
  implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
   integer :: nx,ny,nym,nz,npaire
   real(8),dimension(nx,ny,nz) :: ty
   real(8),dimension(nx,nym,nz) :: uy
   real(8),dimension(nx,nz,ny) :: ry
   real(8),dimension(nx,nz,nym) :: di
   real(8),dimension(ny) :: cfi6y,csi6y,cwi6y
!
   integer :: i,j,k
!
      do k=1,nz
      do j=1,nym
      do i=1,nx
         di(i,k,j)=uy(i,j,k)
      enddo
      enddo
      enddo
! 
      call deciy16(ry,di,cfi6y,csi6y,cwi6y,nx,nym,ny,nz,npaire)
!
      do k=1,nz
      do j=1,ny
      do i=1,nx
         ty(i,j,k)=ry(i,k,j)*yly
      enddo
      enddo
      enddo
!
      return
end subroutine deciy6 
!
!********************************************************************
!
subroutine deciy16(ty,uy,cfi6y,csi6y,cwi6y,nx,nym,ny,nz,npaire) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
  USE dericey6_m 
  USE paramy_m 
  implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
   integer :: nx,ny,nym,nz,npaire
   real(8),dimension(nx,nz,ny) :: ty
   real(8),dimension(nx,nz,nym) :: uy
   real(8),dimension(ny) :: cfi6y,csi6y,cwi6y
!
   integer :: i,j,k
!
         if (npaire==1) then
           do k=1,nz
            do i=1,nx
               ty(i,k,1)=0.            
               ty(i,k,2)=aciy6*(uy(i,k,2)-uy(i,k,1))&
                       +bciy6*(uy(i,k,3)-uy(i,k,1))             
            enddo
            enddo
            do j=3,ny-2
            do k=1,nz
            do i=1,nx
               ty(i,k,j)=aciy6*(uy(i,k,j)-uy(i,k,j-1))&
                       +bciy6*(uy(i,k,j+1)-uy(i,k,j-2))              
            enddo
            enddo
            enddo
            do k=1,nz
            do i=1,nx
               ty(i,k,ny-1)=aciy6*(uy(i,k,ny-1)-uy(i,k,ny-2))&
                           +bciy6*(uy(i,k,ny-1)-uy(i,k,ny-3))         
               ty(i,k,ny)=0.          
            enddo
            enddo
            do j=2,ny
            do k=1,nz
            do i=1,nx
               ty(i,k,j)=ty(i,k,j)-ty(i,k,j-1)*csi6y(j)
            enddo
            enddo
            enddo
            do k=1,nz
            do i=1,nx
               ty(i,k,ny)=ty(i,k,ny)*cwi6y(ny)
            enddo
            enddo
            do j=ny-1,1,-1
            do k=1,nz
            do i=1,nx
               ty(i,k,j)=(ty(i,k,j)-cfi6y(j)*ty(i,k,j+1))*cwi6y(j)
            enddo
            enddo
            enddo 
         endif
!
   return  
end subroutine deciy16
!
!********************************************************************
!
subroutine deci60(tx,ux,rx,sx,cfx6,csx6,cwx6,nx,mx,ny,nz,npaire) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
  USE paramx_m 
  USE dericex6_m 
!...Translated by PSUITE Trans90                  4.3ZH  9:45:39   1/26/ 4  
!...Switches: -yv INDDO=0 -nbejkno
  implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
   integer :: nx,mx,ny,nz,npaire
   real(8),dimension(nx,ny,nz) :: tx
   real(8),dimension(nx,ny,nz) :: ux,rx
   real(8),dimension(ny,nz) :: sx
   real(8),dimension(nx) :: cfx6,csx6,cwx6
!
   integer :: i,j,k
!
   if (nclx==0) then
      do k=1,nz
      do j=1,ny
         tx(1,j,k)=acix6*(ux(2,j,k)-ux(1  ,j,k))&
              +bcix6*(ux(3,j,k)-ux(nx,j,k))
         rx(1,j,k)=-1.
         tx(2,j,k)=acix6*(ux(3,j,k)-ux(2 ,j,k))&
              +bcix6*(ux(4,j,k)-ux(1,j,k))
         rx(2,j,k)=0.
      enddo
      enddo
      do i=3,nx-2
      do k=1,nz
      do j=1,ny
         tx(i,j,k)=acix6*(ux(i+1,j,k)-ux(i,j,k))&
                 +bcix6*(ux(i+2,j,k)-ux(i-1,j,k))
         rx(i,j,k)=0.
      enddo
      enddo
      enddo
      do k=1,nz
      do j=1,ny
         tx(nx-1,j,k)=acix6*(ux(nx,j,k)-ux(nx-1,j,k))&
                    +bcix6*(ux(1 ,j,k)-ux(nx-2,j,k))
         rx(nx-1,j,k)=0.
         tx(nx  ,j,k)=acix6*(ux(1,j,k)-ux(nx,j,k))&
                    +bcix6*(ux(2,j,k)-ux(nx-1,j,k))
         rx(nx  ,j,k)=alcaix6
      enddo
      enddo
      do i=2,nx
      do k=1,nz
      do j=1,ny
         tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*csx6(i)
         rx(i,j,k)=rx(i,j,k)-rx(i-1,j,k)*csx6(i)
      enddo
      enddo
      enddo
      do k=1,nz
      do j=1,ny
         tx(nx,j,k)=tx(nx,j,k)*cwx6(nx)
         rx(nx,j,k)=rx(nx,j,k)*cwx6(nx)
      enddo
      enddo
      do i=nx-1,1,-1
      do k=1,nz
      do j=1,ny
         tx(i,j,k)=(tx(i,j,k)-cfx6(i)*tx(i+1,j,k))*cwx6(i)
         rx(i,j,k)=(rx(i,j,k)-cfx6(i)*rx(i+1,j,k))*cwx6(i)
      enddo
      enddo
      enddo
      do k=1,nz
      do j=1,ny
         sx(j,k)=(tx(1,j,k)-alcaix6*tx(nx,j,k))/&
              (1.+rx(1,j,k)-alcaix6*rx(nx,j,k))
      enddo
      enddo
      do k=1,nz
      do j=1,ny
      do i=1,nx
         tx(i,j,k)=tx(i,j,k)-sx(j,k)*rx(i,j,k)
      enddo
      enddo
      enddo
   endif
!
end subroutine deci60 
!
!********************************************************************
!********************************************************************
!
subroutine interi60(tx,ux,rx,sx,cifx6,cisx6,ciwx6,nx,nxm,mx,ny,nym,nz,npaire) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
  USE paramx_m 
  USE interpol6_m 
!...Translated by PSUITE Trans90                  4.3ZH  9:45:39   1/26/ 4  
!...Switches: -yv INDDO=0 -nbejkno
  implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
   integer :: nx,nxm,mx,ny,nym,nz,npaire
   real(8),dimension(nx,nym,nz) :: tx
   real(8),dimension(nx,nym,nz) :: ux
   real(8),dimension(nx,ny,nz) :: rx
   real(8),dimension(ny,nz) :: sx
   real(8),dimension(nx) :: cifx6,cisx6,ciwx6
!
   integer :: i,j,k
!
   if (nclx==0) then
      do k=1,nz
      do j=1,ny
         tx(1,j,k)=aicix6*(ux(2,j,k)+ux(1  ,j,k))&
              +bicix6*(ux(3,j,k)+ux(nx,j,k))
         rx(1,j,k)=-1.
         tx(2,j,k)=aicix6*(ux(3,j,k)+ux(2 ,j,k))&
              +bicix6*(ux(4,j,k)+ux(1,j,k))
         rx(2,j,k)=0.
      enddo
      enddo
      do i=3,nx-2
      do k=1,nz
      do j=1,ny
         tx(i,j,k)=aicix6*(ux(i+1,j,k)+ux(i,j,k))&
              +bicix6*(ux(i+2,j,k)+ux(i-1,j,k))
         rx(i,j,k)=0.
      enddo
      enddo
      enddo
      do k=1,nz
      do j=1,ny
         tx(nx-1,j,k)=aicix6*(ux(nx,j,k)+ux(nx-1,j,k))&
                    +bicix6*(ux(1 ,j,k)+ux(nx-2,j,k))
         rx(nx-1,j,k)=0.
         tx(nx  ,j,k)=aicix6*(ux(1,j,k)+ux(nx,j,k))&
              +bicix6*(ux(2,j,k)+ux(nx-1,j,k))
         rx(nx  ,j,k)=ailcaix6
      enddo
      enddo
      do i=2,nx
      do k=1,nz
      do j=1,ny
         tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*cisx6(i)
         rx(i,j,k)=rx(i,j,k)-rx(i-1,j,k)*cisx6(i)
      enddo
      enddo
      enddo
      do k=1,nz
      do j=1,ny
         tx(nx,j,k)=tx(nx,j,k)*ciwx6(nx)
         rx(nx,j,k)=rx(nx,j,k)*ciwx6(nx)
      enddo
      enddo
      do i=nx-1,1,-1
      do k=1,nz
      do j=1,ny
         tx(i,j,k)=(tx(i,j,k)-cifx6(i)*tx(i+1,j,k))*ciwx6(i)
         rx(i,j,k)=(rx(i,j,k)-cifx6(i)*rx(i+1,j,k))*ciwx6(i)
      enddo
      enddo
      enddo
      do k=1,nz
      do j=1,ny
         sx(j,k)=(tx(1,j,k)-ailcaix6*tx(nx,j,k))/&
                (1.+rx(1,j,k)-ailcaix6*rx(nx,j,k))
      enddo
      enddo
      do k=1,nz
      do j=1,ny
      do i=1,nx
         tx(i,j,k)=tx(i,j,k)-sx(j,k)*rx(i,j,k)
      enddo
      enddo
      enddo
   endif
   return
end subroutine interi60
!********************************************************************
!
!********************************************************************
!
subroutine deciy60(ty,uy,ry,di,sy,cfy6,csy6,cwy6,nx,ny,nym,my,nz,npaire) 
!...Translated by PSUITE Trans90                  4.3ZH  9:45:39   1/26/ 4  
!...Switches: -yv INDDO=0 -nbejkno
  implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
   integer :: nx,ny,my,nz,npaire,nym
   real(8),dimension(nx,ny,nz) :: ty
   real(8),dimension(nx,nym,nz) :: uy
   real(8),dimension(nx,nz,ny) :: ry
   real(8),dimension(nx,nz,nym) :: di
   real(8),dimension(nx,nz) ::sy
   real(8),dimension(ny) :: cfy6,csy6,cwy6
!
   integer :: i,j,k
!
      do k=1,nz
      do j=1,nym
      do i=1,nx
         di(i,k,j)=uy(i,j,k)
      enddo
      enddo
      enddo
! 
      call deciy160(ry,di,ty,sy,cfy6,csy6,cwy6,nx,ny,my,nz,npaire)
!
      do k=1,nz
      do j=1,ny
      do i=1,nx
         ty(i,j,k)=ry(i,k,j)
      enddo
      enddo
      enddo
!
      return
end subroutine deciy60 
!
!********************************************************************
!
subroutine deciy160(ty,uy,ry,sy,cfy6,csy6,cwy6,nx,ny,nym,my,nz,npaire) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
  USE dericey6_m 
  USE paramy_m 
!...Translated by PSUITE Trans90                  4.3ZH  9:45:39   1/26/ 4  
!...Switches: -yv INDDO=0 -nbejkno
  implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
   integer :: nx,ny,my,nz,npaire,nym
   real(8),dimension(nx,nz,nym) :: ty
   real(8),dimension(nx,nz,ny) :: uy
   real(8),dimension(nx,nz,ny) :: ry
   real(8),dimension(nx,nz) :: sy
   real(8),dimension(nym) :: cfy6,csy6,cwy6
!
   integer :: i,j,k
!   
   if (ncly==0) then
      do k=1,nz 
      do i=1,nx 
         ty(i,k,1)=aciy6*(uy(i,k,1)-uy(i,k,ny))&
                  +bciy6*(uy(i,k,2)-uy(i,k,ny-1)) 
         ry(i,k,1)=-1.
         ty(i,k,2)=aciy6*(uy(i,k,2)-uy(i,k,1))&
                  +bciy6*(uy(i,k,3)-uy(i,k,ny)) 
         ry(i,k,2)=0. 
      enddo
      enddo 
      do j=3,ny-2 
      do k=1,nz 
      do i=1,nx 
         ty(i,k,j)=aciy6*(uy(i,k,j)-uy(i,k,j-1))&
                  +bciy6*(uy(i,k,j+1)-uy(i,k,j-2)) 
         ry(i,k,j)=0. 
      enddo
      enddo 
      enddo
      do k=1,nz 
      do i=1,nx 
         ty(i,k,ny-1)=aciy6*(uy(i,k,ny-1)-uy(i,k,ny-2))&
                     +bciy6*(uy(i,k,ny)-uy(i,k,ny-3)) 
         ry(i,k,ny-1)=0. 
         ty(i,k,ny)=aciy6*(uy(i,k,ny)-uy(i,k,ny-1))&
                   +bciy6*(uy(i,k,1)-uy(i,k,ny-2)) 
         ry(i,k,ny)=alcaiy6
      enddo 
      enddo 
      do j=2,ny 
      do k=1,nz 
      do i=1,nx 
         ty(i,k,j)=ty(i,k,j)-ty(i,k,j-1)*csy6(j) 
         ry(i,k,j)=ry(i,k,j)-ry(i,k,j-1)*csy6(j) 
      enddo 
      enddo 
      enddo
      do k=1,nz 
      do i=1,nx 
         ty(i,k,ny)=ty(i,k,ny)*cwy6(ny) 
         ry(i,k,ny)=ry(i,k,ny)*cwy6(ny) 
      enddo 
      enddo 
      do j=ny-1,1,-1 
      do k=1,nz 
      do i=1,nx 
         ty(i,k,j)=(ty(i,k,j)-cfy6(j)*ty(i,k,j+1))*cwy6(j) 
         ry(i,k,j)=(ry(i,k,j)-cfy6(j)*ry(i,k,j+1))*cwy6(j) 
      enddo 
      enddo 
      enddo
      do k=1,nz 
      do i=1,nx 
         sy(i,k)=(ty(i,k,1)-alcaiy6*ty(i,k,ny))&
                /(1.+ry(i,k,1)-alcaiy6*ry(i,k,ny)) 
      enddo 
      enddo 
      do j=1,ny 
      do k=1,nz 
      do i=1,nx
         ty(i,k,j)=ty(i,k,j)-sy(i,k)*ry(i,k,j) 
      enddo 
      enddo
      enddo 
  endif
!
   return
end subroutine deciy160
!
!********************************************************************
!
!********************************************************************
subroutine interiy60(ty,uy,ry,di,sy,cify6,cisy6,ciwy6,nx,nxm,ny,nym,my,nz,npaire) 
!...Translated by PSUITE Trans90                  4.3ZH  9:45:39   1/26/ 4  
!...Switches: -yv INDDO=0 -nbejkno
  implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
   integer :: nx,nxm,ny,my,nz,npaire,nym
   real(8),dimension(nxm,ny,nz) :: ty
   real(8),dimension(nxm,nym,nz) :: uy
   real(8),dimension(nxm,nz,ny) :: ry
   real(8),dimension(nxm,nz,nym) :: di
   real(8),dimension(nx,nz) :: sy
   real(8),dimension(ny) :: cify6,cisy6,ciwy6
!
   integer :: i,j,k
!
      do k=1,nz
      do j=1,nym
      do i=1,nxm
         di(i,k,j)=uy(i,j,k)
      enddo
      enddo
      enddo
! 
      call interiy160(ry,di,ty,sy,cify6,cisy6,ciwy6,nx,nxm,ny,nym,my,nz,npaire)
!
      do k=1,nz
      do j=1,ny
      do i=1,nxm
         ty(i,j,k)=ry(i,k,j)
      enddo
      enddo
      enddo
!
      return 
end subroutine interiy60
!
!********************************************************************
!
subroutine interiy160(ty,uy,ry,sy,cify6,cisy6,ciwy6,nx,nxm,ny,nym,my,nz,npaire) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
  USE interpoly6_m 
  USE paramy_m 
!...Translated by PSUITE Trans90                  4.3ZH  9:45:39   1/26/ 4  
!...Switches: -yv INDDO=0 -nbejkno
  implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
   integer :: nx,nxm,ny,my,nz,npaire,nym
   real(8),dimension(nxm,nz,ny) :: ty
   real(8),dimension(nxm,nz,nym) :: uy
   real(8),dimension(nxm,nz,ny) :: ry
   real(8),dimension(nx,nz) :: sy
   real(8),dimension(ny) :: cify6,cisy6,ciwy6
!
   integer :: i,j,k
!  
   if (ncly==0) then
      do k=1,nz
      do i=1,nxm
         ty(i,k,1)=aiciy6*(uy(i,k,1)+uy(i,k,ny))&
              +biciy6*(uy(i,k,2)+uy(i,k,ny-1))
         ry(i,k,1)=-1.
         ty(i,k,2)=aiciy6*(uy(i,k,2)+uy(i,k,1))&
              +biciy6*(uy(i,k,3)+uy(i,k,ny))
         ry(i,k,2)=0.
      enddo
      enddo
      do j=3,ny-2
      do k=1,nz
      do i=1,nxm
         ty(i,k,j)=aiciy6*(uy(i,k,j)+uy(i,k,j-1))&
              +biciy6*(uy(i,k,j+1)+uy(i,k,j-2))
         ry(i,k,j)=0.
      enddo
      enddo
      enddo
      do k=1,nz
      do i=1,nxm
         ty(i,k,ny-1)=aiciy6*(uy(i,k,ny-1)+uy(i,k,ny-2))&
                    +biciy6*(uy(i,k,ny)+uy(i,k,ny-3))
         ry(i,k,ny-1)=0.
         ty(i,k,ny)=aiciy6*(uy(i,k,ny)+uy(i,k,ny-1))&
              +biciy6*(uy(i,k,1)+uy(i,k,ny-2))
         ry(i,k,ny)=ailcaiy6
      enddo
      enddo
      do j=2,ny
      do k=1,nz
      do i=1,nxm
         ty(i,k,j)=ty(i,k,j)-ty(i,k,j-1)*cisy6(j)
         ry(i,k,j)=ry(i,k,j)-ry(i,k,j-1)*cisy6(j)
      enddo
      enddo
      enddo
      do k=1,nz
      do i=1,nxm
         ty(i,k,ny)=ty(i,k,ny)*ciwy6(ny)
         ry(i,k,ny)=ry(i,k,ny)*ciwy6(ny)
      enddo
      enddo
      do j=ny-1,1,-1
      do k=1,nz
      do i=1,nxm
         ty(i,k,j)=(ty(i,k,j)-cify6(j)*ty(i,k,j+1))*ciwy6(j)
         ry(i,k,j)=(ry(i,k,j)-cify6(j)*ry(i,k,j+1))*ciwy6(j)
      enddo
      enddo
      enddo
      do k=1,nz
      do i=1,nxm
         sy(i,k)=(ty(i,k,1)-ailcaiy6*ty(i,k,ny))/&
                (1.+ry(i,k,1)-ailcaiy6*ry(i,k,ny))
      enddo
      enddo
      do k=1,nz
      do j=1,ny
      do i=1,nxm
         ty(i,k,j)=ty(i,k,j)-sy(i,k)*ry(i,k,j)
      enddo
      enddo
      enddo
      
   endif
!
   return
end subroutine interiy160
!
!********************************************************************
!
subroutine derx_r(tx,ux,ffx,fsx,fwx,nx,ny,nz,npaire) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
  USE paramx_m
  USE derivex_m
  implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  integer , intent(in) :: nx 
  integer , intent(in) :: ny 
  integer , intent(in) :: nz 
  integer , intent(in) :: npaire 
  real , intent(inout) :: tx(nx,ny,nz) 
  real(8) , intent(in) :: ux(nx,ny,nz) 
  real(8) , intent(in) :: ffx(nx) 
  real(8) , intent(in) :: fsx(nx) 
  real(8) , intent(in) :: fwx(nx) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer :: k, j, i 
!-----------------------------------------------
!
!********************************************************************
!
   if (nclx==2) then 
      do k=1,nz 
      do j=1,ny 
         tx(1,j,k)=af1x*ux(1,j,k)+bf1x*ux(2,j,k)+cf1x*ux(3,j,k) 
         tx(2,j,k)=af2x*(ux(3,j,k)-ux(1,j,k)) 
      enddo 
      enddo 
      do i=3,nx-2 
      do k=1,nz 
      do j=1,ny 
         tx(i,j,k)=afix*(ux(i+1,j,k)-ux(i-1,j,k))&
                  +bfix*(ux(i+2,j,k)-ux(i-2,j,k)) 
      enddo 
      enddo 
      enddo 
      do k=1,nz 
      do j=1,ny 
         tx(nx-1,j,k)=afmx*(ux(nx,j,k)-ux(nx-2,j,k)) 
         tx(nx,j,k)=(-afnx*ux(nx,j,k))-bfnx*ux(nx-1,j,k)-cfnx*ux(nx-2,j,k) 
      enddo
      enddo
      do i=2,nx 
      do k=1,nz 
      do j=1,ny 
         tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*fsx(i) 
      enddo
      enddo
      enddo
      do k=1,nz 
      do j=1,ny 
         tx(nx,j,k)=tx(nx,j,k)*fwx(nx) 
      enddo 
      enddo 
      do i=nx-1,1,-1 
      do k=1,nz 
      do j=1,ny 
         tx(i,j,k)=(tx(i,j,k)-ffx(i)*tx(i+1,j,k))*fwx(i) 
      enddo 
      enddo 
      enddo 
   endif 
!
end subroutine derx_r 
!
!********************************************************************
!
subroutine dery_r(ty,uy,ffy,fsy,fwy,nx,ny,nz,npaire) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
  USE paramy_m 
  USE derivey_m 
  implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  integer , intent(in) :: nx 
  integer , intent(in) :: ny 
  integer , intent(in) :: nz 
  integer , intent(in) :: npaire 
  real(8) , intent(inout) :: ty(nx,nz,ny) 
  real(8) , intent(in) :: uy(nx,nz,ny) 
  real(8) , intent(in) :: ffy(ny) 
  real(8) , intent(in) :: fsy(ny) 
  real(8) , intent(in) :: fwy(ny) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer :: k, j, i 
!-----------------------------------------------
!
!********************************************************************
!
   if (ncly==1) then 
      if (npaire==1) then 
         do k=1,nz 
         do i=1,nx 
            ty(i,k,1)=0. 
            ty(i,k,2)=afjy*(uy(i,k,3)-uy(i,k,1))&
                     +bfjy*(uy(i,k,4)-uy(i,k,2)) 
         enddo 
         enddo 
         do k=1,nz 
         do i=1,nx 
         do j=3,ny-2 
            ty(i,k,j)=afjy*(uy(i,k,j+1)-uy(i,k,j-1))&
                     +bfjy*(uy(i,k,j+2)-uy(i,k,j-2)) 
         enddo 
         enddo 
         enddo 
         do k=1,nz 
         do i=1,nx 
            ty(i,k,ny-1)=afjy*(uy(i,k,ny)-uy(i,k,ny-2))&
                        +bfjy*(uy(i,k,ny-1)-uy(i,k,ny-3)) 
            ty(i,k,ny)=0. 
         enddo 
         enddo 
         do j=2,ny 
         do k=1,nz 
         do i=1,nx 
            ty(i,k,j)=ty(i,k,j)-ty(i,k,j-1)*fsy(j) 
         enddo 
         enddo 
         enddo 
         do k=1,nz 
         do i=1,nx 
            ty(i,k,ny)=ty(i,k,ny)*fwy(ny) 
         enddo 
         enddo 
         do j=ny-1,1,-1 
         do k=1,nz 
         do i=1,nx 
            ty(i,k,j)=(ty(i,k,j)-ffy(j)*ty(i,k,j+1))*fwy(j) 
         enddo 
         enddo 
         enddo 
      endif
   endif
!
!********************STRET STRET************************
   do k=1,nz 
   do j=1,ny 
   do i=1,nx 
      ty(i,j,k)=ty(i,k,j)*ylyr
   enddo 
   enddo
   enddo 
!*******************************************************
!
   return  
end subroutine dery_r
!
!********************************************************************
!
subroutine derx_rg(tx,ux,ffx,fsx,fwx,nx,ny,nz,npaire) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
  USE paramx_m
  USE derivex_m
  implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  integer , intent(in) :: nx 
  integer , intent(in) :: ny 
  integer , intent(in) :: nz 
  integer , intent(in) :: npaire 
  real , intent(inout) :: tx(nx,ny,nz) 
  real(8) , intent(in) :: ux(nx,ny,nz) 
  real(8) , intent(in) :: ffx(nx) 
  real(8) , intent(in) :: fsx(nx) 
  real(8) , intent(in) :: fwx(nx) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer :: k, j, i 
!-----------------------------------------------
!
!********************************************************************
!
   if (nclx==2) then 
      do k=1,nz 
      do j=1,ny 
         tx(1,j,k)=af1x_rg*ux(1,j,k)+bf1x_rg*ux(2,j,k)+cf1x_rg*ux(3,j,k) 
         tx(2,j,k)=af2x_rg*(ux(3,j,k)-ux(1,j,k)) 
      enddo 
      enddo 
      do i=3,nx-2 
      do k=1,nz 
      do j=1,ny 
         tx(i,j,k)=afix_rg*(ux(i+1,j,k)-ux(i-1,j,k))&
                  +bfix_rg*(ux(i+2,j,k)-ux(i-2,j,k)) 
      enddo 
      enddo 
      enddo 
      do k=1,nz 
      do j=1,ny 
         tx(nx-1,j,k)=afmx_rg*(ux(nx,j,k)-ux(nx-2,j,k)) 
         tx(nx,j,k)=(-afnx_rg*ux(nx,j,k))-bfnx_rg*ux(nx-1,j,k)-cfnx_rg*ux(nx-2,j,k) 
      enddo
      enddo
      do i=2,nx 
      do k=1,nz 
      do j=1,ny 
         tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*fsx(i) 
      enddo
      enddo
      enddo
      do k=1,nz 
      do j=1,ny 
         tx(nx,j,k)=tx(nx,j,k)*fwx(nx) 
      enddo 
      enddo 
      do i=nx-1,1,-1 
      do k=1,nz 
      do j=1,ny 
         tx(i,j,k)=(tx(i,j,k)-ffx(i)*tx(i+1,j,k))*fwx(i) 
      enddo 
      enddo 
      enddo 
   endif 
!
end subroutine derx_rg 
!
!********************************************************************
!
subroutine dery_rg(ty,uy,ffy,fsy,fwy,nx,ny,nz,npaire) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
  USE paramy_m 
  USE derivey_m 
  implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  integer , intent(in) :: nx 
  integer , intent(in) :: ny 
  integer , intent(in) :: nz 
  integer , intent(in) :: npaire 
  real(8) , intent(inout) :: ty(nx,nz,ny) 
  real(8) , intent(in) :: uy(nx,nz,ny) 
  real(8) , intent(in) :: ffy(ny) 
  real(8) , intent(in) :: fsy(ny) 
  real(8) , intent(in) :: fwy(ny) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer :: k, j, i 
!-----------------------------------------------
!
!********************************************************************
!
   if (ncly==1) then 
      if (npaire==1) then 
         do k=1,nz 
         do i=1,nx 
            ty(i,k,1)=0. 
            ty(i,k,2)=afjy_rg*(uy(i,k,3)-uy(i,k,1))&
                     +bfjy_rg*(uy(i,k,4)-uy(i,k,2)) 
         enddo 
         enddo 
         do k=1,nz 
         do i=1,nx 
         do j=3,ny-2 
            ty(i,k,j)=afjy_rg*(uy(i,k,j+1)-uy(i,k,j-1))&
                     +bfjy_rg*(uy(i,k,j+2)-uy(i,k,j-2)) 
         enddo 
         enddo 
         enddo 
         do k=1,nz 
         do i=1,nx 
            ty(i,k,ny-1)=afjy_rg*(uy(i,k,ny)-uy(i,k,ny-2))&
                        +bfjy_rg*(uy(i,k,ny-1)-uy(i,k,ny-3)) 
            ty(i,k,ny)=0. 
         enddo 
         enddo 
         do j=2,ny 
         do k=1,nz 
         do i=1,nx 
            ty(i,k,j)=ty(i,k,j)-ty(i,k,j-1)*fsy(j) 
         enddo 
         enddo 
         enddo 
         do k=1,nz 
         do i=1,nx 
            ty(i,k,ny)=ty(i,k,ny)*fwy(ny) 
         enddo 
         enddo 
         do j=ny-1,1,-1 
         do k=1,nz 
         do i=1,nx 
            ty(i,k,j)=(ty(i,k,j)-ffy(j)*ty(i,k,j+1))*fwy(j) 
         enddo 
         enddo 
         enddo 
      endif
   endif
!
!********************STRET STRET************************
   do k=1,nz 
   do j=1,ny 
   do i=1,nx 
      ty(i,j,k)=ty(i,k,j)*ylyr
   enddo 
   enddo
   enddo 
!*******************************************************
!
   return  
end subroutine dery_rg
!
