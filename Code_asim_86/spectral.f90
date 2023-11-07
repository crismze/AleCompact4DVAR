! ***********************************************************
!
subroutine waves (mx,my,mz,xkx,yky,zkz,xk2,yk2,zk2,fk2,&
     knonul)
!
!***********************************************************
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
  USE derivex_m 
  USE derivey_m 
  USE derivez_m 
  USE paramx_m
  USE paramy_m
  USE paramz_m
!
  implicit none
!
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  integer  :: mx 
  integer  :: my 
  integer  :: mz 
  real(8),dimension(mx) :: xkx,xk2
  real(8),dimension(my) :: yky,yk2 
  real(8),dimension(mz) :: zkz,zk2
  real(8) , intent(inout) :: fk2(my,mz) 
  logical , intent(out) :: knonul(my,mz)   
!
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer :: i,ip,j,jj,k,kk,jp 
  real(8) :: twopi,w,wp 

   10 format(i4,4f13.6)
   20 format(a14,i4,f13.6,i4,f13.6,i4,f13.6)
!   30 format(a14,i4,f13.6)
!
      twopi=2.*acos(-1.)
!
      xkx(:)=0. ; xk2(:)=0. ; yky(:)=0. ; yk2(:)=0.
      zkz(:)=0. ; zk2(:)=0.
!
      if (nclx.eq.1) then
         do i=1,mx-3
!            print *,'trtr',dx,afjx,bfjx
            w=twopi*0.5*(i-1)/(mx-3)
            wp=afix*2.*dx*sin(w)+(bfix*2.*dx)*sin(2.*w)
            wp=wp/(1.+2.*alfaix*cos(w))
!            wp=w
            xkx(i)=(mx-3)*wp/xlx
            xk2(i)=xkx(i)*xkx(i)
            write(*,10) i,w,xkx(i),(mx-3)*w/xlx
         enddo
      endif
      if (nclx.eq.0) then
         do i=1,mx-2,2
            ip=i+1
            w=twopi*0.5*(i-1)/(mx-2)
            wp=afix*2.*dx*sin(w)+bfix*2.*dx*sin(2.*w)
            wp=wp/(1.+2.*alfaix*cos(w))
            xkx(i)=(mx-2)*wp/xlx
            xkx(ip)=xkx(i)
            xk2(i)=xkx(i)*xkx(i)
            xk2(ip)=xk2(i) 
            write(*,10) i,w,xkx(i),(mx-2)*w/xlx
         enddo
         xkx(mx-1)=0.
         xk2(mx-1)=0.
         xkx(mx  )=0.
         xk2(mx  )=0.
         if (mz.gt.1) then
            do j=1,my-3
               jj=j-1
               if (j.gt.(my-3)/2) jj=jj-(my-3)
               if (j.eq.(my-3)/2+1) jj=0
               w=twopi*jj/(my-3)
               wp=afjy*2.*dy*sin(w)+bfjy*2.*dy*sin(2.*w)
               wp=wp/(1.+2.*alfajy*cos(w))
               yky(j)=(my-3)*wp/yly
               yk2(j)=yky(j)*yky(j)
               write(*,10) j,w,yky(j),(my-3)*w/yly
            enddo
            yky(my-2)=0.
            yk2(my-2)=0.
            yky(my-1)=0.
            yk2(my-1)=0.
            yky(my  )=0.
            yk2(my  )=0.
            if (nclz.eq.0) then
               do k=1,mz-2
                  kk=k-1
                  if (k.gt.(mz-2)/2) kk=kk-(mz-2)
                  if (k.eq.(mz-2)/2+1) kk=0
                  w=twopi*kk/(mz-2)
                  wp=afkz*2.*dz*sin(w)+(bfkz*2.*dz)*sin(2.*w)
                  wp=wp/(1.+2.*alfakz*cos(w))
                  zkz(k)=(mz-2)*wp/zlz
                  zk2(k)=zkz(k)*zkz(k)
                  write(*,10) k,w,zkz(k),(mz-2)*w/zlz
               enddo
            endif
            if (nclz.eq.1) then
               do k=1,mz-2
                  w=twopi*0.5*(k-1)/(mz-2)
                  wp=afkz*2.*dz*sin(w)+bfkz*2.*dz*sin(2.*w)
                  wp=wp/(1.+2.*alfakz*cos(w))
                  zkz(k)=(mz-2)*wp/zlz
                  zk2(k)=zkz(k)*zkz(k)
                  write(*,10) k,w,zkz(k),(mz-2)*w/zlz
               enddo
            endif
            zkz(mz-1)=0.
            zk2(mz-1)=0.
            zkz(mz  )=0.
            zk2(mz  )=0.
         else
            if (ncly.eq.0) then
               do j=1,my-3
                  jj=j-1
                  if (j.gt.(my-3)/2) jj=jj-(my-3)
                  if (j.eq.(my-3)/2+1) jj=0
                  w=twopi*jj/(my-3)
                  wp=afjy*2.*dy*sin(w)+(bfjy*2.*dy)*sin(2.*w)
                  wp=wp/(1.+2.*alfajy*cos(w))
                  yky(j)=(my-3)*wp/yly
                  yk2(j)=yky(j)*yky(j)
                  write(*,10) j,w,yky(j),(my-3)*w/yly
               enddo
               yky(my-2)=0.
               yk2(my-2)=0.
            endif
!***********par ici la bonne soupe*************** 
            if (ncly.eq.1) then
               do j=1,my-3
!                  w=twopi*0.5*(k-1)/(mz-2)
!                  wp=afkz*2.*dz*sin(w)+bfkz*2.*dz*sin(2.*w)
!                  wp=wp/(1.+2.*alfakz*cos(w))
!                  zkz(k)=(mz-2)*wp/zlz
!                  zk2(k)=zkz(k)*zkz(k)

                  w=twopi*0.5*(j-1)/(my-3)
                  wp=afjy*2.*dy*sin(w)+(bfjy*2.*dy)*sin(2.*w)
                  wp=wp/(1.+2.*alfajy*cos(w))
!                  wp=w
                  yky(j)=(my-3)*wp/yly
                  yk2(j)=yky(j)*yky(j)
                  write(*,10) j,w,yky(j),(my-3)*w/yly
               enddo
            endif
            yky(my-1)=0.
            yk2(my-1)=0.
            yky(my  )=0.
            yk2(my  )=0.
         endif
      else
         if (mz.gt.1) then
            do j=1,my-2,2
               jp=j+1
               w=twopi*0.5*(j-1)/(my-2)
               wp=afjy*2.*dy*sin(w)+bfjy*2.*dy*sin(2.*w)
               wp=wp/(1.+2.*alfajy*cos(w))
               yky(j)=(my-2)*wp/yly
               yky(jp)=yky(j)
               yk2(j)=yky(j)*yky(j)
               yk2(jp)=yk2(j)
               write(*,10) j,w,yky(j),(my-2)*w/yly
            enddo
            yky(my-1)=0.
            yk2(my-1)=0.
            yky(my  )=0.
            yk2(my  )=0.
            if (nclz.eq.0) then
               print *,'nclz.eq.0'
               do k=1,mz-3
                  kk=k-1
                  if (k.gt.(mz-3)/2) kk=kk-(mz-3)
                  if (k.eq.(mz-3)/2+1) kk=0
                  w=twopi*kk/(mz-3)
                  wp=afkz*2.*dz*sin(w)+(bfkz*2.*dz)*sin(2.*w)
                  wp=wp/(1.+2.*alfakz*cos(w))
                  zkz(k)=(mz-3)*wp/zlz
                  zk2(k)=zkz(k)*zkz(k)
               enddo
               zkz(mz-2)=0.
               zk2(mz-2)=0.
            endif
            if (nclz.eq.1) then
               do k=1,mz-3
                  w=twopi*0.5*(k-1)/(mz-3)
                  wp=afkz*2.*dz*sin(w)+(bfkz*2.*dz)*sin(2.*w)
                  wp=wp/(1.+2.*alfakz*cos(w))
                  zkz(k)=(mz-3)*wp/zlz
                  zk2(k)=zkz(k)*zkz(k)
                  write(*,10) k,w,zkz(k),(mz-3)*w/zlz
               enddo
               zkz(mz-2)=0.
               zk2(mz-2)=0.
            endif
            zkz(mz-1)=0.
            zk2(mz-1)=0.
            zkz(mz  )=0.
            zk2(mz  )=0.
         else
            if (ncly.eq.0) then
               do j=1,my-2,2
                  jp=j+1
                  w=twopi*0.5*(j-1)/(my-2)
                  wp=afjy*2.*dy*sin(w)+bfjy*2.*dy*sin(2.*w)
                  wp=wp/(1.+2.*alfajy*cos(w))
                  yky(j)=(my-2)*wp/yly
                  yky(jp)=yky(j)
                  yk2(j)=yky(j)*yky(j)
                  yk2(jp)=yk2(j) 
                  write(*,10) j,w,yky(j),(my-2)*w/yly
               enddo
            endif
            if (ncly.eq.1) then
               do j=1,my-2
                  w=twopi*0.5*(j-1)/(my-2)
                  wp=afjy*2.*dy*sin(w)+bfjy*2.*dy*sin(2.*w)
                  wp=wp/(1.+2.*alfajy*cos(w))
                  yky(j)=(my-2)*wp/yly
                  yk2(j)=yky(j)*yky(j)
                  write(*,10) j,w,yky(j),(my-2)*w/yly
               enddo
            endif
            yky(my-1)=0.
            yk2(my-1)=0.
            yky(my  )=0.
            yk2(my  )=0.
         endif 
      endif
!
      if (mz.gt.1) then
         do k=1,mz
         do j=1,my
            fk2(j,k)=-(yk2(j)+zk2(k))
            if (fk2(j,k).ne.0.) then
               knonul(j,k)=.true.
            else
               knonul(j,k)=.false.
               write(*,20) 'Mode exclu : ',j,yky(j),k,zkz(k)
            endif
         enddo
         enddo
      else
         do j=1,my
            fk2(j,1)=-yk2(j)
            if (fk2(j,1).ne.0.) then
               knonul(j,1)=.true.
            else
               knonul(j,1)=.false.
               write(*,20) 'Mode exclu : ',j,yky(j)
            endif
         enddo
      endif
!
      return
      end
