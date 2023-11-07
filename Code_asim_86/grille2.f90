! **************************************************************
!
   program grille
!
! **************************************************************
!
     implicit none
!
     integer,parameter ::nx=128,ny=128,nz=1
     integer,parameter :: my=ny+1
!
     real,dimension(nx) :: xx(nx)
     real,dimension(ny) :: yy(ny)
     real,dimension(nz) :: zz(nz)
     integer :: i,j,k,ncly
     real :: xlx,yly,zlz,dx,dy,dz,x,y,z
     real,dimension(ny) :: ppy,pp2y,pp3y,pp4y,fux1,d1,d,hprime
     real,dimension(my) :: yp,yeta
     real :: yinf,beta,den,xnum,alpha,xcx,den1,den2,den3,den4,xnum1,cst
     real :: pi,fourpi
!
     xlx=36.
     yly=12.
     zlz=6.
!
     dx=xlx/(nx-1)
     dy=yly/ny
     dz=zlz/nz
!
     x=0.
     do i=1,nx
        xx(i)=x
        x=x+dx 
     enddo
!
     y=0.
     do i=1,ny
        yy(i)=y
        y=y+dy 
     enddo
!
     z=0.
     do i=1,nz
        zz(i)=z
        z=z+dz 
     enddo
!
     open (10,file='avs2.x',form='unformatted')
     write (10) (xx(i),i=1,nx)
     close (10)
!
!*************************************************************************
     ncly=1
   pi=acos(-1.)
   fourpi=4.*pi
!     yinf cest la borne inf du domaine que l'on veut obtenir
   yinf=-yly/2.
      
      beta=3.5
!      beta=2.
!      beta=0.5
      den=2.*beta*yinf
      xnum=-yinf-sqrt(pi*pi*beta*beta+yinf*yinf)
      alpha=abs(xnum/den)
      xcx=1./beta/alpha
      
!
!     c'est ma fonction de stretching
      if (alpha.ne.0.) then 
      do j=1,ny
      if (ncly.eq.0) yeta(j)=(j-0.5)*(1./ny)
      if (ncly.eq.1) yeta(j)=(j-0.5)*(1./(ny-1.))
      den1=sqrt(alpha*beta+1.)
      xnum=den1/sqrt(alpha/pi)/sqrt(beta)/sqrt(pi)     
      den=2.*sqrt(alpha/pi)*sqrt(beta)*pi*sqrt(pi)
      den3=((sin(pi*yeta(j)))*(sin(pi*yeta(j)))/beta/pi)+alpha/pi
      den4=2.*alpha*beta-cos(2.*pi*yeta(j))+1.
      xnum1=(atan(xnum*tan(pi*yeta(j))))*den4/den1/den3/den
      cst=sqrt(beta)*pi/(2.*sqrt(alpha)*sqrt(alpha*beta+1.))
      if (yeta(j).lt.0.5) yp(j)=xnum1-cst-yinf
      if (yeta(j).eq.0.5) yp(j)=0.-yinf
      if (yeta(j).gt.0.5) yp(j)=xnum1+cst-yinf
      enddo
      endif
      if (alpha.eq.0.) then
         yp(1)=-1.e10
         do j=2,ny
            yeta(j)=(j-1.)*(1./ny)
            yp(j)=-beta*cos(pi*yeta(j))/sin(yeta(j)*pi)
         enddo
      endif
      do j=2,ny
         print *,'yy',yp(j),'dif',yp(j)-yp(j-1)
      enddo
      print *, 'min',yp(ny/2+1)-yp(ny/2),'=?',dx
!      pause'yy'
      open (10,file='avs2.y',form='unformatted')
      write (10) (yp(j),j=1,ny)
      close (10)     
!*****************************************************************
!
     open (10,file='avs2.z',form='unformatted')
     write (10) (zz(i),i=1,nz)
     close (10)
!
   end program grille
