!*****************DEBUT DE LA BOUCLE EN XXXXXXXXXXXX******************
       do i=1,nxm
!*********************************************************************
!*****************petit tablo d onde dans b******************** 
      do j=1,ny/2
         d(j)=yky(2*j-1)*transy(i)
         d1(j)=yky(2*j)*transy(i)
      enddo 
      
      ps2d1=0. ; ps2d1=ps2d

         xps=ps2d(i,ny)
         e=0. ; c=0.
         do j=1,ny/2
            e(j)=ps2d(i,2*j-1)
            c(j)=ps2d(i,2*j)
!            print *,e(j),c(j)
         enddo
!         pause
!***************Construction de la matrice a**************************
      xvar=alpha+1./2./beta
      xmo=alpha+3./4./beta
      tw=twopi*twopi
      beta1=beta*beta
      xvar1=xvar*xvar
      a(1,m1+1)=-(d5(i)*d5(i)*transx(1)*transx(1))-d(1)*d(1)*xvar1/tw*4.&
           -(1./4./tw/beta1)*(d(1)*d(1+1))
      do j=2,ny/2-1
         a(j,m1+1)=-(d5(i)*d5(i)*transx(2*j-1)*transx(2*j-1))-d(j)*d(j)*xvar1/tw*4.&
              -(1./4./tw/beta1)*(d(j)*d(j-1)+d(j)*d(j+1))
      enddo
      a(ny/2,m1+1)=-(d5(i)*d5(i)*transx(ny-1)*transx(ny-1))-d(ny/2)*d(ny/2)*xvar1/tw*4.&
           -(1./4./tw/beta1)*(d(ny/2)*d(ny/2-1))
!****************element du cote des m1*******************************
      a(1,m1)=0.
      do j=2,ny/2
         a(j,m1)=xvar/tw/beta*(d(j)*d(j-1)+d(j-1)*d(j-1))
      enddo
!
      a(1,m1-1)=0.
      a(2,m1-1)=0.
      do j=3,ny/2
         a(j,m1-1)=-d(j-1)*d(j-2)/4./tw/beta1
      enddo
!****************element du cote des m2*******************************
      a(1,4)=xvar/tw/beta*(d(1)*d(2)+d(2)*d(2))
      do j=2,ny/2-1
         a(j,m1+m2)=xvar/tw/beta*(d(j)*d(j+1)+d(j+1)*d(j+1))
      enddo
      a(1,m1+m2)=a(1,m1+m2)*2.
!      a(ny/2-1,m1+m2)=xvar/tw/beta*(d(ny/2-1)*d(ny/2))
!      a(ny/2-1,m1+m2)=a(ny/2-1,m1+m2)+xmo/tw/beta*(d(ny/2)*d(ny/2))
      a(ny/2,m1+m2)=0.
!        
      do j=1,ny/2-2
         a(j,m1+m2+1)=-d(j+1)*d(j+2)/4./tw/beta1
      enddo
      a(1,m1+m2+1)=a(1,m1+m2+1)*2.
      a(ny/2-1,m1+m2+1)=0.
      a(ny/2,m1+m2+1)=0.
!********************************************************************
      if (d5(i).eq.0.) then 
         a(1,m1+1)=1.
         a(1,m1+m2)=0.
         a(1,m1+m2+1)=0.
!         a(ny/2,m1+1)=1.
!         a(ny/2,m1+m2)=0.
!         a(ny/2,m1+m2+1)=0.
      endif
!*********************************************************************
!*******construction de la matrice a2*********************************
      xvar=alpha+1./2./beta
      xmo=alpha+3./4./beta
      xmo1=alpha+1./4./beta
      tw=twopi*twopi
      beta1=beta*beta
      xvar1=xvar*xvar
      a2(1,m1+1)=-(d5(i)*d5(i))-d1(1)*d1(1)*xmo*xmo1/tw*4.&
           -(1./4./tw/beta1)*(d1(1)*d1(2))
      do j=2,ny/2-1
         a2(j,m1+1)=-(d5(i)*d5(i))-d1(j)*d1(j)*xvar1/tw*4.&
              -(1./4./tw/beta1)*(d1(j)*d1(j-1)+d1(j)*d1(j+1))
      enddo
      a2(ny/2,m1+1)=-(d5(i)*d5(i))-d1(ny/2)*d1(ny/2)*xmo*xmo/tw*4.&
           -(1./4./tw/beta1)*(d1(ny/2)*d1(ny/2-1))
!***************element du cote des m1********************************
      a2(1,m1)=0.
      a2(2,m1)=xvar/tw/beta*(d1(2)*d1(1))&
           +xmo/tw/beta*(d1(1)*d1(1))
      do j=3,ny/2-1
         a2(j,m1)=xvar/tw/beta*(d1(j)*d1(j-1)+d1(j-1)*d1(j-1))
      enddo
      a2(ny/2,m1)=xmo/tw/beta*(d1(ny/2)*d1(ny/2-1))&
           +xvar/tw/beta*(d1(ny/2-1)*d1(ny/2-1))
      a2(1,m1-1)=0.
      a2(2,m1-1)=0.
      do j=3,ny/2
         a2(j,m1-1)=-d1(j-1)*d1(j-2)/4./tw/beta1
      enddo
!****************element du cote des m2*********************************
      a2(1,m1+m2)=xmo1/tw/beta*(d1(1)*d1(2))&
           +xvar/tw/beta*(d1(2)*d1(2))
      do j=2,ny/2-2
         a2(j,m1+m2)=xvar/tw/beta*(d1(j)*d1(j+1)+d1(j+1)*d1(j+1))
      enddo
      a2(ny/2-1,m1+m2)=xvar/tw/beta*(d1(ny/2-1)*d1(ny/2))&
           +xmo/tw/beta*(d1(ny/2)*d1(ny/2))
      a2(ny/2,m1+m2)=0.
      do j=1,ny/2-2
         a2(j,m1+m2+1)=-d1(j+1)*d1(j+2)/4./tw/beta1
      enddo
      a2(ny/2-1,m1+m2+1)=0.
      a2(ny/2,m1+m2+1)=0.
!***********************************************************************
       call bandec(a,ny,my,m1,m2,al,indx,xd)
       call bandec(a2,ny,my,m1,m2,al1,indx1,xd1)
       call banbks(a,ny,my,m1,m2,al,indx,e)
       call banbks(a2,ny,my,m1,m2,al1,indx1,c)
!********************on reconstruit les a+ib****************************
       do j=1,ny-1,2
          us2d(i,j)=e((j+1)/2)
       enddo
       do j=2,ny,2
          us2d(i,j)=c(j/2)
       enddo 
!       do j=1,my
!          print *,j,ps2d1(i,j)
!       enddo
!***********************************************************************
       if (d5(i).ne.0.) then
          xtoto=xps+e(ny/2-1)*d(ny/2)*d(ny/2-1)/2./beta1/tw&
               -e(ny/2)*d(ny/2)*d(ny/2)*2.*xvar/tw/beta
          xtoto=-xtoto/d5(i)/d5(i)
       endif     
       if (d5(i).eq.0) then
          us2d(i,1)=0.
          us2d(i,ny)=0.
       endif
       
!!     us2d(i,ny)=xtoto
!       print *,xps
!       pause
!       if (i==1) us2d(1,ny)=-xps
!*******FIN BOUCLE EN XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXxxxxx*************
       enddo
!***********************************************************************
