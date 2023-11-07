! ***************************************************************
!
      subroutine fftid (a,mx,my)
!
! ***************************************************************
!
!     inverse 2d fourier transform (x,y)
!     spectral ----> physical space
!     isign=1
!      
      include 'precision.inc'
!
      dimension a(mx,my)
!
      nx=mx-2
      ny=my-3
      nxhp=nx/2+1
!
      call fft2d (a,nxhp,mx,my,1)
      call fft1d (a,mx,my,1)
!
      return
      end
!
! ***************************************************************
!
      subroutine fftdd (a,mx,my)
!
! ***************************************************************
!
!     direct 2d fourier transform (x,y)
!     physical ----> spectral space
!     isign=-1
!      
      include 'precision.inc'
!
      dimension a(mx,my)
!
      nx=mx-2
      ny=my-3
      nxhp=nx/2+1
!
      call fft1d (a,mx,my,-1)
      call fft2d (a,nxhp,mx,my,-1)
!
      return
      end
!
! ***************************************************************
!
      subroutine cfftid (a,mx,my)
!
! ***************************************************************
!
!     inverse cosine 2d fourier transform (x,y)
!     spectral ----> physical space
!     isign=1
!      
      include 'precision.inc'
!
      dimension a(mx,my),b(my,mx)
!
      nx=my-2
      ny=mx-3
      nxhp=ny/2+1
!
      do j=1,my
      do i=1,mx
         b(j,i)=a(i,j)
      enddo
      enddo
      call fft1d  (b,my,mx,1)
      call cfft2d (b,my,mx,1)
      do j=1,my
      do i=1,mx
         a(i,j)=b(j,i)
      enddo
      enddo
!
       nx=mx-3
       ny=my-2
        do j=1,my
!          a(nx+1,j)=0.
!          a(nx+2,j)=0.
!          a(mx,j)=0.
       enddo
       do i=1,mx
!          a(i,ny+1)=0.
!          a(i,ny+2)=0.
!          a(i,my)=0.
       enddo
      return
      end
!
! ***************************************************************
!
      subroutine cfftdd (a,mx,my)
!
! ***************************************************************
!
!     direct cosine 2d fourier transform (x,y)
!     physical ----> spectral space
!     isign=-1
!      
      include 'precision.inc'
!
      dimension a(mx,my),b(my,mx)
!
      nx=my-2
      ny=mx-3
      nxhp=ny/2+1
!
      do j=1,my
      do i=1,mx
         b(j,i)=a(i,j)
      enddo
      enddo
      call cfft2d (b,my,mx,-1)
      call fft1d  (b,my,mx,-1)
      do j=1,my
      do i=1,mx
         a(i,j)=b(j,i)
      enddo
      enddo
       nx=mx-3
       ny=my-2
        do j=1,my
!          a(nx+1,j)=0.
!          a(nx+2,j)=0.
!          a(mx,j)=0.
       enddo
       do i=1,mx
!          a(i,ny+1)=0.
!          a(i,ny+2)=0.
!          a(i,my)=0.
       enddo
!
      return
      end
!
! ************************************************************ 
!
      subroutine mulmi (b,nwdsc)
!
! ************************************************************
!
      include 'precision.inc'
!
      complex b(nwdsc)
!                                  
      do i=1,nwdsc
         b(i)=(0.,-1.)*b(i)
      enddo
!
      return
      end
! 
! ************************************************************
!
      subroutine mulpi (b,nwdsc)
!
! ************************************************************
!
      include 'precision.inc'
!
      complex b(nwdsc)
!                                  
      do i=1,nwdsc
         b(i)=(0.,1.)*b(i)
      enddo
!
      return
      end
!
! ************************************************************************
!
      subroutine fft1d (b,my,mx,isign)
!
! ************************************************************************
!
!     fft reel-->complexe (2D)
!
      include 'precision.inc'
!
      integer,parameter :: n1=4*(n-1)
      integer, parameter :: n2=4+n1
      real, dimension(my,mx) :: a,b
      real, dimension(n2) :: table
      real, dimension(n1) :: work
! 
      ndim=my-2 
!
      if (isign == 1 ) then
         CALL SCFFT(0, ndim, 1., a(1,1), b(1,1), TABLE, WORK, 0)
         CALL SCFFT(1, ndim, 1./ndim, a(1,1),b(1,1), TABLE, WORK, 0)
      endif
      if (isign == -1 ) then
         CALL CSFFT(-1, ndim, 1., b(1,1), a(1,1), TABLE, WORK, 0)
      endif
!
      return
      end 
!
! ************************************************************************
!
      subroutine fft2d (a,nxhp,mx,my,isign)
!
! ***********************************************************************
!
!     fft complexe-->complexe (2D)
!
      include 'precision.inc'
!
      real a(2,nxhp,my) 
!
      common /fty/ifaxy(1)
      common /exy/trigsy(1)
      common /wrk/work(1) 
!
      inc=mx
      nft=nxhp
      jump=2
      ndim=my-3
!
      call cfft999 (a(1,1,1),a(2,1,1),work(1),
     1              trigsy,ifaxy,inc,jump,ndim,nft,isign)   
!   
      return
      end 
!
! ************************************************************************
!
      subroutine cfft2d (b,my,mx,isign)
!
! ***********************************************************************
!               
!     if isign.eq.-1: real   ---> cosine
!     if isign.eq. 1: cosine ---> real
!
      include 'precision.inc'
!
      real b(my,mx) 
!
      inc=my
      nft=my
      jump=1
      ndim=mx-3 
!
      if (isign == 1) then
         init = 'i'
     call c06hbf(1, ndim, b(1,1), init, table, work, ifail)
     do i=1,n+1
        a(i)=a(i)/(0.5*(n-1))
     enddo
      endif
!
      call fct99 (b(1,1),work(1),
     1            trigsx,trigsx2,workx1,ifaxx,
     2           inc,jump,ndim,nft,isign)   
!   
      return
      end  
!
