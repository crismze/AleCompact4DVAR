!******************************************************************
!
subroutine maxmin(ux,uy,uz,nx,ny,nz) 
!...Translated by PSUITE Trans90                  4.3ZH 15:01:47   1/27/ 4  
!...Switches: -nbejknoc -rl -xl   
  implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  integer , intent(in) :: nx 
  integer , intent(in) :: ny 
  integer , intent(in) :: nz 
  real , intent(in) :: ux(nx,ny,nz) 
  real , intent(in) :: uy(nx,ny,nz) 
  real , intent(in) :: uz(nx,ny,nz) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer :: k, j, i, imax, jmax, kmax, imin, jmin, kmin 
  real :: uxmax, uymax, uzmax, uxmin, uymin, uzmin 
!-----------------------------------------------
!
!******************************************************************
!
   uxmax=ux(1,1,1) 
   uymax=uy(1,1,1) 
   uzmax=uz(1,1,1) 
   uxmin=ux(1,1,1) 
   uymin=uy(1,1,1) 
   uzmin=uz(1,1,1) 
!
   do k=1,nz 
   do j=1,ny 
   do i=1,nx 
      if (ux(i,j,k) > uxmax) then 
         uxmax=ux(i,j,k) 
         imax=i 
         jmax=j 
         kmax=k 
      endif
      if (uy(i,j,k) > uzmax) then 
         uymax=uy(i,j,k)  
      endif
      if (uz(i,j,k) > uzmax) then 
         uzmax=uz(i,j,k)  
      endif 
      if (ux(i,j,k) < uxmin) then 
         uxmin=ux(i,j,k) 
         imin=i 
         jmin=j 
         kmin=k 
      endif
      if (uy(i,j,k) < uymin) then 
         uymin=uy(i,j,k) 
      endif
      if (uz(i,j,k) < uzmin) then 
         uzmin=uz(i,j,k) 
      endif
   enddo
   enddo 
   enddo 
!
   print *,'Ux,Uy,Uz max',uxmax,uymax,uzmax 
   print *,'i ,j ,k umax',imax,jmax,kmax 
   print *,'Ux,Uy,Uz min',uxmin,uymin,uzmin 
   print *,'i ,j ,k umin',imin,jmin,kmin 
! 
   return  
end subroutine maxmin
!
!******************************************************************
!
subroutine annule(a,n) 
!...Translated by PSUITE Trans90                  4.3ZH 15:01:47   1/27/ 4  
!...Switches: -nbejknoc -rl -xl   
  implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  integer , intent(in) :: n 
  real , intent(out) :: a(n) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer :: i 
!-----------------------------------------------
!
!******************************************************************
!
   a=0. 
!
   return  
end subroutine annule
!
!******************************************************************
!
real function alea (r) 
!...Translated by PSUITE Trans90                  4.3ZH 15:01:47   1/27/ 4  
!...Switches: -nbejknoc -rl -xl   
  implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  real , intent(in) :: r 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer :: ia1, ia0, ia1ma0, ic, ix1, ix0, iy0, iy1 
!-----------------------------------------------
!
!******************************************************************
!
   data ia1,ia0,ia1ma0/1536,1029,507/ 
   data ic/1731/ 
   data ix1,ix0/0,0/ 
!
!      PRINT *,'alea !!!'
!
   if (r >= 0.) then 
      if (r > 0.) go to 20 
!
!     A*X = 2**22*IA1*IX1 + 2**11*(IA1*IX1 + (IA1-IA0)*(IX0-IX1)
!         + IA0*IX0) + IA0*IX0
!
      iy0=ia0*ix0 
      iy1=ia1*ix1 + ia1ma0*(ix0 - ix1) + iy0 
      iy0=iy0 + ic 
      ix0=mod(iy0,2048) 
      iy1=iy1 + (iy0 - ix0)/2048 
      ix1=mod(iy1,2048) 
!
   endif
10 continue 
   alea=ix1*2048 + ix0 
   alea=alea/4194304. 
   return  
!
20 continue 
   ix1=amod(r,1.)*4194304. + 0.5 
   ix0=mod(ix1,2048) 
   ix1=(ix1 - ix0)/2048 
   go to 10 
!
end function alea
!
!******************************************************************
!
subroutine numcar(num,car) 
!...Translated by PSUITE Trans90                  4.3ZH 15:01:47   1/27/ 4  
!...Switches: -nbejknoc -rl -xl   
  implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  integer , intent(in) :: num 
  character , intent(out) :: car*3 
!-----------------------------------------------
!
!******************************************************************

!
  if (num >= 100) then 
     write (car,1) num 
1    format(i3) 
  else 
     if (num >= 10) then 
        write (car,2) num 
2       format('0',i2) 
     else 
        write (car,3) num 
3       format('00',i1) 
     endif
  endif
!
  return  
end subroutine numcar
!
!********************************************************************
!
real function epmach (dum) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
  USE value_m 
!...Translated by PSUITE Trans90                  4.3ZH 15:01:47   1/27/ 4  
!...Switches: -nbejknoc -rl -xl   
  implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  real  :: dum 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  real :: eps 
!-----------------------------------------------
!
!********************************************************************
!
   eps=1. 
   eps=eps/10. 
   call store (eps + 1.) 
   do while(v - 1. > 0.) 
      eps=eps/10. 
      call store (eps + 1.) 
   enddo
102 continue 
   epmach=100.*eps 
!
   return  
 end function epmach
!
!********************************************************************
!
 subroutine store(x) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
   USE value_m 
!...Translated by PSUITE Trans90                  4.3ZH 15:01:47   1/27/ 4  
!...Switches: -nbejknoc -rl -xl   
   implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
   real , intent(in) :: x 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!-----------------------------------------------
!
!********************************************************************
!
   v=x 
!
   return  
end subroutine store
!
!****************************************************
!
real function rand2 (idum) 
!
! ****************************************************
!
!
! =========================================================
!
!  "Long" period ( > 2.E+18 ) random number generator
!  of L'Ecuyer with Bays-Durham shuffle and added safeguards.
!  Returns a uniform random deviate between 0. and 1.
!  (exclusive of the endpoint values).
!
!  Call with IDUM, a negative integer to initialize;
!  thereafter, do not alter IDUM between succesive deviates
!  in a sequence.  RNMAX should approximate the largest
!  floating value that  is less than 1.
!
!  REF: NUMERICAL RECEPIES -  pp 271
!
! =========================================================
!
!...Translated by PSUITE Trans90                  4.3ZH 15:01:47   1/27/ 4  
!...Switches: -nbejknoc -rl -xl   
  implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  integer , intent(inout) :: idum 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
  integer, parameter :: im1 = 2147483563 
  integer, parameter :: im2 = 2147483399 
  integer, parameter :: imm1 = im1 - 1 
  integer, parameter :: ia1 = 40014 
  integer, parameter :: ia2 = 40692 
  integer, parameter :: iq1 = 53668 
  integer, parameter :: iq2 = 52774 
  integer, parameter :: ir1 = 12211 
  integer, parameter :: ir2 = 3791 
  integer, parameter :: ntab = 32 
  integer, parameter :: ndiv = 1 + imm1/ntab 
  real, parameter :: am = 1./im1 
  real, parameter :: eps = 1.2E-7 
  real, parameter :: rnmx = 1. - eps 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer :: idum2, j, k 
  integer , dimension(ntab) :: iv 
  integer :: iy 
!
  save idum2, iv, iy 
!-----------------------------------------------
!
!
!
   data idum2/123456789/ 
   data iv/ntab*0/ 
   data iy/0/ 
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
   if (idum <= 0) then                        !Initialize. 
      idum=max((-idum),1) 
      idum2=idum 
      do j=ntab + 8,1,-1 
         k=idum/iq1 
         idum=ia1*(idum - k*iq1) - k*ir1 
         if (idum < 0) idum=idum + im1 
         if (j > ntab) cycle  
         iv(j)=idum 
      enddo
   endif
!
   k=idum/iq1 
   idum=ia1*(idum - k*iq1) - k*ir1 
   if (idum < 0) idum=idum + im1 
!
   k=idum2/iq2 
   idum2=ia2*(idum2 - k*iq2) - k*ir2 
   if (idum2 < 0) idum2=idum2 + im2 
!
   j=1 + iy/ndiv 
   iy=iv(j) - idum2 
   iv(j)=idum 
   if (iy < 1) iy=iy + imm1 
!===========================
   rand2=min(am*iy,rnmx) 
!===========================
!
   return  
end function rand2
