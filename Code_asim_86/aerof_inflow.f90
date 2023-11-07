!***********************************************************************
!
subroutine aerof_inflow(ux,uy,pp,ffxp,fsxp,fwxp,ffyp,fsyp,fwyp,&
             ux_ref,uy_ref,pp_ref,xDrag,yLift,nx,ny,nz,&
             i_coe1,i_coe2,j_coe1,j_coe2,itime,nlock)
!
!***********************************************************************
!
!  Calcula drag, lift a partir de campos ux, uy, pp disponibles
!  solo en inflow (cara CD del volumen de control) 
!
   USE paramx_m
   USE paramy_m
   USE paramz_m
   USE param1_m
   USE barreau_m 
!
   implicit none
!
   real(8),dimension(nx,ny,nz) :: ux,uy,gx,gy,tx,ty
   real(8),dimension(nx,ny,nz) :: pp
   real(8),dimension(nx) :: ffxp,fsxp,fwxp
   real(8),dimension(ny) :: ffyp,fsyp,fwyp
   integer :: i,j,k,nx,ny,nz,itime,nlock
   integer :: i_coe1,i_coe2,j_coe1,j_coe2
   integer :: icvrt,jcvlw,jcvup
   real(8) :: ux_ref,uy_ref,pp_ref
   real(8) :: xmom,ymom,f2x,f2y
   real(8) :: xDrag,yLift,xcv
   integer :: longueur1
   character(len=120) dir
!
!  Directorio para guardar momentum and force terms
   dir = '/home/ale/Asim_86/aerof/'
   longueur1=index(dir,' ')-1
!
!  Definicion de la cara "inflow" del volumen de control
   xcv=0.5
   icvrt = i_coe1+nint(xcv*(i_coe2-i_coe1))
   jcvlw = j_coe1
   jcvup = j_coe2
!  
!  Calculation of the momentum terms
   call rmomentum(ux,uy,xmom,ymom,nx,ny,nz,ux_ref,uy_ref,icvrt,jcvlw,jcvup,&
          itime,nlock,dir,longueur1)
!
!  Calculation of forces on the external CV surface
   call sforce(ux,uy,pp,f2x,f2y,nx,ny,nz,pp_ref,icvrt,jcvlw,jcvup,&
          ffxp,fsxp,fwxp,ffyp,fsyp,fwyp,itime,nlock,dir,longueur1)
!
!  Calculation of the aerodynamic force components
   yLift = f2y-ymom
   xDrag = f2x-xmom
   yLift = 2.*yLift
   xDrag = 2.*xDrag
!       
   return
end subroutine aerof_inflow
!
!***********************************************************************
!
subroutine rmomentum(ux,uy,xmom,ymom,nx,ny,nz,ux_ref,uy_ref,icvrt,jcvlw,jcvup,&
             itime,nlock,dir,longueur1)
!             
!***********************************************************************
!
   USE paramx_m
   USE paramy_m
   USE paramz_m
   USE param1_m
   USE barreau_m 
!
   implicit none
!
   real(8),dimension(nx,ny,nz) :: ux,uy
   integer :: i,j,k,nx,ny,nz,itime,nlock
   real(8) :: xmom,ymom
   real(8) :: uxm,uym,ux_ref,uy_ref
   real(8) :: fxab,fyab,pxab,pyab
   real(8) :: fxbc,fybc,pxbc,pybc
   real(8) :: fxcd,fycd,pxcd,pycd
   real(8) :: fxda,fyda,pxda,pyda
   integer :: icvrt,jcvlw,jcvup
   integer :: longueur1
   character(len=120) dir
!
!  The surface momentum fluxes
!
   fxab=0.0
   fyab=0.0
   do j=jcvlw,jcvup-1
      uxm = ux_ref
      uym = uy_ref
      pxab=-uxm*uxm*dy
      fxab= fxab+pxab
      pyab=-uxm*uym*dy
      fyab= fyab+pyab
   enddo
!
   fxbc=0.0
   fybc=0.0
!
   fxcd=0.0
   fycd=0.0
   i=icvrt
   do j=jcvlw,jcvup-1
      uxm = (ux(i,j,1)+ux(i,j+1,1))/2.0
      uym = (uy(i,j,1)+uy(i,j+1,1))/2.0
      pxcd=+uxm*uxm*dy
      fxcd= fxcd+pxcd
      pycd=+uxm*uym*dy
      fycd= fycd+pycd
   enddo
!
   fxda=0.0
   fyda=0.0
!
!  Components of the total momentum 
   xmom=fxab+fxbc+fxcd+fxda
   ymom=fyab+fybc+fycd+fyda
!
   if (nlock.eq.1) then
     open(80,file=dir(1:longueur1)//'rmomentum_ini',access='append')
     write(80,1100) itime,fxab,fxcd,fyab,fycd
     close(80)
   endif
   if (nlock.eq.2) then
     open(90,file=dir(1:longueur1)//'rmomentum_fin',access='append')
     write(90,1100) itime,fxab,fxcd,fyab,fycd
     close(90)
   endif
 1100 format(I5,4F16.10) 
!
   return
end subroutine rmomentum
!
!***********************************************************************
!
subroutine sforce(ux,uy,pp,f2x,f2y,nx,ny,nz,pp_ref,icvrt,jcvlw,jcvup,&
             ffxp,fsxp,fwxp,ffyp,fsyp,fwyp,itime,nlock,dir,longueur1)
!             
!***********************************************************************
!
   USE paramx_m
   USE paramy_m
   USE paramz_m
   USE param1_m
!
   implicit none
!
   real(8),dimension(nx,ny,nz) :: ux,uy,gx,gy,tx,ty
   real(8),dimension(nx,ny,nz) :: pp
   real(8),dimension(nx) :: ffxp,fsxp,fwxp
   real(8),dimension(ny) :: ffyp,fsyp,fwyp
   integer :: i,j,k,nx,ny,nz,itime,nlock  
   real(8) :: f2x,f2y 
   real(8) :: dudxm,dudym,dvdxm,dvdym
   real(8) :: foxab,foyab,pxab,pyab
   real(8) :: foxbc,foybc,pxbc,pybc
   real(8) :: foxcd,foycd,pxcd,pycd
   real(8) :: foxda,foyda,pxda,pyda
   real(8) :: rho,foxpr,foypr,pab,pbc,pcd,pad,difp,pxp,pyp
   real(8) :: pp_ref
   integer :: icvrt,jcvlw,jcvup
   integer :: longueur1
   character(len=120) dir
!
!  Force calculation
!
   call derx (tx,ux,ffxp,fsxp,fwxp,nx,ny,nz,1)
   call dery (ty,ux,ffyp,fsyp,fwyp,nx,ny,nz,1)
   call derx (gx,uy,ffxp,fsxp,fwxp,nx,ny,nz,1)
   call dery (gy,uy,ffyp,fsyp,fwyp,nx,ny,nz,1)
!
!  Along CV's entrance face AB
   foxab=0.0
   foyab=0.0
!
!  Along CV's upper face BC
   foxbc=0.0
   foybc=0.0
!
!  Along CV's exit face CD
   i=icvrt
   foxcd=0.0
   foycd=0.0
   do j=jcvlw,jcvup-1
      dudxm = (tx(i,j,1)+tx(i,j+1,1))/2.0
      dudym = (ty(i,j,1)+ty(i,j+1,1))/2.0
      dvdxm = (gx(i,j,1)+gx(i,j+1,1))/2.0
      pxcd  =+2.0*xnu*dudxm*dy
      foxcd = foxcd+pxcd
      pycd  =+xnu*(dvdxm+dudym)*dy
      foycd = foycd+pycd
   enddo
!
!  Along CV's lower face DA
   foxda=0.0
   foyda=0.0
!
!  The pressure contribution 
!
!  The pressure force between planes AB and CD
   rho=1.0
   foxpr=0.0
   do j=jcvlw,jcvup-1
      pab   = pp_ref
      i     = icvrt
      pcd   = (pp(i,j,1)+pp(i,j+1,1))/2.0
      difp  = (pab-pcd)/rho
      pxp   = difp*dy
      foxpr = foxpr+pxp
   enddo
!
!  The pressure force between planes AD and BC
   foypr=0.0
!
!  Components of the total force 
   f2x=foxab+foxbc+foxcd+foxda+foxpr
   f2y=foyab+foybc+foycd+foyda+foypr
!
   if (nlock.eq.1) then
     open(100,file=dir(1:longueur1)//'sforce_ini',access='append')
     write(100,1200) itime,foxcd,foycd,foxpr
     close(100)
   endif
   if (nlock.eq.2) then
     open(110,file=dir(1:longueur1)//'sforce_fin',access='append')
     write(110,1200) itime,foxcd,foycd,foxpr
     close(110)
   endif
 1200 format(I5,3F16.10) 
!
   return
end subroutine sforce




