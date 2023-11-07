      module aleatoire_m 
!...Created by PSUITE Trans90                  4.3ZH 10:37:46   1/28/ 4  
      integer :: idum=-67 
      end module aleatoire_m 
!
      module barreau_m 
!...Created by PSUITE Trans90                  4.3ZH 11:19:37   1/22/ 4  
      real(8) :: cex, cey, ra, alphav, betav, xkk 
      end module barreau_m 
!
      module controle_m 
!...Created by PSUITE Trans90                  4.3ZH 11:19:37   1/22/ 4  
      integer :: icontrole 
      real(8) :: rlat, distance, ufreq, uampl, phi, theta1, theta2 
      end module controle_m 
!
      module delta_m 
!...Created by PSUITE Trans90                  4.3ZH 11:19:37   1/22/ 4  
      real(8) :: deltx, delty, deltz, y1, z1 
      end module delta_m
!
      module dericex6_m 
!...Created by PSUITE Trans90                  4.3ZH 11:19:37   1/22/ 4  
      real(8) :: alcaix6, acix6, bcix6 
      end module dericex6_m 
!
      module dericex_m 
!...Created by PSUITE Trans90                  4.3ZH 11:19:37   1/22/ 4  
      real(8) :: alcaix, acix, bcix 
      end module dericex_m 
!
      module dericey6_m 
!...Created by PSUITE Trans90                  4.3ZH 11:19:37   1/22/ 4  
      real(8) :: alcaiy6, aciy6, bciy6 
      end module dericey6_m
!
      module dericey_m 
!...Created by PSUITE Trans90                  4.3ZH 11:19:37   1/22/ 4  
      real(8) :: alcaiy, aciy, bciy 
      end module dericey_m
!
      module derivex_m 
!...Created by PSUITE Trans90                  4.3ZH 11:19:37   1/22/ 4  
      real(8) :: alfa1x, af1x, bf1x, cf1x, df1x, alfa2x, af2x, alfanx, afnx, bfnx&
         , cfnx, dfnx, alfamx, afmx, alfaix, afix, bfix, alsa1x, as1x, bs1x, &
         cs1x, ds1x, alsa2x, as2x, alsanx, asnx, bsnx, csnx, dsnx, alsamx, asmx&
         , alsaix, asix, bsix 
      real(8) :: alfa1x_rg, af1x_rg, bf1x_rg, cf1x_rg, df1x_rg, alfa2x_rg, af2x_rg, alfanx_rg, afnx_rg, bfnx_rg&
         , cfnx_rg, dfnx_rg, alfamx_rg, afmx_rg, alfaix_rg, afix_rg, bfix_rg 

      end module derivex_m
!
      module derivey_m 
!...Created by PSUITE Trans90                  4.3ZH 11:19:37   1/22/ 4  
      real(8) :: alfa1y, af1y, bf1y, cf1y, df1y, alfa2y, af2y, alfany, afny, bfny&
         , cfny, dfny, alfamy, afmy, alfajy, afjy, bfjy, alsa1y, as1y, bs1y, &
         cs1y, ds1y, alsa2y, as2y, alsany, asny, bsny, csny, dsny, alsamy, asmy&
         , alsajy, asjy, bsjy 
      real(8) :: alfa1y_rg, af1y_rg, bf1y_rg, cf1y_rg, df1y_rg, alfa2y_rg, af2y_rg, alfany_rg, afny_rg, bfny_rg&
         , cfny_rg, dfny_rg, alfamy_rg, afmy_rg, alfajy_rg, afjy_rg, bfjy_rg 
      end module derivey_m 
!
      module derivez_m 
!...Created by PSUITE Trans90                  4.3ZH 11:19:37   1/22/ 4  
      real(8) :: alfa1z, af1z, bf1z, cf1z, df1z, alfa2z, af2z, alfanz, afnz, bfnz&
         , cfnz, dfnz, alfamz, afmz, alfakz, afkz, bfkz, alsa1z, as1z, bs1z, &
         cs1z, ds1z, alsa2z, as2z, alsanz, asnz, bsnz, csnz, dsnz, alsamz, asmz&
         , alsakz, askz, bskz 
      end module derivez_m 
!
      module ecoulement_m 
!...Created by PSUITE Trans90                  4.3ZH 11:19:37   1/22/ 4  
      integer :: iecoule, iscalaire, ientree, nschema, idebut, ifin, iles, &
         ifiltre, ivirtuel 
      real(8) :: u1, u2, rayon, re 
      end module ecoulement_m 
!
      module interpol6_m 
!...Created by PSUITE Trans90                  4.3ZH 11:19:37   1/22/ 4  
      real(8) :: ailcaix6, aicix6, bicix6 
      end module interpol6_m 
!
      module interpoly6_m 
!...Created by PSUITE Trans90                  4.3ZH 11:19:37   1/22/ 4  
      real(8) :: ailcaiy6, aiciy6, biciy6 
      end module interpoly6_m
!
      module limit_m 
!...Created by PSUITE Trans90                  4.3ZH 11:19:37   1/22/ 4  
      integer :: nbound 
      end module limit_m 
!
      module module_avs_m 
!...Created by PSUITE Trans90                  4.3ZH 11:19:37   1/22/ 4  
      integer :: iavs, ndsize, idebavs, imodavs, nx_min, nx_max, ny_min, ny_max&
         , nz_min, nz_max 
      real(8) :: val_iso 
      character :: path_network*80, nom_network*80, path_files*80,&
           nom_file*80, nom_script*80, nom_film*80 
      end module module_avs_m 
!
      module ondelette_m 
!...Created by PSUITE Trans90                  4.3ZH 11:19:37   1/22/ 4  
      real(8) :: wavl, was, wacpsi, wa8, wa12, wadx 
      end module ondelette_m
!
      module param1_m 
!...Created by PSUITE Trans90                  4.3ZH 11:19:37   1/22/ 4  
      real(8) :: dt, xnu, xkappa 
      end module param1_m 
!
      module param2_m 
!...Created by PSUITE Trans90                  4.3ZH 11:19:37   1/22/ 4  
      real(8) :: bruit, vdebit, bruitb 
      end module param2_m 
!
      module paramx_m 
!...Created by PSUITE Trans90                  4.3ZH 11:19:37   1/22/ 4  
      integer :: nclx 
      real(8) :: xlx, dx, dx2, dxg, xlxr
      end module paramx_m 
!
      module paramy_m 
!...Created by PSUITE Trans90                  4.3ZH 11:19:37   1/22/ 4  
      integer :: ncly 
      real(8) :: yly, dy, dy2, dyg, ylyr
      end module paramy_m 
!
      module paramz_m 
!...Created by PSUITE Trans90                  4.3ZH 11:19:37   1/22/ 4  
      integer :: nclz 
      real(8) :: zlz, dz, dz2 
      end module paramz_m 
!
      module parfix_m 
!...Created by PSUITE Trans90                  4.3ZH 11:19:37   1/22/ 4  
      real(8) :: fia1x, fib1x, fic1x, fid1x, fie1x, fia2x, fib2x, fic2x, fid2x, &
         fie2x, fia3x, fib3x, fic3x, fid3x, fie3x, fianx, fibnx, ficnx, fidnx, &
         fienx, fiamx, fibmx, ficmx, fidmx, fiemx, fiapx, fibpx, ficpx, fidpx, &
         fiepx, fiaix, fibix, ficix, fidix, fialx, fibex, fih1x, fih2x, fih3x, &
         fih4x 
      end module parfix_m 
!
      module parfiy_m 
!...Created by PSUITE Trans90                  4.3ZH 11:19:37   1/22/ 4  
      real(8) :: fia1y, fib1y, fic1y, fid1y, fie1y, fia2y, fib2y, fic2y, fid2y, &
         fie2y, fia3y, fib3y, fic3y, fid3y, fie3y, fiany, fibny, ficny, fidny, &
         fieny, fiamy, fibmy, ficmy, fidmy, fiemy, fiapy, fibpy, ficpy, fidpy, &
         fiepy, fiaiy, fibiy, ficiy, fidiy, fialy, fibey, fih1y, fih2y, fih3y, &
         fih4y 
      end module parfiy_m 
!
      module parfiz_m 
!...Created by PSUITE Trans90                  4.3ZH 11:19:37   1/22/ 4  
      real(8) :: fia1z, fib1z, fic1z, fid1z, fie1z, fia2z, fib2z, fic2z, fid2z, &
         fie2z, fia3z, fib3z, fic3z, fid3z, fie3z, fianz, fibnz, ficnz, fidnz, &
         fienz, fiamz, fibmz, ficmz, fidmz, fiemz, fiapz, fibpz, ficpz, fidpz, &
         fiepz, fiaiz, fibiz, ficiz, fidiz, fialz, fibez, fih1z, fih2z, fih3z, &
         fih4z 
      end module parfiz_m 
!
      module sauvegarde_m 
!...Created by PSUITE Trans90                  4.3ZH 11:19:37   1/22/ 4  
      integer :: isave, ilit, idebmod, imodulo, idemarre, icommence, irecord, &
         nxboite, istat, iobs, imodulob 
      character :: filecharge*80, filesauve*80, filebruit*80, nchamp*80, nchamp1*80, nchampb*80, &
         nchampbr*80 
      end module sauvegarde_m 
!
      module turbi_m 
!...Created by PSUITE Trans90                  4.3ZH 11:19:37   1/22/ 4  
      integer :: iturb 
      end module turbi_m 
!
     module interpol4_m 
!...Created by PSUITE Trans90                  4.3ZH 11:19:37   1/22/ 4  
      real(8) :: ailcaix, aicix, bicix 
      end module interpol4_m 
!
     module interpoly4_m 
!...Created by PSUITE Trans90                  4.3ZH 11:19:37   1/22/ 4  
      real(8) :: ailcaiy, aiciy, biciy 
      end module interpoly4_m 
!
      module value_m 
!...Created by PSUITE Trans90                  4.3ZH 14:53:24   1/27/ 4  
      real(8) :: v 
      end module value_m 
!
      module pwtcom_m 
!...Created by PSUITE Trans90                  4.3ZH 14:53:24   1/27/ 4  
      integer :: ncof, ioff, joff 
      real(8), dimension(20) :: cc, cr, c20
      end module pwtcom_m 
!
      module q
!...Created by PSUITE Trans90                  4.3ZH 14:53:24   1/27/ 4  
      integer :: iq,idebq,ifinq
      character(len=80) :: lect_q
      real(8) :: avs_q
      end module q
!
      module pdf
!...Created by PSUITE Trans90                  4.3ZH 14:53:24   1/27/ 4  
      integer :: ipdf,idebpdf,ifinpdf,npdf
      character(len=80) :: lect_pdf
      end module pdf
!
      module moyen
!...Created by PSUITE Trans90                  4.3ZH 14:53:24   1/27/ 4  
      integer :: imoyen,idebmoy,ifinmoy
      character(len=80) lect_moy
      end module moyen
!
      module expan
!...Created by PSUITE Trans90                  4.3ZH 14:53:24   1/27/ 4  
      integer :: iexp,idebexp,ifinexp
      character(len=80) :: lect_exp
      end module expan
!
      module avs
!...Created by PSUITE Trans90                  4.3ZH 14:53:24   1/27/ 4  
      integer :: iavs,ndsize,idebavs,imodavs
      character(len=80) :: path_network,nom_network,path_files
      character(len=80) :: nom_file,nom_script,nom_film
      real(8) :: nx_min,nx_max,ny_min,ny_max,nz_min,nz_max,val_iso
      end module avs
