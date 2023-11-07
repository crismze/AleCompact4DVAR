
# Set correct variables and paths
EXEC = /home/ale/Asim_86/iasim

###compiler
FC = gfortran-4.4
OPTFC = -O3 -fdefault-integer-8  -fdefault-real-8


$(EXEC): module_param.o asimila2d.o mesure.o pression.o \
	spectral.o outils.o schemas.o derive.o \
	filtre.o parametre.o controle.o \
	front_virt.o navier.o pression6.o module_avs.o \
	jmccfft2d.o jmccfft3d.o jmccfft.o \
	jmccfftm.o jmccm1d2.o jmccm1d3.o \
	jmccm1d4.o jmccm1d5.o jmccm1d.o jmccm1dp.o jmcctranspcs.o \
	jmcfftmlt.o jmcftfax.o jmcsfft2d.o jmcsfft3d.o jmcsfft.o \
	jmcsfftm.o jmcsm1d.o jmcsm1dxy.o jmdecoup3.o jmdecoup.o \
	jmerreur1.o jmerreur2.o jmerreur3.o jmerreur4.o jmfact.o \
	jmfftfax.o jmgetcode.o jmgeterreur.o jmgeterror.o jmgetmessage.o \
	jmgetnwork.o jmgetsetcode.o jmgetseterreur.o jmgetseterror.o jmgetsetnwork.o \
	jmgetsetstop.o jmgetstop.o jmrfftmlt.o jmscfft2d.o jmscfft3d.o \
	jmscfft.o jmscfftm.o jmscm1d.o jmscm1dxy.o jmsetcode.o jmseterreur.o \
	jmsetnwork.o jmsetstop.o jmtable.o jmtransp.o c06hbf.o \
	cosfftmlt.o sinfftmlt.o jmtableCos.o jmtableSin.o slfft2d.o \
	funcion1_d.o derx_d.o dery_d.o derxx_d.o deryy_d.o \
	deci6_d.o deciy6_d.o deciy16_d.o interi6_d.o interiy6_d.o interiy16_d.o \
	decx6_d.o decy6_d.o decy16_d.o inter6_d.o intery6_d.o intery16_d.o \
	funcion3_d.o funcion4_d.o banbks1_d.o funcion5_d.o funcion7_d.o snorm_d.o reconstruye4_d.o vortz_d.o \
	funcion1_b.o derx_b.o dery_b.o derxx_b.o deryy_b.o \
	deci6_b.o deciy6_b.o deciy16_b.o interi6_b.o interiy6_b.o interiy16_b.o \
	decx6_b.o decy6_b.o decy16_b.o inter6_b.o intery6_b.o intery16_b.o \
	funcion3_b.o funcion4_b.o banbks1_b.o funcion5_b.o funcion7_b.o snorm_b.o reconstruye4_b.o vortz_b.o \
	derx_r_b.o dery_r_b.o vortz_r_b.o \
	lbfgs.o adBuffer.o adStack.o
 

	$(FC) $(OPTFC) -o $(EXEC) \
	module_param.o asimila2d.o mesure.o pression.o \
	spectral.o outils.o schemas.o derive.o \
	filtre.o parametre.o controle.o \
	front_virt.o navier.o pression6.o module_avs.o \
	jmccfft2d.o jmccfft3d.o jmccfft.o \
	jmccfftm.o jmccm1d2.o jmccm1d3.o \
	jmccm1d4.o jmccm1d5.o jmccm1d.o jmccm1dp.o jmcctranspcs.o \
	jmcfftmlt.o jmcftfax.o jmcsfft2d.o jmcsfft3d.o jmcsfft.o \
	jmcsfftm.o jmcsm1d.o jmcsm1dxy.o jmdecoup3.o jmdecoup.o \
	jmerreur1.o jmerreur2.o jmerreur3.o jmerreur4.o jmfact.o \
	jmfftfax.o jmgetcode.o jmgeterreur.o jmgeterror.o jmgetmessage.o \
	jmgetnwork.o jmgetsetcode.o jmgetseterreur.o jmgetseterror.o jmgetsetnwork.o \
	jmgetsetstop.o jmgetstop.o jmrfftmlt.o jmscfft2d.o jmscfft3d.o \
	jmscfft.o jmscfftm.o jmscm1d.o jmscm1dxy.o jmsetcode.o jmseterreur.o \
	jmsetnwork.o jmsetstop.o jmtable.o jmtransp.o c06hbf.o \
	cosfftmlt.o sinfftmlt.o jmtableCos.o jmtableSin.o slfft2d.o \
	funcion1_d.o derx_d.o dery_d.o derxx_d.o deryy_d.o \
	deci6_d.o deciy6_d.o deciy16_d.o interi6_d.o interiy6_d.o interiy16_d.o \
	decx6_d.o decy6_d.o decy16_d.o inter6_d.o intery6_d.o intery16_d.o \
	funcion3_d.o funcion4_d.o banbks1_d.o funcion5_d.o funcion7_d.o snorm_d.o reconstruye4_d.o vortz_d.o \
	funcion1_b.o derx_b.o dery_b.o derxx_b.o deryy_b.o \
	deci6_b.o deciy6_b.o deciy16_b.o interi6_b.o interiy6_b.o interiy16_b.o \
	decx6_b.o decy6_b.o decy16_b.o inter6_b.o intery6_b.o intery16_b.o \
	funcion3_b.o funcion4_b.o banbks1_b.o funcion5_b.o funcion7_b.o snorm_b.o reconstruye4_b.o vortz_b.o \
	derx_r_b.o dery_r_b.o vortz_r_b.o \
	lbfgs.o adBuffer.o \
	-lgfortran adStack.o


asimila2d.o : asimila2d.f90
	$(FC) $(OPTFC) -c asimila2d.f90

mesure.o : mesure.f90
	$(FC) $(OPTFC) -c mesure.f90

pression6.o : pression6.f90
	$(FC) $(OPTFC) -c pression6.f90

pression.o : pression.f90
	$(FC) $(OPTFC) -c pression.f90

spectral.o : spectral.f90
	$(FC) $(OPTFC) -c spectral.f90

schemas.o : schemas.f90
	$(FC) $(OPTFC) -c schemas.f90

derive.o : derive.f90
	$(FC) $(OPTFC) -c derive.f90

outils.o : outils.f90
	$(FC) $(OPTFC) -c outils.f90

navier.o : navier.f90
	$(FC) $(OPTFC) -c navier.f90

filtre.o : filtre.f90
	$(FC) $(OPTFC) -c filtre.f90

parametre.o : parametre.f90
	$(FC) $(OPTFC) -c parametre.f90

controle.o : controle.f90
	$(FC) $(OPTFC) -c controle.f90

front_virt.o : front_virt.f90
	$(FC) $(OPTFC) -c front_virt.f90

module_avs.o : module_avs.f90
	$(FC) $(OPTFC) -c module_avs.f90

jmccfft2d.o : jmccfft2d.f90
	$(FC) $(OPTFC) -c jmccfft2d.f90

jmccfft3d.o : jmccfft3d.f90
	$(FC) $(OPTFC) -c jmccfft3d.f90

jmccfft.o : jmccfft.f90
	$(FC) $(OPTFC) -c jmccfft.f90

jmccfftm.o : jmccfftm.f90
	$(FC) $(OPTFC) -c jmccfftm.f90

jmccm1d2.o : jmccm1d2.f90
	$(FC) $(OPTFC) -c jmccm1d2.f90

jmccm1d3.o : jmccm1d3.f90
	$(FC) $(OPTFC) -c jmccm1d3.f90

jmccm1d4.o : jmccm1d4.f90
	$(FC) $(OPTFC) -c jmccm1d4.f90

jmccm1d5.o : jmccm1d5.f90
	$(FC) $(OPTFC) -c jmccm1d5.f90

jmccm1d.o : jmccm1d.f90
	$(FC) $(OPTFC) -c jmccm1d.f90

jmccm1dp.o : jmccm1dp.f90
	$(FC) $(OPTFC) -c jmccm1dp.f90

jmcctranspcs.o : jmcctranspcs.f90
	$(FC) $(OPTFC) -c jmcctranspcs.f90

jmcfftmlt.o : jmcfftmlt.f90
	$(FC) $(OPTFC) -c jmcfftmlt.f90

jmcftfax.o : jmcftfax.f90
	$(FC) $(OPTFC) -c jmcftfax.f90

jmcsfft2d.o : jmcsfft2d.f90
	$(FC) $(OPTFC) -c jmcsfft2d.f90

jmcsfft3d.o : jmcsfft3d.f90
	$(FC) $(OPTFC) -c jmcsfft3d.f90

jmcsfft.o : jmcsfft.f90
	$(FC) $(OPTFC) -c jmcsfft.f90

jmcsfftm.o : jmcsfftm.f90
	$(FC) $(OPTFC) -c jmcsfftm.f90

jmcsm1d.o : jmcsm1d.f90
	$(FC) $(OPTFC) -c jmcsm1d.f90

jmcsm1dxy.o : jmcsm1dxy.f90
	$(FC) $(OPTFC) -c jmcsm1dxy.f90

jmdecoup3.o : jmdecoup3.f90
	$(FC) $(OPTFC) -c jmdecoup3.f90

jmdecoup.o : jmdecoup.f90
	$(FC) $(OPTFC) -c jmdecoup.f90

jmerreur1.o : jmerreur1.f90
	$(FC) $(OPTFC) -c jmerreur1.f90

jmerreur2.o : jmerreur2.f90
	$(FC) $(OPTFC) -c jmerreur2.f90

jmerreur3.o : jmerreur3.f90
	$(FC) $(OPTFC) -c jmerreur3.f90

jmerreur4.o : jmerreur4.f90
	$(FC) $(OPTFC) -c jmerreur4.f90

jmfact.o : jmfact.f90
	$(FC) $(OPTFC) -c jmfact.f90

jmfftfax.o : jmfftfax.f90
	$(FC) $(OPTFC) -c jmfftfax.f90

jmgetcode.o : jmgetcode.f90
	$(FC) $(OPTFC) -c jmgetcode.f90

jmgeterreur.o : jmgeterreur.f90
	$(FC) $(OPTFC) -c jmgeterreur.f90

jmgeterror.o : jmgeterror.f90
	$(FC) $(OPTFC) -c jmgeterror.f90

jmgetmessage.o : jmgetmessage.f90
	$(FC) $(OPTFC) -c jmgetmessage.f90

jmgetnwork.o : jmgetnwork.f90
	$(FC) $(OPTFC) -c jmgetnwork.f90

jmgetsetcode.o : jmgetsetcode.f90
	$(FC) $(OPTFC) -c jmgetsetcode.f90

jmgetseterreur.o : jmgetseterreur.f90
	$(FC) $(OPTFC) -c jmgetseterreur.f90

jmgetseterror.o : jmgetseterror.f90
	$(FC) $(OPTFC) -c jmgetseterror.f90

jmgetsetmessage.o : jmgetsetmessage.f90
	$(FC) $(OPTFC) -c jmgetsetmessage.f90

jmgetsetnwork.o : jmgetsetnwork.f90
	$(FC) $(OPTFC) -c jmgetsetnwork.f90

jmgetsetstop.o : jmgetsetstop.f90
	$(FC) $(OPTFC) -c jmgetsetstop.f90

jmgetstop.o : jmgetstop.f90
	$(FC) $(OPTFC) -c jmgetstop.f90

jmrfftmlt.o : jmrfftmlt.f90
	$(FC) $(OPTFC) -c jmrfftmlt.f90

jmscfft2d.o : jmscfft2d.f90
	$(FC) $(OPTFC) -c jmscfft2d.f90

jmscfft3d.o : jmscfft3d.f90
	$(FC) $(OPTFC) -c jmscfft3d.f90

jmscfft.o : jmscfft.f90
	$(FC) $(OPTFC) -c jmscfft.f90

jmscfftm.o : jmscfftm.f90
	$(FC) $(OPTFC) -c jmscfftm.f90

jmscm1d.o : jmscm1d.f90
	$(FC) $(OPTFC) -c jmscm1d.f90

jmscm1dxy.o : jmscm1dxy.f90
	$(FC) $(OPTFC) -c jmscm1dxy.f90

jmsetcode.o : jmsetcode.f90
	$(FC) $(OPTFC) -c jmsetcode.f90

jmseterreur.o : jmseterreur.f90
	$(FC) $(OPTFC) -c jmseterreur.f90

jmsetnwork.o : jmsetnwork.f90
	$(FC) $(OPTFC) -c jmsetnwork.f90

jmsetstop.o : jmsetstop.f90
	$(FC) $(OPTFC) -c jmsetstop.f90

jmtable.o : jmtable.f90
	$(FC) $(OPTFC) -c jmtable.f90

jmtransp.o : jmtransp.f90
	$(FC) $(OPTFC) -c jmtransp.f90

c06hbf.o : c06hbf.f90
	$(FC) $(OPTFC) -c c06hbf.f90

cosfftmlt.o : cosfftmlt.f90
	$(FC) $(OPTFC) -c cosfftmlt.f90

sinfftmlt.o : sinfftmlt.f90
	$(FC) $(OPTFC) -c sinfftmlt.f90

jmtableCos.o : jmtableCos.f90
	$(FC) $(OPTFC) -c jmtableCos.f90

jmtableSin.o : jmtableSin.f90
	$(FC) $(OPTFC) -c jmtableSin.f90

slfft2d.o : slfft2d.f90
	$(FC) $(OPTFC) -c slfft2d.f90

module_param.o : module_param.f90
	$(FC) $(OPTFC) -c module_param.f90

funcion1_d.o : funcion1_d.f90
	$(FC) $(OPTFC) -c funcion1_d.f90

derx_d.o : derx_d.f90
	$(FC) $(OPTFC) -c derx_d.f90

dery_d.o : dery_d.f90
	$(FC) $(OPTFC) -c dery_d.f90

derxx_d.o : derxx_d.f90
	$(FC) $(OPTFC) -c derxx_d.f90

deryy_d.o : deryy_d.f90
	$(FC) $(OPTFC) -c deryy_d.f90

deci6_d.o : deci6_d.f90
	$(FC) $(OPTFC) -c deci6_d.f90

deciy6_d.o : deciy6_d.f90
	$(FC) $(OPTFC) -c deciy6_d.f90

deciy16_d.o : deciy16_d.f90
	$(FC) $(OPTFC) -c deciy16_d.f90

interi6_d.o : interi6_d.f90
	$(FC) $(OPTFC) -c interi6_d.f90

interiy6_d.o : interiy6_d.f90
	$(FC) $(OPTFC) -c interiy6_d.f90

interiy16_d.o : interiy16_d.f90
	$(FC) $(OPTFC) -c interiy16_d.f90

decx6_d.o : decx6_d.f90
	$(FC) $(OPTFC) -c decx6_d.f90

decy6_d.o : decy6_d.f90
	$(FC) $(OPTFC) -c decy6_d.f90

decy16_d.o : decy16_d.f90
	$(FC) $(OPTFC) -c decy16_d.f90

inter6_d.o : inter6_d.f90
	$(FC) $(OPTFC) -c inter6_d.f90

intery6_d.o : intery6_d.f90
	$(FC) $(OPTFC) -c intery6_d.f90

intery16_d.o : intery16_d.f90
	$(FC) $(OPTFC) -c intery16_d.f90

funcion3_d.o : funcion3_d.f90
	$(FC) $(OPTFC) -c funcion3_d.f90

funcion4_d.o : funcion4_d.f90
	$(FC) $(OPTFC) -c funcion4_d.f90

banbks1_d.o : banbks1_d.f90
	$(FC) $(OPTFC) -c banbks1_d.f90

funcion5_d.o : funcion5_d.f90
	$(FC) $(OPTFC) -c funcion5_d.f90

funcion7_d.o : funcion7_d.f90
	$(FC) $(OPTFC) -c funcion7_d.f90

snorm_d.o : snorm_d.f90
	$(FC) $(OPTFC) -c snorm_d.f90

reconstruye4_d.o : reconstruye4_d.f90
	$(FC) $(OPTFC) -c reconstruye4_d.f90

vortz_d.o : vortz_d.f90
	$(FC) $(OPTFC) -c vortz_d.f90

funcion1_b.o : funcion1_b.f90
	$(FC) $(OPTFC) -c funcion1_b.f90

derx_b.o : derx_b.f90
	$(FC) $(OPTFC) -c derx_b.f90

dery_b.o : dery_b.f90
	$(FC) $(OPTFC) -c dery_b.f90

derxx_b.o : derxx_b.f90
	$(FC) $(OPTFC) -c derxx_b.f90

deryy_b.o : deryy_b.f90
	$(FC) $(OPTFC) -c deryy_b.f90

deci6_b.o : deci6_b.f90
	$(FC) $(OPTFC) -c deci6_b.f90

deciy6_b.o : deciy6_b.f90
	$(FC) $(OPTFC) -c deciy6_b.f90

deciy16_b.o : deciy16_b.f90
	$(FC) $(OPTFC) -c deciy16_b.f90

interi6_b.o : interi6_b.f90
	$(FC) $(OPTFC) -c interi6_b.f90

interiy6_b.o : interiy6_b.f90
	$(FC) $(OPTFC) -c interiy6_b.f90

interiy16_b.o : interiy16_b.f90
	$(FC) $(OPTFC) -c interiy16_b.f90

decx6_b.o : decx6_b.f90
	$(FC) $(OPTFC) -c decx6_b.f90

decy6_b.o : decy6_b.f90
	$(FC) $(OPTFC) -c decy6_b.f90

decy16_b.o : decy16_b.f90
	$(FC) $(OPTFC) -c decy16_b.f90

inter6_b.o : inter6_b.f90
	$(FC) $(OPTFC) -c inter6_b.f90

intery6_b.o : intery6_b.f90
	$(FC) $(OPTFC) -c intery6_b.f90

intery16_b.o : intery16_b.f90
	$(FC) $(OPTFC) -c intery16_b.f90

funcion3_b.o : funcion3_b.f90
	$(FC) $(OPTFC) -c funcion3_b.f90

funcion4_b.o : funcion4_b.f90
	$(FC) $(OPTFC) -c funcion4_b.f90

banbks1_b.o : banbks1_b.f90
	$(FC) $(OPTFC) -c banbks1_b.f90

funcion5_b.o : funcion5_b.f90
	$(FC) $(OPTFC) -c funcion5_b.f90

funcion7_b.o : funcion7_b.f90
	$(FC) $(OPTFC) -c funcion7_b.f90

snorm_b.o : snorm_b.f90
	$(FC) $(OPTFC) -c snorm_b.f90

reconstruye4_b.o : reconstruye4_b.f90
	$(FC) $(OPTFC) -c reconstruye4_b.f90

vortz_b.o : vortz_b.f90
	$(FC) $(OPTFC) -c vortz_b.f90

derx_r_b.o : derx_r_b.f90
	$(FC) $(OPTFC) -c derx_r_b.f90

dery_r_b.o : dery_r_b.f90
	$(FC) $(OPTFC) -c dery_r_b.f90

vortz_r_b.o : vortz_r_b.f90
	$(FC) $(OPTFC) -c vortz_r_b.f90

lbfgs.o : lbfgs.f
	$(FC) $(OPTFC) -c lbfgs.f

adBuffer.o : adBuffer.f
	$(FC) $(OPTFC) -c adBuffer.f

adStack.o : adStack.c
	cc -c adStack.c


################### TAPENADE ###################

### to differentiate in the tangent mode
#TANGENT = -tangent

### to differentiate in the reverse mode
REVERSE = -reverse

### to specify the top differentiation routine
HEAD = -head 

### for the dependent output variables, as a white-space separated list enclosed by double quotes
OUTVARS = -outvars

### for the independent input variables, as a white-space separated list enclosed by double quotes
VARS = -vars

### files containing routines that can be recursively called from the top differentiation routine
FILES1 = navier.f90 derive.f90 
FILES3 = navier.f90 
FILES4 = navier.f90 pression.f90
FILES8 = outils.f90

### to create an additional HTML output, that can be displayed in your browser
#HTML = -html

tapenade1:
	tapenade $(TANGENT) $(REVERSE) $(HEAD) funcion1 $(OUTVARS) "u8x u8y g8x g8y ps2d1" $(VARS) "ux uy gx gy" $(FILES1) $(HTML)
tapenade3:
	tapenade $(TANGENT) $(REVERSE) $(HEAD) funcion3 $(OUTVARS) "ps2d3" $(VARS) "ps2d2" $(FILES3) $(HTML)
tapenade4:
	tapenade $(TANGENT) $(REVERSE) $(HEAD) funcion4 $(OUTVARS) "us2d" $(VARS) "ps2d" $(FILES4) $(HTML)
tapenade5:
	tapenade $(TANGENT) $(REVERSE) $(HEAD) funcion5 $(OUTVARS) "us2d2" $(VARS) "us2d1" $(FILES3) $(HTML)
tapenade7:
	tapenade $(TANGENT) $(REVERSE) $(HEAD) funcion7 $(OUTVARS) "uxn uyn" $(VARS) "us2d3 u8x u8y" $(FILES1) $(HTML)
tapenade8:
	tapenade $(TANGENT) $(REVERSE) $(HEAD) prueba1 $(OUTVARS) "u" $(VARS) "u0" $(FILES8) $(HTML)
tapenade9:
	tapenade $(TANGENT) $(REVERSE) $(HEAD) snorm $(OUTVARS) "error" $(VARS) "uxtest uytest" $(FILES8) $(HTML)
tapenade10:
	tapenade $(TANGENT) $(REVERSE) $(HEAD) reconstruye4 $(OUTVARS) "ux0 uy0 bx1 by1" $(VARS) "ux0r uy0r bx1r by1r" $(FILES8) $(HTML)

#################################################


clean:
	rm -f *.o $(EXEC)
#	rm -f *_d.*
#	rm -f *_b.*
#	rm -rf tapenadehtml

