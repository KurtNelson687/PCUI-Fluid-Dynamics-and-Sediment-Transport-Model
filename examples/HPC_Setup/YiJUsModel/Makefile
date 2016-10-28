INSTALL_DIR = /opt/openmpi
FLINKER = $(INSTALL_DIR)/bin/mpif77
%LIBS = -lm
COMP = mpif90
FLAG = -c -w
NPROC = 32

OBCA = ns.o mp.o io.o dpth.o\
       init.o conv.o pred.o corr.o trid.o scal.o\
       pres.o smad.o resd.o fctr.o exch.o eddy.o\
       sedi.o random.o update.o lu.o

INCS = size.inc mpi.inc para.inc depth.inc mygrid.inc\
       ns.inc metric.inc sedi.inc eddy.inc

sedi : $(OBCA) $(INCS) mygrid.inc mygrid.o 
	$(COMP) -o sedi $(LIBS) $(OBCA) $(FLIBS) mygrid.o

test : sedi
#	$(INSTALL_DIR)/bin/mpirun -np $(NPROC) ./sedi 
	qsub r620.sh

# delete files
clean :
	/bin/rm -f *.o *~ fort.* test.*

# main program
ns.o : ns.f $(INCS)
	$(COMP) $(FLAG) ns.f 

# message passing
mp.o : mp.f $(INCS)
	$(COMP) $(FLAG) mp.f 

# input / output
io.o : io.f $(INCS)
	$(COMP) $(FLAG) io.f 

# les
eddy.o : eddy.f $(INCS)
	$(COMP) $(FLAG) eddy.f 

# initialization
init.o : init.f $(INCS)
	$(COMP) $(FLAG) init.f
 
random.o : random.f $(INCS)
	$(COMP) $(FLAG) random.f

# convection
conv.o : conv.f $(INCS)
	$(COMP) $(FLAG) conv.f 

# prediction
pred.o : pred.f $(INCS)
	$(COMP) $(FLAG) pred.f 

# correction
corr.o : corr.f $(INCS)
	$(COMP) $(FLAG) corr.f 

# tridiagonal solver
trid.o : trid.f $(INCS)
	$(COMP) $(FLAG) trid.f 

# pressure initialization
pres.o : pres.f $(INCS)
	$(COMP) $(FLAG) pres.f 

# ADI smoother	        
smad.o : smad.f $(INCS)
	$(COMP) $(FLAG) smad.f

# residue
resd.o : resd.f $(INCS)
	$(COMP) $(FLAG) resd.f 

# fine / coarse transfer
fctr.o : fctr.f $(INCS)
	$(COMP) $(FLAG) fctr.f  

# message exchange
exch.o : exch.f $(INCS)
	$(COMP) $(FLAG) exch.f 

#sedi
mygrid.o : mygrid.f ${INCS}
	${COMP} ${FLAG} mygrid.f

# scalar 
scal.o : scal.f $(INCS)
	$(COMP) $(FLAG) scal.f 

# sediment
sedi.o : sedi.f $(INCS)
	$(COMP) $(FLAG) sedi.f

# grid update
update.o : update.f $(INCS)
	$(COMP) $(FLAG) update.f

dpth.o : dpth.f $(INCS)
	$(COMP) $(FLAG) dpth.f

# LU-solver
lu.o : lu.f
	$(COMP) $(FLAG) lu.f

 
