COMP = mpif90
FLAG = -c -w -O3 

OBCA = ns.o mp.o io.o \
       init.o conv.o pred.o corr.o trid.o scal.o\
       pres.o smad.o resd.o fctr.o exch.o eddy.o\
			 indexx.o energy.o

OBCAS =nsst.o mp.o io.o stat.o\
       init.o conv.o pred.o corr.o trid.o scal.o\
       pres.o smad.o resd.o fctr.o exch.o eddy.o

OBJG = grid_init.o mp.o io.o cavity.o

INCS = size.inc mpi.inc para.inc \
       ns.inc metric.inc stat.inc eddy.inc

# 3-D cavity
cav : $(OBCA) $(INCS) cavity.inc cavity.o 	
	$(COMP) -o cavity $(OBCA) cavity.o	

# 3-D grid initialization
grid : $(OBJG) $(INCS) cavity.inc
	$(COMP) -o grid $(OBJG)

# delete files
clean :
	/bin/rm -f *.o
cleano:
	/bin/rm -f output* continue_run* uvw* xyz* rho* fort* qoutput grid cavity

# main program ( with statistics )
nsst.o : nsst.f $(INCS)
	$(COMP) $(FLAG) nsst.f 

# statistic
stat.o : stat.f $(INCS)
	 $(COMP) $(FLAG) stat.f 

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

# isw		
cavity.o : cavity.f $(INCS)
	$(COMP) $(FLAG) cavity.f 	

# grid_init		
grid_init.o : grid_init.f $(INCS)
	$(COMP) $(FLAG) grid_init.f 

# scalar 
scal.o : scal.f $(INCS)
	$(COMP) $(FLAG) scal.f 

# sort
indexx.o : indexx.f $(INCS)
	$(COMP) $(FLAG) indexx.f

# energy
energy.o : energy.f $(INCS)
	$(COMP) $(FLAG) energy.f
