ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine parameter

	include "size.inc"
        include "mpif.h"
        include "mpi.inc"
	include "para.inc"
        include "sedi.inc"
	include "eddy.inc"
	
	integer i

	call MPI_COMM_RANK( MPI_COMM_WORLD, MYID, ierr )

	if ( MYID .eq. 0 ) then

	case = 0
	newrun = 0
	periodic = 1
	periodicZ = 1
	iscalar = 0
	sedi = 1
	bedload = 0
	d_diff = 0
	d_grav = 1
        ieddy = 1
	mg_level = 6

        nstep = 20000
	nsave = 400
	ncheck = 1
        dtime = 0.0015
	maxstep = 120
	do i = 1, 6
	   iterchk(i) = 2
	   maxiter(i) = 20
	   tol(i) = 1.
	   ter(i) = 1.
	   slowiter(i) = 0.6
	enddo
	slowiter(1) = 0.7
	vis = 0.000001D0
	ak  = 0.000001D0
	g   = 9.81
	datum = 0.15
	omg_cyl = 0
	omg_lid = 0

	factor = 1.0e-3

	phi1 = 0
	phi2 = 0
	yphi = 0
	aphi = 1.

        amp = 0.5D0
        omega = 0.7852D0

        openWest = 1
        openEast = 1

	schd = 1.D0
        diam = 0.00022D0
        spwght = 2.65D0
	delta = 0.001
	rho = 1000.D0
	logU = 0
        btmAug = 0
	btmBuffer = 0
	movingGrid = 0
	open_top = 0
	if (sedi .eq. 1) then
c	   ws = -10.D0*vis/diam*(dsqrt(1.D0+(0.01*(spwght-1)*g*diam**3.D0
c     <          /vis**2.D0))-1)
	   ws = -0.0001D0
	endif
	h0 = 0.01
	hf = 0.05
	rp_length = 0.25
	const_b = 0.05
	grow_b = 0.005
	por = 0.4
	resp = 30
	fwc = dexp(5.213*(2.5*diam/amp)**0.194-5.977)
	cflsuper = 0.D0
	cut = 0

	endif

	call MPI_BCAST(case,        1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(newrun,      1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(periodic,    1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(periodicZ,   1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(iscalar,     1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(sedi,        1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(bedload,     1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(btmBuffer,   1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(ieddy,       1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(mg_level,    1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

	call MPI_BCAST(nstep,       1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(nsave,       1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(ncheck,       1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(d_diff,       1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(d_diff,      1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

	call MPI_BCAST(maxstep,     1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(iterchk(1),  6,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(maxiter(1),  6,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        call MPI_BCAST(dtime,       1,MPI_DOUBLE_PRECISION,0,
     <	                            MPI_COMM_WORLD,ierr)

	call MPI_BCAST(vis,         1,MPI_DOUBLE_PRECISION,0,
     <                              MPI_COMM_WORLD,ierr)
	call MPI_BCAST(ak,          1,MPI_DOUBLE_PRECISION,0,
     <                              MPI_COMM_WORLD,ierr)
	call MPI_BCAST(g,           1,MPI_DOUBLE_PRECISION,0,
     <                              MPI_COMM_WORLD,ierr)

	call MPI_BCAST(omg_cyl,     1,MPI_DOUBLE_PRECISION,0,
     <                              MPI_COMM_WORLD,ierr)
	call MPI_BCAST(omg_lid,     1,MPI_DOUBLE_PRECISION,0,
     <                              MPI_COMM_WORLD,ierr)

	call MPI_BCAST(factor,      1,MPI_DOUBLE_PRECISION,0,
     <                              MPI_COMM_WORLD,ierr)
	call MPI_BCAST(tol(1),      6,MPI_DOUBLE_PRECISION,0,
     <                              MPI_COMM_WORLD,ierr)
	call MPI_BCAST(ter(1),      6,MPI_DOUBLE_PRECISION,0,
     <                              MPI_COMM_WORLD,ierr)
	call MPI_BCAST(slowiter(1), 6,MPI_DOUBLE_PRECISION,0,
     <                              MPI_COMM_WORLD,ierr)

	call MPI_BCAST(phi1,        1,MPI_DOUBLE_PRECISION,0,
     <                              MPI_COMM_WORLD,ierr)
	call MPI_BCAST(phi2,        1,MPI_DOUBLE_PRECISION,0,
     <                              MPI_COMM_WORLD,ierr)
	call MPI_BCAST(yphi,        1,MPI_DOUBLE_PRECISION,0,
     <                              MPI_COMM_WORLD,ierr)
	call MPI_BCAST(aphi,        1,MPI_DOUBLE_PRECISION,0,
     <                              MPI_COMM_WORLD,ierr)
        call MPI_BCAST(amp,         1,MPI_DOUBLE_PRECISION,0,
     <	                            MPI_COMM_WORLD,ierr)
        call MPI_BCAST(omega,      1,MPI_DOUBLE_PRECISION,0,
     <                              MPI_COMM_WORLD,ierr)
        call MPI_BCAST(schd,        1,MPI_DOUBLE_PRECISION,0,
     <                              MPI_COMM_WORLD,ierr)
        call MPI_BCAST(diam,        1,MPI_DOUBLE_PRECISION,0,
     <                              MPI_COMM_WORLD,ierr)
        call MPI_BCAST(spwght,      1,MPI_DOUBLE_PRECISION,0,
     <                              MPI_COMM_WORLD,ierr)
	call MPI_BCAST(delta,      1,MPI_DOUBLE_PRECISION,0,
     <                              MPI_COMM_WORLD,ierr)
        call MPI_BCAST(rho,         1,MPI_DOUBLE_PRECISION,0,
     <                              MPI_COMM_WORLD,ierr)
	call MPI_BCAST(logU,        1,MPI_INTEGER,         0, 
     <	                            MPI_COMM_WORLD,ierr)
        call MPI_BCAST(btmAug,      1,MPI_INTEGER,         0, 
     <	                            MPI_COMM_WORLD,ierr)
        call MPI_BCAST(movingGrid,  1,MPI_INTEGER,         0,
     <                              MPI_COMM_WORLD,ierr)
	call MPI_BCAST(open_top,    1,MPI_INTEGER,         0,
     <	                            MPI_COMM_WORLD,ierr)
	call MPI_BCAST(ws,          1,MPI_DOUBLE_PRECISION,0,
     <                              MPI_COMM_WORLD,ierr)
	call MPI_BCAST(h0,          1,MPI_DOUBLE_PRECISION,0,
     <                              MPI_COMM_WORLD,ierr)
	call MPI_BCAST(hf,          1,MPI_DOUBLE_PRECISION,0,
     <                              MPI_COMM_WORLD,ierr)
	call MPI_BCAST(rp_length,   1,MPI_DOUBLE_PRECISION,0,
     <                              MPI_COMM_WORLD,ierr)
	call MPI_BCAST(const_b,     1,MPI_DOUBLE_PRECISION,0,
     <                              MPI_COMM_WORLD,ierr)
	call MPI_BCAST(grow_b,      1,MPI_DOUBLE_PRECISION,0,
     <                              MPI_COMM_WORLD,ierr)
	call MPI_BCAST(fwc,         1,MPI_DOUBLE_PRECISION,0,
     <                              MPI_COMM_WORLD,ierr)
	call MPI_BCAST(por,         1,MPI_DOUBLE_PRECISION,0,
     <                              MPI_COMM_WORLD,ierr)
	call MPI_BCAST(resp,        1,MPI_DOUBLE_PRECISION,0,
     <                              MPI_COMM_WORLD,ierr)
	call MPI_BCAST(cflsuper,    1,MPI_DOUBLE_PRECISION,0,
     <                              MPI_COMM_WORLD,ierr)
	call MPI_BCAST(datum,       1,MPI_DOUBLE_PRECISION,0,
     <                              MPI_COMM_WORLD,ierr)
	call MPI_BCAST(cut,         1,MPI_INTEGER,         0,
     <	                            MPI_COMM_WORLD,ierr)
	omg2 = omg_cyl * 2.D0

	return
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine input

	include "size.inc"
        include "mpif.h"
        include "mpi.inc"
	include "para.inc"
	include "ns.inc"
	include "metric.inc"
        include "sedi.inc"

	character*4 ID 
	write(ID, fmt='(I3)') 100+myid
        open(myid, file='./out180.9-210.9sec-E0.005/out.'//ID,
     <       form='unformatted',status='unknown')

	read(myid) kount
	read(myid) vis, ak
	read(myid) time
	read(myid) u
	read(myid) uxi
	read(myid) uej
	read(myid) uzk
	read(myid) p
	read(myid) xp
	read(myid) depth
	if ( iscalar .eq. 1 ) read(myid) phi
        if ( sedi    .eq. 1 ) read(myid) sedc

	close(unit = myid)
	
	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine storeData
        
        include "size.inc"
        include "mpif.h"
        include "mpi.inc"
	include "para.inc"
	include "ns.inc"
	include "metric.inc"
        include "sedi.inc"
	
	character*4 ID 
	
	write(ID, fmt='(I3)') 100+myid
        open(myid, file='./out.'//ID,
     <       form='unformatted',status='unknown')

      
	write(myid) kount
	write(myid) vis, ak
	write(myid) time
	write(myid) u
	write(myid) uxi
	write(myid) uej
	write(myid) uzk
	write(myid) p
	write(myid) xp
	write(myid) depth
	if ( iscalar .eq. 1 ) write(myid) phi
        if ( sedi    .eq. 1 ) write(myid) sedc

        close(unit = myid)
	return
	end
     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine output

	include "size.inc"
        include "mpif.h"
        include "mpi.inc"
	include "para.inc"
	include "ns.inc"
	include "metric.inc"
	include "eddy.inc"
	include "sedi.inc"

	character*4 ID

c	write(ID, fmt='(I3)') 100+myid

	if(istep .eq. nsave) then
	   open(500+myid, form='unformatted',
     <          status='unknown')
	   open(100+myid, form='unformatted',
     <          status='unknown')
	endif

	if(istep .eq. nsave) then
	   write(500+myid) xp
	endif

        write(100+myid) u
	write(100+myid) p
c	write(100+myid) sai
c	write(100+myid) aksd



	if (iscalar .eq. 1) write(100+myid) phi
	if (sedi .eq. 1) then
           write(100+myid) sedc
c	   write(100+myid) sedv
	endif
	if (istep .eq. nstep) then
	   close(unit = 100+myid)
	   close(unit = 500+myid)
	endif

	return
 	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine cfl_check

	include "size.inc"
        include "mpif.h"
	include "mpi.inc"
	include "para.inc"
	include "ns.inc"
	include "metric.inc"

	integer i, j, k
	double precision cflmax, cfl, temp

	cflmax = 0.D0
	cfl = 0.D0
	
	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni
           temp = dabs( uxi(i-1,j,k) + uxi(i,j,k) 
     <                 -kxi(i-1,j,k) - kxi(i,j,k)) 
     <          + dabs( uej(i,j-1,k) + uej(i,j,k)
     <                 -kej(i,j-1,k) - kej(i,j,k)
     <                 +wej(i,j-1,k) + wej(i,j,k))
     <          + dabs( uzk(i,j,k-1) + uzk(i,j,k) 
     <                 -kzk(i,j,k-1) - kzk(i,j,k))
c	   temp = dabs( uxi(i-1,j,k) + uxi(i,j,k))
c     <          + dabs( uej(i,j-1,k) + uej(i,j,k)) 
c     <          + dabs( uzk(i,j,k-1) + uzk(i,j,k)) 
	   temp = temp * jac(i,j,k)
           cfl = max(cfl, temp)
	enddo
	enddo
	enddo
	
	call MPI_REDUCE(cfl,cflmax,1,MPI_DOUBLE_PRECISION,
     <	                MPI_MAX,0,comm3d,ierr)

	if ( myid .eq. 0 ) cflmax = 0.5D0 * dtime * cflmax

	call MPI_BCAST(cflmax,1,MPI_DOUBLE_PRECISION,0,comm3d,ierr)
c	if (istep .eq. 1) then
c	   dtime = 0.2/cflmax*dtime
c	else
c	   if (istep .eq. 2) then
c	      dtime = 0.4/cflmax*dtime
c	   else
c	      dtime = 0.6/cflmax*dtime
c	   endif
c	endif

	if (istep .eq. 1) then
	   cflsuper = cflmax
	else
	   
	   cflsuper = max(cflmax, cflsuper)
	endif
	if ( cflmax .gt. 0.90 ) then
	   if ( MYID .eq. 0 )
     <	      write(*,*) ' CFL = ', cflmax, ' > 0.9, Stop! '
	   call storeData
	   write(*,*) 'Done with writing data'
	   stop
	else
	   if ( MYID .eq. 0 ) 
     <	      write(*,1) cflmax, cflsuper, dtime, time
	endif

    1	format(' CFL = ', f5.3, ' CFL_max = ', f5.3, 
     <        '  dtime = ', e9.3, '  time = ', e9.3)

	return
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine diff_cfl_check
        
	include "size.inc"
        include "mpif.h"
	include "mpi.inc"
	include "para.inc"
	include "ns.inc"
	include "metric.inc"
	include "eddy.inc"

	integer i, j, k
	double precision diff_cfl_max, diff_cfl, temp 

	diff_cfl_max = 0.D0
	diff_cfl = 0.D0

	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni
           temp =4.D0*( 1.D0/(xp(i+1,j,k,1)-xp(i-1,j,k,1))**2 
     <                + 1.D0/(xp(i,j+1,k,2)-xp(i,j-1,k,2))**2 
     <                + 1.D0/(xp(i,j,k+1,3)-xp(i,j,k-1,3))**2)
	   temp = temp*(vis + vst(i,j,k))
           diff_cfl = max(diff_cfl, temp)
	enddo
	enddo
	enddo
        diff_cfl = diff_cfl * dtime
	call MPI_REDUCE(diff_cfl,diff_cfl_max,1,MPI_DOUBLE_PRECISION,
     <	                MPI_MAX,0,comm3d,ierr)

	call MPI_BCAST(diff_cfl_max,1,MPI_DOUBLE_PRECISION,0,comm3d,ierr)
  

        if ( diff_cfl_max .gt. 0.5 ) then
	   if ( MYID .eq. 0 )
     <	      write(*,*) ' S = ', diff_cfl_max, ' > 0.5, Stop! '
	   stop
	else
	   if ( MYID .eq. 0 ) 
     <	      write(*,1) diff_cfl_max, dtime, time
	endif
    1	format(' S   = ', f5.3, '  dtime = ', e9.3, '  time = ', e9.3)
  
	return
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
