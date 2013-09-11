cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine parameter

	include "size.inc"
        include "mpif.h"
        include "mpi.inc"
	include "para.inc"

	integer :: i, nP
	double precision :: PI, omega, period

	call MPI_COMM_RANK( MPI_COMM_WORLD, MYID, ierr )

	if ( MYID .eq. 0 ) then

	dtime    = 0.004D0
	case     = 0
	newrun   = 1
	periodic = 0
	iscalar  = 1
	ieddy    = 1
	mg_level = 5
	nstep    = 20000
	nsave    = 40
	maxstep  = 10

	do i = 1, 5
           iterchk(i)  = 2
           maxiter(i)  = 20
           tol(i)      = 1.D0
           ter(i)      = 1.D0
           slowiter(i) = 0.6D0
        enddo

        slowiter(1) = 0.7D0
	maxiter(5)  = 30
        vis         = 10e-6
        ak          = 10e-6
        g           = 9.81D0
        omg_cyl     = 0
        omg_lid     = 0
        factor      = 1.0e-4
        phi1        = 0
        phi2        = 0
        yphi        = 0
        aphi        = 1.

	endif

	call MPI_BCAST(case,        1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(newrun,      1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(periodic,    1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(iscalar,     1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(ieddy,       1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(mg_level,    1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

	call MPI_BCAST(nstep,       1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(nsave,       1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

	call MPI_BCAST(maxstep,     1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(iterchk(1),  5,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(maxiter(1),  5,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

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
	call MPI_BCAST(tol(1),      5,MPI_DOUBLE_PRECISION,0,
     <                              MPI_COMM_WORLD,ierr)
	call MPI_BCAST(ter(1),      5,MPI_DOUBLE_PRECISION,0,
     <                              MPI_COMM_WORLD,ierr)
	call MPI_BCAST(slowiter(1), 5,MPI_DOUBLE_PRECISION,0,
     <                              MPI_COMM_WORLD,ierr)

	call MPI_BCAST(phi1,        1,MPI_DOUBLE_PRECISION,0,
     <                              MPI_COMM_WORLD,ierr)
	call MPI_BCAST(phi2,        1,MPI_DOUBLE_PRECISION,0,
     <                              MPI_COMM_WORLD,ierr)
	call MPI_BCAST(yphi,        1,MPI_DOUBLE_PRECISION,0,
     <                              MPI_COMM_WORLD,ierr)
	call MPI_BCAST(aphi,        1,MPI_DOUBLE_PRECISION,0,
     <                              MPI_COMM_WORLD,ierr)

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

	read(10+myid) kount
	read(10+myid) vis, ak
	read(10+myid) time
	read(10+myid) u
	read(10+myid) uxi
	read(10+myid) uej
	read(10+myid) uzk
	read(10+myid) p
	if ( iscalar .eq. 1 ) read(10+myid) phi

	return
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine store_X_params
        
        include "size.inc"
        include "mpif.h"
        include "mpi.inc"
	include "para.inc"
	include "ns.inc"
	include "metric.inc"
	include "cavity.inc"

	character*4 ID 

	if (myid .eq. 0) then
	   print *, 'Recording parameters ...'
	   open(123,file='params',form='unformatted',status='unknown')
	   write(123) bx
	   write(123) by
	   write(123) bz
	   write(123) vis
	   write(123) ak
	   write(123) dtime
	   write(123) nstep
	   write(123) nsave
	   write(123) ni
	   write(123) nj
	   write(123) nk
	   write(123) px
	   write(123) py
	   write(123) pz
	   close(unit = 123)
	end if

	write(ID, fmt='(I3)') 600+myid
	open(500+myid,file='xyz.'//ID,form='unformatted',status='unknown')
	write(500+myid) xp
        close(unit = 500+myid)
	
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
	
	character*4 :: ID

	write(ID, fmt='(I3)') 500+myid

	if (myid .eq. 0) print *, 'Recording to output files ...'

	if (istep.gt.1) then
	   open(50+myid, file='output_S.'//ID, form='unformatted',
     >          status='old',position='append')
	   open(200+myid, file='output_UVW.'//ID, form='unformatted',
     >          status='old',position='append')
	else
	   open(50+myid, file='output_S.'//ID, form='unformatted',
     >          status='unknown')	   
	   open(200+myid, file='output_UVW.'//ID, form='unformatted',
     >          status='unknown')	   
	end if

	write(200+myid) u
	write(50+myid) phi
	close(unit = 50+myid)
	close(unit = 200+myid)

	return
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine outputbackup

	include "size.inc"
        include "mpif.h"
        include "mpi.inc"
	include "para.inc"
	include "ns.inc"
	include "metric.inc"
	include "eddy.inc"

	write(50+myid) kount
	write(50+myid) vis, ak
	write(50+myid) time
	write(50+myid) u
	write(50+myid) uxi
	write(50+myid) uej
	write(50+myid) uzk
	write(50+myid) p

	if ( iscalar .eq. 1 ) write(50+myid) phi

	if ( ieddy .gt. 0 ) then
	   write(50+myid) vst
	   write(50+myid) akst
	   write(50+myid) sab
	   write(50+myid) jac
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

	integer*4 today(3), now(3)
	integer i, j, k
	double precision cflmax, cfl, temp 

	cflmax = 0.D0
	cfl = 0.D0

	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni
           temp = dabs( uxi(i-1,j,k) + uxi(i,j,k) ) 
     <          + dabs( uej(i,j-1,k) + uej(i,j,k) )
     <          + dabs( uzk(i,j,k-1) + uzk(i,j,k) )
	   temp = temp * jac(i,j,k)
           cfl = max(cfl, temp)
	enddo
	enddo
	enddo

	call MPI_REDUCE(cfl,cflmax,1,MPI_DOUBLE_PRECISION,
     <	                MPI_MAX,0,comm3d,ierr)

	if ( myid .eq. 0 ) cflmax = 0.5D0 * dtime * cflmax

	call MPI_BCAST(cflmax,1,MPI_DOUBLE_PRECISION,0,comm3d,ierr)

	if ( cflmax .gt. 0.9 .and. istep .gt. 5) then
	   if ( MYID .eq. 0 )
     <	      write(*,*) ' CFL = ', cflmax, ' > 0.9, Stop! '
	   stop
	else
	   if ( MYID .eq. 0 ) then
	      print *, ''
	      print *, '************************************************' 
     	      write(*,1) cflmax, dtime, time
	      print *, '************************************************'
	      print *, ''
	   end if
	endif
    1	format(' CFL = ', e15.9, '  dtime = ', e9.3, '  time = ', e9.3)

	if ( MYID .eq. 0) then
	   call idate(today)	! today(1)=day, (2)=month, (3)=year
	   call itime(now)	! now(1)=hour, (2)=minute, (3)=second
	   write ( *, 1000 )  today(2), today(1), today(3), now
	end if
 1000	format ( 'Date ', i2.2, '/', i2.2, '/', i4.4, '; time ',
     &	               i2.2, ':', i2.2, ':', i2.2 )

	   if (MYID .eq. 0) then
	      if(istep.eq.1) then
		 open(123, file='qoutput',status='replace')
	      else
		 open(123, file='qoutput',position='append')
	      endif

	      write(123,105) istep, nstep, time, cflmax	   
	      call idate(today)	! today(1)=day, (2)=month, (3)=year
	      call itime(now)	! now(1)=hour, (2)=minute, (3)=second
 105          format('Time step = ',i7,'/',i7', time = ',e10.3, 
     &                ', cfl = ',e10.3)
	      close(unit=123)
	   endif

	return
	end
