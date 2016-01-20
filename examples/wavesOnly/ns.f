cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	program ns

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "para.inc"

	double precision ta, tt, t0, t1, t2, t3, t4, t5, t6
	real, parameter :: PI=3.1415926535897932
	call MPI_INIT( ierr ) !initializes the MPI execution environment

	call MPI_Barrier(MPI_COMM_WORLD, ierr) !stops program until all processes have reached this routine. MPI_COMM_WORLD is the type of communication used 
	ta = MPI_Wtime() !Time in seconds since some arbitrary time

C	Time variables used to track how much time each major component of the code is taking
	t0 = 0.D0
	t1 = 0.D0
	t2 = 0.D0
	t3 = 0.D0
	t4 = 0.D0
	t5 = 0.D0
	t6 = 0.D0

	call parameter !defines computational parameters and broadcasts them
	call mpi_initial !this sets up the partitioning and creates a map for the processors.
	call grid !This sets up the grid and calculates all the grid variables needed including the inverse Jacobian and the  mesh skewness tensor.
	call output_xyz !This writes the x, y, and z coordinates of the created grid.
	call initial !This initializes velocities, density field, and turbulence properties

	do istep = 1, nstep !time stepping of simulation

C          if ( mod(istep, nsave) .eq. 0 .and. MYID .EQ. 0 )
	   
	   if ( MYID .EQ. 0 .AND. MOD(istep,10) .EQ. 0) then
	      write(*,*) ' istep = ', istep, ' kount  = ', kount
	   end if

C	Compute pressure gradient from waves if presents
	   if ( waves .eq. 1 ) then
	      dpdxWave = -2*PI/Twave*waveMag*
     \               sin(2*PI*(istep-1)*dtime/Twave)
	   end if

	   if(mod(istep,nsave) .eq. 0 .or. istep .eq. 1) then
	      call MPI_Barrier(MPI_COMM_WORLD, ierr)
	      tt =  MPI_Wtime()
	      call output !writes density, and velocity field
	      call MPI_Barrier(MPI_COMM_WORLD, ierr)
	      t6 = t6 + MPI_Wtime() - tt
	   end if

C	Compute subgrid scale stress if turbulence is turned on (i.e. ieddy = 1 in io.f)
	   call MPI_Barrier(MPI_COMM_WORLD, ierr)
	   tt = MPI_Wtime()
	   if ( ieddy .eq. 1 ) call eddy
	   call MPI_Barrier(MPI_COMM_WORLD, ierr)
	   t0 = t0 + MPI_Wtime() - tt
	  
C	Compute  
	   call MPI_Barrier(MPI_COMM_WORLD, ierr)
	   tt = MPI_Wtime()
	   if ( iscalar .eq. 1 ) call scalar_rhs
	   call MPI_Barrier(MPI_COMM_WORLD, ierr)
	   t5 = t5 + MPI_Wtime() - tt

C	This setups the source term for u_star, then applies approximate factorization to solve for u_star. The tridiagonal matrix systems are solved with the tridiagonal matrix algorithm	   
	   call MPI_Barrier(MPI_COMM_WORLD, ierr)
	   tt = MPI_Wtime()
	   call predictor
	   call MPI_Barrier(MPI_COMM_WORLD, ierr)
	   t1 = t1 + MPI_Wtime() - tt
	   
	   call MPI_Barrier(MPI_COMM_WORLD, ierr)
	   tt = MPI_Wtime()
C	   call pressure
	   call MPI_Barrier(MPI_COMM_WORLD, ierr)
	   t2 = t2 + MPI_Wtime() - tt
	   
	   call MPI_Barrier(MPI_COMM_WORLD, ierr)
	   tt = MPI_Wtime()
	   call corrector
	   call MPI_Barrier(MPI_COMM_WORLD, ierr)
	   t3 = t3 + MPI_Wtime() - tt
	   
	   call MPI_Barrier(MPI_COMM_WORLD, ierr)
	   tt = MPI_Wtime()
	   if ( iscalar .eq. 1 ) call scalar_solve
	   call MPI_Barrier(MPI_COMM_WORLD, ierr)
	   t4 = t4 + MPI_Wtime() - tt
	   
	   kount = kount + 1
	   time = time + dtime

	   call cfl_check

c     call MPI_Barrier(MPI_COMM_WORLD, ierr)
c     call output_continue_run
c     call MPI_Barrier(MPI_COMM_WORLD, ierr)

	enddo

	call MPI_Barrier(MPI_COMM_WORLD, ierr)
	ta = MPI_Wtime() - ta

	call MPI_Barrier(MPI_COMM_WORLD, ierr)
	tt = MPI_Wtime()
	call output
	call MPI_Barrier(MPI_COMM_WORLD, ierr)
	t6 = t6 + MPI_Wtime() - tt

	if ( MYID .eq. 0 ) then
     	   WRITE(*,*) ' Total CPU time = ', ta, ta/nstep		
     	   WRITE(*,*) ' eddy:      ', t0/nstep, int(t0/ta*100.D0), '%'		
     	   WRITE(*,*) ' Predictor: ', t1/nstep, int(t1/ta*100.D0), '%'		
     	   WRITE(*,*) ' Pressure:  ', t2/nstep, int(t2/ta*100.D0), '%'		
     	   WRITE(*,*) ' Corrector: ', t3/nstep, int(t3/ta*100.D0), '%'		
     	   WRITE(*,*) ' Scalarp:   ', t5/nstep, int(t5/ta*100.D0), '%' 		
     	   WRITE(*,*) ' Scalars:   ', t4/nstep, int(t4/ta*100.D0), '%' 		
     	   WRITE(*,*) ' Output:    ', t6/nstep, int(t6/ta*100.D0), '%' 		
	endif

	call MPI_FINALIZE( ierr )

	stop
	end
