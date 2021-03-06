cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	program ns

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "para.inc"
	include "sed.inc"

	double precision ta, tt, t0, t1, t2, t3, t4, t5, t6, t7, t8
	double precision :: PI=3.1415926535897932
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
	t7 = 0.D0
	t8 = 0.D0

	call parameter !defines computational parameters and broadcasts them
	call mpi_initial !this sets up the partitioning and creates a map for the processors.
	call grid !This sets up the grid and calculates all the grid variables needed including the inverse Jacobian and the  mesh skewness tensor.
	call output_xyz !This writes the x, y, and z coordinates of the created grid.
C	call init_pSteady !initialize steady pressure gradient
	call getUtheo !This computes the steady state profile from the constant pressure gardient
	call initial !This initializes velocities, density field, and turbulence properties
	call readMeanProfiles !Reads in mean velocity profiles
	call output_Utheo !Outputs the theoretical log profile
	if(iTKE .eq. 1) call computeMeanAndPrimes !Compute mean and perturbation values if tracking volume integrated TKE

	do istep = 1, nstep !time stepping of simulation
	   if ( MYID .EQ. 0 .AND. MOD(istep,200) .EQ. 0) then
	      write(*,*) ' istep = ', istep, ' kount  = ', kount
	   end if

	   if ( pAdjust .eq. 0 ) then
	   call MPI_Barrier(MPI_COMM_WORLD, ierr)
	   call adjustPressure(PI)
	   call MPI_Barrier(MPI_COMM_WORLD, ierr)
	   end if


	   if(mod(kount-1,nsavePro) .eq. 0) then
	   call MPI_Barrier(MPI_COMM_WORLD, ierr)
	   tt =  MPI_Wtime()
	   call output_profiles(PI)
	   call MPI_Barrier(MPI_COMM_WORLD, ierr)
	   t6 = t6 + MPI_Wtime() - tt
	   end if

C	   if(mod(kount-2,nsavePro) .eq. 0) then
C	   call MPI_Barrier(MPI_COMM_WORLD, 
C	   tt =  MPI_Wtime()
C	   call output_profilesTwo(PI)
C	   call MPI_Barrier(MPI_COMM_WORLD, ierr)
C	   t6 = t6 + MPI_Wtime() - tt 

	   if(mod(kount-1,nsave) .eq. 0) then
	      call MPI_Barrier(MPI_COMM_WORLD, ierr)
	      tt =  MPI_Wtime()
	      call output(PI) !writes density, and velocity field
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

C	Compute right hand side of sediment transport model  
	   call MPI_Barrier(MPI_COMM_WORLD, ierr)
	   tt = MPI_Wtime()
	   if ( ised .eq. 1 ) call Csed_rhs
	   call MPI_Barrier(MPI_COMM_WORLD, ierr)
	   t7 = t7 + MPI_Wtime() - tt

C	This setups the source term for u_star, then applies approximate factorization to solve for u_star.   
	   call MPI_Barrier(MPI_COMM_WORLD, ierr)
	   tt = MPI_Wtime()
	   call predictor
	   call MPI_Barrier(MPI_COMM_WORLD, ierr)
	   t1 = t1 + MPI_Wtime() - tt
	   
	   call MPI_Barrier(MPI_COMM_WORLD, ierr)
	   tt = MPI_Wtime()
	   call pressure
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
	   
C          Solve for new sediment concentration
	   call MPI_Barrier(MPI_COMM_WORLD, ierr)
	   tt = MPI_Wtime()
	   if ( ised .eq. 1 ) call Csed_solve
	   call MPI_Barrier(MPI_COMM_WORLD, ierr)
	   t8 = t8 + MPI_Wtime() - tt
	   
	   if (iTKE .eq. 1) then
	   call MPI_Barrier(MPI_COMM_WORLD, ierr)
	   tt =  MPI_Wtime()
	   call getProAndDis
	   t6 = t6 + MPI_Wtime() - tt
	   call MPI_Barrier(MPI_COMM_WORLD, ierr)
	   endif

	   call MPI_Barrier(MPI_COMM_WORLD, ierr)
	   call eqstate
	   call MPI_Barrier(MPI_COMM_WORLD, ierr)
	  
	   kount = kount + 1
	   time = time + dtime

	   call cfl_check
 
	   if(mod(istep,nsave) .eq. 0) then
           call MPI_Barrier(MPI_COMM_WORLD, ierr)
           call output_continue_run
           call MPI_Barrier(MPI_COMM_WORLD, ierr)
           endif
	enddo

	call MPI_Barrier(MPI_COMM_WORLD, ierr)
	call output_continue_run
	call MPI_Barrier(MPI_COMM_WORLD, ierr)

	ta = MPI_Wtime() - ta


	if ( MYID .eq. 0 ) then
     	   WRITE(*,*) ' Total CPU time = ', ta, ta/nstep		
     	   WRITE(*,*) ' eddy:      ', t0/nstep, int(t0/ta*100.D0), '%'		
     	   WRITE(*,*) ' Predictor: ', t1/nstep, int(t1/ta*100.D0), '%'		
     	   WRITE(*,*) ' Pressure:  ', t2/nstep, int(t2/ta*100.D0), '%'		
     	   WRITE(*,*) ' Corrector: ', t3/nstep, int(t3/ta*100.D0), '%'		
     	   WRITE(*,*) ' Scalarp:   ', t5/nstep, int(t5/ta*100.D0), '%' 		
     	   WRITE(*,*) ' Scalars:   ', t4/nstep, int(t4/ta*100.D0), '%' 	
	   WRITE(*,*) ' SedimentRHS', t7/nstep, int(t7/ta*100.D0), '%'
	   WRITE(*,*) ' SedimentSol', t8/nstep, int(t8/ta*100.D0), '%'
     	   WRITE(*,*) ' Output:    ', t6/nstep, int(t6/ta*100.D0), '%' 		
	endif

	call MPI_FINALIZE( ierr )

	stop
	end
