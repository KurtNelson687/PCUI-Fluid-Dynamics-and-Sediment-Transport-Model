cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	program ns

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "para.inc"

        double precision ta, tt, t0, t1, t2, t3, t4, t5, t6

	call MPI_INIT( ierr )!this intializes the mpi environment

	call MPI_Barrier(MPI_COMM_WORLD, ierr)
	ta = MPI_Wtime()

	t0 = 0.D0
	t1 = 0.D0
	t2 = 0.D0
	t3 = 0.D0
	t4 = 0.D0
	t5 = 0.D0
	t6 = 0.D0

	call parameter

	call mpi_initial

	call grid

	call initial

	CALL STAT_INIT

	do istep = 1, nstep

C	   if ( mod(istep, nsave) .eq. 0 .and. MYID .EQ. 0 )
	   if ( MYID .EQ. 0 )
     <	      write(*,*) ' istep = ', istep, ' kount  = ', kount

	   call MPI_Barrier(MPI_COMM_WORLD, ierr)
           tt = MPI_Wtime()
           if ( ieddy .eq. 1 ) call eddy
	   call MPI_Barrier(MPI_COMM_WORLD, ierr)
	   t0 = t0 + MPI_Wtime() - tt

	   call MPI_Barrier(MPI_COMM_WORLD, ierr)
	   tt = MPI_Wtime()
           if ( iscalar .eq. 1 ) call scalar_rhs
	   call MPI_Barrier(MPI_COMM_WORLD, ierr)
	   t5 = t5 + MPI_Wtime() - tt

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
	
	   kount = kount + 1
	   time = time + dtime

	   call cfl_check

	   CALL STAT_COLL

	enddo

	CALL STAT_POST

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
