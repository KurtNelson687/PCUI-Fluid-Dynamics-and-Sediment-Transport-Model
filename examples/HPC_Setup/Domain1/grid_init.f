cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	program grid_init

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "para.inc"

	call MPI_INIT( ierr )

	call parameter
	call mpi_initial
	call grid
	call output_xyz

	if ( MYID .eq. 0 ) WRITE(*,*) ' Wrote the grid.'

	call MPI_FINALIZE( ierr )

	stop
	end
