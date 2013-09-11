cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine mpi_initial

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "para.inc"

	integer nprocs, ndims
	integer dims(3), coords(3)
	logical reorder 

	call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
	call MPI_COMM_SIZE( MPI_COMM_WORLD, nprocs, ierr )

	ndims = 3
C	call MPI_DIMS_CREATE( nprocs, ndims, dims, ierr )
	dims(1) = px
	dims(2) = py
	dims(3) = pz

	periods(1) = .false.
	periods(2) = .false.
	periods(3) = .false.
	if ( periodic .eq. 1 ) periods(3) = .true.
	reorder    = .true.
	call MPI_CART_CREATE( MPI_COMM_WORLD, ndims, dims, 
     <	                      periods, reorder, comm3d, ierr )

	coords(1) = 0
	coords(2) = 0
	coords(3) = 0
	call MPI_CART_COORDS( comm3d, myid, 3, coords, ierr )
	npx = coords(1)
	npy = coords(2)
	npz = coords(3)
	
c	WRITE(*,*) myid, "Coordinates ", npx, npy, npz

	call MPI_Cart_shift( comm3d, 0, 1, n_west, n_east, ierr )
	call MPI_Cart_shift( comm3d, 1, 1, n_suth, n_nrth, ierr )
	call MPI_Cart_shift( comm3d, 2, 1, n_back, n_frnt, ierr )

	
	WRITE(*,420) myid, "Neighbor (W E S N B F)", 
     <	           n_west, n_east, n_suth, n_nrth, n_back, n_frnt

 420	format (i6,a,i6,i6,i6,i6,i6,i6)

	if ( n_west .ne. MPI_PROC_NULL ) then
	   ius = 0
	else
	   ius = 1
	endif

	if ( n_east .ne. MPI_PROC_NULL ) then
	   iue = nni
	else
           iue = nni-1
	endif

	if ( n_suth .ne. MPI_PROC_NULL ) then
           jus = 0
	else
           jus = 1
	endif

	if ( n_nrth .ne. MPI_PROC_NULL ) then
	   jue = nnj
	else
           jue = nnj-1
	endif

	if ( n_back .ne. MPI_PROC_NULL ) then
           kus = 0
	else
           kus = 1
	endif

	if ( n_frnt .ne. MPI_PROC_NULL ) then
	   kue = nnk
	else
           kue = nnk-1
	endif

	return
	end
