cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine mpi_initial

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "para.inc"

	integer nprocs, ndims
	integer dims(3), coords(3)
	logical reorder 

	call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr ) !This tells each processor who they are - gives an identification to each processor
	call MPI_COMM_SIZE( MPI_COMM_WORLD, nprocs, ierr ) !This tells each processor how many others there are

	ndims = 3
C	call MPI_DIMS_CREATE( nprocs, ndims, dims, ierr )
	dims(1) = px !number of processors in x direction
	dims(2) = py
	dims(3) = pz

C	periods is a logical array specifying which dimensions the grid is periodic in
	periods(1) = .false.
	periods(2) = .false.
	periods(3) = .false.
!	if ( periodic .eq. 1 ) periods(3) = .true. !Commented by Kurt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Added by Kurt
	if ( periodic .eq. 1) then
	   periods(1) = .true.
	   periods(3) = .true.
	endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	reorder    = .true.
	call MPI_CART_CREATE( MPI_COMM_WORLD, ndims, dims, 
     <	                      periods, reorder, comm3d, ierr ) !This reorders the topology and creates a communicator storing the topology information. note: comm3d and ierr are outputs here. comm3d is the handle

	coords(1) = 0
	coords(2) = 0
	coords(3) = 0
	call MPI_CART_COORDS( comm3d, myid, 3, coords, ierr )!creates mapping of ranks to Cartesian coordinates. For a 3-d simulation, coords is an array of length 3, unique to each processor, containing its cartesian coordinates
	npx = coords(1) !unique coordinates for each processor.
	npy = coords(2)
	npz = coords(3)
	
c	WRITE(*,*) myid, "Coordinates ", npx, npy, npz

!These lines tell each processor who is around them. n_west gives the rank (i.e. myid) of who is 1 processor behind in the cartesian grid.
! If no one is there, then n_west = MPI_PROC_NULL = -1. If the boundaries are perioidic, then n_west equals the rank of the processor on the far left. 
	call MPI_Cart_shift( comm3d, 0, 1, n_west, n_east, ierr )
	call MPI_Cart_shift( comm3d, 1, 1, n_suth, n_nrth, ierr )
	call MPI_Cart_shift( comm3d, 2, 1, n_back, n_frnt, ierr )

	
c	WRITE(*,420) myid, "Neighbor (W E S N B F)", 
c     <	           n_west, n_east, n_suth, n_nrth, n_back, n_frnt

c 420	format (i6,a,i6,i6,i6,i6,i6,i6)

C	all these lines do is assign indices for convection and correction stepping.
	if ( n_west .ne. MPI_PROC_NULL ) then
	   ius = 0
	else
	   ius = 1
	endif

	if ( n_east .ne. MPI_PROC_NULL ) then
	   iue = nni !nni is the number of points on each processor 
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
