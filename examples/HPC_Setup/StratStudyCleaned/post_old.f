
	subroutine horizontalAverage(inArray, outArray, numGhost)
C       This subroutine horizontally Averages an array. It is written for cartesian grids with uniform hosizontal spacing

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"

	integer, intent(in) :: numGhost 
	double precision, intent(in) :: inArray(1-numGhost:nni+numGhost,
     <                1-numGhost:nnj+numGhost,1-numGhost:nnk+numGhost)
	double precision, intent(out) :: outArray(1:nnj)
	
	double precision isum(1:nnj,1:nnk),myArray(1:nnj)

	isum = sum(inArray(1:nni,1:nnj,1:nnk), DIM = 1)
	myArray = sum(isum(1:nnj,1:nnk), DIM = 2)
	myArray = myArray/(nni*nnk)
	call MPI_REDUCE(myArray, outArray, nnj, MPI_DOUBLE_PRECISION,
     <               MPI_SUM, 0, hor_comm, ierr)
	
	outArray = outArray/(px*pz)
	call MPI_BCAST(outArray,nnj,MPI_DOUBLE_PRECISION,0,hor_comm,ierr)
	return
	end
