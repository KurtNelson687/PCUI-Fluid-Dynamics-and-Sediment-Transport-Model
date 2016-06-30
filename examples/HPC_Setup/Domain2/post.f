
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
	



	subroutine get_pertVel(u, uMean, vMean, wMean, velPrimes)
C       This subroutine routine computes velociy purturbations

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"

	double precision, intent(in) :: u(-1:nni+2,-1:nnj+2,-1:nnk+2,3),
     <        uMean(1:nnj), vMean(1:nnj), wMean(1:nnj)
	double precision, intent(out) :: velPrimes(1:nni, 1:nnj, 1:nnk,3)
	
	integer i, j, k
	
	do i = 1, nni
	   do j = 1, nnj
	      do k = 1, nnk
	      velPrimes(i,j,k,1) = u(i,j,k,1)-uMean(j)
	      velPrimes(i,j,k,2) = u(i,j,k,2)-vMean(j)
	      velPrimes(i,j,k,3) = u(i,j,k,3)-wMean(j)
	      enddo
	   enddo
	enddo
	
	return
	end


	subroutine get_turbIntensity(velPrimes, uTurb, vTurb, wTurb)
C       This subroutine routine computes mean turbulent intensity profiles

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"

	double precision, intent(in) :: velPrimes(1:nni, 1:nnj, 1:nnk,3)
	double precision, intent(out) :: uTurb(1:nnj), vTurb(1:nnj),
     <      wTurb(1:nnj)
	double precision velPrimesS(1:nni, 1:nnj, 1:nnk,3)
	integer i, j, k
	
	do i = 1, nni
	   do j = 1, nnj
	      do k = 1, nnk
	      velPrimesS(i,j,k,1) = velPrimes(i,j,k,1)**2
	      velPrimesS(i,j,k,2) = velPrimes(i,j,k,2)**2
	      velPrimesS(i,j,k,3) = velPrimes(i,j,k,3)**2
	      enddo
	   enddo
	enddo
	
	call MPI_Barrier(MPI_COMM_WORLD, ierr)
	call horizontalAverage(velPrimesS(:,:,:,1),uTurb,0) 
	call horizontalAverage(velPrimesS(:,:,:,2),vTurb,0) 
	call horizontalAverage(velPrimesS(:,:,:,3),wTurb,0) 
	
	do j = 1, nnj
	   uTurb(j) = SQRT(uTurb(j))
	   vTurb(j) = SQRT(vTurb(j))
	   wTurb(j) = SQRT(wTurb(j))
	enddo

	return
	end



	subroutine get_reynoldsStress(velPrimes, uvPrime, uwPrime, vwPrime)
C       This subroutine routine computes velociy purturbations

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"

	double precision, intent(in) :: velPrimes(1:nni, 1:nnj, 1:nnk,3)
	double precision, intent(out) :: uvPrime(1:nnj), uwPrime(1:nnj),
     <      vwPrime(1:nnj)
	double precision reynoldsStresses(1:nni, 1:nnj, 1:nnk,3)
	integer i, j, k
	
	do i = 1, nni
	   do j = 1, nnj
	      do k = 1, nnk
	      reynoldsStresses(i,j,k,1) = velPrimes(i,j,k,1)*velPrimes(i,j,k,2)
	      reynoldsStresses(i,j,k,2) = velPrimes(i,j,k,1)*velPrimes(i,j,k,3)
	      reynoldsStresses(i,j,k,3) = velPrimes(i,j,k,3)*velPrimes(i,j,k,2)
	      enddo
	   enddo
	enddo

	call MPI_Barrier(MPI_COMM_WORLD, ierr)
	call horizontalAverage(reynoldsStresses(:,:,:,1),uvPrime,0) 
	call horizontalAverage(reynoldsStresses(:,:,:,2),uwPrime,0) 
	call horizontalAverage(reynoldsStresses(:,:,:,3),vwPrime,0) 
	
	uvPrime = -1*uvPrime
	uwPrime = -1*uwPrime
	vwPrime = -1*vwPrime
	return
	end
