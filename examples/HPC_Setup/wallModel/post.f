cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine depthAverage(inArray, outValue)
C       This subroutine depth averages an array that has been horizontally averaged
	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "metric.inc"
	include "cavity.inc"
	double precision, intent(in) :: inArray(1:nnj)
	double precision, intent(out) :: outValue
	double precision H, mySum
	double precision integrationArray(1:nnj)
	integer j
	
	
	H = yAll(nj)+(yAll(nj)-yAll(nj-1))/2
	do j = 1, nnj
	    integrationArray(j) = 1/H*inArray(j)*0.5*
     <         (xp(1,j+1,1,2)-xp(1,j-1,1,2))
	enddo
	
	mySum = sum(integrationArray)

	call MPI_REDUCE(mySum, outValue, 1, MPI_DOUBLE_PRECISION,
     <               MPI_SUM, 0, vert_comm, ierr)
	call MPI_BCAST(outValue,1,MPI_DOUBLE_PRECISION,0,vert_comm,ierr)
	return
	end
	

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine volumeAve(inArray, outValue,numGhost)
C       This subroutine volume averagesover the entire domain
	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "metric.inc"
	include "cavity.inc"
	

	double precision, intent(out) :: outValue
	integer, intent(in) :: numGhost 
	double precision, intent(in) :: inArray(1-numGhost:nni+numGhost,
     <                1-numGhost:nnj+numGhost,1-numGhost:nnk+numGhost)
	double precision myIntegration
	integer i,j,k
	
	do i = 1, nni 
	   do j = 1, nnj 
	      do k = 1, nnk
	        myIntegration = myIntegration
     <                 +inArray(i,j,k)/jac(i,j,k)/domainVol
	      enddo
	    enddo
	enddo
	
	call MPI_REDUCE(myIntegration, outValue, 1, MPI_DOUBLE_PRECISION,
     <               MPI_SUM, 0, MPI_COMM_WORLD, ierr)
	call MPI_BCAST(outValue,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
	return
	end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


	subroutine getDomainVol(inArray, outValue,numGhost)
C       This subroutine sums over the domain
	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "metric.inc"
	include "cavity.inc"
	

	double precision, intent(out) :: outValue
	integer, intent(in) :: numGhost 
	double precision, intent(in) :: inArray(1-numGhost:nni+numGhost,
     <                1-numGhost:nnj+numGhost,1-numGhost:nnk+numGhost)
	double precision mySum
	integer i,j,k

	mySum = 0.D0
	do i = 1, nni 
	   do j = 1, nnj 
	      do k = 1, nnk
	        mySum = mySum+1/inArray(i,j,k)
	      enddo
	    enddo
	enddo
	
	call MPI_REDUCE(mySum, outValue, 1, MPI_DOUBLE_PRECISION,
     <               MPI_SUM, 0, MPI_COMM_WORLD, ierr)
	call MPI_BCAST(outValue,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
	return
	end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
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


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	subroutine get_pertVel(u, uMean, vMean, wMean, velPrimes)
C       This subroutine routine computes velociy purturbations

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"

	double precision, intent(in) :: u(-1:nni+2,-1:nnj+2,-1:nnk+2,3),
     <        uMean(0:nnj+1), vMean(0:nnj+1), wMean(0:nnj+1)
	double precision, intent(out) :: velPrimes(0:nni+1, 0:nnj+1, 0:nnk+1,3)
	
	integer i, j, k
	
	do i = 0, nni+1
	   do j = 0, nnj+1
	      do k = 0, nnk+1
	      velPrimes(i,j,k,1) = u(i,j,k,1)-uMean(j)
	      velPrimes(i,j,k,2) = u(i,j,k,2)-vMean(j)
	      velPrimes(i,j,k,3) = u(i,j,k,3)-wMean(j)
	      enddo
	   enddo
	enddo
	
	return
	end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine get_turbIntensity(velPrimes, uTurb, vTurb, wTurb)
C       This subroutine routine computes mean turbulent intensity profiles

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"

	double precision, intent(in) :: velPrimes(0:nni+1, 0:nnj+1, 0:nnk+1,3)
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

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine get_reynoldsStress(velPrimes, uvPrime, uwPrime, vwPrime)
C       This subroutine routine computes velociy purturbations

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"

	double precision, intent(in) :: velPrimes(0:nni+1, 0:nnj+1, 0:nnk+1,3)
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


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine compute_totalkinetic(input,output)
C       This subroutine routine computes velociy purturbations

	include "size.inc"

	double precision, intent(in) :: input(-1:nni+2,-1:nnj+2,-1:nnk+2,3)
	double precision, intent(out) :: output(-1:nni+2,-1:nnj+2,-1:nnk+2)
	integer i, j, k
	
	do i = -1, nni+2
	   do j = -1, nnj+2
	      do k = -1, nnk+2
	      output(i,j,k) = 0.5*(input(i,j,k,1)**2
     <          + input(i,j,k,2)**2+input(i,j,k,3)**2) 
	      enddo
	   enddo
	enddo
	return
	end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine compute_dissipation(velPrimes,dissipationMean)
C       This subroutine routine computes a mean dissipation profile

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "para.inc"
	include "metric.inc"
	double precision, intent(in) :: velPrimes(0:nni+1, 0:nnj+1, 0:nnk+1,3)
	double precision, intent(out) :: dissipationMean(1:nnj)
	double precision dissipationAll(1:nni, 1:nnj, 1:nnk)
	integer i, j, k
	
	do i = 1, nni
	   do j = 1, nnj
	      do k = 1, nnk
	      dissipationAll(i,j,k) =vis*(
     <          ((velPrimes(i+1,j,k,1)-velPrimes(i-1,j,k,1))
     <               /(xp(i+1,j,k,1)-xp(i-1,j,k,1)))**2+
     <          ((velPrimes(i,j+1,k,1)-velPrimes(i,j-1,k,1))
     <               /(xp(i,j+1,k,2)-xp(i,j-1,k,2)))**2+
     <          ((velPrimes(i,j,k+1,1)-velPrimes(i,j,k-1,1))
     <               /(xp(i,j,k+1,3)-xp(i,j,k-1,3)))**2+

     <          ((velPrimes(i+1,j,k,2)-velPrimes(i-1,j,k,2))
     <               /(xp(i+1,j,k,1)-xp(i-1,j,k,1)))**2+
     <          ((velPrimes(i,j+1,k,2)-velPrimes(i,j-1,k,2))
     <               /(xp(i,j+1,k,2)-xp(i,j-1,k,2)))**2+
     <          ((velPrimes(i,j,k+1,2)-velPrimes(i,j,k-1,2))
     <               /(xp(i,j,k+1,3)-xp(i,j,k-1,3)))**2+
     
     <          ((velPrimes(i+1,j,k,3)-velPrimes(i-1,j,k,3))
     <               /(xp(i+1,j,k,1)-xp(i-1,j,k,1)))**2+
     <          ((velPrimes(i,j+1,k,3)-velPrimes(i,j-1,k,3))
     <               /(xp(i,j+1,k,2)-xp(i,j-1,k,2)))**2+
     <          ((velPrimes(i,j,k+1,3)-velPrimes(i,j,k-1,3))
     <               /(xp(i,j,k+1,3)-xp(i,j,k-1,3)))**2)
	      enddo
	   enddo
	enddo

	call MPI_Barrier(MPI_COMM_WORLD, ierr)
	call horizontalAverage(dissipationAll,dissipationMean,0) 
	return
	end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine getMeanVelPro(inArray, outArray)

C       This subroutine horizontally Averages an array. It is written for cartesian grids with uniform hosizontal spacing

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"

	double precision, intent(in) :: inArray(-1:nni+2,
     <                -1:nnj+2,-1:nnk+2)
	double precision, intent(out) :: outArray(0:nnj+1)
	
	double precision isum(0:nnj+1,1:nnk),myArray(0:nnj+1)

	isum = sum(inArray(1:nni,0:nnj+1,1:nnk), DIM = 1)
	myArray = sum(isum(0:nnj+1,1:nnk), DIM = 2)
	myArray = myArray/(nni*nnk)
	call MPI_REDUCE(myArray, outArray, nnj+2, MPI_DOUBLE_PRECISION,
     <               MPI_SUM, 0, hor_comm, ierr)
	
	outArray = outArray/(px*pz)
	call MPI_BCAST(outArray,nnj+2,MPI_DOUBLE_PRECISION,0,hor_comm,ierr)
	return
	end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
