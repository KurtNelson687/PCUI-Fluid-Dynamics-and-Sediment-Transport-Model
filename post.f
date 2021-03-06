cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine depthAverage(inArray, outValue)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Name: depthAverage
      !
      ! Author: Kurt Nelson
      !
      ! Purpose: This subroutine depth averages a profile specified by
      ! "inArray" and returns the value of the depth averaging as a
      ! double precision in "outValue".  
      !
      ! Inputs: inArray (double precision)
      !
      ! Outputs: outValue (double precision)
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "metric.inc"
	include "cavity.inc"
	double precision, intent(in) :: inArray(1:nnj)
	double precision, intent(out) :: outValue
	double precision H, mySum
	double precision integrationArray(1:nnj)
	integer j, jj0

	!gives bottom reference index in j direction for processor at
        !height npy
	jj0 = npy * nnj  
	
        H = yAll(nj)+(yAll(nj)-yAll(nj-1))/2 !domain depth
	do j = 1, nnj !summation for average
	    integrationArray(j) = 1/H*inArray(j)*dyAll(jj0+j)  
	enddo
	
        !sum cells for processor
	mySum = sum(integrationArray)

        !sum over processors in the vertical direction.
	call MPI_REDUCE(mySum, outValue, 1, MPI_DOUBLE_PRECISION,
     <               MPI_SUM, 0, vert_comm, ierr)
	call MPI_BCAST(outValue,1,MPI_DOUBLE_PRECISION,0,vert_comm,ierr)
	return
	end
	

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	subroutine horizontalAverage(inArray, outArray, numGhost)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Name: horizontalAverage 
      !
      ! Author: Kurt Nelson
      !
      ! Purpose: This subroutine horizontally averages the 3d array specified by
      ! "inArray" and returns the horizontally averaged profile as a double
      ! precision array in "outValue". When averaging arrays with ghost
      ! points, specify the number of ghost points wiht "numGhost". Ghost
      ! points are excluded from the averaging. 
      !
      ! Inputs: inArray (double precision)
      !         numGhost (double precision)
      !
      ! Outputs: outValue (double precision)
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

	subroutine get_pertVel(u, uAve, velPrimes)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Name: get_pertVel
      !
      ! Author: Kurt Nelson
      !
      ! Purpose: This subroutine computes perturbation velocities. It
      ! takes in the 4d velocity array specified by "u" and subtracts
      ! out the horzontal average of the streamwise velocity component
      ! u(:,:,:,1) specified by "uAve". As written, this subroutine only
      ! works for square cartesian grids, and it assumes the horizontal
      ! mean of the spanwise u(:,:,:,3) and vertical u(:,:,:,2) velocity
      ! field is zero. This however can be easily generalized by modifying
      ! lines within the do loop.    
      !
      ! Inputs: u (double precision)
      !         uAve (double precision)
      !
      ! Outputs: velPrimes (double precision)
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "para.inc"

	double precision, intent(in) :: u(-1:nni+2,-1:nnj+2,-1:nnk+2,3),
     <      uAve(0:nnj+1)
	double precision, intent(out) :: 
     <      velPrimes(0:nni+1, 0:nnj+1, 0:nnk+1,3)
	
	integer i, j, k

	do i = 0, nni+1
	  do j = 0, nnj+1
	    do k = 0, nnk+1
	      velPrimes(i,j,k,1) = u(i,j,k,1)-uAve(j)
	      velPrimes(i,j,k,2) = u(i,j,k,2)
	      velPrimes(i,j,k,3) = u(i,j,k,3)
	    enddo
	  enddo
	enddo

	return
	end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine get_turbIntensity(uTurb, vTurb, wTurb)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Name: 
      !
      ! Author: Kurt Nelson
      !
      ! Purpose: This subroutine computes turbulent intensities.
      ! Because "velPrimes" is stored within a block in "padjust.inc" it
      ! is not passed.    
      !
      ! Inputs: inArray
      !
      ! Outputs: outValue
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C       This subroutine routine computes mean turbulent intensity profiles

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "padjust.inc"

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

	subroutine get_reynoldsStress(uvPrime, uwPrime, vwPrime)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Name: 
      !
      ! Author: Kurt Nelson
      !
      ! Purpose: 
      !
      ! Inputs: inArray
      !
      ! Outputs: outValue
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C       This subroutine routine computes velociy purturbations

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "padjust.inc"

	double precision, intent(out) :: uvPrime(1:nnj), uwPrime(1:nnj),
     <      vwPrime(1:nnj)
	double precision reynoldsStresses(1:nni, 1:nnj, 1:nnk,3)
	integer i, j, k
	
	do i = 1, nni
	   do j = 1, nnj
	      do k = 1, nnk
	      reynoldsStresses(i,j,k,1) = 
     <           velPrimes(i,j,k,1)*velPrimes(i,j,k,2)
	      reynoldsStresses(i,j,k,2) = 
     <           velPrimes(i,j,k,1)*velPrimes(i,j,k,3)
	      reynoldsStresses(i,j,k,3) = 
     <           velPrimes(i,j,k,3)*velPrimes(i,j,k,2)
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
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Name: 
      !
      ! Author: Kurt Nelson
      !
      ! Purpose: 
      !
      ! Inputs: inArray
      !
      ! Outputs: outValue
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C       This subroutine routine computes velociy purturbations

	include "size.inc"

	double precision, intent(in) ::
     <       input(-1:nni+2,-1:nnj+2,-1:nnk+2,3)
	double precision, intent(out) ::
     <       output(-1:nni+2,-1:nnj+2,-1:nnk+2)
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

	subroutine compute_dissipation(dissipationMean)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Name: 
      !
      ! Author: Kurt Nelson
      !
      ! Purpose: 
      !
      ! Inputs: inArray
      !
      ! Outputs: outValue
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C       This subroutine routine computes a mean dissipation profile

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "para.inc"
	include "metric.inc"
	include "padjust.inc"
	double precision, intent(out) :: dissipationMean(1:nnj)
	double precision dissipationAll(1:nni, 1:nnj, 1:nnk)
	integer i, j, k
	
	do i = 1, nni
	   do j = 1, nnj
	      do k = 1, nnk
	      dissipationAll(i,j,k) =vis*(
     <          ((velPrimes(i,j,k,1)-velPrimes(i-1,j,k,1))
     <               /(xp(i,j,k,1)-xp(i-1,j,k,1)))**2+
     <          ((velPrimes(i,j,k,1)-velPrimes(i,j-1,k,1))
     <               /(xp(i,j,k,2)-xp(i,j-1,k,2)))**2+
     <          ((velPrimes(i,j,k,1)-velPrimes(i,j,k-1,1))
     <               /(xp(i,j,k,3)-xp(i,j,k-1,3)))**2+

     <          ((velPrimes(i,j,k,2)-velPrimes(i-1,j,k,2))
     <               /(xp(i,j,k,1)-xp(i-1,j,k,1)))**2+
     <          ((velPrimes(i,j,k,2)-velPrimes(i,j-1,k,2))
     <               /(xp(i,j,k,2)-xp(i,j-1,k,2)))**2+
     <          ((velPrimes(i,j,k,2)-velPrimes(i,j,k-1,2))
     <               /(xp(i,j,k,3)-xp(i,j,k-1,3)))**2+
     
     <          ((velPrimes(i,j,k,3)-velPrimes(i-1,j,k,3))
     <               /(xp(i,j,k,1)-xp(i-1,j,k,1)))**2+
     <          ((velPrimes(i,j,k,3)-velPrimes(i,j-1,k,3))
     <               /(xp(i,j,k,2)-xp(i,j-1,k,2)))**2+
     <          ((velPrimes(i,j,k,3)-velPrimes(i,j,k-1,3))
     <               /(xp(i,j,k,3)-xp(i,j,k-1,3)))**2)
	      enddo
	   enddo
	enddo

	call MPI_Barrier(MPI_COMM_WORLD, ierr)
	call horizontalAverage(dissipationAll,dissipationMean,0) 
	return
	end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine getMeanVelPro(inArray, outArray)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Name: 
      !
      ! Author: Kurt Nelson
      !
      ! Purpose: 
      !
      ! Inputs: inArray
      !
      ! Outputs: outValue
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
	call MPI_BCAST(outArray,nnj+2,
     <      MPI_DOUBLE_PRECISION,0,hor_comm,ierr)
	return
	end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	subroutine computeMeanAndPrimes
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Name: 
      !
      ! Author: Kurt Nelson
      !
      ! Purpose: 
      !
      ! Inputs: inArray
      !
      ! Outputs: outValue
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "ns.inc"
	include "padjust.inc"
	include "para.inc"

	character*4 :: ID
	double precision, dimension(0:nnj+1) :: uMean, vMean, wMean
	double precision uDepth
	call getMeanVelPro(u(:,:,:,1), uMean)
	call getMeanVelPro(u(:,:,:,2), vMean)
	call getMeanVelPro(u(:,:,:,3), wMean)
	call MPI_Barrier(MPI_COMM_WORLD, ierr)
	call depthAverage(uMean(1:nnj), uDepth)

	if(iTKE .eq. 1) then
	  call get_pertVel(u, uMean,velPrimes)
	else
	  if(numPhase .eq. phasePerT) then
	    numPhase = 1
	  else
	    numPhase = numPhase + 1 
	  endif

	  call get_pertVel(u, uMeanAve(:,numPhase), velPrimes)
	endif

	if (100+vert_id .gt. 999) then
	  write(ID, fmt='(I4)') 100+vert_id
	else
	  write(ID, fmt='(I3)') 100+vert_id
	endif
	

	if (myid .eq. 0) then
	   if (istep .gt.1) then !appends to existing file - use this if you want continue run to create new files
	     open(50+myid, file='outputpVal_udepth.'//ID,
     >         form='unformatted',status='old',position='append')
	     write(50+myid) uDepth
	     close(unit = 50+myid)
	   else
	     open(50+myid, file='outputpVal_udepth.'//ID, 
     >          form='unformatted',status='unknown')
	     write(50+myid) uDepth
	     close(unit = 50+myid)
	   endif
	endif

	if (hor_id .eq. 0) then

	   if (istep .gt.1) then !appends to existing file - use this if you want continue run to create new files
	     open(50+myid, file='outputp_umean.'//ID, form='unformatted',
     >          status='old',position='append')
	     write(50+myid) uMean(1:nnj)
	     close(unit = 50+myid)

	     open(50+myid, file='outputp_vmean.'//ID, form='unformatted',
     >          status='old',position='append')
	     write(50+myid) vMean(1:nnj)
	     close(unit = 50+myid)

	     open(50+myid, file='outputp_wmean.'//ID, form='unformatted',
     >          status='old',position='append')
	     write(50+myid) wMean(1:nnj)
	     close(unit = 50+myid)
	   else
	     open(50+myid, file='outputp_umean.'//ID, form='unformatted',
     >          status='unknown')
	     write(50+myid) uMean(1:nnj)
	     close(unit = 50+myid)

	     open(50+myid, file='outputp_vmean.'//ID, form='unformatted',
     >          status='unknown')
	     write(50+myid) vMean(1:nnj)
	     close(unit = 50+myid)

	     open(50+myid, file='outputp_wmean.'//ID, form='unformatted',
     >          status='unknown')
	     write(50+myid) wMean(1:nnj)
	     close(unit = 50+myid)
	   endif
	endif

	return
	end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	subroutine getProAndDis
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Name: 
      !
      ! Author: Kurt Nelson
      !
      ! Purpose: 
      !
      ! Inputs: inArray
      !
      ! Outputs: outValue
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C       This subroutine Saves Cn and puts it into the Production matrix for computing global production

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "metric.inc"
	include "padjust.inc"
	include "para.inc"
	include "eddy.inc"
	include "ns.inc"


	integer i,j,k,m
	double precision, dimension(1:nnj) ::
     <     ProMean, temp, DisMean
	double precision ProDepth, DisDepth
	double precision, dimension(1:nni,1:nnj,1:nnk,1:3) ::
     <     OldPrimes
	double precision, dimension(1:nni,1:nnj,1:nnk) ::
     <     Produc, Dissip

CCCCCCCCCCCCC Add u^n parts to production and dissipation calculation

	do m = 1, 3
	call horizontalAverage(ConNew(:,:,:,m), temp, 0)
	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni
           ConNew(i,j,k,m) = ConNew(i,j,k,m)-temp(j)
	enddo
	enddo
	enddo 
	enddo

	do m = 1, 3
	call horizontalAverage(DisNew(:,:,:,m), temp, 0)
	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni
           DisNew(i,j,k,m) = DisNew(i,j,k,m)-temp(j)
	enddo
	enddo
	enddo 
	enddo


	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni
           Produc(i,j,k) =
     <       (1.5D0*ConNew(i,j,k,1)-0.5D0*ConOld(i,j,k,1))
     <           *velPrimes(i,j,k,1)+
     <       (1.5D0*ConNew(i,j,k,2)-0.5D0*ConOld(i,j,k,2))
     <           *velPrimes(i,j,k,2)+
     <       (1.5D0*ConNew(i,j,k,3)-0.5D0*ConOld(i,j,k,3))
     <           *velPrimes(i,j,k,3)
	enddo
	enddo
	enddo 

	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni
           Dissip(i,j,k) =
     <       (1.5D0*DisNew(i,j,k,1)-0.5D0*DisOld(i,j,k,1)+
     <         0.5D0*DisINew(i,j,k,1))*velPrimes(i,j,k,1)+
     <       (1.5D0*DisNew(i,j,k,2)-0.5D0*DisOld(i,j,k,2)+
     <         0.5D0*DisINew(i,j,k,2))*velPrimes(i,j,k,2)+
     <       (1.5D0*DisNew(i,j,k,3)-0.5D0*DisOld(i,j,k,3)+
     <         0.5D0*DisINew(i,j,k,3))*velPrimes(i,j,k,3)
	enddo
	enddo
	enddo 

	do m = 1,3
	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni
           OldPrimes(i,j,k,m) = velPrimes(i,j,k,m)
	enddo
	enddo
	enddo
	enddo

	call computeMeanAndPrimes

CCCCCCCCCCCCC Add u^(n+1) parts to Production calculation

	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni
           Produc(i,j,k) = Produc(i,j,k)+
     <       (1.5D0*ConNew(i,j,k,1)-0.5D0*ConOld(i,j,k,1))
     <           *velPrimes(i,j,k,1)+
     <       (1.5D0*ConNew(i,j,k,2)-0.5D0*ConOld(i,j,k,2))
     <           *velPrimes(i,j,k,2)+
     <       (1.5D0*ConNew(i,j,k,3)-0.5D0*ConOld(i,j,k,3))
     <           *velPrimes(i,j,k,3)

	   Produc(i,j,k) = 0.5D0*jac(i,j,k)*Produc(i,j,k)
	enddo
	enddo
	enddo 

	call horizontalAverage(Produc, ProMean, 0)
	call depthAverage(ProMean, ProDepth)


	do m = 1, 3
	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni
           ConOld(i,j,k,m) = ConNew(i,j,k,m)
	enddo
	enddo
	enddo 
	enddo


CCCCCCCCCCCCC Add u^(n+1) parts to Dissipation calculation

	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni
           Dissip(i,j,k) = Dissip(i,j,k)+
     <       (1.5D0*DisNew(i,j,k,1)-0.5D0*DisOld(i,j,k,1)+
     <         0.5D0*DisINew(i,j,k,1))*velPrimes(i,j,k,1)+
     <       (1.5D0*DisNew(i,j,k,2)-0.5D0*DisOld(i,j,k,2)+
     <         0.5D0*DisINew(i,j,k,2))*velPrimes(i,j,k,2)+
     <       (1.5D0*DisNew(i,j,k,3)-0.5D0*DisOld(i,j,k,3)+
     <         0.5D0*DisINew(i,j,k,3))*velPrimes(i,j,k,3)
	enddo
	enddo
	enddo 

	do m = 1, 3
	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni
	   DisINew(i,j,k,m) = 
     <         ( vis + 0.5D0*(vst(i,j,k) + vst(i+1,j,k)) ) * 
     <		g11(i,  j,  k  ) * ( u(i+1,j,  k,  m) - u(i,j,k,m) )  
     <        + ( vis + 0.5D0*(vst(i,j,k) + vst(i-1,j,k)) ) * 
     <		g11(i-1,j,  k  ) * ( u(i-1,j,  k,  m) - u(i,j,k,m) ) 
     <        + ( vis + 0.5D0*(vst(i,j,k) + vst(i,j+1,k)) ) * 
     <		g22(i,  j,  k  ) * ( u(i,  j+1,k,  m) - u(i,j,k,m) )  
     <        + ( vis + 0.5D0*(vst(i,j,k) + vst(i,j-1,k)) ) * 
     <		g22(i,  j-1,k  ) * ( u(i,  j-1,k,  m) - u(i,j,k,m) )  
     <        + ( vis + 0.5D0*(vst(i,j,k) + vst(i,j,k+1)) ) * 
     <		g33(i,  j,  k  ) * ( u(i,  j,  k+1,m) - u(i,j,k,m) )  
     <        + ( vis + 0.5D0*(vst(i,j,k) + vst(i,j,k-1)) ) * 
     <		g33(i,  j,  k-1) * ( u(i,  j,  k-1,m) - u(i,j,k,m) )  
	enddo
	enddo
	enddo 
	enddo 

	do m = 1, 3
	call horizontalAverage(DisINew(:,:,:,m), temp, 0)
	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni
           DisINew(i,j,k,m) = DisINew(i,j,k,m)-temp(j)
	enddo
	enddo
	enddo 
	enddo

	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni
           Dissip(i,j,k) = Dissip(i,j,k)+
     <         0.5D0*DisINew(i,j,k,1)*velPrimes(i,j,k,1)+
     <         0.5D0*DisINew(i,j,k,1)*OldPrimes(i,j,k,1)+
     <         0.5D0*DisINew(i,j,k,2)*velPrimes(i,j,k,2)+
     <         0.5D0*DisINew(i,j,k,2)*OldPrimes(i,j,k,2)+
     <         0.5D0*DisINew(i,j,k,3)*velPrimes(i,j,k,3)+
     <         0.5D0*DisINew(i,j,k,3)*OldPrimes(i,j,k,3)

	   Dissip(i,j,k) = 0.5D0*jac(i,j,k)*Dissip(i,j,k)
	enddo
	enddo
	enddo 

	call horizontalAverage(Dissip, DisMean, 0)
	call depthAverage(DisMean, DisDepth)

	if ( MYID .eq. 0 ) then

	   if (istep .gt.1) then !appends to existing file - use this if you want continue run to create new files
	     open(50+myid, file='outputpval_DisDepth.100',
     >          form='unformatted',status='old',position='append')
	     write(50+myid) DisDepth
	     close(unit = 50+myid)

	     open(50+myid, file='outputpval_ProDepth.100', 
     >          form='unformatted', status='old',position='append')
	     write(50+myid) ProDepth
	     close(unit = 50+myid)

	   else
	     open(50+myid, file='outputpval_DisDepth.100',
     >          form='unformatted',status='unknown')
	     write(50+myid) DisDepth
	     close(unit = 50+myid)

	     open(50+myid, file='outputpval_ProDepth.100', 
     >          form='unformatted',status='unknown')
	     write(50+myid) ProDepth
	     close(unit = 50+myid)

	   endif
	endif


	do m = 1, 3
	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni
           DisOld(i,j,k,m) = DisNew(i,j,k,m)
	enddo
	enddo
	enddo 
	enddo

	return
	end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	subroutine sedMass1(inArray, outValue,numGhost)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Name: 
      !
      ! Author: Kurt Nelson
      !
      ! Purpose: 
      !
      ! Inputs: inArray
      !
      ! Outputs: outValue
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

	myIntegration=0.D0
	
	do i = 1, nni 
	   do j = 1, nnj 
	      do k = 1, nnk
	        myIntegration = myIntegration
     <                 +inArray(i,j,k)/jac(i,j,k)
	      enddo
	    enddo
	enddo
	
	call MPI_REDUCE(myIntegration, outValue, 1, MPI_DOUBLE_PRECISION,
     <               MPI_SUM, 0, MPI_COMM_WORLD, ierr)
	call MPI_BCAST(outValue,1,MPI_DOUBLE_PRECISION,
     <       0,MPI_COMM_WORLD,ierr)
	return
	end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	subroutine sedMass2(inArray, outValue,numGhost)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Name: 
      !
      ! Author: Kurt Nelson
      !
      ! Purpose: 
      !
      ! Inputs: inArray
      !
      ! Outputs: outValue
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

	myIntegration=0.D0
	
	do i = 1, nni 
	   do j = 1, nnj 
	      do k = 1, nnk
	        myIntegration = myIntegration +inArray(i,j,k)*dyAll(j)
     <            *(xp(3,1,1,1)-xp(2,1,1,1))*(xp(1,1,3,3)-xp(1,1,2,3))
	      enddo
	    enddo
	enddo
	
	call MPI_REDUCE(myIntegration, outValue, 1, MPI_DOUBLE_PRECISION,
     <               MPI_SUM, 0, MPI_COMM_WORLD, ierr)
	call MPI_BCAST(outValue,1,MPI_DOUBLE_PRECISION,
     <       0,MPI_COMM_WORLD,ierr)
	return
	end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	subroutine get_sedTurbFlux(sedMean,vCsed)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Name: 
      !
      ! Author: Kurt Nelson
      !
      ! Purpose: 
      !
      ! Inputs: inArray
      !
      ! Outputs: outValue
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C       This subroutine routine computes velociy purturbations

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "padjust.inc"
	include "sed.inc"

	double precision, intent(in) :: sedMean(1:nnj)
	double precision, intent(out) :: vCsed(1:nnj)
	double precision sedFlux(1:nni, 1:nnj, 1:nnk)
	double precision cPrime
	integer i, j, k
	
	do i = 1, nni
	   do j = 1, nnj
	      do k = 1, nnk
	       cPrime = Csed(i,j,k)-sedMean(j)
	       sedFlux(i,j,k) = 
     <           cPrime*velPrimes(i,j,k,2)
	      enddo
	   enddo
	enddo

	call MPI_Barrier(MPI_COMM_WORLD, ierr)
	call horizontalAverage(sedFlux,vCsed,0) 
	
	vCsed = -1*vCsed
	return
	end



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	subroutine get_BruntN(den,BruntNMean,rhoMean,rhoSqrMean,rhoPrimeMean)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Name: 
      !
      ! Author: Kurt Nelson
      !
      ! Purpose: 
      !
      ! Inputs: inArray
      !
      ! Outputs: outValue
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C       This subroutine computes the local gradient Richardson Number
C       note: this is defined at cell edges

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "metric.inc"
	include "para.inc"
	
	double precision, intent(in) :: den(-1:nni+2,-1:nnj+2,-1:nnk+2)
	double precision, intent(out) :: BruntNMean(1:nnj)
	double precision, intent(out) :: rhoMean(1:nnj)
	double precision, intent(out) :: rhoSqrMean(1:nnj)
	double precision, intent(out) :: rhoPrimeMean(1:nnj)
        double precision rhoSqr(1:nni,1:nnj,1:nnk)
        double precision rhoSqrPrime(1:nni,1:nnj,1:nnk)
	integer i, j, k


	call horizontalAverage(den,rhoMean,2)
	do i = 1, nni
	   do j = 1, nnj
	      do k = 1, nnk
	       rhoSqr(i,j,k) = den(i,j,k)**2
	       rhoSqrPrime(i,j,k) = (den(i,j,k)-rhoMean(j))**2
	      enddo
	   enddo
	enddo
	call horizontalAverage(rhoSqr,rhoSqrMean,0)
	call horizontalAverage(rhoSqrPrime,rhoPrimeMean,0)


	do j = 1, nnj
	   rhoPrimeMean(j) = SQRT(rhoPrimeMean(j))
	enddo

        do j = 1, nnj
                if(j .eq. nnj) then
	  BruntNMean(j) = -g/(0.5*(rhoMean(j)+rhoMean(j-1)))
     <     *(rhoMean(j)-rhoMean(j-1))/(xp(i,j,k,2)-xp(i,j-1,k,2))
                else
	  BruntNMean(j) = -g/(0.5*(rhoMean(j)+rhoMean(j+1)))
     <     *(rhoMean(j+1)-rhoMean(j))/(xp(i,j+1,k,2)-xp(i,j,k,2))
                endif
        enddo

	call MPI_Barrier(MPI_COMM_WORLD, ierr)
	return
	end









