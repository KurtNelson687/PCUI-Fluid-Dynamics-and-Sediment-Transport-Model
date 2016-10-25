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
	integer j, jj0

	
	jj0 = npy * nnj
	
	H = yAll(nj)+(yAll(nj)-yAll(nj-1))/2
	do j = 1, nnj

C	    integrationArray(j) = 1/H*inArray(j)*dyAll(jj0+j)
	    integrationArray(j) = 1/H*inArray(j)*0.5*
     <         (xp(1,j+1,1,2)-xp(1,j-1,1,2))
	enddo
	
	mySum = sum(integrationArray)

	call MPI_REDUCE(mySum, outValue, 1, MPI_DOUBLE_PRECISION,
     <               MPI_SUM, 0, vert_comm, ierr)
	call MPI_BCAST(outValue,1,MPI_DOUBLE_PRECISION,0,vert_comm,ierr)
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
	double precision, intent(out) :: 
     <      velPrimes(0:nni+1, 0:nnj+1, 0:nnk+1,3)
	
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

	subroutine get_turbIntensity(uTurb, vTurb, wTurb)
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
	call get_pertVel(u, uMean,vMean,wMean, velPrimes)
	call depthAverage(uMean(1:nnj), uDepth)

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

C       This subroutine Saves Cn and puts it into the Production matrix for computing global production

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "ns.inc"
	include "metric.inc"
	include "padjust.inc"
	include "para.inc"
	include "eddy.inc"

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

