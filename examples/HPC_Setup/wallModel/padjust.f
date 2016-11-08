cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine init_pSteady

	include "size.inc"
	include "para.inc"
	include "padjust.inc"

	integer j
	
	if ( waves .eq. 1 ) then
	   do j = 1, nnj
	       steadyPall(j) = waveMag
	   enddo
	else
	   do j = 1, nnj
	       steadyPall(j) = dpdxSteady
	   enddo
	endif
	return
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine adjustPressure

	include "size.inc"
	include "para.inc"
	include "ns.inc"
	include "padjust.inc"
	include "mpi.inc"
	include "metric.inc"

	double precision, dimension(1:nnj) :: uMean,vertStressHor
	double precision uDepth, vertStressMean, dpdxSource
	double precision vertStress(-1:nni+2,-1:nnj+2,-1:nnk+2)
	integer i, j, k
	real, parameter :: PI=3.1415926535897932
	character*4 :: ID

C	call horizontalAverage(u(:,:,:,1), uMean, 2)
C	call depthAverage(uMean, uDepth)
	
	call volumeAve(u(:,:,:,1), uDepth,2)
	
	if ( waves .eq. 1 ) then
	   do j = 1,nnj
	       steadyPall(j)= waveMag*cos(2*PI*time/Twave)
     <          +dpdxSteady
	   enddo
	 else
	   do j = 1, nnj
	         steadyPall(j) = 1/dtime*(Ubulk-uDepth)
	   enddo
	endif

	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine getUtheo
	
	include "size.inc"
	include "padjust.inc"
	include "cavity.inc"
	include "mpi.inc"
	include "para.inc"
	integer j, jj0
	double precision u_fric, karman, zo, H
	
	H = yAll(nj)+(yAll(nj)-yAll(nj-1))*0.5
	karman = 0.41
	u_fric = SQRT(dpdxSteady*H/rhoWater)
	zo = vis/(9*u_fric)

C	Ubulk = u_fric/karman*(log(H/zo)+zo/H-1)
	Ubulk = uTarget
	if (myid .eq. 0) write(*,*) "Ubulk = ", Ubulk
	jj0 = npy * nnj

	do j = 1, nnj
	   if (yAll(jj0+j)<11*vis/u_fric) then
	       uTheo(j) = 99999
	   else
	       uTheo(j) = u_fric/karman*LOG(yAll(jj0+j)/zo) 
	   endif
	enddo
	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine adjust_u

	include "size.inc"
	include "mpif.h"
	include "para.inc"
	include "ns.inc"
	include "padjust.inc"
	include "mpi.inc"
	include "metric.inc"

	integer i, j, k



	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni
	   u(i,j,k,1) = u(i,j,k,1)+dtime*steadyPall(j)
	enddo
	enddo
	enddo


	return
	end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine adjust_uxi

	include "size.inc"
	include "mpif.h"
	include "para.inc"
	include "ns.inc"
	include "padjust.inc"
	include "mpi.inc"
	include "metric.inc"

	integer i, j, k

	do k = 1, nnk
	do j = 1, nnj
	do i = ius, iue
	   uxi(i,j,k) = uxi(i,j,k) + dtime*xix(i,j,k)*steadyPall(j)

	enddo
	enddo
	enddo
	
	return
	end