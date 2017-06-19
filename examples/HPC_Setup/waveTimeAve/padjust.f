cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine adjustPressure(PI)

	include "size.inc"
	include "para.inc"
	include "ns.inc"
	include "padjust.inc"
	include "mpi.inc"
	include "metric.inc"

	double precision, intent(in) :: PI
	double precision, dimension(1:nnj) :: uMean
	double precision uDepth, vertStressMean
	integer i, j, k

	if ( waves .eq. 1 ) then
	   do j = 1,nnj
	       steadyPall(j)= dpdxWave*cos(2*PI*time/Twave)
	   enddo
	 else
	   do j = 1, nnj
	         steadyPall(j) = dpdxSteady
	   enddo
	endif

	return
	end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine adjustS

	include "size.inc"
	include "para.inc"
	include "ns.inc"
	include "padjust.inc"
	include "mpi.inc"
	include "metric.inc"

	double precision, dimension(1:nnj) :: uMean
	double precision uDepth
	integer i, j, k

	call horizontalAverage(u(:,:,:,1), uMean, 2)
	call depthAverage(uMean, uDepth)
        uRun(MOD(kount-1,1500))=uDepth
        if(kount .eq. 1) then
                Sn = 1/dtime*(Ubulk-uDepth)
        elseif(kount .le. 1500) then
                Sn = kount/dtime*(Ubulk-SUM(uRun(0:kount-1))/kount)
     <          -(kount-1)* SUM(SnRun(0:kount-2))/(kount-1)     
        else
                Sn = kount/dtime*(Ubulk-SUM(uRun)/1500)
     <          -(1500-1)* SUM(SnRun)/(1500-1)
        endif     
        SnRun(MOD(kount-1,1499))=Sn
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
C	u_fric = SQRT(dpdxSteady*H/rhoWater)
C	zo = vis/(9*u_fric)

	Ubulk = 0.0956
	if (myid .eq. 0) write(*,*) "Ubulk = ", Ubulk
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
	   u(i,j,k,1) = u(i,j,k,1)+dtime*Sn
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
	   uxi(i,j,k) = uxi(i,j,k) + dtime*xix(i,j,k)*Sn

	enddo
	enddo
	enddo
	
	return
	end
