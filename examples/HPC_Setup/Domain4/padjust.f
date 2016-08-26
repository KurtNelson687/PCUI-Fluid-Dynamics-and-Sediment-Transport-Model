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

	double precision uMean(1:nnj)
	double precision uDepth
	integer j
	real, parameter :: PI=3.1415926535897932

	call horizontalAverage(u(:,:,:,1), uMean, 2)
	call depthAverage(uMean, uDepth)

	if ( waves .eq. 1 ) then
	   do j = 1,nnj
	       steadyPall(j)= waveMag*cos(2*PI*time/Twave)
     <          +dpdxSteady
	   enddo
	   else
	      do j = 1, nnj
C	         if (uMean(j)/uTheo(j)>1.02) then
C	            steadyPall(j) = 0
C	         else
C	             steadyPall(j) = dpdxSteady
C	         endif

	       steadyPall(j) = steadyPall(j)+driveFac*dtime*(Ubulk-uDepth)
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

	Ubulk = u_fric/karman*(log(H/zo)+zo/H-1)
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
