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
	       steadyPall(j) = dpdxLaminar
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

	call horizontalAverage(u(:,:,:,1), uMean, 2)
	call depthAverage(uMean, uDepth)

	if ( waves .eq. 1 ) then
	   do j = 1,nnj
	       steadyPall(j)= waveMag*cos(2*PI*time/Twave)
     <          +dpdxSteady
	   enddo
	 else
	  do i = 1,nni
	   do j = 1,nnj
	    do k = 1,nnk
	       vertStress(i,j,k)=vis*rhoWater
     <           /(0.5*(xp(i,j+1,k,2)-xp(i,j-1,k,2)))*(
     <           (u(i,j+1,k,1)-u(i,j,k,1))/(xp(i,j+1,k,2)-xp(i,j,k,2)) 
     <          -(u(i,j,k,1)-u(i,j-1,k,1))/(xp(i,j,k,2)-xp(i,j-1,k,2))) 
	    enddo
	   enddo
	  enddo

	   call horizontalAverage(vertStress, vertStressHor, 2)
	   call depthAverage(vertStressHor,vertStressMean)
	   dpdxSource = -dpdxSteady-vertStressMean

	if(mod(istep,nsave) .eq. 0 .or. istep .eq. 1) then
	if (myid .eq. 0) then
	   write(ID, fmt='(I3)') 100+vert_id
	   if (istep .gt.1) then !appends to existing file - use this if you want continue run to create new files
	     open(50+myid, file='outputpVal_Source.'//ID, form='unformatted',
     >          status='old',position='append')
	     write(50+myid) dpdxSource
	     close(unit = 50+myid)
	   else
	     open(50+myid, file='outputpVal_Source.'//ID, form='unformatted',
     >          status='unknown')
	     write(50+myid) dpdxSource
	     close(unit = 50+myid)
	   endif
	endif
	endif

	   do j = 1, nnj
C	     if (uMean(j)/uTheo(j)>1.02) then
C	         steadyPall(j) = 0
C	      else
C	         steadyPall(j) = dpdxSteady
C            endif

C This is what Ive been using
C	       steadyPall(j) = steadyPall(j)+driveFac*dtime*(Ubulk-uDepth)
	  steadyPall(j) = dpdxSteady+dpdxSource+rhoWater/dtime*(Ubulk-uDepth)
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
