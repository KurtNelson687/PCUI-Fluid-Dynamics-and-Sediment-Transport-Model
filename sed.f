cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine Csed_rhs
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Name: Csed_rhs
      !
      ! Author: Kurt Nelson
      !
      ! Purpose: This subroutine constructs the right hand side of the
      ! sediment transport equation outlined in Nelson et al (2018,
      ! JGR-Oceans).
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "ns.inc"
	include "metric.inc"
	include "para.inc"
	include "eddy.inc"
	include "sed.inc"

	double precision coef, eps
	double precision delf, Csede, Settle
	
	double precision, dimension(0:nni+1,0:nnj+1,0:nnk+1,1:3) :: 
     <		cf, Csedf
	double precision  deposition(nni,nnk), erosion(nni,nnk)
	integer i, j, k

	double precision Csedmin, Csedmax, pmin, pmax, debug

C...... Compute erosion and depostion
	call getBedFluxes(deposition, erosion)
        coef = 1.5d0

C......	Take an Euler step on the first step

        if ( istep .eq. 1 ) then
           coef = 1.d0
           do k = 0, nnk+1
           do j = 0, nnj+1
           do i = 0, nni+1
              hbCsed(i,j,k) = 0.D0
	   enddo
	   enddo 
	   enddo 
        endif
 
C......	First put in the part of Adams Bashforth from step n-2

	do k = 0, nnk+1
	do j = 0, nnj+1
	do i = 0, nni+1
           suCsed(i,j,k) = - 0.5D0 * hbCsed(i,j,k)
	enddo
	enddo 
	enddo 

C......	SHARP advection is added in all three directions independently

	eps = TINY(1.D0)

C...... Subtract out the settling volume flux from the contravariant volume flux 
	do k = 0, nnk+1
	do j = 0, nnj+1
	do i = 0, nni+1
	   uej(i,j,k)=uej(i,j,k)-ws*ety(i,j,k)
	enddo
	enddo
	enddo

C...... For I direction
	
	do k = 0, nnk
	do j = 0, nnj
	do i = 0, nni
	if ( uxi(i,j,k) .ge. 0.D0 ) then
	   delf = Csed(i+1,j,k) - Csed(i-1,j,k)  
	   if( dabs(delf) .lt. 1.D-5 ) then !Is this the correct cutoff for sediment?
	      cf(i,j,k,1) = 0.125D0 !This is just QUICK interpolation
	   else
              Csede = ( Csed(i,j,k)-Csed(i-1,j,k) ) / delf !locally normalized variable
	      if ( Csede .le. -1.D0 .or. Csede .ge. 1.5D0 ) then
	         cf(i,j,k,1) = 0.125D0
	      else
	         if ( Csede .le. 0.25D0 ) then
	            if ( Csede .le. 0.D0 ) then
	               cf(i,j,k,1) = 0.5D0 + 0.375D0 * Csede
	            else
	               cf(i,j,k,1) = 0.5D0 - 0.625D0 * dsqrt(Csede)
	            endif
	         else
	            if ( Csede .le. 1.D0 ) then
	               cf(i,j,k,1) = 0.25D0 * ( 1.D0 - Csede )
	            else
	               cf(i,j,k,1) =-0.25d0 * ( 1.D0 - Csede )
	            endif
	         endif
	      endif
    	   endif
           Csedf(i,j,k,1) = 0.5D0 * ( Csed(i,j,k) + Csed(i+1,j,k) )
     <	                 - cf(i,j,k,1) * ( Csed(i-1,j,k) 
     <	                                 - Csed(i,j,k)
     <	                                 - Csed(i,j,k)
     <	                                 + Csed(i+1,j,k) )
	else
	   delf = Csed(i,j,k) - Csed(i+2,j,k)
	   if( dabs(delf) .lt. 1.D-5 ) then
	      cf(i,j,k,1) = 0.125D0
	   else
              Csede = ( Csed(i+1,j,k)-Csed(i+2,j,k) ) / delf
	      if ( Csede .le. -1.D0 .or. Csede .ge. 1.5D0 ) then
	         cf(i,j,k,1) = 0.125D0
	      else
	         if ( Csede .le. 0.25D0 ) then
	            if ( Csede .le. 0.D0 ) then
	               cf(i,j,k,1) = 0.5D0 + 0.375D0 * Csede
	            else
	               cf(i,j,k,1) = 0.5D0 - 0.625D0 * dsqrt(Csede)
	            endif
	         else
	            if ( Csede .le. 1.D0 ) then
	               cf(i,j,k,1) = 0.25D0 * ( 1.D0 - Csede )
	            else
	               cf(i,j,k,1) =-0.25d0 * ( 1.D0 - Csede )
	            endif
	         endif
	      endif
    	   endif
           Csedf(i,j,k,1) = 0.5D0 * ( Csed(i,j,k) + Csed(i+1,j,k) )
     <	                 - cf(i,j,k,1) * ( Csed(i,j,k) 
     <	                                 - Csed(i+1,j,k)
     <	                                 - Csed(i+1,j,k)
     <	                                 + Csed(i+2,j,k) )
        endif
	enddo
	enddo 
	enddo 

C...... For J direction
	
	do k = 0, nnk
	do j = 0, nnj
	do i = 0, nni
	if ( uej(i,j,k) .ge. 0.D0 ) then
	   delf = Csed(i,j+1,k) - Csed(i,j-1,k)
	   if( dabs(delf) .lt. 1.D-5 ) then
	      cf(i,j,k,2) = 0.125D0 !This is just quick
	   else
              Csede = ( Csed(i,j,k) - Csed(i,j-1,k) ) / delf
	      if ( Csede .le. -1.D0 .or. Csede .ge. 1.5D0 ) then
	         cf(i,j,k,2) = 0.125D0
	      else
	         if ( Csede .le. 0.25D0 ) then
	            if ( Csede .le. 0.D0 ) then
	               cf(i,j,k,2) = 0.5D0 + 0.375D0 * Csede
	            else
	               cf(i,j,k,2) = 0.5D0 - 0.625D0 * dsqrt(Csede)
	            endif
	         else
	            if ( Csede .le. 1.D0 ) then
	               cf(i,j,k,2) = 0.25D0 * ( 1.D0 - Csede )
	            else
	               cf(i,j,k,2) =-0.25d0 * ( 1.D0 - Csede )
	            endif
	         endif
	      endif
	   endif
           Csedf(i,j,k,2) = 0.5D0 * ( Csed(i,j,k) + Csed(i,j+1,k) )
     <	                 - cf(i,j,k,2) * ( Csed(i,j-1,k) 
     <	                                 - Csed(i,j,k)
     <	                                 - Csed(i,j,k)
     <	                                 + Csed(i,j+1,k) )
	else
	   delf = Csed(i,j,k) - Csed(i,j+2,k)
	   if( dabs(delf) .lt. 1.D-5 ) then
	      cf(i,j,k,2) = 0.125D0
	   else
              Csede = ( Csed(i,j+1,k) - Csed(i,j+2,k) ) / delf
	      if ( Csede .le. -1.D0 .or. Csede .ge. 1.5D0 ) then
	         cf(i,j,k,2) = 0.125D0
	      else
	         if ( Csede .le. 0.25D0 ) then
	            if ( Csede .le. 0.D0 ) then
	               cf(i,j,k,2) = 0.5D0 + 0.375D0 * Csede
	            else
	               cf(i,j,k,2) = 0.5D0 - 0.625D0 * dsqrt(Csede)
	            endif
	         else
	            if ( Csede .le. 1.D0 ) then
	               cf(i,j,k,2) = 0.25D0 * ( 1.D0 - Csede )
	            else
	               cf(i,j,k,2) =-0.25d0 * ( 1.D0 - Csede )
	            endif
	         endif
	      endif
	   endif
           Csedf(i,j,k,2) = 0.5D0 * ( Csed(i,j,k) + Csed(i,j+1,k) )
     <	                 - cf(i,j,k,2) * ( Csed(i,j,k) 
     <	                                 - Csed(i,j+1,k)
     <	                                 - Csed(i,j+1,k)
     <	                                 + Csed(i,j+2,k) )
	endif
	enddo
	enddo 
	enddo 

C...... For K direction
	
	do k = 0, nnk
	do j = 0, nnj
	do i = 0, nni
	if ( uzk(i,j,k) .ge. 0.D0 ) then
	   delf = Csed(i,j,k+1) - Csed(i,j,k-1)
	   if( dabs(delf) .lt. 1.D-5 ) then
	      cf(i,j,k,3) = 0.125D0
	   else
              Csede = ( Csed(i,j,k) - Csed(i,j,k-1) ) / delf
	      if ( Csede .le. -1.D0 .or. Csede .ge. 1.5D0 ) then
	         cf(i,j,k,3) = 0.125D0
	      else
	         if ( Csede .le. 0.25D0 ) then
	            if ( Csede .le. 0.D0 ) then
	               cf(i,j,k,3) = 0.5D0 + 0.375D0 * Csede
	            else
	               cf(i,j,k,3) = 0.5D0 - 0.625D0 * dsqrt(Csede)
	            endif
	         else
	            if ( Csede .le. 1.D0 ) then
	               cf(i,j,k,3) = 0.25D0 * ( 1.D0 - Csede )
	            else
	               cf(i,j,k,3) =-0.25d0 * ( 1.D0 - Csede )
	            endif
	         endif
	      endif
	   endif
           Csedf(i,j,k,3) = 0.5D0 * ( Csed(i,j,k) + Csed(i,j,k+1) )
     <	                 - cf(i,j,k,3) * ( Csed(i,j,k-1) 
     <	                                 - Csed(i,j,k)
     <	                                 - Csed(i,j,k)
     <	                                 + Csed(i,j,k+1) )
	else
	   delf = Csed(i,j,k) - Csed(i,j,k+2)
	   if( dabs(delf) .lt. 1.D-5 ) then
	      cf(i,j,k,3) = 0.125D0
	   else
              Csede = ( Csed(i,j,k+1) - Csed(i,j,k+2) ) / delf
	      if ( Csede .le. -1.D0 .or. Csede .ge. 1.5D0 ) then
	         cf(i,j,k,3) = 0.125D0
	      else
	         if ( Csede .le. 0.25D0 ) then
	            if ( Csede .le. 0.D0 ) then
	               cf(i,j,k,3) = 0.5D0 + 0.375D0 * Csede
	            else
	               cf(i,j,k,3) = 0.5D0 - 0.625D0 * dsqrt(Csede)
	            endif
	         else
	            if ( Csede .le. 1.D0 ) then
	               cf(i,j,k,3) = 0.25D0 * ( 1.D0 - Csede )
	            else
	               cf(i,j,k,3) =-0.25d0 * ( 1.D0 - Csede )
	            endif
	         endif
	      endif
	   endif
           Csedf(i,j,k,3) = 0.5D0 * ( Csed(i,j,k) + Csed(i,j,k+1) )
     <	                 - cf(i,j,k,3) * ( Csed(i,j,k) 
     <	                                 - Csed(i,j,k+1)
     <	                                 - Csed(i,j,k+1)
     <	                                 + Csed(i,j,k+2) )
	endif
	enddo
	enddo 
	enddo 


C......	Convective terms (explicit)

	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni
	   hbCsed(i,j,k) = - uxi(i,  j,k) * Csedf(i,  j,k,1) 
     <	   + uxi(i-1,j,k) * Csedf(i-1,j,k,1)
     <	   - uej(i,j  ,k) * Csedf(i,j,  k,2)
     <	   + uej(i,j-1,k) * Csedf(i,j-1,k,2)
     <	   - uzk(i,j,k  ) * Csedf(i,j,k,  3)
     <	   + uzk(i,j,k-1) * Csedf(i,j,k-1,3)

	   hbCsed(I,J,K) = hbCsed(I,J,K) + Csed(I,J,K) *
     <	      (UXI(I,J,K) - UXI(I-1,J,K)+
     <        UEJ(I,J,K)- UEJ(I,J-1,K)+
     <	      UZK(I,J,K) - UZK(I,J,K-1))
	enddo
	enddo
	enddo

C......	LES self-similarity term
	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni
	   hbCsed(i,j,k) = hbCsed(i,j,k) + rr(i,j,k,4) !I need to change rr for bottom cell
	enddo
	enddo
	enddo

C......	Cross viscous terms at step n-1 from Crank-Nicolson  
!!!!!!!!!!!!! Note:  akst for sediment model should be checked !!!!!!!!!!!!!!!!!
	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni

	   hbCsed(i,j,k) = hbCsed(i,j,k) 
     <        + ( ak + 0.5D0*(akst(i,j,k) + akst(i+1,j,k)) ) *
     <		( g12(i,  j,k) * ( Csed(i,  j+1,k) - Csed(i,  j-1,k) 
     <		                 + Csed(i+1,j+1,k) - Csed(i+1,j-1,k) )
     <		+ g13(i,  j,k) * ( Csed(i,  j,k+1) - Csed(i,  j,k-1)
     <		                 + Csed(i+1,j,k+1) - Csed(i+1,j,k-1) ) ) 
     <        - ( ak + 0.5D0*(akst(i,j,k) + akst(i-1,j,k)) ) *
     <		( g12(i-1,j,k) * ( Csed(i,  j+1,k) - Csed(i,  j-1,k)
     <		                 + Csed(i-1,j+1,k) - Csed(i-1,j-1,k) )
     <		+ g13(i-1,j,k) * ( Csed(i,  j,k+1) - Csed(i,  j,k-1)
     <		                 + Csed(i-1,j,k+1) - Csed(i-1,j,k-1) ) )

	   hbCsed(i,j,k) = hbCsed(i,j,k) 
     <        + ( ak + 0.5D0*(akst(i,j,k) + akst(i,j+1,k)) ) *
     <		( g23(i,j,  k) * ( Csed(i,j,  k+1) - Csed(i,j,  k-1)
     <		                 + Csed(i,j+1,k+1) - Csed(i,j+1,k-1) )
     <		+ g21(i,j,  k) * ( Csed(i+1,j,  k) - Csed(i-1,j,  k)
     <		                 + Csed(i+1,j+1,k) - Csed(i-1,j+1,k) ) ) 
     <        - ( ak + 0.5D0*(akst(i,j,k) + akst(i,j-1,k)) ) *
     <		( g23(i,j-1,k) * ( Csed(i,j,  k+1) - Csed(i,j,  k-1)
     <		                 + Csed(i,j-1,k+1) - Csed(i,j-1,k-1) ) 
     <		+ g21(i,j-1,k) * ( Csed(i+1,j,  k) - Csed(i-1,j,  k)
     <		                 + Csed(i+1,j-1,k) - Csed(i-1,j-1,k) ) )

	   hbCsed(i,j,k) = hbCsed(i,j,k)  
     <        + ( ak + 0.5D0*(akst(i,j,k) + akst(i,j,k+1)) ) *
     <		( g31(i,j,k  ) * ( Csed(i+1,j,k  ) - Csed(i-1,j,k  )
     <		                 + Csed(i+1,j,k+1) - Csed(i-1,j,k+1) )
     <		+ g32(i,j,k  ) * ( Csed(i,j+1,k  ) - Csed(i,j-1,k  )
     <		                 + Csed(i,j+1,k+1) - Csed(i,j-1,k+1) ) ) 
     <        - ( ak + 0.5D0*(akst(i,j,k) + akst(i,j,k-1)) ) *
     <		( g31(i,j,k-1) * ( Csed(i+1,j,k  ) - Csed(i-1,j,k  )
     <		                 + Csed(i+1,j,k-1) - Csed(i-1,j,k-1) )
     <		+ g32(i,j,k-1) * ( Csed(i,j+1,k  ) - Csed(i,j-1,k  )
     <		                 + Csed(i,j+1,k-1) - Csed(i,j-1,k-1) ) )

	enddo
	enddo
	enddo

C...... bottom wall BC
	if ( n_suth .eq. MPI_PROC_NULL ) then
	   do i = 1, nni
	   do k = 1, nnk

C        Add back diffusion, and advection from bottom faces for bottom cells
C        and add erosion and deposition

	   hbCsed(i,1,k) = hbCsed(i,1,k)
C        Diffussion terms
     <		 + ( ak + 0.5D0*(akst(i,1,k) + akst(i,0,k)) ) *
     <		( g23(i,0,k) * ( Csed(i,1,  k+1) - Csed(i,1,  k-1)
     <		                 + Csed(i,0,k+1) - Csed(i,0,k-1) ) 
     <		+ g21(i,0,k) * ( Csed(i+1,1,  k) - Csed(i-1,1,  k)
     <		                 + Csed(i+1,0,k) - Csed(i-1,0,k) ) )
     <		- uej(i,0,k) * Csedf(i,0,k,2)
     <		-deposition(i,k)
     <		+erosion(i,k)

	   enddo
	   enddo
	endif

C...... top  wall BC
	
	if ( n_nrth .eq. MPI_PROC_NULL ) then
	   do i = 1, nni
	   do k = 1, nnk
	   hbCsed(i,nnj,k) = hbCsed(i,nnj,k)
C        Diffussion terms
	   
     <        - ( ak + 0.5D0*(akst(i,nnj,k) + akst(i,nnj+1,k)) ) *
     <		( g23(i,nnj,  k) * ( Csed(i,nnj,  k+1) - Csed(i,nnj,  k-1)
     <		                 + Csed(i,nnj+1,k+1) - Csed(i,nnj+1,k-1) )
     <		+ g21(i,nnj,  k) * ( Csed(i+1,nnj,  k) - Csed(i-1,nnj,  k)
     <		                 + Csed(i+1,nnj+1,k) - Csed(i-1,nnj+1,k) ) ) 
C	Advection term
     <		+ uej(i,nnj  ,k) * Csedf(i,nnj,  k,2)
	   enddo
	   enddo
	endif


C...... Add back the settling velocity flux to the contravariant volume flux
	do k = 0, nnk+1
	do j = 0, nnj+1
	do i = 0, nni+1
	   uej(i,j,k)=uej(i,j,k)+ws*ety(i,j,k)

	enddo
	enddo
	enddo

C......	Diagonal viscous terms at step n-1 from Crank-Nicolson 

	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni
	   suCsed(i,j,k) = suCsed(i,j,k) + coef * hbCsed(i,j,k)
	enddo
	enddo 
	enddo 
	
	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni
	   suCsed(i,j,k) = suCsed(i,j,k) 
     <        + ( ak + 0.5D0*(akst(i,j,k) + akst(i+1,j,k)) ) *
     <		g11(i,  j,  k  ) * ( Csed(i+1,j,  k  ) - Csed(i,j,k) )  
     <        + ( ak + 0.5D0*(akst(i,j,k) + akst(i-1,j,k)) ) *
     <		g11(i-1,j,  k  ) * ( Csed(i-1,j,  k  ) - Csed(i,j,k) )   
     <        + ( ak + 0.5D0*(akst(i,j,k) + akst(i,j+1,k)) ) *
     <		g22(i,  j,  k  ) * ( Csed(i,  j+1,k  ) - Csed(i,j,k) )   
     <        + ( ak + 0.5D0*(akst(i,j,k) + akst(i,j-1,k)) ) *
     <		g22(i,  j-1,k  ) * ( Csed(i,  j-1,k  ) - Csed(i,j,k) )  
     <        + ( ak + 0.5D0*(akst(i,j,k) + akst(i,j,k+1)) ) *
     <		g33(i,  j,  k  ) * ( Csed(i,  j,  k+1) - Csed(i,j,k) )  
     <        + ( ak + 0.5D0*(akst(i,j,k) + akst(i,j,k-1)) ) *
     <		g33(i,  j,  k-1) * ( Csed(i,  j,  k-1) - Csed(i,j,k) ) 
	enddo
	enddo 
	enddo

C...... add difussion back at bottom wall
	if ( n_suth .eq. MPI_PROC_NULL ) then
	   do i = 1, nni
	   do k = 1, nnk
	   suCsed(i,1,k) = suCsed(i,1,k) 
     <        - ( ak + 0.5D0*(akst(i,0,k) + akst(i,1,k)) ) *
     <		g22(i,  0,k  ) * ( Csed(i,0,k  ) - Csed(i,1,k) )  
	   enddo
	   enddo
	endif


C...... add difussion back at top wall
	if ( n_nrth .eq. MPI_PROC_NULL ) then
	   do i = 1, nni
	   do k = 1, nnk
	   suCsed(i,nnj,k) = suCsed(i,nnj,k) 
     <        - ( ak + 0.5D0*(akst(i,nnj,k) + akst(i,nnj+1,k)) ) *
     <		g22(i,  nnj,k  ) * ( Csed(i,nnj+1,k  ) - Csed(i,nnj,k) )  
	   enddo
	   enddo
	endif

	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni
	   suCsed(i,j,k) = dtime * jac(i,j,k) * suCsed(i,j,k)
	enddo
	enddo 
	enddo 

	return
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine Csed_solve
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Name: Csed_solve
      !
      ! Author: Kurt Nelson
      !
      ! Purpose: This subroutine uses approximate factorization and 
      ! the Thomas algorithm to solve the sediment concentration at step
      ! n+1  
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "ns.inc"
	include "metric.inc"
	include "para.inc"
	include "eddy.inc"
	include "sed.inc"

	double precision, dimension(nnj,0:nni+1) ::
     <		ax, bx, cx, fx
	double precision, dimension(nni,0:nnj+1) ::
     <		ay, by, cy, fy
	double precision, dimension(nni,0:nnk+1) ::
     <		az, bz, cz, fz

	integer i, j, k
	double precision coef
	
	coef = 0.5D0 * dtime

C...... solve for I-direction

	do k = 1, nnk

	do j = 1, nnj
	do i = 1, nni
	   ax(j,i) = -( ak + 0.5D0*(akst(i,j,k) + akst(i-1,j,k)) ) * 
     <               coef * jac(i,j,k) * g11(i-1,j,k)
	   cx(j,i) = -( ak + 0.5D0*(akst(i,j,k) + akst(i+1,j,k)) ) * 
     <               coef * jac(i,j,k) * g11(i,  j,k)
	   bx(j,i) = 1.D0 - ax(j,i) - cx(j,i)
	   fx(j,i) = suCsed(i,j,k)
	enddo
	enddo

	if ( n_west .eq. MPI_PROC_NULL ) then
	   do j = 1, nnj
	      hbCsed(0,j,k)=( g12(0,j,k) * ( Csed(0,j+1,k) - Csed(0,j-1,k)
     <	                              + Csed(1,j+1,k) - Csed(1,j-1,k) )
     <	                 + g13(0,j,k) * ( Csed(0,j,k+1) - Csed(0,j,k-1)
     <	                              + Csed(1,j,k+1) - Csed(1,j,k-1) ))
     <	                 / g11(0,j,k)
	      ax(j,0) =  0.D0
	      bx(j,0) =  1.D0
	      cx(j,0) = -1.D0
	      fx(j,0) =  hbCsed(0,j,k)
	   enddo
	endif
	
	if ( n_east .eq. MPI_PROC_NULL ) then
	   do j = 1, nnj
	      hbCsed(nni+1,j,k)=( g12(nni,j,k) 
     <	                     * ( Csed(nni,  j+1,k) - Csed(nni,  j-1,k)
     <	                       + Csed(nni+1,j+1,k) - Csed(nni+1,j-1,k) )
     <	                     + g13(nni,j,k) 
     <	                     * ( Csed(nni,  j,k+1) - Csed(nni,  j,k-1)
     <	                       + Csed(nni+1,j,k+1) - Csed(nni+1,j,k-1) ))
     <	                     / g11(nni,j,k)
	      ax(j,nni+1) =  1.D0
	      bx(j,nni+1) = -1.D0
	      cx(j,nni+1) =  0.D0
	      fx(j,nni+1) =  hbCsed(nni+1,j,k)
	   enddo
	endif

	if ( periodic .eq. 1 ) then
	call trip( ax, bx, cx, fx, nnj, 1, nni, n_west, n_east )
	else
	call trid( ax, bx, cx, fx, nnj, 1, nni, n_west, n_east )
	endif


!	call trid( ax, bx, cx, fx, nnj, 1, nni, n_west, n_east ) !Commented by Kurt

	do j = 1, nnj
	do i = 1, nni
	   suCsed(i,j,k) = fx(j,i)
	enddo
	enddo

	enddo

C...... solve for J-direction

	do k = 1, nnk

	do j = 1, nnj
	do i = 1, nni
	   ay(i,j) = -( ak + 0.5D0*(akst(i,j,k) + akst(i,j-1,k)) ) * 
     <               coef * jac(i,j,k) * g22(i,j-1,k)
	   cy(i,j) = -( ak + 0.5D0*(akst(i,j,k) + akst(i,j+1,k)) ) * 
     <               coef * jac(i,j,k) * g22(i,j,  k)
	   by(i,j) = 1.D0 - ay(i,j) - cy(i,j)
	   fy(i,j) = suCsed(i,j,k)
	enddo
	enddo

	if ( n_suth .eq. MPI_PROC_NULL ) then
	   do i = 1, nni
C	      hbCsed(i,0,k)=( g23(i,0,k) * ( Csed(i,0,k+1) - Csed(i,0,k-1)
C     <	                              + Csed(i,1,k+1) - Csed(i,1,k-1) )
C     <	               + g21(i,0,k) * ( Csed(i+1,0,k) - Csed(i-1,0,k)
C     <	                              + Csed(i+1,1,k) - Csed(i-1,1,k) ))
C     <	               / g22(i,0,k)
C	      fy(i,0) =  hbCsed(i,0,k)
	      ay(i,0) =  0.D0
	      by(i,0) =  1.D0
	      cy(i,0) = -1.D0
	      fy(i,0) =  0.D0

C	      ay(i,1) =  0.D0
C	      cy(i,1) = -( ak + 0.5D0*(akst(i,1,k) + akst(i,2,k)) ) * 
C     <               coef * jac(i,1,k) * g22(i,1,k)
C	      by(i,1) = 1.D0 - cy(i,1)

	   enddo
	endif
	
	if ( n_nrth .eq. MPI_PROC_NULL ) then
	   do i = 1, nni
C	      hbCsed(i,nnj+1,k)=( g23(i,nnj,k) 
C     <	                   * ( Csed(i,nnj,  k+1) - Csed(i,nnj,  k-1)
C     <	                     + Csed(i,nnj+1,k+1) - Csed(i,nnj+1,k-1) )
C     <                    + g21(i,nnj,k) 
C     <	                   * ( Csed(i+1,nnj,  k) - Csed(i-1,nnj,  k)
C     <	                     + Csed(i+1,nnj+1,k) - Csed(i-1,nnj+1,k) ) )
C     <	                  / g22(i,nnj,k)
C	      fy(i,nnj+1) =  hbCsed(i,nnj+1,k)
	      ay(i,nnj+1) =  1.D0
	      by(i,nnj+1) = -1.D0
	      cy(i,nnj+1) =  0.D0
	      fy(i,nnj+1) =  0.D0

C	      ay(i,nnj) = -( ak + 0.5D0*(akst(i,nnj,k) + akst(i,nnj-1,k)) ) * 
C     <               coef * jac(i,nnj,k) * g22(i,nnj-1,k)
C	      cy(i,nnj) =  0.D0
C	      by(i,nnj) = 1.D0 - ay(i,nnj)

	   enddo
	endif
	
	call trid( ay, by, cy, fy, nni, 1, nnj, n_suth, n_nrth )

	do j = 1, nnj
	do i = 1, nni
	   suCsed(i,j,k) = fy(i,j)
	enddo
	enddo

	enddo

C...... solve for K-direction

	do j = 1, nnj

	do k = 1, nnk
	do i = 1, nni
	   az(i,k) = -( ak + 0.5D0*(akst(i,j,k) + akst(i,j,k-1)) ) * 
     <               coef * jac(i,j,k) * g33(i,j,k-1)
	   cz(i,k) = -( ak + 0.5D0*(akst(i,j,k) + akst(i,j,k+1)) ) * 
     <               coef * jac(i,j,k) * g33(i,j,k  )
	   bz(i,k) = 1.D0 - az(i,k) - cz(i,k)
	   fz(i,k) = suCsed(i,j,k)
	enddo
	enddo

	if ( n_back .eq. MPI_PROC_NULL ) then
	   do i = 1, nni
	      hbCsed(i,j,0)=( g31(i,j,0) * ( Csed(i+1,j,0) - Csed(i-1,j,0)
     <	                              + Csed(i+1,j,1) - Csed(i-1,j,1) )
     <	                 + g32(i,j,0) * ( Csed(i,j+1,0) - Csed(i,j-1,0)
     <	                              + Csed(i,j+1,1) - Csed(i,j-1,1) ))
     <	                 / g33(i,j,0)
	      az(i,0) =  0.D0
	      bz(i,0) =  1.D0
	      cz(i,0) = -1.D0
	      fz(i,0) =  hbCsed(i,j,0)
	   enddo
	endif
	
	if ( n_frnt .eq. MPI_PROC_NULL ) then
	   do i = 1, nni
	      hbCsed(i,j,nnk+1)=( g31(i,j,nnk) 
     <	                     * ( Csed(i+1,j,nnk  ) - Csed(i-1,j,nnk  )
     <	                      + Csed(i+1,j,nnk+1) - Csed(i-1,j,nnk+1) )
     <	                     + g32(i,j,nnk) 
     <	                     * ( Csed(i,j+1,nnk  ) - Csed(i,j-1,nnk  )
     <	                     + Csed(i,j+1,nnk+1) - Csed(i,j-1,nnk+1) ) )
     <	                     / g33(i,j,nnk)
	      az(i,nnk+1) =  1.D0
	      bz(i,nnk+1) = -1.D0
	      cz(i,nnk+1) =  0.D0
	      fz(i,nnk+1) =  hbCsed(i,j,nnk+1)
	   enddo
	endif
	
	if ( periodic .eq. 1 ) then
	call trip( az, bz, cz, fz, nni, 1, nnk, n_back, n_frnt )
	else
	call trid( az, bz, cz, fz, nni, 1, nnk, n_back, n_frnt )
	endif

	do k = 1, nnk
	do i = 1, nni
	   suCsed(i,j,k) = fz(i,k)
	enddo
	enddo

	enddo
	
C...... update sediment concentration (Csed)

	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni
	   Csed(i,j,k) = Csed(i,j,k) + suCsed(i,j,k)
	enddo
	enddo
	enddo

C...... Exchange information for ghost cells
	call Csed_bc
	call Csed_exchange

C...... Calculate energy
C      call energy

	return
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine Csed_bc
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Name: Csed_bc
      !
      ! Author: Kurt Nelson
      !
      ! Purpose: This subroutine simply fills ghost cells for outer 
      ! processors. It uses linear interpolation at each wall.
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "ns.inc"
	include "sed.inc"
	
	integer i, j, k

	if ( n_west .eq. MPI_PROC_NULL ) then
	   do k = -1, nnk+2
	   do j = -1, nnj+2
	      Csed( 0,j,k) = Csed(1,j,k) + hbCsed(0,j,k)
	   enddo
	   enddo
	   do k = -1, nnk+2
	   do j = -1, nnj+2
	      Csed(-1,j,k) = 3.D0 * ( Csed(0,j,k) - Csed(1,j,k) )
     <	                  + Csed(2,j,k)
	   enddo
	   enddo
	endif

	if ( n_east .eq. MPI_PROC_NULL ) then
	   do k = -1, nnk+2
	   do j = -1, nnj+2
	      Csed(nni+1,j,k) = Csed(nni,j,k) - hbCsed(nni+1,j,k)
	   enddo
	   enddo
	   do k = -1, nnk+2
	   do j = -1, nnj+2
	      Csed(nni+2,j,k) = 3.D0 * ( Csed(nni+1,j,k) - Csed(nni,j,k) )
     <	                     + Csed(nni-1,j,k)
	   enddo
	   enddo
	endif

	if ( n_suth .eq. MPI_PROC_NULL ) then
	   do k = -1, nnk+2
	   do i = -1, nni+2
	      Csed(i, 0,k) = Csed(i,1,k)
	   enddo
	   enddo
	   do k = -1, nnk+2
	   do i = -1, nni+2
	      Csed(i,-1,k) = Csed(i,1,k) 
	   enddo
	   enddo
	endif

	if ( n_nrth .eq. MPI_PROC_NULL ) then
	   do k = -1, nnk+2
	   do i = -1, nni+2
	      Csed(i,nnj+1,k) = Csed(i,nnj,k)
	   enddo
	   enddo
	   do k = -1, nnk+2
	   do i = -1, nni+2
	      Csed(i,nnj+2,k) = Csed(i,nnj,k)
	   enddo
	   enddo
	endif

	if ( n_back .eq. MPI_PROC_NULL ) then
	   do j = -1, nnj+2
	   do i = -1, nni+2
	      Csed(i,j, 0) = Csed(i,j,1) + hbCsed(i,j,0)
	   enddo
	   enddo
	   do j = -1, nnj+2
	   do i = -1, nni+2
	      Csed(i,j,-1) = 3.D0 * ( Csed(i,j,0) - Csed(i,j,1) )
     <	                  + Csed(i,j,2)
	   enddo
	   enddo
	endif

	if ( n_frnt .eq. MPI_PROC_NULL ) then
	   do j = -1, nnj+2
	   do i = -1, nni+2
	      Csed(i,j,nnk+1) = Csed(i,j,nnk) - hbCsed(i,j,nnk+1)
	   enddo
	   enddo
	   do j = -1, nnj+2
	   do i = -1, nni+2
	      Csed(i,j,nnk+2) = 3.D0 * ( Csed(i,j,nnk+1) - Csed(i,j,nnk) )
     <	                     + Csed(i,j,nnk-1)
	   enddo
	   enddo
	endif

	if ( n_west .eq. MPI_PROC_NULL ) then
	   if ( n_suth .eq. MPI_PROC_NULL ) then
	      do k = 1, nnk
C******	0 - 0 - *
   	Csed(0,0,k) = 0.5D0 * ( Csed(1,0,k) + Csed(0,1,k) )
	      enddo
	   endif
	   if ( n_nrth .eq. MPI_PROC_NULL ) then
	      do k = 1, nnk
C******	0 - Y - *
 	Csed(0,nnj+1,k) = 0.5D0 * ( Csed(1,nnj+1,k) + Csed(0,nnj,k) )
	      enddo
	   endif
	   if ( n_back .eq. MPI_PROC_NULL ) then
	      do j = 1, nnj
C******	0 - * - 0
   	Csed(0,j,0) = 0.5D0 * ( Csed(1,j,0) + Csed(0,j,1) )
	      enddo
	   endif
	   if ( n_frnt .eq. MPI_PROC_NULL ) then
	      do j = 1, nnj
C******	0 - * - Z
   	Csed(0,j,nnk+1) = 0.5D0 * ( Csed(1,j,nnk+1) + Csed(0,j,nnk) )
	      enddo
	   endif
	endif	      

	if ( n_east .eq. MPI_PROC_NULL ) then
	   if ( n_suth .eq. MPI_PROC_NULL ) then
	      do k = 1, nnk
C******	X - 0 - *
   	Csed(nni+1,0,k) = 0.5D0 * ( Csed(nni,0,k) + Csed(nni+1,1,k) )
	      enddo
	   endif
	   if ( n_nrth .eq. MPI_PROC_NULL ) then
	      do k = 1, nnk
C******	X - Y - *
  	Csed(nni+1,nnj+1,k) = 0.5D0 * ( Csed(nni,nnj+1,k) 
     <	                             + Csed(nni+1,nnj,k) )
	      enddo
	   endif
	   if ( n_back .eq. MPI_PROC_NULL ) then
	      do j = 1, nnj
C******	X - * - 0
   	Csed(nni+1,j,0) = 0.5D0 * ( Csed(nni,j,0) + Csed(nni+1,j,1) )
	      enddo
	   endif
	   if ( n_frnt .eq. MPI_PROC_NULL ) then
	      do j = 1, nnj
C******	X - * - Z
  	Csed(nni+1,j,nnk+1) = 0.5D0 * ( Csed(nni,j,nnk+1) 
     <	                             + Csed(nni+1,j,nnk) )
	      enddo
	   endif
	endif	
	
	if ( n_suth .eq. MPI_PROC_NULL ) then
	   if ( n_back .eq. MPI_PROC_NULL ) then
	      do i = 1, nni
C******	* - 0 - 0
  	Csed(i,0,0) = 0.5D0 * ( Csed(i,1,0) + Csed(i,0,1) )
	      enddo
	   endif
	   if ( n_frnt .eq. MPI_PROC_NULL ) then
	      do i = 1, nni
C******	* - 0 - Z
  	Csed(i,0,nnk+1) = 0.5D0 * ( Csed(i,1,nnk+1) + Csed(i,0,nnk) )
	      enddo
	   endif
	endif      

	if ( n_nrth .eq. MPI_PROC_NULL ) then
	   if ( n_back .eq. MPI_PROC_NULL ) then
	      do i = 1, nni
C******	* - Y - 0
  	Csed(i,nnj+1,0) = 0.5D0 * ( Csed(i,nnj,0) + Csed(i,nnj+1,1) )
	      enddo
	   endif
	   if ( n_frnt .eq. MPI_PROC_NULL ) then
	      do i = 1, nni
C******	* - Y - Z
  	Csed(i,nnj+1,nnk+1) = 0.5D0 * ( Csed(i,nnj,nnk+1)
     <	                             + Csed(i,nnj+1,nnk) )
	      enddo
	   endif
	endif  

	return
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine getBedFluxes(deposition, erosion)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Name: getBedFluxes
      !
      ! Author: Kurt Nelson
      !
      ! Purpose: This subroutine computes erosion and depostion at the
      ! bottom wall.
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "ns.inc"
	include "metric.inc"
	include "para.inc"
	include "sed.inc"

	double precision, intent(out) :: deposition(nni,nnk),
     <          erosion(nni,nnk)
	double precision bedShear
	integer i,k

	   do i = 1, nni
	   do k = 1, nnk

C        Compute deposition
	   deposition(i,k) = ws*ety(i,0,k)*(Csed(i,1,k)+(0-xp(i,1,k,2))
     <       /(xp(i,2,k,2)-xp(i,1,k,2))*(Csed(i,2,k)-Csed(i,1,k)))
C	Compute erosion
	   bedShear = rhoWater*vis
     <         *DSQRT(u(i,1,k,1)**2+u(i,1,k,3)**2)/xp(i,1,k,2)
	   if(bedShear .ge. tauCrit) then
C       note: this is based on Jones and Jaffe (2013). The 0.01 is because E is in cm/s
C	     erosion(i,k) = 0.01*ety(i,0,k)*dryBulk*Ased*bedShear**nsed
C      note: This is based on Brand (2015) linear model. Note the M they report must be converted to kg
	     erosion(i,k) = ety(i,0,k)*Ased*(bedShear-tauCrit)
	   else
	     erosion(i,k) = 0
	   endif  

	   enddo
	   enddo

	return
	end
	
