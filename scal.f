cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine scalar_rhs

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "ns.inc"
	include "metric.inc"
	include "para.inc"
	include "eddy.inc"

	double precision coef, eps
	double precision delf, phie
	
	double precision, dimension(0:nni+1,0:nnj+1,0:nnk+1,1:3) :: 
     <		cf, phif

	integer i, j, k

	double precision phimin, phimax, pmin, pmax
	
        coef = 1.5d0

C......	Take an Euler step on the first step

        if ( istep .eq. 1 ) then
           coef = 1.d0
           do k = 0, nnk+1
           do j = 0, nnj+1
           do i = 0, nni+1
              hbs(i,j,k) = 0.D0
	   enddo
	   enddo 
	   enddo 
        endif

C......	First put in the part of Adams Bashforth from step n-2

	do k = 0, nnk+1
	do j = 0, nnj+1
	do i = 0, nni+1
           sus(i,j,k) = - 0.5D0 * hbs(i,j,k)
	enddo
	enddo 
	enddo 

C......	SHARP formula

	eps = TINY(1.D0)

C...... I direction
	
	do k = 0, nnk
	do j = 0, nnj
	do i = 0, nni
	if ( uxi(i,j,k) .ge. 0.D0 ) then
	   delf = phi(i+1,j,k) - phi(i-1,j,k)  
	   if( dabs(delf) .lt. 1.D-5 ) then
	      cf(i,j,k,1) = 0.125D0
	   else
              phie = ( phi(i,j,k)-phi(i-1,j,k) ) / delf
	      if ( phie .le. -1.D0 .or. phie .ge. 1.5D0 ) then
	         cf(i,j,k,1) = 0.125D0
	      else
	         if ( phie .le. 0.25D0 ) then
	            if ( phie .le. 0.D0 ) then
	               cf(i,j,k,1) = 0.5D0 + 0.375D0 * phie
	            else
	               cf(i,j,k,1) = 0.5D0 - 0.625D0 * dsqrt(phie)
	            endif
	         else
	            if ( phie .le. 1.D0 ) then
	               cf(i,j,k,1) = 0.25D0 * ( 1.D0 - phie )
	            else
	               cf(i,j,k,1) =-0.25d0 * ( 1.D0 - phie )
	            endif
	         endif
	      endif
    	   endif
           phif(i,j,k,1) = 0.5D0 * ( phi(i,j,k) + phi(i+1,j,k) )
     <	                 - cf(i,j,k,1) * ( phi(i-1,j,k) 
     <	                                 - phi(i,j,k)
     <	                                 - phi(i,j,k)
     <	                                 + phi(i+1,j,k) )
	else
	   delf = phi(i,j,k) - phi(i+2,j,k)
	   if( dabs(delf) .lt. 1.D-5 ) then
	      cf(i,j,k,1) = 0.125D0
	   else
              phie = ( phi(i+1,j,k)-phi(i+2,j,k) ) / delf
	      if ( phie .le. -1.D0 .or. phie .ge. 1.5D0 ) then
	         cf(i,j,k,1) = 0.125D0
	      else
	         if ( phie .le. 0.25D0 ) then
	            if ( phie .le. 0.D0 ) then
	               cf(i,j,k,1) = 0.5D0 + 0.375D0 * phie
	            else
	               cf(i,j,k,1) = 0.5D0 - 0.625D0 * dsqrt(phie)
	            endif
	         else
	            if ( phie .le. 1.D0 ) then
	               cf(i,j,k,1) = 0.25D0 * ( 1.D0 - phie )
	            else
	               cf(i,j,k,1) =-0.25d0 * ( 1.D0 - phie )
	            endif
	         endif
	      endif
    	   endif
           phif(i,j,k,1) = 0.5D0 * ( phi(i,j,k) + phi(i+1,j,k) )
     <	                 - cf(i,j,k,1) * ( phi(i,j,k) 
     <	                                 - phi(i+1,j,k)
     <	                                 - phi(i+1,j,k)
     <	                                 + phi(i+2,j,k) )
        endif
	enddo
	enddo 
	enddo 

C...... J direction
	
	do k = 0, nnk
	do j = 0, nnj
	do i = 0, nni
	if ( uej(i,j,k) .ge. 0.D0 ) then
	   delf = phi(i,j+1,k) - phi(i,j-1,k)
	   if( dabs(delf) .lt. 1.D-5 ) then
	      cf(i,j,k,2) = 0.125D0
	   else
              phie = ( phi(i,j,k) - phi(i,j-1,k) ) / delf
	      if ( phie .le. -1.D0 .or. phie .ge. 1.5D0 ) then
	         cf(i,j,k,2) = 0.125D0
	      else
	         if ( phie .le. 0.25D0 ) then
	            if ( phie .le. 0.D0 ) then
	               cf(i,j,k,2) = 0.5D0 + 0.375D0 * phie
	            else
	               cf(i,j,k,2) = 0.5D0 - 0.625D0 * dsqrt(phie)
	            endif
	         else
	            if ( phie .le. 1.D0 ) then
	               cf(i,j,k,2) = 0.25D0 * ( 1.D0 - phie )
	            else
	               cf(i,j,k,2) =-0.25d0 * ( 1.D0 - phie )
	            endif
	         endif
	      endif
	   endif
           phif(i,j,k,2) = 0.5D0 * ( phi(i,j,k) + phi(i,j+1,k) )
     <	                 - cf(i,j,k,2) * ( phi(i,j-1,k) 
     <	                                 - phi(i,j,k)
     <	                                 - phi(i,j,k)
     <	                                 + phi(i,j+1,k) )
	else
	   delf = phi(i,j,k) - phi(i,j+2,k)
	   if( dabs(delf) .lt. 1.D-5 ) then
	      cf(i,j,k,2) = 0.125D0
	   else
              phie = ( phi(i,j+1,k) - phi(i,j+2,k) ) / delf
	      if ( phie .le. -1.D0 .or. phie .ge. 1.5D0 ) then
	         cf(i,j,k,2) = 0.125D0
	      else
	         if ( phie .le. 0.25D0 ) then
	            if ( phie .le. 0.D0 ) then
	               cf(i,j,k,2) = 0.5D0 + 0.375D0 * phie
	            else
	               cf(i,j,k,2) = 0.5D0 - 0.625D0 * dsqrt(phie)
	            endif
	         else
	            if ( phie .le. 1.D0 ) then
	               cf(i,j,k,2) = 0.25D0 * ( 1.D0 - phie )
	            else
	               cf(i,j,k,2) =-0.25d0 * ( 1.D0 - phie )
	            endif
	         endif
	      endif
	   endif
           phif(i,j,k,2) = 0.5D0 * ( phi(i,j,k) + phi(i,j+1,k) )
     <	                 - cf(i,j,k,2) * ( phi(i,j,k) 
     <	                                 - phi(i,j+1,k)
     <	                                 - phi(i,j+1,k)
     <	                                 + phi(i,j+2,k) )
	endif
	enddo
	enddo 
	enddo 

C...... K direction
	
	do k = 0, nnk
	do j = 0, nnj
	do i = 0, nni
	if ( uzk(i,j,k) .ge. 0.D0 ) then
	   delf = phi(i,j,k+1) - phi(i,j,k-1)
	   if( dabs(delf) .lt. 1.D-5 ) then
	      cf(i,j,k,3) = 0.125D0
	   else
              phie = ( phi(i,j,k) - phi(i,j,k-1) ) / delf
	      if ( phie .le. -1.D0 .or. phie .ge. 1.5D0 ) then
	         cf(i,j,k,3) = 0.125D0
	      else
	         if ( phie .le. 0.25D0 ) then
	            if ( phie .le. 0.D0 ) then
	               cf(i,j,k,3) = 0.5D0 + 0.375D0 * phie
	            else
	               cf(i,j,k,3) = 0.5D0 - 0.625D0 * dsqrt(phie)
	            endif
	         else
	            if ( phie .le. 1.D0 ) then
	               cf(i,j,k,3) = 0.25D0 * ( 1.D0 - phie )
	            else
	               cf(i,j,k,3) =-0.25d0 * ( 1.D0 - phie )
	            endif
	         endif
	      endif
	   endif
           phif(i,j,k,3) = 0.5D0 * ( phi(i,j,k) + phi(i,j,k+1) )
     <	                 - cf(i,j,k,3) * ( phi(i,j,k-1) 
     <	                                 - phi(i,j,k)
     <	                                 - phi(i,j,k)
     <	                                 + phi(i,j,k+1) )
	else
	   delf = phi(i,j,k) - phi(i,j,k+2)
	   if( dabs(delf) .lt. 1.D-5 ) then
	      cf(i,j,k,3) = 0.125D0
	   else
              phie = ( phi(i,j,k+1) - phi(i,j,k+2) ) / delf
	      if ( phie .le. -1.D0 .or. phie .ge. 1.5D0 ) then
	         cf(i,j,k,3) = 0.125D0
	      else
	         if ( phie .le. 0.25D0 ) then
	            if ( phie .le. 0.D0 ) then
	               cf(i,j,k,3) = 0.5D0 + 0.375D0 * phie
	            else
	               cf(i,j,k,3) = 0.5D0 - 0.625D0 * dsqrt(phie)
	            endif
	         else
	            if ( phie .le. 1.D0 ) then
	               cf(i,j,k,3) = 0.25D0 * ( 1.D0 - phie )
	            else
	               cf(i,j,k,3) =-0.25d0 * ( 1.D0 - phie )
	            endif
	         endif
	      endif
	   endif
           phif(i,j,k,3) = 0.5D0 * ( phi(i,j,k) + phi(i,j,k+1) )
     <	                 - cf(i,j,k,3) * ( phi(i,j,k) 
     <	                                 - phi(i,j,k+1)
     <	                                 - phi(i,j,k+1)
     <	                                 + phi(i,j,k+2) )
	endif
	enddo
	enddo 
	enddo 


C......	Convective terms (explicit)

	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni
	   hbs(i,j,k) = - uxi(i,  j,k) * phif(i,  j,k,1) 
     <	                + uxi(i-1,j,k) * phif(i-1,j,k,1)
     <	                - uej(i,j  ,k) * phif(i,j,  k,2)
     <	                + uej(i,j-1,k) * phif(i,j-1,k,2)
     <	                - uzk(i,j,k  ) * phif(i,j,k,  3)
     <	                + uzk(i,j,k-1) * phif(i,j,k-1,3)
	   HBS(I,J,K) = HBS(I,J,K) + PHI(I,J,K) *
     <	              ( UXI(I,J,K) - UXI(I-1,J,K)
     <	              + UEJ(I,J,K) - UEJ(I,J-1,K)
     <	              + UZK(I,J,K) - UZK(I,J,K-1) )
	enddo
	enddo
	enddo

C......	LES self-similarity term

	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni
	   hbs(i,j,k) = hbs(i,j,k) + rr(i,j,k,4)
	enddo
	enddo
	enddo

C......	Cross viscous terms at step n-1 from Crank-Nicolson  

	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni

	   hbs(i,j,k) = hbs(i,j,k) 
     <        + ( ak + 0.5D0*(akst(i,j,k) + akst(i+1,j,k)) ) *
     <		( g12(i,  j,k) * ( phi(i,  j+1,k) - phi(i,  j-1,k) 
     <		                 + phi(i+1,j+1,k) - phi(i+1,j-1,k) )
     <		+ g13(i,  j,k) * ( phi(i,  j,k+1) - phi(i,  j,k-1)
     <		                 + phi(i+1,j,k+1) - phi(i+1,j,k-1) ) ) 
     <        - ( ak + 0.5D0*(akst(i,j,k) + akst(i-1,j,k)) ) *
     <		( g12(i-1,j,k) * ( phi(i,  j+1,k) - phi(i,  j-1,k)
     <		                 + phi(i-1,j+1,k) - phi(i-1,j-1,k) )
     <		+ g13(i-1,j,k) * ( phi(i,  j,k+1) - phi(i,  j,k-1)
     <		                 + phi(i-1,j,k+1) - phi(i-1,j,k-1) ) )

	   hbs(i,j,k) = hbs(i,j,k) 
     <        + ( ak + 0.5D0*(akst(i,j,k) + akst(i,j+1,k)) ) *
     <		( g23(i,j,  k) * ( phi(i,j,  k+1) - phi(i,j,  k-1)
     <		                 + phi(i,j+1,k+1) - phi(i,j+1,k-1) )
     <		+ g21(i,j,  k) * ( phi(i+1,j,  k) - phi(i-1,j,  k)
     <		                 + phi(i+1,j+1,k) - phi(i-1,j+1,k) ) ) 
     <        - ( ak + 0.5D0*(akst(i,j,k) + akst(i,j-1,k)) ) *
     <		( g23(i,j-1,k) * ( phi(i,j,  k+1) - phi(i,j,  k-1)
     <		                 + phi(i,j-1,k+1) - phi(i,j-1,k-1) ) 
     <		+ g21(i,j-1,k) * ( phi(i+1,j,  k) - phi(i-1,j,  k)
     <		                 + phi(i+1,j-1,k) - phi(i-1,j-1,k) ) )

	   hbs(i,j,k) = hbs(i,j,k)  
     <        + ( ak + 0.5D0*(akst(i,j,k) + akst(i,j,k+1)) ) *
     <		( g31(i,j,k  ) * ( phi(i+1,j,k  ) - phi(i-1,j,k  )
     <		                 + phi(i+1,j,k+1) - phi(i-1,j,k+1) )
     <		+ g32(i,j,k  ) * ( phi(i,j+1,k  ) - phi(i,j-1,k  )
     <		                 + phi(i,j+1,k+1) - phi(i,j-1,k+1) ) ) 
     <        - ( ak + 0.5D0*(akst(i,j,k) + akst(i,j,k-1)) ) *
     <		( g31(i,j,k-1) * ( phi(i+1,j,k  ) - phi(i-1,j,k  )
     <		                 + phi(i+1,j,k-1) - phi(i-1,j,k-1) )
     <		+ g32(i,j,k-1) * ( phi(i,j+1,k  ) - phi(i,j-1,k  )
     <		                 + phi(i,j+1,k-1) - phi(i,j-1,k-1) ) )

	enddo
	enddo
	enddo

C......	Diagonal viscous terms at step n-1 from Crank-Nicolson 

	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni
	   sus(i,j,k) = sus(i,j,k) + coef * hbs(i,j,k)
	enddo
	enddo 
	enddo 

	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni
	   sus(i,j,k) = sus(i,j,k) 
     <        + ( ak + 0.5D0*(akst(i,j,k) + akst(i+1,j,k)) ) *
     <		g11(i,  j,  k  ) * ( phi(i+1,j,  k  ) - phi(i,j,k) )  
     <        + ( ak + 0.5D0*(akst(i,j,k) + akst(i-1,j,k)) ) *
     <		g11(i-1,j,  k  ) * ( phi(i-1,j,  k  ) - phi(i,j,k) )   
     <        + ( ak + 0.5D0*(akst(i,j,k) + akst(i,j+1,k)) ) *
     <		g22(i,  j,  k  ) * ( phi(i,  j+1,k  ) - phi(i,j,k) )   
     <        + ( ak + 0.5D0*(akst(i,j,k) + akst(i,j-1,k)) ) *
     <		g22(i,  j-1,k  ) * ( phi(i,  j-1,k  ) - phi(i,j,k) )  
     <        + ( ak + 0.5D0*(akst(i,j,k) + akst(i,j,k+1)) ) *
     <		g33(i,  j,  k  ) * ( phi(i,  j,  k+1) - phi(i,j,k) )  
     <        + ( ak + 0.5D0*(akst(i,j,k) + akst(i,j,k-1)) ) *
     <		g33(i,  j,  k-1) * ( phi(i,  j,  k-1) - phi(i,j,k) ) 
	enddo
	enddo 
	enddo 

	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni
	   sus(i,j,k) = dtime * jac(i,j,k) * sus(i,j,k)
	enddo
	enddo 
	enddo 

	return
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine scalar_solve

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "ns.inc"
	include "metric.inc"
	include "para.inc"
	include "eddy.inc"

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
	   fx(j,i) = sus(i,j,k)
	enddo
	enddo

	if ( n_west .eq. MPI_PROC_NULL ) then
	   do j = 1, nnj
	      hbs(0,j,k)=( g12(0,j,k) * ( phi(0,j+1,k) - phi(0,j-1,k)
     <	                                + phi(1,j+1,k) - phi(1,j-1,k) )
     <	                 + g13(0,j,k) * ( phi(0,j,k+1) - phi(0,j,k-1)
     <	                                + phi(1,j,k+1) - phi(1,j,k-1) ))
     <	                 / g11(0,j,k)
	      ax(j,0) =  0.D0
	      bx(j,0) =  1.D0
	      cx(j,0) = -1.D0
	      fx(j,0) =  hbs(0,j,k)
	   enddo
	endif
	
	if ( n_east .eq. MPI_PROC_NULL ) then
	   do j = 1, nnj
	      hbs(nni+1,j,k)=( g12(nni,j,k) 
     <	                     * ( phi(nni,  j+1,k) - phi(nni,  j-1,k)
     <	                       + phi(nni+1,j+1,k) - phi(nni+1,j-1,k) )
     <	                     + g13(nni,j,k) 
     <	                     * ( phi(nni,  j,k+1) - phi(nni,  j,k-1)
     <	                       + phi(nni+1,j,k+1) - phi(nni+1,j,k-1) ))
     <	                     / g11(nni,j,k)
	      ax(j,nni+1) =  1.D0
	      bx(j,nni+1) = -1.D0
	      cx(j,nni+1) =  0.D0
	      fx(j,nni+1) =  hbs(nni+1,j,k)
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
	   sus(i,j,k) = fx(j,i)
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
	   fy(i,j) = sus(i,j,k)
	enddo
	enddo

	if ( n_suth .eq. MPI_PROC_NULL ) then
	   do i = 1, nni
	      hbs(i,0,k)=( g23(i,0,k) * ( phi(i,0,k+1) - phi(i,0,k-1)
     <	                                + phi(i,1,k+1) - phi(i,1,k-1) )
     <	                 + g21(i,0,k) * ( phi(i+1,0,k) - phi(i-1,0,k)
     <	                                + phi(i+1,1,k) - phi(i-1,1,k) ))
     <	                 / g22(i,0,k)
	      ay(i,0) =  0.D0
	      by(i,0) =  1.D0
	      cy(i,0) = -1.D0
	      fy(i,0) =  hbs(i,0,k)
	   enddo
	endif
	
	if ( n_nrth .eq. MPI_PROC_NULL ) then
	   do i = 1, nni
	      hbs(i,nnj+1,k)=( g23(i,nnj,k) 
     <	                     * ( phi(i,nnj,  k+1) - phi(i,nnj,  k-1)
     <	                       + phi(i,nnj+1,k+1) - phi(i,nnj+1,k-1) )
     <	                     + g21(i,nnj,k) 
     <	                     * ( phi(i+1,nnj,  k) - phi(i-1,nnj,  k)
     <	                       + phi(i+1,nnj+1,k) - phi(i-1,nnj+1,k) ) )
     <	                   / g22(i,nnj,k)
	      ay(i,nnj+1) =  1.D0
	      by(i,nnj+1) = -1.D0
	      cy(i,nnj+1) =  0.D0
	      fy(i,nnj+1) =  hbs(i,nnj+1,k)
	   enddo
	endif
	
	call trid( ay, by, cy, fy, nni, 1, nnj, n_suth, n_nrth )

	do j = 1, nnj
	do i = 1, nni
	   sus(i,j,k) = fy(i,j)
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
	   fz(i,k) = sus(i,j,k)
	enddo
	enddo

	if ( n_back .eq. MPI_PROC_NULL ) then
	   do i = 1, nni
	      hbs(i,j,0)=( g31(i,j,0) * ( phi(i+1,j,0) - phi(i-1,j,0)
     <	                                + phi(i+1,j,1) - phi(i-1,j,1) )
     <	                 + g32(i,j,0) * ( phi(i,j+1,0) - phi(i,j-1,0)
     <	                                + phi(i,j+1,1) - phi(i,j-1,1) ))
     <	                 / g33(i,j,0)
	      az(i,0) =  0.D0
	      bz(i,0) =  1.D0
	      cz(i,0) = -1.D0
	      fz(i,0) =  hbs(i,j,0)
	   enddo
	endif
	
	if ( n_frnt .eq. MPI_PROC_NULL ) then
	   do i = 1, nni
	      hbs(i,j,nnk+1)=( g31(i,j,nnk) 
     <	                     * ( phi(i+1,j,nnk  ) - phi(i-1,j,nnk  )
     <	                       + phi(i+1,j,nnk+1) - phi(i-1,j,nnk+1) )
     <	                     + g32(i,j,nnk) 
     <	                     * ( phi(i,j+1,nnk  ) - phi(i,j-1,nnk  )
     <	                       + phi(i,j+1,nnk+1) - phi(i,j-1,nnk+1) ) )
     <	                     / g33(i,j,nnk)
	      az(i,nnk+1) =  1.D0
	      bz(i,nnk+1) = -1.D0
	      cz(i,nnk+1) =  0.D0
	      fz(i,nnk+1) =  hbs(i,j,nnk+1)
	   enddo
	endif
	
	if ( periodic .eq. 1 ) then
	call trip( az, bz, cz, fz, nni, 1, nnk, n_back, n_frnt )
	else
	call trid( az, bz, cz, fz, nni, 1, nnk, n_back, n_frnt )
	endif

	do k = 1, nnk
	do i = 1, nni
	   sus(i,j,k) = fz(i,k)
	enddo
	enddo

	enddo
	
C...... update scalar field (phi)

	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni
	   phi(i,j,k) = phi(i,j,k) + sus(i,j,k)
	enddo
	enddo
	enddo

C...... Change information at boundaries

	call phi_bc
	call phi_exchange

C...... Calculate energy
      call energy

	return
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine phi_bc

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "ns.inc"
	
	integer i, j, k

	if ( n_west .eq. MPI_PROC_NULL ) then
	   do k = -1, nnk+2
	   do j = -1, nnj+2
	      phi( 0,j,k) = phi(1,j,k) + hbs(0,j,k)
	   enddo
	   enddo
	   do k = -1, nnk+2
	   do j = -1, nnj+2
	      phi(-1,j,k) = 3.D0 * ( phi(0,j,k) - phi(1,j,k) )
     <	                  + phi(2,j,k)
	   enddo
	   enddo
	endif

	if ( n_east .eq. MPI_PROC_NULL ) then
	   do k = -1, nnk+2
	   do j = -1, nnj+2
	      phi(nni+1,j,k) = phi(nni,j,k) - hbs(nni+1,j,k)
	   enddo
	   enddo
	   do k = -1, nnk+2
	   do j = -1, nnj+2
	      phi(nni+2,j,k) = 3.D0 * ( phi(nni+1,j,k) - phi(nni,j,k) )
     <	                     + phi(nni-1,j,k)
	   enddo
	   enddo
	endif

	if ( n_suth .eq. MPI_PROC_NULL ) then
	   do k = -1, nnk+2
	   do i = -1, nni+2
	      phi(i, 0,k) = phi(i,1,k) + hbs(i,0,k)
	   enddo
	   enddo
	   do k = -1, nnk+2
	   do i = -1, nni+2
	      phi(i,-1,k) = 3.D0 * ( phi(i,0,k) - phi(i,1,k) )
     <	                  + phi(i,2,k)
	   enddo
	   enddo
	endif

	if ( n_nrth .eq. MPI_PROC_NULL ) then
	   do k = -1, nnk+2
	   do i = -1, nni+2
	      phi(i,nnj+1,k) = phi(i,nnj,k) - hbs(i,nnj+1,k)
	   enddo
	   enddo
	   do k = -1, nnk+2
	   do i = -1, nni+2
	      phi(i,nnj+2,k) = 3.D0 * ( phi(i,nnj+1,k) - phi(i,nnj,k) )
     <	                     + phi(i,nnj-1,k)
	   enddo
	   enddo
	endif

	if ( n_back .eq. MPI_PROC_NULL ) then
	   do j = -1, nnj+2
	   do i = -1, nni+2
	      phi(i,j, 0) = phi(i,j,1) + hbs(i,j,0)
	   enddo
	   enddo
	   do j = -1, nnj+2
	   do i = -1, nni+2
	      phi(i,j,-1) = 3.D0 * ( phi(i,j,0) - phi(i,j,1) )
     <	                  + phi(i,j,2)
	   enddo
	   enddo
	endif

	if ( n_frnt .eq. MPI_PROC_NULL ) then
	   do j = -1, nnj+2
	   do i = -1, nni+2
	      phi(i,j,nnk+1) = phi(i,j,nnk) - hbs(i,j,nnk+1)
	   enddo
	   enddo
	   do j = -1, nnj+2
	   do i = -1, nni+2
	      phi(i,j,nnk+2) = 3.D0 * ( phi(i,j,nnk+1) - phi(i,j,nnk) )
     <	                     + phi(i,j,nnk-1)
	   enddo
	   enddo
	endif

	if ( n_west .eq. MPI_PROC_NULL ) then
	   if ( n_suth .eq. MPI_PROC_NULL ) then
	      do k = 1, nnk
C******	0 - 0 - *
   	phi(0,0,k) = 0.5D0 * ( phi(1,0,k) + phi(0,1,k) )
	      enddo
	   endif
	   if ( n_nrth .eq. MPI_PROC_NULL ) then
	      do k = 1, nnk
C******	0 - Y - *
 	phi(0,nnj+1,k) = 0.5D0 * ( phi(1,nnj+1,k) + phi(0,nnj,k) )
	      enddo
	   endif
	   if ( n_back .eq. MPI_PROC_NULL ) then
	      do j = 1, nnj
C******	0 - * - 0
   	phi(0,j,0) = 0.5D0 * ( phi(1,j,0) + phi(0,j,1) )
	      enddo
	   endif
	   if ( n_frnt .eq. MPI_PROC_NULL ) then
	      do j = 1, nnj
C******	0 - * - Z
   	phi(0,j,nnk+1) = 0.5D0 * ( phi(1,j,nnk+1) + phi(0,j,nnk) )
	      enddo
	   endif
	endif	      

	if ( n_east .eq. MPI_PROC_NULL ) then
	   if ( n_suth .eq. MPI_PROC_NULL ) then
	      do k = 1, nnk
C******	X - 0 - *
   	phi(nni+1,0,k) = 0.5D0 * ( phi(nni,0,k) + phi(nni+1,1,k) )
	      enddo
	   endif
	   if ( n_nrth .eq. MPI_PROC_NULL ) then
	      do k = 1, nnk
C******	X - Y - *
  	phi(nni+1,nnj+1,k) = 0.5D0 * ( phi(nni,nnj+1,k) 
     <	                             + phi(nni+1,nnj,k) )
	      enddo
	   endif
	   if ( n_back .eq. MPI_PROC_NULL ) then
	      do j = 1, nnj
C******	X - * - 0
   	phi(nni+1,j,0) = 0.5D0 * ( phi(nni,j,0) + phi(nni+1,j,1) )
	      enddo
	   endif
	   if ( n_frnt .eq. MPI_PROC_NULL ) then
	      do j = 1, nnj
C******	X - * - Z
  	phi(nni+1,j,nnk+1) = 0.5D0 * ( phi(nni,j,nnk+1) 
     <	                             + phi(nni+1,j,nnk) )
	      enddo
	   endif
	endif	
	
	if ( n_suth .eq. MPI_PROC_NULL ) then
	   if ( n_back .eq. MPI_PROC_NULL ) then
	      do i = 1, nni
C******	* - 0 - 0
  	phi(i,0,0) = 0.5D0 * ( phi(i,1,0) + phi(i,0,1) )
	      enddo
	   endif
	   if ( n_frnt .eq. MPI_PROC_NULL ) then
	      do i = 1, nni
C******	* - 0 - Z
  	phi(i,0,nnk+1) = 0.5D0 * ( phi(i,1,nnk+1) + phi(i,0,nnk) )
	      enddo
	   endif
	endif      

	if ( n_nrth .eq. MPI_PROC_NULL ) then
	   if ( n_back .eq. MPI_PROC_NULL ) then
	      do i = 1, nni
C******	* - Y - 0
  	phi(i,nnj+1,0) = 0.5D0 * ( phi(i,nnj,0) + phi(i,nnj+1,1) )
	      enddo
	   endif
	   if ( n_frnt .eq. MPI_PROC_NULL ) then
	      do i = 1, nni
C******	* - Y - Z
  	phi(i,nnj+1,nnk+1) = 0.5D0 * ( phi(i,nnj,nnk+1)
     <	                             + phi(i,nnj+1,nnk) )
	      enddo
	   endif
	endif  

	return
	end

