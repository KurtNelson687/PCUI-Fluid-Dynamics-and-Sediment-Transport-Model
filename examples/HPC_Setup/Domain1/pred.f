cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine predictor

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "ns.inc"
	include "metric.inc"
	include "para.inc"
	include "eddy.inc"

	double precision, dimension(0:nni+1,0:nnj+1,0:nnk+1,3) :: su

	double precision, dimension(nnj,  0:nni+1) :: ax, bx, cx
	double precision, dimension(nnj,3,0:nni+1) :: fx
	double precision, dimension(nni,  0:nnj+1) :: ay, by, cy
	double precision, dimension(nni,3,0:nnj+1) :: fy
	double precision, dimension(nni,  0:nnk+1) :: az, bz, cz
	double precision, dimension(nni,3,0:nnk+1) :: fz

	integer i, j, k, m, L
	double precision coef, temp, dpdxi, dpdet, dpdzt
	
        coef = 1.5d0

C......	Take an Euler step on the first step

        if ( istep .eq. 1 ) then
           coef = 1.d0
           do m = 1, 3
           do k = 0, nnk+1
           do j = 0, nnj+1
           do i = 0, nni+1
              hb(i,j,k,m) = 0.D0
	   enddo
	   enddo
	   enddo 
	   enddo 
        endif

C......	First put in the part of Adams Bashforth from step n-2 (step 1 in notes)

	do m = 1, 3
	do k = 0, nnk+1
	do j = 0, nnj+1
	do i = 0, nni+1
           su(i,j,k,m) = - 0.5D0 * hb(i,j,k,m)
	enddo
	enddo
	enddo 
	enddo 

C......	Convective terms (explicit) (step 2 in notes)

        call convection

C......	LES self-similarity term (step 3 in notes)

	do m = 1, 3
	do k = 0, nnk+1
	do j = 0, nnj+1
	do i = 0, nni+1
           hb(i,j,k,m) = hb(i,j,k,m) + rr(i,j,k,m)
	enddo
	enddo
	enddo 
	enddo 

C......	Cross viscous terms at step n-1 from Crank-Nicolson (step 4 in notes - I think this should be from Adams Bashforth, not Crank-Nicolson) 

	do m = 1, 3
	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni

	   hb(i,j,k,m) = hb(i,j,k,m) 
     <        + ( vis + 0.5D0*(vst(i,j,k)+vst(i+1,j,k)) ) *
     <		( g12(i,  j,k) * ( u(i,  j+1,k,m) - u(i,  j-1,k,m) 
     <		                 + u(i+1,j+1,k,m) - u(i+1,j-1,k,m) )
     <		+ g13(i,  j,k) * ( u(i,  j,k+1,m) - u(i,  j,k-1,m)
     <		                 + u(i+1,j,k+1,m) - u(i+1,j,k-1,m) ) ) 
     <        - ( vis + 0.5D0*(vst(i,j,k)+vst(i-1,j,k)) ) *
     <		( g12(i-1,j,k) * ( u(i,  j+1,k,m) - u(i,  j-1,k,m)
     <		                 + u(i-1,j+1,k,m) - u(i-1,j-1,k,m) )
     <		+ g13(i-1,j,k) * ( u(i,  j,k+1,m) - u(i,  j,k-1,m)
     <		                 + u(i-1,j,k+1,m) - u(i-1,j,k-1,m) ) )

	   hb(i,j,k,m) = hb(i,j,k,m) 
     <        + ( vis + 0.5D0*(vst(i,j,k)+vst(i,j+1,k)) ) *
     <		( g23(i,j,  k) * ( u(i,j,  k+1,m) - u(i,j,  k-1,m)
     <		                 + u(i,j+1,k+1,m) - u(i,j+1,k-1,m) )
     <		+ g21(i,j,  k) * ( u(i+1,j,  k,m) - u(i-1,j,  k,m)
     <		                 + u(i+1,j+1,k,m) - u(i-1,j+1,k,m) ) ) 
     <        - ( vis + 0.5D0*(vst(i,j,k)+vst(i,j-1,k)) ) *
     <		( g23(i,j-1,k) * ( u(i,j,  k+1,m) - u(i,j,  k-1,m)
     <		                 + u(i,j-1,k+1,m) - u(i,j-1,k-1,m) ) 
     <		+ g21(i,j-1,k) * ( u(i+1,j,  k,m) - u(i-1,j,  k,m)
     <		                 + u(i+1,j-1,k,m) - u(i-1,j-1,k,m) ) )

	   hb(i,j,k,m) = hb(i,j,k,m)  
     <        + ( vis + 0.5D0*(vst(i,j,k)+vst(i,j,k+1)) ) *
     <		( g31(i,j,k  ) * ( u(i+1,j,k,  m) - u(i-1,j,k,  m)
     <		                 + u(i+1,j,k+1,m) - u(i-1,j,k+1,m) )
     <		+ g32(i,j,k  ) * ( u(i,j+1,k,  m) - u(i,j-1,k,  m)
     <		                 + u(i,j+1,k+1,m) - u(i,j-1,k+1,m) ) ) 
     <        - ( vis + 0.5D0*(vst(i,j,k)+vst(i,j,k-1)) ) *
     <		( g31(i,j,k-1) * ( u(i+1,j,k,  m) - u(i-1,j,k,  m)
     <		                 + u(i+1,j,k-1,m) - u(i-1,j,k-1,m) )
     <		+ g32(i,j,k-1) * ( u(i,j+1,k,  m) - u(i,j-1,k,  m)
     <		                 + u(i,j+1,k-1,m) - u(i,j-1,k-1,m) ) )

	enddo
	enddo
	enddo
	enddo

C......	Coriolis and bouyance force terms

	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni
	   temp = 1.D0 / jac(i,j,k)
	   hb(i,j,k,1) = hb(i,j,k,1) - omg2 * u(i,j,k,3) * temp
     <              + dpdxSteady*temp/rhoWater+dpdxWave*temp/rhoWater !Started with1.5D-6 and halved it each consecutive run 
	   hb(i,j,k,3) = hb(i,j,k,3) + omg2 * u(i,j,k,1) * temp
	   hb(i,j,k,2) = hb(i,j,k,2) ! - g * (  rho(i,j,k)
!     <                                      - rhoWater)/rhoWater * temp
	enddo
	enddo
	enddo

C......	Add to the source terms

	do m = 1, 3
	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni
	   su(i,j,k,m) = su(i,j,k,m) + coef * hb(i,j,k,m)
	enddo
	enddo
	enddo 
	enddo 

C......	Diagonal viscous terms at step n-1 from Crank-Nicolson 

	do m = 1, 3
	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni
	   su(i,j,k,m) = su(i,j,k,m)  
     <        + ( vis + 0.5D0*(vst(i,j,k) + vst(i+1,j,k)) ) * 
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

C......	Multiply dt

	do m = 1, 3
	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni
	   su(i,j,k,m) = dtime * jac(i,j,k) * su(i,j,k,m)
	enddo
	enddo
	enddo 
	enddo

C......	........................................

	coef = 0.5D0 * dtime

C...... solve for I-direction (This is setting up the approximate factorization in the x direction. See Olivers notes in lecture 11 section four, or Zang's 1994 paper.

	do k = 1, nnk

	do j = 1, nnj
	do i = 1, nni
	   ax(j,i) = -( vis + 0.5D0*(vst(i,j,k) + vst(i-1,j,k)) ) *
     <               coef * jac(i,j,k) * g11(i-1,j,k)
	   cx(j,i) = -( vis + 0.5D0*(vst(i,j,k) + vst(i+1,j,k)) ) * 
     <               coef * jac(i,j,k) * g11(i,  j,k)
	   bx(j,i) = 1.D0 - ax(j,i) - cx(j,i)
	enddo
	enddo

	do m = 1, 3
	do j = 1, nnj
	do i = 1, nni
	   fx(j,m,i) = su(i,j,k,m) !This is the right hand side for the above ax, bx, cx values for each velocity component
	enddo
	enddo
	enddo

	if ( n_west .eq. MPI_PROC_NULL ) then

	   do j = 1, nnj
	      ax(j,0) = 0.D0
	      bx(j,0) = -1.D0
	      cx(j,0) = 1.D0
	   enddo

CBCBCBC	BCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBC

	   do j = 1, nnj

	      dpdxi = p(1,j,  k) - p(0,j  ,k)
	      dpdet = p(0,j+1,k) - p(0,j-1,k) 
     <              + p(1,j+1,k) - p(1,j-1,k)
	      dpdzt = p(0,j,k+1) - p(0,j,k-1)
     <              + p(1,j,k+1) - p(1,j,k-1) 

	      fx(j,1,0) =   2.D0 * dtime * jac(1,j,k)  *      
     <          (                          xix(0,j,k)  * dpdxi
     <           + 0.125D0 * (etx(1,j-1,k)+etx(1,j,k)) * dpdet          
     <           + 0.125D0 * (ztx(1,j,k-1)+ztx(1,j,k)) * dpdzt )  
	      fx(j,2,0) =   2.D0 * dtime * jac(1,j,k)  *      
     <          (                          xiy(0,j,k)  * dpdxi
     <           + 0.125D0 * (ety(1,j-1,k)+ety(1,j,k)) * dpdet          
     <           + 0.125D0 * (zty(1,j,k-1)+zty(1,j,k)) * dpdzt )  
	      fx(j,3,0) =   2.D0 * dtime * jac(1,j,k)  *      
     <          (                          xiz(0,j,k)  * dpdxi
     <           + 0.125D0 * (etz(1,j-1,k)+etz(1,j,k)) * dpdet          
     <           + 0.125D0 * (ztz(1,j,k-1)+ztz(1,j,k)) * dpdzt )  

	   enddo

CBCBCBC	BCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBC

	endif
	
	if ( n_east .eq. MPI_PROC_NULL ) then

	   do j = 1, nnj
	      ax(j,nni+1) = 1.D0
	      bx(j,nni+1) = 1.D0
	      cx(j,nni+1) = 0.D0
	   enddo

CBCBCBC	BCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBC

	   do j = 1, nnj

	      dpdxi = p(nni+1,j  ,k) - p(nni  ,j  ,k)
	      dpdet = p(nni  ,j+1,k) - p(nni  ,j-1,k) 
     <              + p(nni+1,j+1,k) - p(nni+1,j-1,k) 
	      dpdzt = p(nni  ,j,k+1) - p(nni  ,j,k-1)
     <              + p(nni+1,j,k+1) - p(nni+1,j,k-1) 

	      fx(j,1,nni+1) = 2.D0 * dtime * jac(nni,j,k)  *      
     <          (                            xix(nni,j,k)  * dpdxi
     <           + 0.125D0 * (etx(nni,j-1,k)+etx(nni,j,k)) * dpdet          
     <           + 0.125D0 * (ztx(nni,j,k-1)+ztx(nni,j,k)) * dpdzt )  

	      fx(j,2,nni+1) = 2.D0 * dtime * jac(nni,j,k)  *      
     <          (                            xiy(nni,j,k)  * dpdxi
     <           + 0.125D0 * (ety(nni,j-1,k)+ety(nni,j,k)) * dpdet          
     <		 + 0.125D0 * (zty(nni,j,k-1)+zty(nni,j,k)) * dpdzt )  

	      fx(j,3,nni+1) = 2.D0 * dtime * jac(nni,j,k)  *      
     <          (                            xiz(nni,j,k)  * dpdxi
     <           + 0.125D0 * (etz(nni,j-1,k)+etz(nni,j,k)) * dpdet          
     <           + 0.125D0 * (ztz(nni,j,k-1)+ztz(nni,j,k)) * dpdzt )  

	   enddo

CBCBCBC	BCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBC

	endif
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Added by Kurt
	if ( periodic .eq. 1 ) then
	call trip( ax, bx, cx, fx, nnj, 3, nni, n_west, n_east )
	else 
	call trid( ax, bx, cx, fx, nnj, 3, nni, n_west, n_east )
	endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	call trid( ax, bx, cx, fx, nnj, 3, nni, n_west, n_east ) !Commented by Kurt

	do m = 1, 3
	do j = 1, nnj
	do i = 1, nni
	   su(i,j,k,m) = fx(j,m,i)
	enddo
	enddo
	enddo

	enddo

C...... solve for J-direction

	do k = 1, nnk

	do j = 1, nnj
	do i = 1, nni
	   ay(i,j) = -( vis + 0.5D0*(vst(i,j,k) + vst(i,j-1,k)) ) * 
     <	             coef * jac(i,j,k) * g22(i,j-1,k)
	   cy(i,j) = -( vis + 0.5D0*(vst(i,j,k) + vst(i,j+1,k)) ) * 
     <	             coef * jac(i,j,k) * g22(i,j,  k)
	   by(i,j) = 1.D0 - ay(i,j) - cy(i,j)
	enddo
	enddo

	do m = 1, 3
	do j = 1, nnj
	do i = 1, nni
	   fy(i,m,j) = su(i,j,k,m)
	enddo
	enddo
	enddo

	if ( n_suth .eq. MPI_PROC_NULL ) then

	   do i = 1, nni
	      ay(i,0) = 0.D0
	      by(i,0) = 1.D0
	      cy(i,0) = 1.D0
	   enddo

CBCBCBC	BCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBC

	   do i = 1, nni

	      dpdet = p(i  ,1,k) - p(i  ,0,k)
	      dpdxi = p(i+1,0,k) - p(i-1,0,k) 
     <              + p(i+1,1,k) - p(i-1,1,k)
	      dpdzt = p(i,0,k+1) - p(i,0,k-1)
     <		    + p(i,1,k+1) - p(i,1,k-1)

	      fy(i,1,0) =  2.D0 * dtime * jac(i,1,k)  * 
     <          ( 0.125D0 * (xix(i-1,1,k)+xix(i,1,k)) * dpdxi
     <          +                         etx(i,0,k)  * dpdet
     <          + 0.125D0 * (ztx(i,1,k-1)+ztx(i,1,k)) * dpdzt  )

	      fy(i,2,0) =  2.D0 * dtime * jac(i,1,k)  * 
     <          ( 0.125D0 * (xiy(i-1,1,k)+xiy(i,1,k)) * dpdxi
     <          +                         ety(i,0,k)  * dpdet
     <          + 0.125D0 * (zty(i,1,k-1)+zty(i,1,k)) * dpdzt  )

	      fy(i,3,0) =  2.D0 * dtime * jac(i,1,k)  * 
     <          ( 0.125D0 * (xiz(i-1,1,k)+xiz(i,1,k)) * dpdxi
     <          +                         etz(i,0,k)  * dpdet
     <          + 0.125D0 * (ztz(i,1,k-1)+ztz(i,1,k)) * dpdzt  )

	   enddo

CBCBCBC	BCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBC

	endif
	
	if ( n_nrth .eq. MPI_PROC_NULL ) then

	   do i = 1, nni
	      ay(i,nnj+1) = 1.D0
	      by(i,nnj+1) = -1.D0
	      cy(i,nnj+1) = 0.D0
	   enddo

CBCBCBC	BCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBC

	   do i = 1, nni

	      dpdet = p(i  ,nnj+1,k) - p(i  ,nnj,k)
	      dpdxi = p(i+1,nnj  ,k) - p(i-1,nnj  ,k) 
     <              + p(i+1,nnj+1,k) - p(i-1,nnj+1,k)
	      dpdzt = p(i,nnj  ,k+1) - p(i,nnj  ,k-1)
     <		    + p(i,nnj+1,k+1) - p(i,nnj+1,k-1)

	      fy(i,1,nnj+1) = 2.D0 * dtime * jac(i,nnj,k)  * 
     <          (  0.125D0 * (xix(i-1,nnj,k)+xix(i,nnj,k)) * dpdxi
     <           +                           etx(i,nnj,k)  * dpdet
     <           + 0.125D0 * (ztx(i,nnj,k-1)+ztx(i,nnj,k)) * dpdzt  )

	      fy(i,2,nnj+1) = 2.D0 * dtime * jac(i,nnj,k)  * 
     <          (  0.125D0 * (xiy(i-1,nnj,k)+xiy(i,nnj,k)) * dpdxi
     <           +                           ety(i,nnj,k)  * dpdet
     <           + 0.125D0 * (zty(i,nnj,k-1)+zty(i,nnj,k)) * dpdzt  )

	      fy(i,3,nnj+1) = 2.D0 * dtime * jac(i,nnj,k)  * 
     <          (  0.125D0 * (xiz(i-1,nnj,k)+xiz(i,nnj,k)) * dpdxi
     <           +                           etz(i,nnj,k)  * dpdet
     <           + 0.125D0 * (ztz(i,nnj,k-1)+ztz(i,nnj,k)) * dpdzt  )

	   enddo

CBCBCBC	BCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBC

	endif
	
	call trid( ay, by, cy, fy, nni, 3, nnj, n_suth, n_nrth )

	do m = 1, 3
	do j = 1, nnj
	do i = 1, nni
	  if (m .ne. 2) then
	   su(i,j,k,m) = fy(i,m,j)
	  endif
	enddo
	enddo
	enddo

	enddo


C...... solve for J-direction (only obtain the vertical component)
	do k = 1, nnk

	do j = 1, nnj
	do i = 1, nni
	   ay(i,j) = -( vis + 0.5D0*(vst(i,j,k) + vst(i,j-1,k)) ) * 
     <	             coef * jac(i,j,k) * g22(i,j-1,k)
	   cy(i,j) = -( vis + 0.5D0*(vst(i,j,k) + vst(i,j+1,k)) ) * 
     <	             coef * jac(i,j,k) * g22(i,j,  k)
	   by(i,j) = 1.D0 - ay(i,j) - cy(i,j)
           
	enddo
	enddo
	

	do j = 1, nnj
	do i = 1, nni
	   fy(i,2,j) = su(i,j,k,2)
	   
	enddo
	enddo
	

	if ( n_suth .eq. MPI_PROC_NULL ) then
       	      do i = 1, nni
	         ay(i,0) = 0.D0
	         by(i,0) = 1.D0
	         cy(i,0) = 1.D0
	      enddo
CBCBCBC	BCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBC

	   do i = 1, nni

	      dpdet = p(i  ,1,k) - p(i  ,0,k)
	      dpdxi = p(i+1,0,k) - p(i-1,0,k) 
     <              + p(i+1,1,k) - p(i-1,1,k)
	      dpdzt = p(i,0,k+1) - p(i,0,k-1)
     <		    + p(i,1,k+1) - p(i,1,k-1)

	      fy(i,1,0) =  2.D0 * dtime * jac(i,1,k)  * 
     <          ( 0.125D0 * (xix(i-1,1,k)+xix(i,1,k)) * dpdxi
     <          +                         etx(i,0,k)  * dpdet
     <          + 0.125D0 * (ztx(i,1,k-1)+ztx(i,1,k)) * dpdzt  )

	      fy(i,2,0) =  2.D0 * dtime * jac(i,1,k)  * 
     <          ( 0.125D0 * (xiy(i-1,1,k)+xiy(i,1,k)) * dpdxi
     <          +                         ety(i,0,k)  * dpdet
     <          + 0.125D0 * (zty(i,1,k-1)+zty(i,1,k)) * dpdzt  )
	      
	      fy(i,3,0) =  2.D0 * dtime * jac(i,1,k)  * 
     <          ( 0.125D0 * (xiz(i-1,1,k)+xiz(i,1,k)) * dpdxi
     <          +                         etz(i,0,k)  * dpdet
     <          + 0.125D0 * (ztz(i,1,k-1)+ztz(i,1,k)) * dpdzt  )

	   enddo

CBCBCBC	BCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBC

	endif
	
	if ( n_nrth .eq. MPI_PROC_NULL ) then

	   do i = 1, nni
	      ay(i,nnj+1) = 1.D0
	      by(i,nnj+1) = 1.D0
	      cy(i,nnj+1) = 0.D0
	   enddo

CBCBCBC	BCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBC

	   do i = 1, nni

	      dpdet = p(i  ,nnj+1,k) - p(i  ,nnj,k)
	      dpdxi = p(i+1,nnj  ,k) - p(i-1,nnj  ,k) 
     <              + p(i+1,nnj+1,k) - p(i-1,nnj+1,k)
	      dpdzt = p(i,nnj  ,k+1) - p(i,nnj  ,k-1)
     <		    + p(i,nnj+1,k+1) - p(i,nnj+1,k-1)

	      fy(i,1,nnj+1) = 2.D0 * dtime * jac(i,nnj,k)  * 
     <          (  0.125D0 * (xix(i-1,nnj,k)+xix(i,nnj,k)) * dpdxi
     <           +                           etx(i,nnj,k)  * dpdet
     <           + 0.125D0 * (ztx(i,nnj,k-1)+ztx(i,nnj,k)) * dpdzt  )

	      fy(i,2,nnj+1) = 2.D0 * dtime * jac(i,nnj,k)  * 
     <          (  0.125D0 * (xiy(i-1,nnj,k)+xiy(i,nnj,k)) * dpdxi
     <           +                           ety(i,nnj,k)  * dpdet
     <           + 0.125D0 * (zty(i,nnj,k-1)+zty(i,nnj,k)) * dpdzt  )

	      fy(i,3,nnj+1) = 2.D0 * dtime * jac(i,nnj,k)  * 
     <          (  0.125D0 * (xiz(i-1,nnj,k)+xiz(i,nnj,k)) * dpdxi
     <           +                           etz(i,nnj,k)  * dpdet
     <           + 0.125D0 * (ztz(i,nnj,k-1)+ztz(i,nnj,k)) * dpdzt  )

	   enddo

CBCBCBC	BCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBC

	endif
	
	
        call trid1( ay, by, cy, fy, nni, 2, 2,nnj, n_suth, n_nrth ) !This was written by YiJu	
c.....so, here we only obtain the result for the second (jth) component
	do j = 1, nnj
	do i = 1, nni
	   su(i,j,k,2) = fy(i,2,j) 
	enddo
	enddo

	enddo

C...... solve for K-direction
	
	do j = 1, nnj

	do k = 1, nnk
	do i = 1, nni
	   az(i,k) = -( vis + 0.5D0*(vst(i,j,k)+vst(i,j,k-1)) ) *
     <               coef * jac(i,j,k) * g33(i,j,k-1)
	   cz(i,k) = -( vis + 0.5D0*(vst(i,j,k)+vst(i,j,k+1)) ) * 
     <               coef * jac(i,j,k) * g33(i,j,k  )
	   bz(i,k) = 1.D0 - az(i,k) - cz(i,k)
	enddo
	enddo

	do m = 1, 3
	do k = 1, nnk
	do i = 1, nni
	   fz(i,m,k) = su(i,j,k,m)
	enddo
	enddo
	enddo

	if ( n_back .eq. MPI_PROC_NULL ) then

	   do i = 1, nni
	      az(i,0) = 0.D0
	      bz(i,0) = 1.D0
	      cz(i,0) = 1.D0
	   enddo

CBCBCBC	BCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBC

	   do i = 1, nni

	      dpdzt = p(i  ,j,1) - p(i  ,j,0) 
	      dpdxi = p(i+1,j,0) - p(i-1,j,0)
     <              + p(i+1,j,1) - p(i-1,j,1)
	      dpdet = p(i,j+1,0) - P(i,j-1,0) 
     <              + p(i,j+1,1) - p(i,j-1,1) 

	      fz(i,1,0) =   2.D0 * dtime * jac(i,j,1)  *
     <          (  0.125D0 * (xix(i-1,j,1)+xix(i,j,1)) * dpdxi
     <           + 0.125D0 * (etx(i,j-1,1)+etx(i,j,1)) * dpdet
     <           +                         ztx(i,j,0)  * dpdzt )

	      fz(i,2,0) =   2.D0 * dtime * jac(i,j,1)  *
     <          (  0.125D0 * (xiy(i-1,j,1)+xiy(i,j,1)) * dpdxi
     <           + 0.125D0 * (ety(i,j-1,1)+ety(i,j,1)) * dpdet
     <           +                         zty(i,j,0)  * dpdzt )

	      fz(i,3,0) =   2.D0 * dtime * jac(i,j,1)  *
     <          (  0.125D0 * (xiz(i-1,j,1)+xiz(i,j,1)) * dpdxi
     <           + 0.125D0 * (etz(i,j-1,1)+etz(i,j,1)) * dpdet
     <           +                         ztz(i,j,0)  * dpdzt )

	   enddo

CBCBCBC	BCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBC

	endif
	
	if ( n_frnt .eq. MPI_PROC_NULL ) then

	   do i = 1, nni
	      az(i,nnk+1) = 1.D0
	      bz(i,nnk+1) = 1.D0
	      cz(i,nnk+1) = 0.D0
	   enddo

CBCBCBC	BCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBC

	   do i = 1, nni

	      dpdzt = p(i  ,j,nnk+1) - p(i  ,j,nnk  ) 
	      dpdxi = p(i+1,j,nnk  ) - p(i-1,j,nnk  )
     <              + p(i+1,j,nnk+1) - p(i-1,j,nnk+1)
	      dpdet = p(i,j+1,nnk  ) - P(i,j-1,nnk  ) 
     <              + p(i,j+1,nnk+1) - p(i,j-1,nnk+1)

	      fz(i,1,nnk+1) = 2.D0 * dtime * jac(i,j,nnk)  *
     <          (  0.125D0 * (xix(i-1,j,nnk)+xix(i,j,nnk)) * dpdxi
     <           + 0.125D0 * (etx(i,j-1,nnk)+etx(i,j,nnk)) * dpdet
     <           +                           ztx(i,j,nnk)  * dpdzt )

	      fz(i,2,nnk+1) = 2.D0 * dtime * jac(i,j,nnk)  *
     <          (  0.125D0 * (xiy(i-1,j,nnk)+xiy(i,j,nnk)) * dpdxi
     <           + 0.125D0 * (ety(i,j-1,nnk)+ety(i,j,nnk)) * dpdet
     <           +                           zty(i,j,nnk)  * dpdzt )

	      fz(i,3,nnk+1) = 2.D0 * dtime * jac(i,j,nnk)  *
     <          (  0.125D0 * (xiz(i-1,j,nnk)+xiz(i,j,nnk)) * dpdxi
     <           + 0.125D0 * (etz(i,j-1,nnk)+etz(i,j,nnk)) * dpdet
     <           +                           ztz(i,j,nnk)  * dpdzt )
	   enddo

CBCBCBC	BCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBC

	endif

	if ( periodic .eq. 1 ) then
	call trip( az, bz, cz, fz, nni, 3, nnk, n_back, n_frnt )
	else 
	call trid( az, bz, cz, fz, nni, 3, nnk, n_back, n_frnt )
	endif

	do m = 1, 3
	do k = 1, nnk
	do i = 1, nni
	   su(i,j,k,m) = fz(i,m,k)
	enddo
	enddo
	enddo

	enddo

C...... Update intermediate velocity (u*)

	do m = 1, 3
	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni
	   u(i,j,k,m) = u(i,j,k,m) + su(i,j,k,m)
	enddo
	enddo
	enddo
	enddo

CBCBCBC	BCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBC
	call ustar_bc
CBCBCBC	BCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBC

	call u_exchange
	
	return
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine ustar_bc
C All this subroutine does is extropolate from interal values to fill ghost cells.
	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "ns.inc"
	
	integer i, j, k, L
	
	if ( n_west .eq. MPI_PROC_NULL ) then
	   do L = 1, 3
	   do k = -1, nnk+2
	   do j = -1, nnj+2
	      u( 0,j,k,L) = 3.D0*( u(1,j,k,L)-u(2,j,k,L) ) + u(3,j,k,L)
	   enddo
	   enddo
	   enddo
	   do L = 1, 3
	   do k = -1, nnk+2
	   do j = -1, nnj+2
	      u(-1,j,k,L) = 3.D0*( u(0,j,k,L)-u(1,j,k,L) ) + u(2,j,k,L)
	   enddo
	   enddo
	   enddo
	endif

	if ( n_east .eq. MPI_PROC_NULL ) then
	   do L = 1, 3
	   do k = -1, nnk+2
	   do j = -1, nnj+2
	      u(nni+1,j,k,L) = 3.D0*( u(nni,j,k,L)-u(nni-1,j,k,L) ) 
     <                       + u(nni-2,j,k,L)
	   enddo
	   enddo
	   enddo
	   do L = 1, 3
	   do k = -1, nnk+2
	   do j = -1, nnj+2
	      u(nni+2,j,k,L) = 3.D0*( u(nni+1,j,k,L)-u(nni,j,k,L) ) 
     <                       + u(nni-1,j,k,L)
	   enddo
	   enddo
	   enddo
	endif

	if ( n_suth .eq. MPI_PROC_NULL ) then
	   do L = 1, 3
	   do k = -1, nnk+2
	   do i = -1, nni+2
	      u(i, 0,k,L) = 3.D0*( u(i,1,k,L)-u(i,2,K,L) ) + u(i,3,k,L)
	   enddo
	   enddo
	   enddo
	   do L = 1, 3
	   do k = -1, nnk+2
	   do i = -1, nni+2
	      u(i,-1,k,L) = 3.D0*( u(i,0,k,L)-u(i,1,K,L) ) + u(i,2,k,L)
	   enddo
	   enddo
	   enddo
	endif

	if ( n_nrth .eq. MPI_PROC_NULL ) then
	   do L = 1, 3
	   do k = -1, nnk+2
	   do i = -1, nni+2
	      u(i,nnj+1,k,L) = 3.D0 * ( u(i,nnj,k,L)-u(i,nnj-1,k,L) ) 
     <                       + u(i,nnj-2,k,L)
	   enddo
	   enddo
	   enddo
	   do L = 1, 3
	   do k = -1, nnk+2
	   do i = -1, nni+2
	      u(i,nnj+2,k,L) = 3.D0 * ( u(i,nnj+1,k,L)-u(i,nnj,k,L) ) 
     <                       + u(i,nnj-1,k,L)
	   enddo
	   enddo
	   enddo
	endif

	if ( n_back .eq. MPI_PROC_NULL ) then
	   do L = 1, 3
	   do j = -1, nnj+2
	   do i = -1, nni+2
	      u(i,j, 0,L) = 3.D0*( u(i,j,1,L)-u(i,j,2,L) ) + u(i,j,3,L)
	   enddo
	   enddo
	   enddo
	   do L = 1, 3
	   do j = -1, nnj+2
	   do i = -1, nni+2
	      u(i,j,-1,L) = 3.D0*( u(i,j,0,L)-u(i,j,1,L) ) + u(i,j,2,L)
	   enddo
	   enddo
	   enddo
	endif

	if ( n_frnt .eq. MPI_PROC_NULL ) then
	   do L = 1, 3
	   do j = -1, nnj+2
	   do i = -1, nni+2
	      u(i,j,nnk+1,L) = 3.D0*( u(i,j,nnk,L)-u(i,j,nnk-1,L) ) 
     <                       + u(i,j,nnk-2,L)
	   enddo
	   enddo
	   enddo
	   do L = 1, 3
	   do j = -1, nnj+2
	   do i = -1, nni+2
	      u(i,j,nnk+2,L) = 3.D0*( u(i,j,nnk+1,L)-u(i,j,nnk,L) )
     <                       + u(i,j,nnk-1,L)

           enddo
	   enddo
	   enddo
	endif

	return
	end
