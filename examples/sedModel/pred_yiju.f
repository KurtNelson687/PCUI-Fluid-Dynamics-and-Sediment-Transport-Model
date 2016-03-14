cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine predictor

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "ns.inc"
	include "metric.inc"
	include "para.inc"
	include "eddy.inc"
	include "sedi.inc"

	double precision su(0:nni+1,0:nnj+1,0:nnk+1,3)

	double precision ax(nnj,  0:nni+1), bx(nnj,  0:nni+1), 
     &                   cx(nnj,  0:nni+1)
	double precision fx(nnj,3,0:nni+1)
	double precision ay(nni,  0:nnj+1), 
     &                   by(nni,  0:nnj+1), cy(nni,  0:nnj+1)
	double precision fy(nni,3,0:nnj+1)
	double precision az(nni,  0:nnk+1), 
     &                   bz(nni,  0:nnk+1), cz(nni,  0:nnk+1)      
	double precision fz(nni,3,0:nnk+1)

	integer i, j, k, m, L
	double precision coef, temp, dpdxi, dpdet, dpdzt, dz0
	
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
           
           do k = 0, nnk+1
           do j = 0, nnj+1
           do i = 0, nni+1
              drive(i,j,k) = 0.D0
	   enddo
	   enddo
	   enddo 
        endif

C......	First put in the part of Adams Bashforth from step n-2

	do m = 1, 3
	do k = 0, nnk+1
	do j = 0, nnj+1
	do i = 0, nni+1
           su(i,j,k,m) = - 0.5D0 * hb(i,j,k,m)
	enddo
	enddo
	enddo 
	enddo 

C...... Driving pressure gradient
        call driveP(drive)
        
	if( btmAug .eq. 1 ) then
	   call btmAugmentation
           do m = 1, 3
	   do k = 1, nnk
	   do j = 1, nnj
	   do i = 1, nni
	      su(i,j,k,m) = su(i,j,k,m) + btau(i,j,k,m)
c              write(*,*) btau(i,j,k,m)
	   enddo
	   enddo
	   enddo
	   enddo
	endif

	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni
	   su(i,j,k,1) = su(i,j,k,1) + drive(i,j,k)
	enddo
	enddo
	enddo
	call convection

	if (movingGrid .eq. 1) then
	   if (mod(istep, ncheck) .eq. 0) then
	      if (myid .eq. 0)
     <	      write(*,*) '....Updating bed elevation and grid......'
	      call oldGrid
	      call btm_metric
	      call updateDepth2
	      call updateGrid
	      call bedGeometry
	      call getEtadot
	      call update_wbc
	      call movingGrid_convection
	      if (myid .eq. 0)
     <	      write(*,*) '....Done with updating grid..............'
	   endif
	endif
C......	Convective terms (explicit)
c	if (movingGrid .eq. 1) then
c        call convection
c	call movingGrid_convection
c	else
c        call convection
c	endif

C......	LES self-similarity term

	do m = 1, 3
	do k = 0, nnk+1
	do j = 0, nnj+1
	do i = 0, nni+1
           hb(i,j,k,m) = hb(i,j,k,m) - rr(i,j,k,m)
	enddo
	enddo
	enddo 
	enddo 

C......	Cross viscous terms at step n-1 from Crank-Nicolson  

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

C......	Coriolis and bouyancy force terms

	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni
	   temp = 1.D0 / jac(i,j,k)
c	   hb(i,j,k,1) = hb(i,j,k,1) - omg2 * u(i,j,k,3) * temp 
c	   hb(i,j,k,3) = hb(i,j,k,3) + omg2 * u(i,j,k,1) * temp
	   if (iscalar .eq. 1) then
	      hb(i,j,k,2) = hb(i,j,k,2) - g * (  phi(i,j,k)
     <                                      - phi_init(i,j,k) ) * temp
	   endif
	   if (sedi .eq. 1) then
	      hb(i,j,k,2) = hb(i,j,k,2)-g*(1.D0-1.D0/spwght)
     <                     *spwght/rho_tot(i,j,k)*sedc(i,j,k)*temp
	   endif
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

	if (movingGrid .eq. 1) then
	   if (mod(istep, ncheck) .eq. 0) then
	    
	      temp = 1.D0/dtime
	      do m = 1, 3
	      do k = 1, nnk
	      do j = 1, nnj
	      do i = 1, nni

	         su(i,j,k,m) = su(i,j,k,m) - jac_diff(i,j,k)*u(i,j,k,m)
c     <                        (1.D0/jac_old(i,j,k)-1.D0/jac(i,j,k))
c     <                        *u(i,j,k,m)
	      enddo
	      enddo
	      enddo
	      enddo
	   endif
	endif


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

C...... solve for I-direction

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
	   fx(j,m,i) = su(i,j,k,m)
	enddo
	enddo
	enddo

	if ( n_west .eq. MPI_PROC_NULL ) then

	   do j = 1, nnj
	      ax(j,0) = 0.D0
	      bx(j,0) = 1.D0
c	      cx(j,0) = 1.D0
              cx(j,0) = -1.D0
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
c	      ax(j,nni+1) = 1.D0
              ax(j,nni+1) = -1.D0
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

	if ( periodic .eq. 1 ) then
	call trip( ax, bx, cx, fx, nnj, 3, nni, n_west, n_east )
	else 
	call trid( ax, bx, cx, fx, nnj, 3, nni, n_west, n_east )
	endif
	
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
	

c	if (movingGrid .eq. 1) then
c	   if ( n_suth .eq. MPI_PROC_NULL ) then
c	      do i = 1, nni
c		 fy(i,1,1) = fy(i,1,1)-2.0*ay(i,1)*u_bed(i,k,1)
c		 fy(i,2,1) = fy(i,2,1)-2.0*ay(i,1)*u_bed(i,k,2)
c		 fy(i,3,1) = fy(i,3,1)-2.0*ay(i,1)*u_bed(i,k,3)
c                 by(i,1) = by(i,1)-ay(i,1)
c		 ay(i,1) = 0.D0
c	      enddo
c	   endif
c	else
c	   if (logU .eq. 1) then
c	      do i = 1, nni
c		 dz0 = xp(i, 1, k, 2) - xp(i, 0, k, 2)
c		 fct(i,k) = Cdxi(i, k)*utan(i,k)*dz0/(vis + vst(i, 1, k))
c		 by(i, 1) = ay(i, 1)*(1.D0-fct(i,k)) + by(i, 1)
c		 ay(i, 1) = 0.D0
c	      enddo          
c	   endif
c	endif

	if ( n_suth .eq. MPI_PROC_NULL ) then
	   if (logU .eq. 1) then
	      do i = 1, nni
		 dz0 = xp(i, 1, k, 2) - xp(i, 0, k, 2)
		 fct(i,k) = Cdxi(i, k)*utan(i,k)*dz0/(vis + vst(i, 1, k))
		 by(i, 1) = ay(i, 1)*(1.D0-fct(i,k)) + by(i, 1)
		 ay(i, 1) = 0.D0
	      enddo          
	   endif	   
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
	
	
        call trid1( ay, by, cy, fy, nni, 2, 2,nnj, n_suth, n_nrth )	
c.....so, here we only obtain the result for the second (jth) component
	do j = 1, nnj
	do i = 1, nni
	   su(i,j,k,2) = fy(i,2,j) 
	enddo
	enddo

	enddo

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
	      cz(i,0) = -1.D0
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
	      bz(i,nnk+1) = -1.D0
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

	if ( periodicZ .eq. 1 ) then
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
        subroutine driveP(pres)

        include "size.inc"
        include "para.inc"
	include "metric.inc"
	include "ns.inc"
        
        integer i, j, k
        double precision pres(0:nni+1, 0:nnj+1, 0:nnk+1)
        double precision p0, tmp

        p0 = amp*omega*omega
        
        do i=0,nni+1
	do j=0,nnj+1
        do k=0,nnk+1   
            pres(i, j, k)=Pgrd/jac(i,j,k)
c           pres(i, j, k)=p0*(dcos(omega*time))/jac(i,j,k)
c	   pres(i,j,k) = 0.0;
	enddo
	enddo
	enddo
c        write(*,*) 'time = ', time
c        write(*,*) 'pressure = ', pres(10, 10, 10)
        return
	end 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine btmAugmentation
        
	include "size.inc"
        include "para.inc"
	include "metric.inc"
	include "ns.inc"
	include "eddy.inc"
	include "sedi.inc"

        double precision hc, C, absU, az, dz, temp, sum, bz, ddz, dh
	integer i, j, k, jj
	
        hc = 2.D0*(xp(2, 1, 1, 1) - xp(1, 1, 1, 1))
        C = 0.5D0

        do k = 0, nnk+1
	do j = 0, nnj+1
	do i = 0, nni+1
c           dz = xp(i,j,k,2) - xp(i,0,k,2)
           dz = xp(i,j,k,2) + depth(i,k)
	   dh = (xp(i,j+1,k,2) + xp(i,j,k,2))/2.D0 + depth(i,k)
	   ddz = 0.01D0
	    absU = dsqrt(  u(i,j,k,1)**2 + u(i,j,k,2)**2
     <                     + u(i,j,k,3)**2 )
	   if (dz .le. hc) then
	     
              az = (dcos(pi*dz/(2.D0*hc)))**2
	      bz = ((dcos(pi*dh/(2.D0*hc)))**2)
c                  *dtanh((dh)/ddz)
              btau(i,j,k,1) = -C*az*absU*u(i,j,k,1)
	      btau(i,j,k,2) = -C*az*absU*u(i,j,k,2) 
c	      btau(i,j,k,2) = -C*az*absU*u(i,j,k,2)
	      btau(i,j,k,3) = -C*az*absU*u(i,j,k,3)
	      
c	      if (j .eq. 0) then
c		 bchi(i,j,k)   =  C*bz*absU*sedc(i,j,k)
c     <                        *(xp(i,j+1,k,2)-xp(i,j,k,2))
c		 bchi(i,j,k)   = 0.D0
c	      else
c	 	 bchi(i,j,k)   =  C*bz*dabs(u(i,j,k,2))*sedc(i,j,k)
c     <                          *(xp(i,j+1,k,2)-xp(i,j,k,2))
c     <                          + bchi(i,j-1,k) 
c	      endif
c	      akeq(i,j,k) = (az*vst(i,j,k)/0.5D0
c	      write(*,*) akeq(i,j,k)
	      
	      
	   else
	      btau(i,j,k,1) = 0.D0
	      btau(i,j,k,2) = 0.D0
	      btau(i,j,k,3) = 0.D0
	      bchi(i,j,k)   = 0.D0
	   endif
c	   absU = dsqrt(  uxi(i,j,k)**2.D0 + uej(i,j,k)**2.D0
c     <                     + uzk(i,j,k)**2.D0 )
	   if (j .eq. 0) then
	      bchi(i,j,k) = 0.D0
	   else
	      bchi(i,j,k) = 0.15D0*dexp(-(xp(i,j,k,2)+depth(i,k))
     <	   /(1.5D0*hc))*ushr(i,k)*sedc(i,j,k)
c     <                      *(xp(i,j+1,k,2)-xp(i,j,k,2))
c     <                      + bchi(i,j-1,k)
	   endif
	   
           temp = 1.D0/jac(i,j,k)
           btau(i,j,k,1) = btau(i,j,k,1)*temp
           btau(i,j,k,2) = btau(i,j,k,2)*temp
	   btau(i,j,k,3) = btau(i,j,k,3)*temp
	   bchi(i,j,k)   = bchi(i,j,k)*ety(i,j,k)
c	   if (k .eq. 8) then
c	   write(*,*) i,j,bchi(i,j,k), sedc(i,j,k)
c	   endif
	enddo
	enddo
	enddo

c        do k = 1, nnk
c	do j = 1, nnj
c	do i = 1, nni
c	   taus(i,j,k,4) = taus(i,j,k,4)
c     <                   + btau(i,j,k,1)*0.5D0
c     <                   *(xp(i,j+1,k,2)-xp(i,j-1,k,2))
c	enddo
c	enddo
c	enddo
	
	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine ustar_bc

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "ns.inc"
	include "sedi.inc"
	include "para.inc"
	integer i, j, k, L
	
	if ( n_west .eq. MPI_PROC_NULL ) then
	   do L = 1, 3
	   do k = -1, nnk+2
	   do j = -1, nnj+2
c	      u( 0,j,k,L) = 3.D0*( u(1,j,k,L)-u(2,j,k,L) ) + u(3,j,k,L)
	      u( 0,j,k,L) = u(1,j,k,L)
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
c	      u(nni+1,j,k,L) = 3.D0*( u(nni,j,k,L)-u(nni-1,j,k,L) ) 
c     <                       + u(nni-2,j,k,L)
	      u(nni+1,j,k,L) = u(nni,j,k,L) 
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
c	      if (logU .eq. 1) then
c		u(i, 0,k,L) = u(i, 1,k,L) - fct(i, k)*u(i, 1,k,L)
c	      else
  	        u(i, 0,k,L) = 3.D0*( u(i,1,k,L)-u(i,2,K,L) ) + u(i,3,k,L)
c		u(i, 0,k,L) = u(i,1,k,L)
c	      endif      
c	      write(*,*) fct(i,k)
	   enddo
	   enddo
	   enddo

	   if (movingGrid .eq. 1) then
	   do k = -1, nnk+2
	   do i = -1, nni+2
c	      if (logU .eq. 1) then
c		u(i, 0,k,L) = u(i, 1,k,L) - fct(i, k)*u(i, 1,k,L)
c	      else
c  	        u(i, 0,k,2) = w_bed(i,k)
c		u(i, 0,k,L) = u(i,1,k,L)
c	      endif      
c	      write(*,*) fct(i,k)
	   enddo
	   enddo
	   endif

	   do L = 1, 3
	   do k = -1, nnk+2
	   do i = -1, nni+2
c	      u(i,-1,k,L) = 3.D0*( u(i,0,k,L)-u(i,1,K,L) ) + u(i,2,k,L)
	      u(i,-1,k,L) = u(i,0,k,L)
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
ccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine bottom_mean_stress

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "ns.inc"
	include "metric.inc"
	include "para.inc"
	include "eddy.inc"
	include "sedi.inc"

	integer i,j,k
	double precision dy, utmp, tau_id
	
	tau_id = 0.0
	j = 1
	dy = depth(1,1) + xp(1,1,1,2)
	tau_mean = 0.0
	
	if (n_suth .eq. MPI_PROC_NULL) then
	   do k = 1, nnk
	   do i = 1, nni
	      utmp = dsqrt(u(i,1,k,1)**2.0+u(i,1,k,3)**2.0)
	      tau_id = tau_id + utmp
	   enddo
	   enddo
	   tau_id = tau_id/dble(nni*nnk)
	   tau_id = tau_id/dy*vis*1000.0
	else
	   tau_id = 0.0
	endif
	

	call MPI_REDUCE(tau_id, tau_mean, 1, MPI_DOUBLE_PRECISION,
     <                  MPI_SUM, 0, comm3d, ierr)

	call MPI_BCAST(tau_mean, 1, MPI_DOUBLE_PRECISION,0,
     <                  comm3d,ierr)
	
	tau_mean = tau_mean/dble(px*pz)
	return
	end
