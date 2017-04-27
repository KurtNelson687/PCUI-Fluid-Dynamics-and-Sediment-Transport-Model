cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine corrector

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "para.inc"
	include "ns.inc"
	include "metric.inc"

	integer i, j, k, m, L
	double precision pe, pw, pn, ps, pf, pb, tt, tp

	tt = dtime * 0.5D0

	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni

	   tp  = tt * jac(i,j,k) 
	   pe = p(i+1,j,k) + p(i,j,k)
	   pw = p(i-1,j,k) + p(i,j,k)
	   pn = p(i,j+1,k) + p(i,j,k)
	   ps = p(i,j-1,k) + p(i,j,k)
	   pf = p(i,j,k+1) + p(i,j,k)
	   pb = p(i,j,k-1) + p(i,j,k)

	   u(i,j,k,1) = u(i,j,k,1) 
     <	              - tp * ( xix(i,j,k) * pe - xix(i-1,j,k) * pw
     <	                     + etx(i,j,k) * pn - etx(i,j-1,k) * ps 
     <	                     + ztx(i,j,k) * pf - ztx(i,j,k-1) * pb )

	   u(i,j,k,2) = u(i,j,k,2) 
     <	              - tp * ( xiy(i,j,k) * pe - xiy(i-1,j,k) * pw 
     <	                     + ety(i,j,k) * pn - ety(i,j-1,k) * ps 
     <	                     + zty(i,j,k) * pf - zty(i,j,k-1) * pb )

	   u(i,j,k,3) = u(i,j,k,3) 
     <	              - tp * ( xiz(i,j,k) * pe - xiz(i-1,j,k) * pw 
     <	                     + etz(i,j,k) * pn - etz(i,j-1,k) * ps 
     <	                     + ztz(i,j,k) * pf - ztz(i,j,k-1) * pb )

	enddo
	enddo
	enddo

	if (pAdjust .eq. 1) then
	call adjustS
	call adjust_u
	endif
	call u_bc
	call u_exchange

	do k = 1, nnk
	do j = 1, nnj
	do i = ius, iue
	   uxi(i,j,k) = uxi(i,j,k) - dtime *
     <	      ( g11(i,j,k) * ( p(i+1,j,  k  ) - p(i,  j,  k  ) ) 
     <	      + g12(i,j,k) * ( p(i,  j+1,k  ) - p(i,  j-1,k  )
     <	                     + p(i+1,j+1,k  ) - p(i+1,j-1,k  ) )
     <	      + g13(i,j,k) * ( p(i,  j,  k+1) - p(i,  j,  k-1)
     <                       + p(i+1,j,  k+1) - p(i+1,j,  k-1) ) )
	enddo
	enddo
	enddo


	if (pAdjust .eq. 1) then
	call adjust_uxi
	endif

	do k = 1, nnk
	do j = jus, jue
	do i = 1, nni
	   uej(i,j,k) = uej(i,j,k) - dtime *
     <	      ( g22(i,j,k) * ( p(i,  j+1,k  ) - p(i,  j,  k  ) ) 
     <	      + g23(i,j,k) * ( p(i,  j,  k+1) - p(i,  j,  k-1)
     <	                     + p(i,  j+1,k+1) - p(i,  j+1,k-1) )
     <	      + g21(i,j,k) * ( p(i+1,j,  k  ) - p(i-1,j,  k  )
     <	                     + p(i+1,j+1,k  ) - p(i-1,j+1,k  ) ) )
	enddo
	enddo
	enddo

	do k = kus, kue
	do j = 1, nnj
	do i = 1, nni
	   uzk(i,j,k) = uzk(i,j,k) - dtime * 
     <	      ( g33(i,j,k) * ( p(i,  j,  k+1) - p(i,  j,  k  ) ) 
     <	      + g31(i,j,k) * ( p(i+1,j,  k  ) - p(i-1,j,  k  )
     <	                     + p(i+1,j,  k+1) - p(i-1,j,  k+1) )
     < 	      + g32(i,j,k) * ( p(i,  j+1,k  ) - p(i,  j-1,k  )
     <	                     + p(i,  j+1,k+1) - p(i,  j-1,k+1) ) )
	enddo
	enddo
	enddo

	return
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine u_bc

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "ns.inc"
	include "metric.inc"
	
	integer i, j, k, L

	if (n_nrth .eq. MPI_PROC_NULL ) then
	   do k = -1, nnk+2
	   do i = -1, nni+2
c             Free-slip
	      u(i,nnj+1,k,1) =   u(i,nnj,k,1)
	      u(i,nnj+1,k,2) =  -u(i,nnj,k,2)
	      u(i,nnj+1,k,3) =   u(i,nnj,k,3)

c             No-slip
C	      u(i,nnj+1,k,1) = -  u(i,nnj,k,1)
C	      u(i,nnj+1,k,2) = - u(i,nnj,k,2)
C	      u(i,nnj+1,k,3) = -  u(i,nnj,k,3)
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

	if ( n_suth .eq. MPI_PROC_NULL ) then
c	   do L = 1, 3
	   do k = -1, nnk+2
	   do i = -1, nni+2
c             Free-slip
c	      u(i,0,k,1) =   u(i,1,k,1)
c	      u(i,0,k,2) = - u(i,1,k,2)
c	      u(i,0,k,3) =   u(i,1,k,3)

c             No-slip
	      u(i, 0,k,1) = - u(i,1,k,1)
	      u(i, 0,k,2) = - u(i,1,k,2) 
	      u(i, 0,k,3) = - u(i,1,k,3)
	   enddo
	   enddo
c	   enddo
	   do L = 1, 3
	   do k = -1, nnk+2
	   do i = -1, nni+2
	      u(i,-1,k,L) = 3.D0*( u(i,0,k,L)-u(i,1,K,L) ) + u(i,2,k,L)
	   enddo
	   enddo
	   enddo
	endif

	if ( n_west .eq. MPI_PROC_NULL ) then
c	   do L = 1, 3
	   do k = -1, nnk+2
	   do j = -1, nnj+2
c             Free-slip
	      u(0,j,k,1) =   u(1,j,k,1)
	      u(0,j,k,2) = - u(1,j,k,2)
	      u(0,j,k,3) =   u(1,j,k,3)

c             No-slip
c	      u( 0,j,k,L) = - u(1,j,k,L)
	   enddo
	   enddo
c	   enddo
	   do L = 1, 3
	   do k = -1, nnk+2
	   do j = -1, nnj+2
	      u(-1,j,k,L) = 3.D0*( u(0,j,k,L)-u(1,j,k,L) ) + u(2,j,k,L)
	   enddo
	   enddo
	   enddo
	endif

	if ( n_east .eq. MPI_PROC_NULL ) then
c	   do L = 1, 3
	   do k = -1, nnk+2
	   do j = -1, nnj+2
c             Free-slip
	      u(nni+1,j,k,1) =   u(nni,j,k,1)
	      u(nni+1,j,k,2) = - u(nni,j,k,2)
	      u(nni+1,j,k,3) =   u(nni,j,k,3)

c             No-slip
c	      u(nni+1,j,k,L) = - u(nni,j,k,L)
	   enddo
	   enddo
c	   enddo
	   do L = 1, 3
	   do k = -1, nnk+2
	   do j = -1, nnj+2
	      u(nni+2,j,k,L) = 3.D0*( u(nni+1,j,k,L)-u(nni,j,k,L) ) 
     <                       + u(nni-1,j,k,L)
	   enddo
	   enddo
	   enddo
	endif

	if ( n_back .eq. MPI_PROC_NULL ) then
c	   do L = 1, 3
	   do j = -1, nnj+2
	   do i = -1, nni+2
c             Free-slip
	      u(i,j,0,1) =   u(i,j,1,1)
	      u(i,j,0,2) = - u(i,j,1,2)
	      u(i,j,0,3) =   u(i,j,1,3)

c             No-slip
c	      u(i,j, 0,L) = - u(i,j,1,L)
	   enddo
	   enddo
c	   enddo
	   do L = 1, 3
	   do j = -1, nnj+2
	   do i = -1, nni+2
	      u(i,j,-1,L) = 3.D0*( u(i,j,0,L)-u(i,j,1,L) ) + u(i,j,2,L)
	   enddo
	   enddo
	   enddo
	endif

	if ( n_frnt .eq. MPI_PROC_NULL ) then
c	   do L = 1, 3
	   do j = -1, nnj+2
	   do i = -1, nni+2
c             Free-slip
	      u(i,j,nnk+1,1) =   u(i,j,nnk,1)
	      u(i,j,nnk+1,2) = - u(i,j,nnk,2)
	      u(i,j,nnk+1,3) =   u(i,j,nnk,3)

c             No-slip
c	      u(i,j,nnk+1,L) = - u(i,j,nnk,L)
	   enddo
	   enddo
c	   enddo
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
	
