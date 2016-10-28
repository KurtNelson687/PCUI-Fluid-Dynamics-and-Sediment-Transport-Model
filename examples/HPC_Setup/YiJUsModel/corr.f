cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine corrector

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "para.inc"
	include "ns.inc"
	include "metric.inc"
	include "sedi.inc"

	integer i, j, k, m, L
	double precision pe, pw, pn, ps, pf, pb, tt, tp
	double precision tmp

	tt = dtime * 0.5D0
	
c	if (open_top .eq. 1) then
c	call openP_bc_nrth
c	endif

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
	   
c	   if (j .eq. 1) then
c	      write(*,*) p(10,j-1,10), p(10,j,10),p(10,j+1,10),p(10,j+2,10)
c	      ps = 3.D0*p(i,j,k)-1.D0*p(i,j+1,k)

c	   endif
	   tmp = 1.D0/(1.D0*(1.D0-sedc(i,j,k) + spwght*sedc(i,j,k)))

	   u(i,j,k,1) = u(i,j,k,1) 
     <	              - tp * ( xix(i,j,k) * pe - xix(i-1,j,k) * pw
     <	                     + etx(i,j,k) * pn - etx(i,j-1,k) * ps 
     <	                     + ztx(i,j,k) * pf - ztx(i,j,k-1) * pb )
     <                *tmp      	   
	   u(i,j,k,2) = u(i,j,k,2) 
     <	              - tp * ( xiy(i,j,k) * pe - xiy(i-1,j,k) * pw 
     <	                     + ety(i,j,k) * pn - ety(i,j-1,k) * ps 
     <	                     + zty(i,j,k) * pf - zty(i,j,k-1) * pb )
     <                *tmp      	   
	   u(i,j,k,3) = u(i,j,k,3) 
     <	              - tp * ( xiz(i,j,k) * pe - xiz(i-1,j,k) * pw 
     <	                     + etz(i,j,k) * pn - etz(i,j-1,k) * ps 
     <	                     + ztz(i,j,k) * pf - ztz(i,j,k-1) * pb )
     <                *tmp      	   	   

	enddo
	enddo
	enddo

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
     <        *pcorr(i,j,k,1)
	enddo
	enddo
	enddo

	do k = 1, nnk
	do j = jus, jue
	do i = 1, nni
	   uej(i,j,k) = uej(i,j,k) - dtime *
     <	      ( g22(i,j,k) * ( p(i,  j+1,k  ) - p(i,  j,  k  ) ) 
     <	      + g23(i,j,k) * ( p(i,  j,  k+1) - p(i,  j,  k-1)
     <	                     + p(i,  j+1,k+1) - p(i,  j+1,k-1) )
     <	      + g21(i,j,k) * ( p(i+1,j,  k  ) - p(i-1,j,  k  )
     <	                     + p(i+1,j+1,k  ) - p(i-1,j+1,k  ) ) )
     <        *pcorr(i,j,k,2)
	enddo
	enddo
	enddo
c	write(*,*) 'jus = ', jus
	do k = kus, kue
	do j = 1, nnj
	do i = 1, nni
	   uzk(i,j,k) = uzk(i,j,k) - dtime * 
     <	      ( g33(i,j,k) * ( p(i,  j,  k+1) - p(i,  j,  k  ) ) 
     <	      + g31(i,j,k) * ( p(i+1,j,  k  ) - p(i-1,j,  k  )
     <	                     + p(i+1,j,  k+1) - p(i-1,j,  k+1) )
     < 	      + g32(i,j,k) * ( p(i,  j+1,k  ) - p(i,  j-1,k  )
     <	                     + p(i,  j+1,k+1) - p(i,  j-1,k+1) ) )
     <        *pcorr(i,j,k,3)
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
	include "para.inc"
        include "sedi.inc"
        include "eddy.inc"

	integer i, j, k, L

	if ( n_nrth .eq. MPI_PROC_NULL ) then
	   do k = -1, nnk+2
	   do i = -1, nni+2
c	      u(i,nnj+1,k,1) = 2.D0 * u_lid(i,k) - u(i,nnj,k,1)
c	      u(i,nnj+1,k,2) =                   - u(i,nnj,k,2)
c	      u(i,nnj+1,k,3) = 2.D0 * w_lid(i,k) - u(i,nnj,k,3)
c              u(i,nnj+1,k,1) =  u(i,nnj,k,1)
              u(i,nnj+1,k,1) =  u(i,nnj,k,1)
	      u(i,nnj+1,k,2) =  u(i,nnj,k,2)
              u(i,nnj+1,k,3) =  u(i,nnj,k,3)
	   enddo
	   enddo
	   do L = 1, 3
	   do k = -1, nnk+2
	   do i = -1, nni+2
	      u(i,nnj+2,k,L) = 3.D0 * ( u(i,nnj+1,k,L)-u(i,nnj,k,L) ) 
     <                       + u(i,nnj-1,k,L)
c              u(i, nnj+1, k, L) = u(i, nnj, k, L)
c              u(i, nnj+2, k, L) = u(i, nnj+1, k, L)
	   enddo
	   enddo
	   enddo
	endif

	if ( n_suth .eq. MPI_PROC_NULL ) then
	   do L = 1, 3
	   do k = -1, nnk+2
	   do i = -1, nni+2
c	      if (logU .eq. 1) then
c                fct(i,k) =  Cd(i,k)*dsqrt(u(i,1,k,1)**2.D0  
c     <	            + u(i,1,k,2)**2.D0 + u(i,1,k,3)**2.D0)
c     <              *(xp(i,1,k,2)-xp(i,0,k,2))/(vis + vst(i,j,k))
c                u(i, 0,k,L) = (1.D0-fct(i,k))*u(i, 1,k,L) 		 
c	      else
	      
c	      if (L .ne. 2) then
c	        u(i, 0,k,L) = u(i,1,k,L)
c	      else
c	        u(i, 0,k,L) = w_bed(i,k)
c	      endif
c	      u(i,0,k,1) =  u(i,1,k,1)
c              u(i,0,k,L) =  -u(i,1,k,L)
c              u(i,0,k,3) =  u(i,1,k,3)
c	      if (L .eq. 2) then
c		 u(i, 0,k,L) = u(i, 0,k,L)
c	      else
c                 u(i, 0,k,L) = 2.D0*u(i,1,k,L) - u(i,2,k,L)
c	      endif
	      if (movingGrid .eq. 1) then
		    u(i,0,k,L) = 2.0*u_bed(i,k,L)-u(i,1,k,L)
	      else
		    u(i,0,k,L) = -u(i,1,k,L)
	      endif
	   enddo
	   enddo
	   enddo
	   do L = 1, 3
	   do k = -1, nnk+2
	   do i = -1, nni+2
c	      u(i,-1,k,L) = 3.D0*( u(i,0,k,L)-u(i,1,K,L) ) + u(i,2,k,L)
c              u(i,-1,k,L) = u(i,0,k,L)
c              if (L .ne. 1) then
c                 if (logU .eq. 1) then
c                   u(i,-1,k,L) = u(i,0,k,L) + u(i,0,k,L)-u(i,1,k,L)
c	         else
                   u(i,-1,k,L) = 3.D0*
     <                           ( u(i,0,k,L)-u(i,1,K,L) ) + u(i,2,k,L)
c                 endif
c	      else
	         
	   enddo
	   enddo
	   enddo
	endif


	if ( n_west .eq. MPI_PROC_NULL ) then
	   do k = -1, nnk+2
	   do j = -1, nnj+2
c	      u( 0,j,k,1) = 2.D0*u_flux(j,k) - u(1,j,k,1)
c	      u( 0,j,k,2) = - u(1,j,k,2)
c	      u( 0,j,k,3) = - u(1,j,k,3)
c              u(0, j, k, 1) = u(1, j, k, 1)
c              u(0, j, k, 1) = 1.D0*amp*omega*dcos(omega*time)
              u(0, j, k, 1) = u(1, j, k, 1)
              u(0, j, k, 2) = u(1, j, k, 2)
              u(0, j, k, 3) = u(1, j, k, 3)
	   enddo
	   enddo
	   do L = 1, 3
	   do k = -1, nnk+2
	   do j = -1, nnj+2
	      u(-1,j,k,L) = 3.D0*( u(0,j,k,L)-u(1,j,k,L) ) + u(2,j,k,L)
c              u(-1, j, k, L) = u(0, j, k, L)
	   enddo
	   enddo
	   enddo
	endif

	if ( n_east .eq. MPI_PROC_NULL ) then
	   do k = -1, nnk+2
	   do j = -1, nnj+2
c             u(nni+1, j, k, 1) = u(nni, j, k, 1)
c             u(nni+1, j, k, 1) = 1.D0*amp*omega*dcos(omega*time) 
              u(nni+1, j, k, 1) = u(nni, j, k, 1)              
              u(nni+1, j, k, 2) = u(nni, j, k, 2)
              u(nni+1, j, k, 3) = u(nni, j, k, 3)
c	      u(nni+1,j,k,1) = 2.D0*u_flux(j,k) - u(nni,j,k,1)
c	      u(nni+1,j,k,2) = - u(nni,j,k,2)
c	      u(nni+1,j,k,3) = - u(nni,j,k,3)
	   enddo
	   enddo
	   do L = 1, 3
	   do k = -1, nnk+2
	   do j = -1, nnj+2
	      u(nni+2,j,k,L) = 3.D0*( u(nni+1,j,k,L)-u(nni,j,k,L) ) 
     <                       + u(nni-1,j,k,L)
c              u(nni+2, j, k, L) = u(nni+1, j, k, L)
	   enddo
	   enddo
	   enddo
	endif

	if ( n_back .eq. MPI_PROC_NULL ) then
	   do L = 1, 3
	   do j = -1, nnj+2
	   do i = -1, nni+2
              u(i,j, 0,L) = u(i,j,1,L)
c	      u(i,j, 0,L) = - u(i,j,1,L)
	   enddo
	   enddo
	   enddo
	   do L = 1, 3
	   do j = -1, nnj+2
	   do i = -1, nni+2
              u(i, j, -1, L) = u(i, j, 0, L)
c	      u(i,j,-1,L) = 3.D0*( u(i,j,0,L)-u(i,j,1,L) ) + u(i,j,2,L)
	   enddo
	   enddo
	   enddo
	endif

	if ( n_frnt .eq. MPI_PROC_NULL ) then
	   do L = 1, 3
	   do j = -1, nnj+2
	   do i = -1, nni+2
              u(i,j,nnk+1,L) = u(i,j,nnk,L)
c	      u(i,j,nnk+1,L) = - u(i,j,nnk,L)
	   enddo
	   enddo
	   enddo
	   do L = 1, 3
	   do j = -1, nnj+2
	   do i = -1, nni+2
              u(i, j, nnk+2, L) = u(i, j, nnk+1, L)
c	      u(i,j,nnk+2,L) = 3.D0*( u(i,j,nnk+1,L)-u(i,j,nnk,L) ) 
c     <                       + u(i,j,nnk-1,L)
	   enddo
	   enddo
	   enddo
	endif

	return
	end
	
