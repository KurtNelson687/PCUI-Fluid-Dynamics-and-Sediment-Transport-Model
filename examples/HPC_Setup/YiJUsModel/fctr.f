cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine thrd_f2c(nif, njf, nkf, nic, njc, nkc, pf, pc)

	INCLUDE "size.inc"
	include "mpif.h"
	include "mpi.inc"

	integer nif, njf, nkf, nic, njc, nkc
	double precision pf(0:nif+1,0:njf+1,0:nkf+1)
	double precision pc(0:nic+1,0:njc+1,0:nkc+1)

	integer i, j, k, ii, jj, kk
	
	integer iib, jjb, kkb
	integer iie, jje, kke

	iib = 0
	jjb = 0
	kkb = 0
	iie = 0
	jje = 0
	kke = 0

	do k = 1, nkc
	   kk = 2 * k - 1
	do j = 1, njc
	   jj = 2 * j - 1
	do i = 1, nic
	   ii = 2 * i - 1
	   pc(i,j,k) = pf(ii,  jj,  kk  ) + pf(ii+1,jj,  kk  )
     <	             + pf(ii,  jj+1,kk  ) + pf(ii+1,jj+1,kk  )
     <	             + pf(ii,  jj,  kk+1) + pf(ii+1,jj,  kk+1)
     <	             + pf(ii,  jj+1,kk+1) + pf(ii+1,jj+1,kk+1)
 	enddo
	enddo
	enddo

	if ( n_west .eq. MPI_PROC_NULL ) then
	   do k = 1, nkc
	      kk = 2 * k - 1
	   do j = 1, njc
	      jj = 2 * j - 1
	      pc(0,j,k) = pf(0, jj,  kk  ) 
     <                  + pf(0, jj+1,kk  ) 
     <	                + pf(0, jj,  kk+1) 
     <                  + pf(0, jj+1,kk+1) 
 	   enddo
	   enddo
	endif

        if ( n_east .eq. MPI_PROC_NULL ) then
	   do k = 1, nkc
	      kk = 2 * k - 1
	   do j = 1, njc
	      jj = 2 * j - 1
	      pc(nic+1,j,k) = pf(nif+1, jj,   kk  ) 
     <	                    + pf(nif+1, jj+1, kk  ) 
     <	                    + pf(nif+1, jj,   kk+1)
     <	                    + pf(nif+1, jj+1, kk+1) 
	   enddo
	   enddo
	endif

        if ( n_suth .eq. MPI_PROC_NULL ) then
	   do k = 1, nkc
	      kk = 2 * k - 1
	   do i = 1, nic
	      ii = 2 * i - 1
	      pc(i,0,k) = pf(ii,   0, kk  ) 
     <	                + pf(ii+1, 0, kk  )
     <	                + pf(ii,   0, kk+1) 
     <	                + pf(ii+1, 0, kk+1)
	   enddo
	   enddo
	endif

        if ( n_nrth .eq. MPI_PROC_NULL ) then
	   do k = 1, nkc
	      kk = 2 * k - 1
	   do i = 1, nic
	      ii = 2 * i - 1
	      pc(i,njc+1,k) = pf(ii,   njf+1, kk  ) 
     <	                    + pf(ii+1, njf+1, kk  )
     <	                    + pf(ii,   njf+1, kk+1) 
     <	                    + pf(ii+1, njf+1, kk+1)
	   enddo
	   enddo
	endif

        if ( n_back .eq. MPI_PROC_NULL ) then
	   do j = 1, njc
	      jj = 2 * j - 1
	   do i = 1, nic
	      ii = 2 * i - 1
	      pc(i,j,0) = pf(ii,  jj,  0) 
     <	                + pf(ii+1,jj,  0)
     <	                + pf(ii,  jj+1,0) 
     <	                + pf(ii+1,jj+1,0) 
	   enddo
	   enddo
	endif
        if ( n_frnt .eq. MPI_PROC_NULL ) then
	   do j = 1, njc
	      jj = 2 * j - 1
	   do i = 1, nic
	      ii = 2 * i - 1
	      pc(i,j,nkc+1) = pf(ii,  jj,  nkf+1) 
     <	                    + pf(ii+1,jj,  nkf+1)
     <	                    + pf(ii,  jj+1,nkf+1) 
     <	                    + pf(ii+1,jj+1,nkf+1) 
	   enddo
	   enddo
	endif

	return
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine thrd_c2f(nic, njc, nkc, nif, njf, nkf, pc, pf)

	implicit none

	integer nic, njc, nkc, nif, njf, nkf
	double precision pc(0:nic+1,0:njc+1,0:nkc+1)
	double precision pf(0:nif+1,0:njf+1,0:nkf+1)

	double precision p000, p100, p010, p001, p011, p101, p110, p111
	integer i, j, k, ii, jj, kk

	call pres_corner(nic, njc, nkc, pc)

	do k = 0, nkc
	   kk = 2 * k
	do j = 0, njc
	   jj = 2 * j
	do i = 0, nic
	   ii = 2 * i

	   p000 = pc(i,  j,  k  )
	   p100 = pc(i+1,j,  k  ) 
	   p010 = pc(i,  j+1,k  ) 
	   p001 = pc(i,  j,  k+1)  
	   p011 = pc(i,  j+1,k+1)  
	   p101 = pc(i+1,j,  k+1)  
	   p110 = pc(i+1,j+1,k  )
	   p111 = pc(i+1,j+1,k+1) 

C......	   000
           pf(ii,  jj,  kk  ) = pf(ii,  jj,  kk  ) + 0.015625D0 * 
     <			( 27.D0 * p000 
     <			+  9.D0 * ( p100 + p010 + p001 ) 
     <			+  3.D0 * ( p011 + p101 + p110 ) 
     <			+         p111 )

C......	   100
           pf(ii+1,jj,  kk  ) = pf(ii+1,jj,  kk  ) + 0.015625D0 * 
     <			( 27.D0 * p100 
     <			+  9.D0 * ( p000 + p110 + p101 ) 
     <			+  3.D0 * ( p111 + p001 + p010 ) 
     <			+         p011 )

C......	   010
	   pf(ii,  jj+1,kk  ) = pf(ii,  jj+1,kk  ) + 0.015625D0 * 
     <			( 27.D0 * p010 
     <			+  9.D0 * ( p110 + p000 + p011 ) 
     <			+  3.D0 * ( p001 + p111 + p100 ) 
     <			+         p101 )

C......	   001
	   pf(ii,  jj,  kk+1) = pf(ii,  jj,  kk+1) + 0.015625D0 *
     <			( 27.D0 * p001 
     <			+  9.D0 * ( p101 + p011	+ p000 ) 
     <			+  3.D0 * ( p010 + p100 + p111 ) 
     <			+         p110 )

C......	   011
           pf(ii,  jj+1,kk+1) = pf(ii,  jj+1,kk+1) + 0.015625D0 * 
     <			( 27.D0 * p011 
     <			+  9.D0 * ( p111 + p001 + p010 ) 
     <			+  3.D0 * ( p000 + p110 + p101 ) 
     <			+         p100 )

C......	   101
	   pf(ii+1,jj,  kk+1) = pf(ii+1,jj,  kk+1) + 0.015625D0 * 
     <			( 27.D0 * p101 
     <			+  9.D0 * ( p001 + p111 + p100 ) 
     <			+  3.D0 * ( p110 + p000 + p011 ) 
     <			+         p010 )

C......	   110
	   pf(ii+1,jj+1,kk  ) = pf(ii+1,jj+1,kk  ) + 0.015625D0 * 
     <			( 27.D0 * p110 
     <			+  9.D0 * ( p010 + p100 + p111 ) 
     <			+  3.D0 * ( p101 + p011 + p000 ) 
     <			+         p001 )

C......	   111
	   pf(ii+1,jj+1,kk+1) = pf(ii+1,jj+1,kk+1) + 0.015625D0 * 
     <			( 27.D0 * p111 
     <			+  9.D0 * ( p011 + p101 + p110 ) 
     <			+  3.D0 * ( p100 + p010 + p001 ) 
     <			+         p000 )

	enddo
	enddo
	enddo

	return
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine pres_corner(ii, jj, kk, r)

	implicit none
	include "mpif.h"
	include "mpi.inc"

	integer ii, jj, kk
	double precision r(0:ii+1,0:jj+1,0:kk+1)

	integer i, j, k

	if ( n_west .eq. MPI_PROC_NULL ) then
	   if ( n_suth .eq. MPI_PROC_NULL ) then
	      if ( n_back .eq. MPI_PROC_NULL ) then
C......	0 - 0 - 0
	r(0,0,0) = r(0,1,1) + r(1,0,1) + r(1,1,0)
     >	         - 2.D0 * r(1,1,1)
	      endif
	      if ( n_frnt .eq. MPI_PROC_NULL ) then
C......	0 - 0 - Z
	r(0,0,kk+1) = r(0,1,kk) + r(1,0,kk) + r(1,1,kk+1) 
     <	            - 2.D0 * r(1,1,kk)
	      endif
	   endif
	   if ( n_nrth .eq. MPI_PROC_NULL ) then
	      if ( n_back .eq. MPI_PROC_NULL ) then
C......	0 - Y - 0
	r(0,jj+1,0) = r(0,jj,1) + r(1,jj+1,1) + r(1,jj,0)
     <	            - 2.D0 * r(1,jj,1)
	      endif
	      if ( n_frnt .eq. MPI_PROC_NULL ) then
C......	0 - Y - Z
	r(0,jj+1,kk+1) = r(0,jj,kk) + r(1,jj+1,kk) + r(1,jj,kk+1)
     <	               - 2.D0 * r(1,jj,kk)
	      endif
	   endif
	endif	      

	if ( n_east .eq. MPI_PROC_NULL ) then
	   if ( n_suth .eq. MPI_PROC_NULL ) then
	      if ( n_back .eq. MPI_PROC_NULL ) then
C......	X - 0 - 0
	r(ii+1,0,0) = r(ii+1,1,1) + r(ii,0,1) + r(ii,1,0)
     <	            - 2.D0 * r(ii,1,1)
	      endif
	      if ( n_frnt .eq. MPI_PROC_NULL ) then
C......	X - 0 - Z
	r(ii+1,0,kk+1) = r(ii+1,1,kk) + r(ii,0,kk) + r(ii,1,kk+1)
     <	               - 2.D0 * r(ii,1,kk)
	      endif
	   endif
	   if ( n_nrth .eq. MPI_PROC_NULL ) then
	      if ( n_back .eq. MPI_PROC_NULL ) then
C......	X - Y - 0
	r(ii+1,jj+1,0) = r(ii+1,jj,1) + r(ii,jj+1,1) + r(ii,jj,0)
     <	               - 2.D0 * r(ii,jj,1)
	      endif
	      if ( n_frnt .eq. MPI_PROC_NULL ) then
C......	X - Y - Z
	r(ii+1,jj+1,kk+1)= r(ii+1,jj,kk) + r(ii,jj+1,kk) + r(ii,jj,kk+1)
     <	                 - 2.D0 * r(ii,jj,kk)
	      endif
	   endif
	endif	
	
	return
	end
