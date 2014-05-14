cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine u_exchange

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "ns.inc"
	include "metric.inc"

	double precision, dimension(-1:nnj+2,1:2,1:3) :: 
     <		rWmsg, sWmsg, rEmsg, sEmsg
	double precision, dimension(-1:nni+2,1:2,1:3) :: 
     <		rSmsg, sSmsg, rNmsg, sNmsg
	double precision, dimension(-1:nni+2,1:2,1:3) :: 
     <		rBmsg, sBmsg, rFmsg, sFmsg

	integer status(MPI_STATUS_SIZE,4), req(4), nreq, len
	
	integer i, j, k, L

C......	i-direction

	len = ( nnj + 4 ) * 6

	do k = -1, nnk+2

	nreq = 0

	if ( n_west .ne. MPI_PROC_NULL ) then
	   nreq = nreq + 1
	   call MPI_IRECV( rWmsg(-1,1,1), len, 
     <	                   MPI_DOUBLE_PRECISION,  
     <	                   n_west, 0, comm3d, req(nreq), ierr )
	   do L = 1, 3
	   do j = -1, nnj+2
	      sWmsg(j,1,L) = u(1,j,k,L)
	      sWmsg(j,2,L) = u(2,j,k,L)
	   enddo
	   enddo
	   nreq = nreq + 1
	   call MPI_ISEND( sWmsg(-1,1,1), len, 
     <	                   MPI_DOUBLE_PRECISION,  
     <	                   n_west, 1, comm3d, req(nreq), ierr )
	endif
	   
	if ( n_east .ne. MPI_PROC_NULL ) then
	   nreq = nreq + 1
	   call MPI_IRECV( rEmsg(-1,1,1), len, 
     <	                   MPI_DOUBLE_PRECISION,  
     <	                   n_east, 1, comm3d, req(nreq), ierr )
	   do L = 1, 3
	   do j = -1, nnj+2
	      sEmsg(j,1,L) = u(nni-1,j,k,L)
	      sEmsg(j,2,L) = u(nni,  j,k,L)
	   enddo
	   enddo
	   nreq = nreq + 1
	   call MPI_ISEND( sEmsg(-1,1,1), len, 
     <	                   MPI_DOUBLE_PRECISION,  
     <	                   n_east, 0, comm3d, req(nreq), ierr )
	endif

	call MPI_WAITALL( nreq, req, status, ierr )

	if ( n_west .ne. MPI_PROC_NULL ) then
	   do L = 1, 3
	   do j = -1, nnj+2
	      u(-1,j,k,L) = rWmsg(j,1,L)
	      u( 0,j,k,L) = rWmsg(j,2,L)
	   enddo
	   enddo
	endif

	if ( n_east .ne. MPI_PROC_NULL ) then
	   do L = 1, 3
	   do j = -1, nnj+2
	      u(nni+1,j,k,L) = rEmsg(j,1,L)
	      u(nni+2,j,k,L) = rEmsg(j,2,L)
	   enddo
	   enddo
	endif

	enddo

C......	j-direction

	len = ( nni + 4 ) * 6

	do k = -1, nnk+2

	nreq = 0

	if ( n_suth .ne. MPI_PROC_NULL ) then
	   nreq = nreq + 1
	   call MPI_IRECV( rSmsg(-1,1,1), len, 
     <	                   MPI_DOUBLE_PRECISION,  
     <	                   n_suth, 0, comm3d, req(nreq), ierr )
	   do L = 1, 3
	   do i = -1, nni+2
	      sSmsg(i,1,L) = u(i,1,k,L)
	      sSmsg(i,2,L) = u(i,2,k,L)
	   enddo
	   enddo
	   nreq = nreq + 1
	   call MPI_ISEND( sSmsg(-1,1,1), len, 
     <	                   MPI_DOUBLE_PRECISION,  
     <	                   n_suth, 1, comm3d, req(nreq), ierr )
	endif
	   
	if ( n_nrth .ne. MPI_PROC_NULL ) then
	   nreq = nreq + 1
	   call MPI_IRECV( rNmsg(-1,1,1), len, 
     <	                   MPI_DOUBLE_PRECISION,  
     <	                   n_nrth, 1, comm3d, req(nreq), ierr )
	   do L = 1, 3
	   do i = -1, nni+2
	      sNmsg(i,1,L) = u(i,nnj-1,k,L)
	      sNmsg(i,2,L) = u(i,nnj,  k,L)
	   enddo
	   enddo
	   nreq = nreq + 1
	   call MPI_ISEND( sNmsg(-1,1,1), len, 
     <	                   MPI_DOUBLE_PRECISION,  
     <	                   n_nrth, 0, comm3d, req(nreq), ierr )
	endif

	call MPI_WAITALL( nreq, req, status, ierr )

	if ( n_suth .ne. MPI_PROC_NULL ) then
	   do L = 1, 3
	   do i = -1, nni+2
	      u(i,-1,k,L) = rSmsg(i,1,L)
	      u(i, 0,k,L) = rSmsg(i,2,L)
	   enddo
	   enddo
	endif

	if ( n_nrth .ne. MPI_PROC_NULL ) then
	   do L = 1, 3
	   do i = -1, nni+2
	      u(i,nnj+1,k,L) = rNmsg(i,1,L)
	      u(i,nnj+2,k,L) = rNmsg(i,2,L)
	   enddo
	   enddo
	endif

	enddo

C......	k-direction

	len = ( nni + 4 ) * 6

	do j = -1, nnj+2

	nreq = 0

	if ( n_back .ne. MPI_PROC_NULL ) then
	   nreq = nreq + 1
	   call MPI_IRECV( rBmsg(-1,1,1), len, 
     <	                   MPI_DOUBLE_PRECISION,  
     <	                   n_back, 0, comm3d, req(nreq), ierr )
	   do L = 1, 3
	   do i = -1, nni+2
	      sBmsg(i,1,L) = u(i,j,1,L)
	      sBmsg(i,2,L) = u(i,j,2,L)
	   enddo
	   enddo
	   nreq = nreq + 1
	   call MPI_ISEND( sBmsg(-1,1,1), len, 
     <	                   MPI_DOUBLE_PRECISION,  
     <	                   n_back, 1, comm3d, req(nreq), ierr )
	endif
	   
	if ( n_frnt .ne. MPI_PROC_NULL ) then
	   nreq = nreq + 1
	   call MPI_IRECV( rFmsg(-1,1,1), len, 
     <	                   MPI_DOUBLE_PRECISION,  
     <	                   n_frnt, 1, comm3d, req(nreq), ierr )
	   do L = 1, 3
	   do i = -1, nni+2
	      sFmsg(i,1,L) = u(i,j,nnk-1,L)
	      sFmsg(i,2,L) = u(i,j,nnk,  L)
	   enddo
	   enddo
	   nreq = nreq + 1
	   call MPI_ISEND( sFmsg(-1,1,1), len, 
     <	                   MPI_DOUBLE_PRECISION,  
     <	                   n_frnt, 0, comm3d, req(nreq), ierr )
	endif

	call MPI_WAITALL( nreq, req, status, ierr )

	if ( n_back .ne. MPI_PROC_NULL ) then
	   do L = 1, 3
	   do i = -1, nni+2
	      u(i,j,-1,L) = rBmsg(i,1,L)
	      u(i,j, 0,L) = rBmsg(i,2,L)
	   enddo
	   enddo
	endif

	if ( n_frnt .ne. MPI_PROC_NULL ) then
	   do L = 1, 3
	   do i = -1, nni+2
	      u(i,j,nnk+1,L) = rFmsg(i,1,L)
	      u(i,j,nnk+2,L) = rFmsg(i,2,L)
	   enddo
	   enddo
	endif

	enddo

	return
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine phi_exchange

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "ns.inc"
	include "metric.inc"

	double precision, dimension(-1:nnj+2,1:2) :: 
     <		rWmsg, sWmsg, rEmsg, sEmsg
	double precision, dimension(-1:nni+2,1:2) :: 
     <		rSmsg, sSmsg, rNmsg, sNmsg
	double precision, dimension(-1:nni+2,1:2) :: 
     <		rBmsg, sBmsg, rFmsg, sFmsg

	integer status(MPI_STATUS_SIZE,4), req(4), nreq, len
	
	integer i, j, k

C......	i-direction

	len = ( nnj + 4 ) * 2

	do k = -1, nnk+2

	nreq = 0

	if ( n_west .ne. MPI_PROC_NULL ) then
	   nreq = nreq + 1
	   call MPI_IRECV( rWmsg(-1,1), len, 
     <	                   MPI_DOUBLE_PRECISION,  
     <	                   n_west, 0, comm3d, req(nreq), ierr )
	   do j = -1, nnj+2
	      sWmsg(j,1) = phi(1,j,k)
	      sWmsg(j,2) = phi(2,j,k)
	   enddo
	   nreq = nreq + 1
	   call MPI_ISEND( sWmsg(-1,1), len, 
     <	                   MPI_DOUBLE_PRECISION,  
     <	                   n_west, 1, comm3d, req(nreq), ierr )
	endif
	   
	if ( n_east .ne. MPI_PROC_NULL ) then
	   nreq = nreq + 1
	   call MPI_IRECV( rEmsg(-1,1), len, 
     <	                   MPI_DOUBLE_PRECISION,  
     <	                   n_east, 1, comm3d, req(nreq), ierr )
	   do j = -1, nnj+2
	      sEmsg(j,1) = phi(nni-1,j,k)
	      sEmsg(j,2) = phi(nni,  j,k)
	   enddo
	   nreq = nreq + 1
	   call MPI_ISEND( sEmsg(-1,1), len, 
     <	                   MPI_DOUBLE_PRECISION,  
     <	                   n_east, 0, comm3d, req(nreq), ierr )
	endif

	call MPI_WAITALL( nreq, req, status, ierr )

	if ( n_west .ne. MPI_PROC_NULL ) then
	   do j = -1, nnj+2
	      phi(-1,j,k) = rWmsg(j,1)
	      phi( 0,j,k) = rWmsg(j,2)
	   enddo
	endif

	if ( n_east .ne. MPI_PROC_NULL ) then
	   do j = -1, nnj+2
	      phi(nni+1,j,k) = rEmsg(j,1)
	      phi(nni+2,j,k) = rEmsg(j,2)
	   enddo
	endif

	enddo

C......	j-direction

	len = ( nni + 4 ) * 2

	do k = -1, nnk+2

	nreq = 0

	if ( n_suth .ne. MPI_PROC_NULL ) then
	   nreq = nreq + 1
	   call MPI_IRECV( rSmsg(-1,1), len, 
     <	                   MPI_DOUBLE_PRECISION,  
     <	                   n_suth, 0, comm3d, req(nreq), ierr )
	   do i = -1, nni+2
	      sSmsg(i,1) = phi(i,1,k)
	      sSmsg(i,2) = phi(i,2,k)
	   enddo
	   nreq = nreq + 1
	   call MPI_ISEND( sSmsg(-1,1), len, 
     <	                   MPI_DOUBLE_PRECISION,  
     <	                   n_suth, 1, comm3d, req(nreq), ierr )
	endif
	   
	if ( n_nrth .ne. MPI_PROC_NULL ) then
	   nreq = nreq + 1
	   call MPI_IRECV( rNmsg(-1,1), len, 
     <	                   MPI_DOUBLE_PRECISION,  
     <	                   n_nrth, 1, comm3d, req(nreq), ierr )
	   do i = -1, nni+2
	      sNmsg(i,1) = phi(i,nnj-1,k)
	      sNmsg(i,2) = phi(i,nnj,  k)
	   enddo
	   nreq = nreq + 1
	   call MPI_ISEND( sNmsg(-1,1), len, 
     <	                   MPI_DOUBLE_PRECISION,  
     <	                   n_nrth, 0, comm3d, req(nreq), ierr )
	endif

	call MPI_WAITALL( nreq, req, status, ierr )

	if ( n_suth .ne. MPI_PROC_NULL ) then
	   do i = -1, nni+2
	      phi(i,-1,k) = rSmsg(i,1)
	      phi(i, 0,k) = rSmsg(i,2)
	   enddo
	endif

	if ( n_nrth .ne. MPI_PROC_NULL ) then
	   do i = -1, nni+2
	      phi(i,nnj+1,k) = rNmsg(i,1)
	      phi(i,nnj+2,k) = rNmsg(i,2)
	   enddo
	endif

	enddo

C......	k-direction

	len = ( nni + 4 ) * 2

	do j = -1, nnj+2

	nreq = 0

	if ( n_back .ne. MPI_PROC_NULL ) then
	   nreq = nreq + 1
	   call MPI_IRECV( rBmsg(-1,1), len, 
     <	                   MPI_DOUBLE_PRECISION,  
     <	                   n_back, 0, comm3d, req(nreq), ierr )
	   do i = -1, nni+2
	      sBmsg(i,1) = phi(i,j,1)
	      sBmsg(i,2) = phi(i,j,2)
	   enddo
	   nreq = nreq + 1
	   call MPI_ISEND( sBmsg(-1,1), len, 
     <	                   MPI_DOUBLE_PRECISION,  
     <	                   n_back, 1, comm3d, req(nreq), ierr )
	endif
	   
	if ( n_frnt .ne. MPI_PROC_NULL ) then
	   nreq = nreq + 1
	   call MPI_IRECV( rFmsg(-1,1), len, 
     <	                   MPI_DOUBLE_PRECISION,  
     <	                   n_frnt, 1, comm3d, req(nreq), ierr )
	   do i = -1, nni+2
	      sFmsg(i,1) = phi(i,j,nnk-1)
	      sFmsg(i,2) = phi(i,j,nnk  )
	   enddo
	   nreq = nreq + 1
	   call MPI_ISEND( sFmsg(-1,1), len, 
     <	                   MPI_DOUBLE_PRECISION,  
     <	                   n_frnt, 0, comm3d, req(nreq), ierr )
	endif

	call MPI_WAITALL( nreq, req, status, ierr )

	if ( n_back .ne. MPI_PROC_NULL ) then
	   do i = -1, nni+2
	      phi(i,j,-1) = rBmsg(i,1)
	      phi(i,j, 0) = rBmsg(i,2)
	   enddo
	endif

	if ( n_frnt .ne. MPI_PROC_NULL ) then
	   do i = -1, nni+2
	      phi(i,j,nnk+1) = rFmsg(i,1)
	      phi(i,j,nnk+2) = rFmsg(i,2)
	   enddo
	endif

	enddo

	return
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine phi2_exchange

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "ns.inc"
	include "metric.inc"

	double precision, dimension(-1:nnj+2,1:2) :: 
     <		rWmsg, sWmsg, rEmsg, sEmsg
	double precision, dimension(-1:nni+2,1:2) :: 
     <		rSmsg, sSmsg, rNmsg, sNmsg
	double precision, dimension(-1:nni+2,1:2) :: 
     <		rBmsg, sBmsg, rFmsg, sFmsg

	integer status(MPI_STATUS_SIZE,4), req(4), nreq, len
	
	integer i, j, k

C......	i-direction

	len = ( nnj + 4 ) * 2

	do k = -1, nnk+2

	nreq = 0

	if ( n_west .ne. MPI_PROC_NULL ) then
	   nreq = nreq + 1
	   call MPI_IRECV( rWmsg(-1,1), len, 
     <	                   MPI_DOUBLE_PRECISION,  
     <	                   n_west, 0, comm3d, req(nreq), ierr )
	   do j = -1, nnj+2
	      sWmsg(j,1) = phi2(1,j,k)
	      sWmsg(j,2) = phi2(2,j,k)
	   enddo
	   nreq = nreq + 1
	   call MPI_ISEND( sWmsg(-1,1), len, 
     <	                   MPI_DOUBLE_PRECISION,  
     <	                   n_west, 1, comm3d, req(nreq), ierr )
	endif
	   
	if ( n_east .ne. MPI_PROC_NULL ) then
	   nreq = nreq + 1
	   call MPI_IRECV( rEmsg(-1,1), len, 
     <	                   MPI_DOUBLE_PRECISION,  
     <	                   n_east, 1, comm3d, req(nreq), ierr )
	   do j = -1, nnj+2
	      sEmsg(j,1) = phi2(nni-1,j,k)
	      sEmsg(j,2) = phi2(nni,  j,k)
	   enddo
	   nreq = nreq + 1
	   call MPI_ISEND( sEmsg(-1,1), len, 
     <	                   MPI_DOUBLE_PRECISION,  
     <	                   n_east, 0, comm3d, req(nreq), ierr )
	endif

	call MPI_WAITALL( nreq, req, status, ierr )

	if ( n_west .ne. MPI_PROC_NULL ) then
	   do j = -1, nnj+2
	      phi2(-1,j,k) = rWmsg(j,1)
	      phi2( 0,j,k) = rWmsg(j,2)
	   enddo
	endif

	if ( n_east .ne. MPI_PROC_NULL ) then
	   do j = -1, nnj+2
	      phi2(nni+1,j,k) = rEmsg(j,1)
	      phi2(nni+2,j,k) = rEmsg(j,2)
	   enddo
	endif

	enddo

C......	j-direction

	len = ( nni + 4 ) * 2

	do k = -1, nnk+2

	nreq = 0

	if ( n_suth .ne. MPI_PROC_NULL ) then
	   nreq = nreq + 1
	   call MPI_IRECV( rSmsg(-1,1), len, 
     <	                   MPI_DOUBLE_PRECISION,  
     <	                   n_suth, 0, comm3d, req(nreq), ierr )
	   do i = -1, nni+2
	      sSmsg(i,1) = phi2(i,1,k)
	      sSmsg(i,2) = phi2(i,2,k)
	   enddo
	   nreq = nreq + 1
	   call MPI_ISEND( sSmsg(-1,1), len, 
     <	                   MPI_DOUBLE_PRECISION,  
     <	                   n_suth, 1, comm3d, req(nreq), ierr )
	endif
	   
	if ( n_nrth .ne. MPI_PROC_NULL ) then
	   nreq = nreq + 1
	   call MPI_IRECV( rNmsg(-1,1), len, 
     <	                   MPI_DOUBLE_PRECISION,  
     <	                   n_nrth, 1, comm3d, req(nreq), ierr )
	   do i = -1, nni+2
	      sNmsg(i,1) = phi2(i,nnj-1,k)
	      sNmsg(i,2) = phi2(i,nnj,  k)
	   enddo
	   nreq = nreq + 1
	   call MPI_ISEND( sNmsg(-1,1), len, 
     <	                   MPI_DOUBLE_PRECISION,  
     <	                   n_nrth, 0, comm3d, req(nreq), ierr )
	endif

	call MPI_WAITALL( nreq, req, status, ierr )

	if ( n_suth .ne. MPI_PROC_NULL ) then
	   do i = -1, nni+2
	      phi2(i,-1,k) = rSmsg(i,1)
	      phi2(i, 0,k) = rSmsg(i,2)
	   enddo
	endif

	if ( n_nrth .ne. MPI_PROC_NULL ) then
	   do i = -1, nni+2
	      phi2(i,nnj+1,k) = rNmsg(i,1)
	      phi2(i,nnj+2,k) = rNmsg(i,2)
	   enddo
	endif

	enddo

C......	k-direction

	len = ( nni + 4 ) * 2

	do j = -1, nnj+2

	nreq = 0

	if ( n_back .ne. MPI_PROC_NULL ) then
	   nreq = nreq + 1
	   call MPI_IRECV( rBmsg(-1,1), len, 
     <	                   MPI_DOUBLE_PRECISION,  
     <	                   n_back, 0, comm3d, req(nreq), ierr )
	   do i = -1, nni+2
	      sBmsg(i,1) = phi2(i,j,1)
	      sBmsg(i,2) = phi2(i,j,2)
	   enddo
	   nreq = nreq + 1
	   call MPI_ISEND( sBmsg(-1,1), len, 
     <	                   MPI_DOUBLE_PRECISION,  
     <	                   n_back, 1, comm3d, req(nreq), ierr )
	endif
	   
	if ( n_frnt .ne. MPI_PROC_NULL ) then
	   nreq = nreq + 1
	   call MPI_IRECV( rFmsg(-1,1), len, 
     <	                   MPI_DOUBLE_PRECISION,  
     <	                   n_frnt, 1, comm3d, req(nreq), ierr )
	   do i = -1, nni+2
	      sFmsg(i,1) = phi2(i,j,nnk-1)
	      sFmsg(i,2) = phi2(i,j,nnk  )
	   enddo
	   nreq = nreq + 1
	   call MPI_ISEND( sFmsg(-1,1), len, 
     <	                   MPI_DOUBLE_PRECISION,  
     <	                   n_frnt, 0, comm3d, req(nreq), ierr )
	endif

	call MPI_WAITALL( nreq, req, status, ierr )

	if ( n_back .ne. MPI_PROC_NULL ) then
	   do i = -1, nni+2
	      phi2(i,j,-1) = rBmsg(i,1)
	      phi2(i,j, 0) = rBmsg(i,2)
	   enddo
	endif

	if ( n_frnt .ne. MPI_PROC_NULL ) then
	   do i = -1, nni+2
	      phi2(i,j,nnk+1) = rFmsg(i,1)
	      phi2(i,j,nnk+2) = rFmsg(i,2)
	   enddo
	endif

	enddo

	return
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine phi3_exchange

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "ns.inc"
	include "metric.inc"

	double precision, dimension(-1:nnj+2,1:2) :: 
     <		rWmsg, sWmsg, rEmsg, sEmsg
	double precision, dimension(-1:nni+2,1:2) :: 
     <		rSmsg, sSmsg, rNmsg, sNmsg
	double precision, dimension(-1:nni+2,1:2) :: 
     <		rBmsg, sBmsg, rFmsg, sFmsg

	integer status(MPI_STATUS_SIZE,4), req(4), nreq, len
	
	integer i, j, k

C......	i-direction

	len = ( nnj + 4 ) * 2

	do k = -1, nnk+2

	nreq = 0

	if ( n_west .ne. MPI_PROC_NULL ) then
	   nreq = nreq + 1
	   call MPI_IRECV( rWmsg(-1,1), len, 
     <	                   MPI_DOUBLE_PRECISION,  
     <	                   n_west, 0, comm3d, req(nreq), ierr )
	   do j = -1, nnj+2
	      sWmsg(j,1) = phi3(1,j,k)
	      sWmsg(j,2) = phi3(2,j,k)
	   enddo
	   nreq = nreq + 1
	   call MPI_ISEND( sWmsg(-1,1), len, 
     <	                   MPI_DOUBLE_PRECISION,  
     <	                   n_west, 1, comm3d, req(nreq), ierr )
	endif
	   
	if ( n_east .ne. MPI_PROC_NULL ) then
	   nreq = nreq + 1
	   call MPI_IRECV( rEmsg(-1,1), len, 
     <	                   MPI_DOUBLE_PRECISION,  
     <	                   n_east, 1, comm3d, req(nreq), ierr )
	   do j = -1, nnj+2
	      sEmsg(j,1) = phi3(nni-1,j,k)
	      sEmsg(j,2) = phi3(nni,  j,k)
	   enddo
	   nreq = nreq + 1
	   call MPI_ISEND( sEmsg(-1,1), len, 
     <	                   MPI_DOUBLE_PRECISION,  
     <	                   n_east, 0, comm3d, req(nreq), ierr )
	endif

	call MPI_WAITALL( nreq, req, status, ierr )

	if ( n_west .ne. MPI_PROC_NULL ) then
	   do j = -1, nnj+2
	      phi3(-1,j,k) = rWmsg(j,1)
	      phi3( 0,j,k) = rWmsg(j,2)
	   enddo
	endif

	if ( n_east .ne. MPI_PROC_NULL ) then
	   do j = -1, nnj+2
	      phi3(nni+1,j,k) = rEmsg(j,1)
	      phi3(nni+2,j,k) = rEmsg(j,2)
	   enddo
	endif

	enddo

C......	j-direction

	len = ( nni + 4 ) * 2

	do k = -1, nnk+2

	nreq = 0

	if ( n_suth .ne. MPI_PROC_NULL ) then
	   nreq = nreq + 1
	   call MPI_IRECV( rSmsg(-1,1), len, 
     <	                   MPI_DOUBLE_PRECISION,  
     <	                   n_suth, 0, comm3d, req(nreq), ierr )
	   do i = -1, nni+2
	      sSmsg(i,1) = phi3(i,1,k)
	      sSmsg(i,2) = phi3(i,2,k)
	   enddo
	   nreq = nreq + 1
	   call MPI_ISEND( sSmsg(-1,1), len, 
     <	                   MPI_DOUBLE_PRECISION,  
     <	                   n_suth, 1, comm3d, req(nreq), ierr )
	endif
	   
	if ( n_nrth .ne. MPI_PROC_NULL ) then
	   nreq = nreq + 1
	   call MPI_IRECV( rNmsg(-1,1), len, 
     <	                   MPI_DOUBLE_PRECISION,  
     <	                   n_nrth, 1, comm3d, req(nreq), ierr )
	   do i = -1, nni+2
	      sNmsg(i,1) = phi3(i,nnj-1,k)
	      sNmsg(i,2) = phi3(i,nnj,  k)
	   enddo
	   nreq = nreq + 1
	   call MPI_ISEND( sNmsg(-1,1), len, 
     <	                   MPI_DOUBLE_PRECISION,  
     <	                   n_nrth, 0, comm3d, req(nreq), ierr )
	endif

	call MPI_WAITALL( nreq, req, status, ierr )

	if ( n_suth .ne. MPI_PROC_NULL ) then
	   do i = -1, nni+2
	      phi3(i,-1,k) = rSmsg(i,1)
	      phi3(i, 0,k) = rSmsg(i,2)
	   enddo
	endif

	if ( n_nrth .ne. MPI_PROC_NULL ) then
	   do i = -1, nni+2
	      phi3(i,nnj+1,k) = rNmsg(i,1)
	      phi3(i,nnj+2,k) = rNmsg(i,2)
	   enddo
	endif

	enddo

C......	k-direction

	len = ( nni + 4 ) * 2

	do j = -1, nnj+2

	nreq = 0

	if ( n_back .ne. MPI_PROC_NULL ) then
	   nreq = nreq + 1
	   call MPI_IRECV( rBmsg(-1,1), len, 
     <	                   MPI_DOUBLE_PRECISION,  
     <	                   n_back, 0, comm3d, req(nreq), ierr )
	   do i = -1, nni+2
	      sBmsg(i,1) = phi3(i,j,1)
	      sBmsg(i,2) = phi3(i,j,2)
	   enddo
	   nreq = nreq + 1
	   call MPI_ISEND( sBmsg(-1,1), len, 
     <	                   MPI_DOUBLE_PRECISION,  
     <	                   n_back, 1, comm3d, req(nreq), ierr )
	endif
	   
	if ( n_frnt .ne. MPI_PROC_NULL ) then
	   nreq = nreq + 1
	   call MPI_IRECV( rFmsg(-1,1), len, 
     <	                   MPI_DOUBLE_PRECISION,  
     <	                   n_frnt, 1, comm3d, req(nreq), ierr )
	   do i = -1, nni+2
	      sFmsg(i,1) = phi3(i,j,nnk-1)
	      sFmsg(i,2) = phi3(i,j,nnk  )
	   enddo
	   nreq = nreq + 1
	   call MPI_ISEND( sFmsg(-1,1), len, 
     <	                   MPI_DOUBLE_PRECISION,  
     <	                   n_frnt, 0, comm3d, req(nreq), ierr )
	endif

	call MPI_WAITALL( nreq, req, status, ierr )

	if ( n_back .ne. MPI_PROC_NULL ) then
	   do i = -1, nni+2
	      phi3(i,j,-1) = rBmsg(i,1)
	      phi3(i,j, 0) = rBmsg(i,2)
	   enddo
	endif

	if ( n_frnt .ne. MPI_PROC_NULL ) then
	   do i = -1, nni+2
	      phi3(i,j,nnk+1) = rFmsg(i,1)
	      phi3(i,j,nnk+2) = rFmsg(i,2)
	   enddo
	endif

	enddo

	return
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine phi4_exchange

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "ns.inc"
	include "metric.inc"

	double precision, dimension(-1:nnj+2,1:2) :: 
     <		rWmsg, sWmsg, rEmsg, sEmsg
	double precision, dimension(-1:nni+2,1:2) :: 
     <		rSmsg, sSmsg, rNmsg, sNmsg
	double precision, dimension(-1:nni+2,1:2) :: 
     <		rBmsg, sBmsg, rFmsg, sFmsg

	integer status(MPI_STATUS_SIZE,4), req(4), nreq, len
	
	integer i, j, k

C......	i-direction

	len = ( nnj + 4 ) * 2

	do k = -1, nnk+2

	nreq = 0

	if ( n_west .ne. MPI_PROC_NULL ) then
	   nreq = nreq + 1
	   call MPI_IRECV( rWmsg(-1,1), len, 
     <	                   MPI_DOUBLE_PRECISION,  
     <	                   n_west, 0, comm3d, req(nreq), ierr )
	   do j = -1, nnj+2
	      sWmsg(j,1) = phi4(1,j,k)
	      sWmsg(j,2) = phi4(2,j,k)
	   enddo
	   nreq = nreq + 1
	   call MPI_ISEND( sWmsg(-1,1), len, 
     <	                   MPI_DOUBLE_PRECISION,  
     <	                   n_west, 1, comm3d, req(nreq), ierr )
	endif
	   
	if ( n_east .ne. MPI_PROC_NULL ) then
	   nreq = nreq + 1
	   call MPI_IRECV( rEmsg(-1,1), len, 
     <	                   MPI_DOUBLE_PRECISION,  
     <	                   n_east, 1, comm3d, req(nreq), ierr )
	   do j = -1, nnj+2
	      sEmsg(j,1) = phi4(nni-1,j,k)
	      sEmsg(j,2) = phi4(nni,  j,k)
	   enddo
	   nreq = nreq + 1
	   call MPI_ISEND( sEmsg(-1,1), len, 
     <	                   MPI_DOUBLE_PRECISION,  
     <	                   n_east, 0, comm3d, req(nreq), ierr )
	endif

	call MPI_WAITALL( nreq, req, status, ierr )

	if ( n_west .ne. MPI_PROC_NULL ) then
	   do j = -1, nnj+2
	      phi4(-1,j,k) = rWmsg(j,1)
	      phi4( 0,j,k) = rWmsg(j,2)
	   enddo
	endif

	if ( n_east .ne. MPI_PROC_NULL ) then
	   do j = -1, nnj+2
	      phi4(nni+1,j,k) = rEmsg(j,1)
	      phi4(nni+2,j,k) = rEmsg(j,2)
	   enddo
	endif

	enddo

C......	j-direction

	len = ( nni + 4 ) * 2

	do k = -1, nnk+2

	nreq = 0

	if ( n_suth .ne. MPI_PROC_NULL ) then
	   nreq = nreq + 1
	   call MPI_IRECV( rSmsg(-1,1), len, 
     <	                   MPI_DOUBLE_PRECISION,  
     <	                   n_suth, 0, comm3d, req(nreq), ierr )
	   do i = -1, nni+2
	      sSmsg(i,1) = phi4(i,1,k)
	      sSmsg(i,2) = phi4(i,2,k)
	   enddo
	   nreq = nreq + 1
	   call MPI_ISEND( sSmsg(-1,1), len, 
     <	                   MPI_DOUBLE_PRECISION,  
     <	                   n_suth, 1, comm3d, req(nreq), ierr )
	endif
	   
	if ( n_nrth .ne. MPI_PROC_NULL ) then
	   nreq = nreq + 1
	   call MPI_IRECV( rNmsg(-1,1), len, 
     <	                   MPI_DOUBLE_PRECISION,  
     <	                   n_nrth, 1, comm3d, req(nreq), ierr )
	   do i = -1, nni+2
	      sNmsg(i,1) = phi4(i,nnj-1,k)
	      sNmsg(i,2) = phi4(i,nnj,  k)
	   enddo
	   nreq = nreq + 1
	   call MPI_ISEND( sNmsg(-1,1), len, 
     <	                   MPI_DOUBLE_PRECISION,  
     <	                   n_nrth, 0, comm3d, req(nreq), ierr )
	endif

	call MPI_WAITALL( nreq, req, status, ierr )

	if ( n_suth .ne. MPI_PROC_NULL ) then
	   do i = -1, nni+2
	      phi4(i,-1,k) = rSmsg(i,1)
	      phi4(i, 0,k) = rSmsg(i,2)
	   enddo
	endif

	if ( n_nrth .ne. MPI_PROC_NULL ) then
	   do i = -1, nni+2
	      phi4(i,nnj+1,k) = rNmsg(i,1)
	      phi4(i,nnj+2,k) = rNmsg(i,2)
	   enddo
	endif

	enddo

C......	k-direction

	len = ( nni + 4 ) * 2

	do j = -1, nnj+2

	nreq = 0

	if ( n_back .ne. MPI_PROC_NULL ) then
	   nreq = nreq + 1
	   call MPI_IRECV( rBmsg(-1,1), len, 
     <	                   MPI_DOUBLE_PRECISION,  
     <	                   n_back, 0, comm3d, req(nreq), ierr )
	   do i = -1, nni+2
	      sBmsg(i,1) = phi4(i,j,1)
	      sBmsg(i,2) = phi4(i,j,2)
	   enddo
	   nreq = nreq + 1
	   call MPI_ISEND( sBmsg(-1,1), len, 
     <	                   MPI_DOUBLE_PRECISION,  
     <	                   n_back, 1, comm3d, req(nreq), ierr )
	endif
	   
	if ( n_frnt .ne. MPI_PROC_NULL ) then
	   nreq = nreq + 1
	   call MPI_IRECV( rFmsg(-1,1), len, 
     <	                   MPI_DOUBLE_PRECISION,  
     <	                   n_frnt, 1, comm3d, req(nreq), ierr )
	   do i = -1, nni+2
	      sFmsg(i,1) = phi4(i,j,nnk-1)
	      sFmsg(i,2) = phi4(i,j,nnk  )
	   enddo
	   nreq = nreq + 1
	   call MPI_ISEND( sFmsg(-1,1), len, 
     <	                   MPI_DOUBLE_PRECISION,  
     <	                   n_frnt, 0, comm3d, req(nreq), ierr )
	endif

	call MPI_WAITALL( nreq, req, status, ierr )

	if ( n_back .ne. MPI_PROC_NULL ) then
	   do i = -1, nni+2
	      phi4(i,j,-1) = rBmsg(i,1)
	      phi4(i,j, 0) = rBmsg(i,2)
	   enddo
	endif

	if ( n_frnt .ne. MPI_PROC_NULL ) then
	   do i = -1, nni+2
	      phi4(i,j,nnk+1) = rFmsg(i,1)
	      phi4(i,j,nnk+2) = rFmsg(i,2)
	   enddo
	endif

	enddo

	return
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine p_exchange( ii, jj, kk, f )

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"

	integer ii, jj, kk
	double precision, dimension(0:ii+1,0:jj+1,0:kk+1) :: f

	double precision, dimension(0:jj+1,0:kk+1) ::
     <		sEmsg, rEmsg, sWmsg, rWmsg

	double precision, dimension(0:ii+1,0:kk+1) ::
     <		sSmsg, rSmsg, sNmsg, rNmsg

	double precision, dimension(0:ii+1,0:jj+1) ::
     <		sBmsg, rBmsg, sFmsg, rFmsg

	integer status(MPI_STATUS_SIZE,4), req(4), nreq, len

	integer i, j, k

	nreq = 0

	len = ( jj + 2 ) * ( kk + 2 )

	if ( n_west .ne. MPI_PROC_NULL ) then
	   nreq = nreq + 1
	   call MPI_IRECV( rWmsg(0,0), len, MPI_DOUBLE_PRECISION,  
     <	                   n_west, 0, comm3d, req(nreq), ierr )
	   do k = 0, kk+1
	   do j = 0, jj+1
	      sWmsg(j,k) = f(1,j,k)
	   enddo
	   enddo
	   nreq = nreq + 1
	   call MPI_ISEND( sWmsg(0,0), len, MPI_DOUBLE_PRECISION,  
     <	                   n_west, 1, comm3d, req(nreq), ierr )
	endif

	if ( n_east .ne. MPI_PROC_NULL ) then
	   nreq = nreq + 1
	   call MPI_IRECV( rEmsg(0,0), len, MPI_DOUBLE_PRECISION,  
     <	                   n_east, 1, comm3d, req(nreq), ierr )
	   do k = 0, kk+1
	   do j = 0, jj+1
	      sEmsg(j,k) = f(ii,j,k)
	   enddo
	   enddo
	   nreq = nreq + 1
	   call MPI_ISEND( sEmsg(0,0), len, MPI_DOUBLE_PRECISION,  
     <	                   n_east, 0, comm3d, req(nreq), ierr )
	endif

	call MPI_WAITALL( nreq, req, status, ierr )

	if ( n_west .ne. MPI_PROC_NULL ) then
	   do k = 0, kk+1
	   do j = 0, jj+1
	      f(0,j,k) = rWmsg(j,k)
	   enddo
	   enddo
	endif

	if ( n_east .ne. MPI_PROC_NULL ) then
	   do k = 0, kk+1
	   do j = 0, jj+1
	      f(ii+1,j,k) = rEmsg(j,k)
	   enddo
	   enddo
	endif

	nreq = 0

	len = ( ii + 2 ) * ( kk + 2 )

	if ( n_suth .ne. MPI_PROC_NULL ) then
	   nreq = nreq + 1
	   call MPI_IRECV( rSmsg(0,0), len, MPI_DOUBLE_PRECISION,  
     <	                   n_suth, 0, comm3d, req(nreq), ierr )
	   do k = 0, kk+1
	   do i = 0, ii+1
	      sSmsg(i,k) = f(i,1,k)
	   enddo
	   enddo
	   nreq = nreq + 1
	   call MPI_ISEND( sSmsg(0,0), len, MPI_DOUBLE_PRECISION,  
     <	                   n_suth, 1, comm3d, req(nreq), ierr )
	endif

	if ( n_nrth .ne. MPI_PROC_NULL ) then
	   nreq = nreq + 1
	   call MPI_IRECV( rNmsg(0,0), len, MPI_DOUBLE_PRECISION,  
     <	                   n_nrth, 1, comm3d, req(nreq), ierr )
	   do k = 0, kk+1
	   do i = 0, ii+1
	      sNmsg(i,k) = f(i,jj,k)
	   enddo
	   enddo
	   nreq = nreq + 1
	   call MPI_ISEND( sNmsg(0,0), len, MPI_DOUBLE_PRECISION,  
     <	                   n_nrth, 0, comm3d, req(nreq), ierr )
	endif

	call MPI_WAITALL( nreq, req, status, ierr )

	if ( n_suth .ne. MPI_PROC_NULL ) then
	   do k = 0, kk+1
	   do i = 0, ii+1
	      f(i,0,k) = rSmsg(i,k)
	   enddo
	   enddo
	endif

	if ( n_nrth .ne. MPI_PROC_NULL ) then
	   do k = 0, kk+1
	   do i = 0, ii+1
	      f(i,jj+1,k) = rNmsg(i,k)
	   enddo
	   enddo
	endif

	nreq = 0

	len = ( ii + 2 ) * ( jj + 2 )

	if ( n_back .ne. MPI_PROC_NULL ) then
	   nreq = nreq + 1
	   call MPI_IRECV( rBmsg(0,0), len, MPI_DOUBLE_PRECISION,  
     <	                   n_back, 0, comm3d, req(nreq), ierr )
	   do j = 0, jj+1
	   do i = 0, ii+1
	      sBmsg(i,j) = f(i,j,1)
	   enddo
	   enddo
	   nreq = nreq + 1
	   call MPI_ISEND( sBmsg(0,0), len, MPI_DOUBLE_PRECISION,  
     <	                   n_back, 1, comm3d, req(nreq), ierr )
	endif

	if ( n_frnt .ne. MPI_PROC_NULL ) then
	   nreq = nreq + 1
	   call MPI_IRECV( rFmsg(0,0), len, MPI_DOUBLE_PRECISION,  
     <	                   n_frnt, 1, comm3d, req(nreq), ierr )
	   do j = 0, jj+1
	   do i = 0, ii+1
	      sFmsg(i,j) = f(i,j,kk)
	   enddo
	   enddo
	   nreq = nreq + 1
	   call MPI_ISEND( sFmsg(0,0), len, MPI_DOUBLE_PRECISION,  
     <	                   n_frnt, 0, comm3d, req(nreq), ierr )
	endif

	call MPI_WAITALL( nreq, req, status, ierr )

	if ( n_back .ne. MPI_PROC_NULL ) then
	   do j = 0, jj+1
	   do i = 0, ii+1
	      f(i,j,0) = rBmsg(i,j)
	   enddo
	   enddo
	endif

	if ( n_frnt .ne. MPI_PROC_NULL ) then
	   do j = 0, jj+1
	   do i = 0, ii+1
	      f(i,j,kk+1) = rFmsg(i,j)
	   enddo
	   enddo
	endif

	return
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine pi_exchange( ii, jj, kk, f, js, je, ji, ks, ke, ki,
     <		is, ir, n_send, n_recv )

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"

	integer ii, jj, kk
	double precision, dimension(0:ii+1,0:jj+1,0:kk+1) :: f

	integer js, je, ji, ks, ke, ki
	double precision, dimension((je-js)/ji+1,(ke-ks)/ki+1) ::
     <		smsg, rmsg

	integer is, ir, n_send, n_recv

	integer status(MPI_STATUS_SIZE,2), req(2), nreq, len

	integer j, k, jd, kd

	nreq = 0

	len = ((je-js)/ji+1) * ((ke-ks)/ki+1)

	if ( n_recv .ne. MPI_PROC_NULL ) then
	   nreq = nreq + 1
	   call MPI_IRECV( rmsg(1,1), len, MPI_DOUBLE_PRECISION,  
     <	                   n_recv, 0, comm3d, req(nreq), ierr )
     	endif
     
	if ( n_send .ne. MPI_PROC_NULL ) then
	   kd = 0
	   do k = ks, ke, ki
	      kd = kd + 1
	      jd = 0
	      do j = js, je, ji
	         jd = jd + 1
	         smsg(jd,kd) = f(is,j,k)
	      enddo
	   enddo
	   nreq = nreq + 1
	   call MPI_ISEND( smsg(1,1), len, MPI_DOUBLE_PRECISION,  
     <	                   n_send, 0, comm3d, req(nreq), ierr )
	endif

	call MPI_WAITALL( nreq, req, status, ierr )

	if ( n_recv .ne. MPI_PROC_NULL ) then
	   kd = 0
	   do k = ks, ke, ki
	      kd = kd + 1
	      jd = 0
	      do j = js, je, ji
	         jd = jd + 1
	         f(ir,j,k) = rmsg(jd,kd)
	      enddo
	   enddo
	endif

	return
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine pj_exchange( ii, jj, kk, f, is, ie, ic, ks, ke, ki,
     <		js, jr, n_send, n_recv )

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"

	integer ii, jj, kk
	double precision, dimension(0:ii+1,0:jj+1,0:kk+1) :: f

	integer is, ie, ic, ks, ke, ki
	double precision, dimension((ie-is)/ic+1,(ke-ks)/ki+1) ::
     <		smsg, rmsg

	integer js, jr, n_send, n_recv

	integer status(MPI_STATUS_SIZE,2), req(2), nreq, len

	integer i, k, id, kd

	nreq = 0

	len = ((ie-is)/ic+1) * ((ke-ks)/ki+1)

	if ( n_recv .ne. MPI_PROC_NULL ) then
	   nreq = nreq + 1
	   call MPI_IRECV( rmsg(1,1), len, MPI_DOUBLE_PRECISION,  
     <	                   n_recv, 0, comm3d, req(nreq), ierr )
	endif
	
	if ( n_send .ne. MPI_PROC_NULL ) then
	   kd = 0
	   do k = ks, ke, ki
	      kd = kd + 1
	      id = 0
	      do i = is, ie, ic
	         id = id + 1
	         smsg(id,kd) = f(i,js,k)
	      enddo
	   enddo
	   nreq = nreq + 1
	   call MPI_ISEND( smsg(1,1), len, MPI_DOUBLE_PRECISION,  
     <	                   n_send, 0, comm3d, req(nreq), ierr )
	endif

	call MPI_WAITALL( nreq, req, status, ierr )

	if ( n_recv .ne. MPI_PROC_NULL ) then
	   kd = 0
	   do k = ks, ke, ki
	      kd = kd + 1
	      id = 0
	      do i = is, ie, ic
	         id = id + 1
	         f(i,jr,k) = rmsg(id,kd)
	      enddo
	   enddo
	endif

	return
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine pk_exchange( ii, jj, kk, f, is, ie, ic, js, je, ji,
     <		ks, kr, n_send, n_recv )

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"

	integer ii, jj, kk
	double precision, dimension(0:ii+1,0:jj+1,0:kk+1) :: f

	integer is, ie, ic, js, je, ji
	double precision, dimension((ie-is)/ic+1,(je-js)/ji+1) ::
     <		smsg, rmsg

	integer ks, kr, n_send, n_recv

	integer status(MPI_STATUS_SIZE,2), req(2), nreq, len

	integer i, j, id, jd

	nreq = 0

	len = ((ie-is)/ic+1) * ((je-js)/ji+1)

	if ( n_recv .ne. MPI_PROC_NULL ) then
	   nreq = nreq + 1
	   call MPI_IRECV( rmsg(1,1), len, MPI_DOUBLE_PRECISION,  
     <	                   n_recv, 0, comm3d, req(nreq), ierr )
     	endif
     	
	if ( n_send .ne. MPI_PROC_NULL ) then
	   jd = 0
	   do j = js, je, ji
	      jd = jd + 1
	      id = 0
	      do i = is, ie, ic
	         id = id + 1
	         smsg(id,jd) = f(i,j,ks)
	      enddo
	   enddo
	   nreq = nreq + 1
	   call MPI_ISEND( smsg(1,1), len, MPI_DOUBLE_PRECISION,  
     <	                   n_send, 0, comm3d, req(nreq), ierr )
	endif

	call MPI_WAITALL( nreq, req, status, ierr )

	if ( n_recv .ne. MPI_PROC_NULL ) then
	   jd = 0
	   do j = js, je, ji
	      jd = jd + 1
	      id = 0
	      do i = is, ie, ic
	         id = id + 1
	         f(i,j,kr) = rmsg(id,jd)
	      enddo
	   enddo
	endif

	return
	end
