cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine trid( a, b, c, f, n, L, m, nbrlower, nbrupper )

	implicit none
	include "mpif.h"
	include "mpi.inc" 

	integer n, L, m
	double precision a(n,  0:m+1), b(n,  0:m+1), c(n,  0:m+1)
	double precision f(n,L,0:m+1)
	integer nbrlower, nbrupper

	double precision rmsg(n,L+1), smsg(n,L+1)
	integer status(MPI_STATUS_SIZE), req
	integer i, K, j
	
	if ( nbrlower .ne. MPI_PROC_NULL ) then
	   call MPI_IRECV( rmsg(1,1), n*(L+1), MPI_DOUBLE_PRECISION,  
     <	                   nbrlower, 0, comm3d, req, ierr )
	   call MPI_WAIT( req, status, ierr )
	   do K = 1, L
	   do i = 1, n
	      f(i,K,0) = rmsg(i,K)
	   enddo
	   enddo
	   do i = 1, n
	      c(i,0) = rmsg(i,L+1)
	   enddo
	else
	   do K = 1, L
	   do i = 1, n
	      f(i,K,0) = f(i,K,0) / b(i,0)
	   enddo
           enddo
	   do i = 1, n
	      c(i,0) = c(i,0) / b(i,0)
           enddo
	   
	endif	

        do j = 1, m
	do i = 1, n
	   b(i,j) = 1.D0 / ( b(i,j) - a(i,j) * c(i,j-1) )
	   c(i,j) = c(i,j) * b(i,j)
	enddo
	enddo

        do j = 1, m
	do K = 1, L
	do i = 1, n
	   f(i,K,j) = ( f(i,K,j) - a(i,j) * f(i,K,j-1) ) * b(i,j)
	enddo
	enddo
	enddo

	if ( nbrupper .ne. MPI_PROC_NULL ) then
	   do K = 1, L
	   do i = 1, n
	      smsg(i,K) = f(i,K,m)
	   enddo
	   enddo
	   do i = 1, n
	      smsg(i,L+1) = c(i,m)
	   enddo
	   call MPI_ISEND( smsg(1,1), n*(L+1), MPI_DOUBLE_PRECISION, 
     <	                   nbrupper, 0, comm3d, req, ierr )
	   call MPI_WAIT( req, status, ierr )
	else
	   do K = 1, L
	   do i = 1, n
	      f(i,K,m+1) = ( f(i,K,m+1) - a(i,m+1) * f(i,K,m) )
     <	                 / ( b(i,  m+1) - a(i,m+1) * c(i,  m) )
	   enddo
	   enddo
	endif

	if ( nbrupper .ne. MPI_PROC_NULL ) then
	   call MPI_IRECV( rmsg(1,1), n*L, MPI_DOUBLE_PRECISION,  
     <	                   nbrupper, 1, comm3d, req, ierr )
	   call MPI_WAIT( req, status, ierr )
	   do K = 1, L
	   do i = 1, n
	      f(i,K,m+1) = rmsg(i,K)
	   enddo
	   enddo
	endif	

        do j = m, 1, -1
	do K = 1, L
	do i = 1, n
	   
	   f(i,K,j) = f(i,K,j) - c(i,j) * f(i,K,j+1)
	enddo
	enddo
	enddo

	if ( nbrlower .ne. MPI_PROC_NULL ) then
	   do K = 1, L
	   do i = 1, n
	      smsg(i,K) = f(i,K,1)
	   enddo
	   enddo
	   call MPI_ISEND( smsg(1,1), n*L, MPI_DOUBLE_PRECISION, 
     <	                   nbrlower, 1, comm3d, req, ierr )
	   call MPI_WAIT( req, status, ierr )
	else
	   do K = 1, L
	   do i = 1, n
	      f(i,K,0) = f(i,K,0) - c(i,0) * f(i,K,1)
	      
	   enddo
	   enddo
	endif

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine trip( a, b, c, f, n, L, m, nbrlower, nbrupper )

	implicit none
	include "mpif.h"
	include "mpi.inc" 

	integer n, L, m
	double precision a(n,    0:m+1), b(n,    0:m+1), c(n,    0:m+1)
	double precision f(n,L,  0:m+1)
	double precision g(n,L+1,0:m+1)
	integer nbrlower, nbrupper, nbrl, nbru

	double precision rmsg(n,L+1), smsg(n,L+1)
	integer status(MPI_STATUS_SIZE), req
	integer i, K, j

	do j = 0, m+1
	do K = 1, L
	do i = 1, n
	   g(i,K,j) = f(i,K,j)
	enddo
	enddo
	enddo
	do j = 0, m+1
	do i = 1, n
	   g(i,L+1,j) = 0.D0
	enddo
	enddo
	
	nbrl = nbrlower
	if ( nbrlower .ge. myid ) then
	   nbrl = MPI_PROC_NULL
	   do K = 1, L+1
	   do i = 1, n
	      g(i,K,0) = 0.D0
	   enddo
	   enddo
	   do i = 1, n
	      a(i,0) = 0.D0
	      b(i,0) = 1.D0
	      c(i,0) = 0.D0
	      g(i,L+1,1) = a(i,1)
	      b(i,1) = b(i,1) - a(i,1)
	      a(i,1) = 0.D0
	   enddo
	endif

	nbru = nbrupper
	if ( nbrupper .le. myid ) then
	   nbru = MPI_PROC_NULL
	   do K = 1, L+1
	   do i = 1, n
	      g(i,K,m+1) = 0.D0
	   enddo
	   enddo
	   do i = 1, n
	      a(i,m+1) = 0.D0
	      b(i,m+1) = 1.D0
	      c(i,m+1) = 0.D0
	      g(i,L+1,m) = c(i,m)
	      b(i,m) = b(i,m) - c(i,m)
	      c(i,m) = 0.D0
	   enddo
	endif

	call trid( a, b, c, g, n, L+1, m, nbrl, nbru )

	if ( nbrlower .eq. myid ) then
	   do K = 1, L+1
	   do i = 1, n
	      rmsg(i,K) = g(i,K,m) + g(i,K,1)
	   enddo
	   enddo
	   do K = 1, L
	   do i = 1, n
	      smsg(i,K) = rmsg(i,K) / ( 1.D0 + rmsg(i,L+1) )
	   enddo
	   enddo
	endif
      	   
	if ( nbrlower .gt. myid ) then
	   call MPI_IRECV( rmsg(1,1), n*(L+1), MPI_DOUBLE_PRECISION,  
     <	                   nbrlower, 0, comm3d, req, ierr )
	   call MPI_WAIT( req, status, ierr )
	   do K = 1, L+1
	   do i = 1, n
	      rmsg(i,K) = rmsg(i,K) + g(i,K,1)
	   enddo
	   enddo
	   do K = 1, L
	   do i = 1, n
	      smsg(i,K) = rmsg(i,K) / ( 1.D0 + rmsg(i,L+1) )
	   enddo
	   enddo
	endif

	if ( nbrupper .lt. myid ) then
	   do K = 1, L+1
	   do i = 1, n
	      smsg(i,K) = g(i,K,m)
	   enddo
	   enddo
	   call MPI_ISEND( smsg(1,1), n*(L+1), MPI_DOUBLE_PRECISION, 
     <	                   nbrupper, 0, comm3d, req, ierr )
	   call MPI_WAIT( req, status, ierr )
	endif

	if ( nbrlower .lt. myid ) then
	   call MPI_IRECV( rmsg(1,1), n*L, MPI_DOUBLE_PRECISION,  
     <	                   nbrlower, 0, comm3d, req, ierr )
	   call MPI_WAIT( req, status, ierr )
	   do K = 1, L
	   do i = 1, n
	      smsg(i,K) = rmsg(i,K)
	   enddo
	   enddo
	endif

	if ( nbrupper .gt. myid ) then
	   call MPI_ISEND( smsg(1,1), n*L, MPI_DOUBLE_PRECISION, 
     <	                   nbrupper, 0, comm3d, req, ierr )
	   call MPI_WAIT( req, status, ierr )
	endif

	do j = 0, m+1
	do K = 1, L
	do i = 1, n
	   f(i,K,j) = g(i,K,j) - smsg(i,K) * g(i,L+1,j)
	enddo
	enddo
	enddo

	return
	end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc




c.......This function solve the system matrix as same as "trid" but instead
c       of solving three direction (L = 3) as "trid", this solves for
c       the spcified directions L1~L2. 
	subroutine trid1( a, b, c, f, n, L1, L2, m, nbrlower, nbrupper )

	implicit none
	include "mpif.h"
	include "mpi.inc" 

	integer n, L1, L2, m
	double precision a(n,  0:m+1), b(n,  0:m+1), c(n,  0:m+1)
	double precision f(n,3,0:m+1)
	integer nbrlower, nbrupper

	double precision rmsg(n,L2-L1+1+1), smsg(n,L2-L1+1+1)
	integer status(MPI_STATUS_SIZE), req
	integer i, K, j, L
	L = L2-L1+1
	
	if ( nbrlower .ne. MPI_PROC_NULL ) then
	   call MPI_IRECV( rmsg(1,1), n*(L+1), MPI_DOUBLE_PRECISION,  
     <	                   nbrlower, 0, comm3d, req, ierr )
	   call MPI_WAIT( req, status, ierr )
	   do K = L1, L2
	   do i = 1, n
	      f(i,K,0) = rmsg(i,K)
	   enddo
	   enddo
	   do i = 1, n
	      c(i,0) = rmsg(i,L+1)
	   enddo
	else
	   do K = L1, L2
	   do i = 1, n
	      f(i,K,0) = f(i,K,0) / b(i,0)
	   enddo
           enddo
	   do i = 1, n
	      c(i,0) = c(i,0) / b(i,0)
           enddo
	endif	

        do j = 1, m
	do i = 1, n
	   b(i,j) = 1.D0 / ( b(i,j) - a(i,j) * c(i,j-1) )
	   c(i,j) = c(i,j) * b(i,j)
	enddo
	enddo

        do j = 1, m
	do K = L1, L2
	do i = 1, n
	   f(i,K,j) = ( f(i,K,j) - a(i,j) * f(i,K,j-1) ) * b(i,j)
	enddo
	enddo
	enddo

	if ( nbrupper .ne. MPI_PROC_NULL ) then
	   do K = L1, L2
	   do i = 1, n
	      smsg(i,K) = f(i,K,m)
	   enddo
	   enddo
	   do i = 1, n
	      smsg(i,L+1) = c(i,m)
	   enddo
	   call MPI_ISEND( smsg(1,1), n*(L+1), MPI_DOUBLE_PRECISION, 
     <	                   nbrupper, 0, comm3d, req, ierr )
	   call MPI_WAIT( req, status, ierr )
	else
	   do K = L1, L2
	   do i = 1, n
	      f(i,K,m+1) = ( f(i,K,m+1) - a(i,m+1) * f(i,K,m) )
     <	                 / ( b(i,  m+1) - a(i,m+1) * c(i,  m) )
	   enddo
	   enddo
	endif

	if ( nbrupper .ne. MPI_PROC_NULL ) then
	   call MPI_IRECV( rmsg(1,1), n*L, MPI_DOUBLE_PRECISION,  
     <	                   nbrupper, 1, comm3d, req, ierr )
	   call MPI_WAIT( req, status, ierr )
	   do K = L1, L2
	   do i = 1, n
	      f(i,K,m+1) = rmsg(i,K)
	   enddo
	   enddo
	endif	

        do j = m, 1, -1
	do K = L1, L2
	do i = 1, n
	   f(i,K,j) = f(i,K,j) - c(i,j) * f(i,K,j+1)
	enddo
	enddo
	enddo

	if ( nbrlower .ne. MPI_PROC_NULL ) then
	   do K = L1, L2
	   do i = 1, n
	      smsg(i,K) = f(i,K,1)
	   enddo
	   enddo
	   call MPI_ISEND( smsg(1,1), n*L, MPI_DOUBLE_PRECISION, 
     <	                   nbrlower, 1, comm3d, req, ierr )
	   call MPI_WAIT( req, status, ierr )
	else
	   do K = L1, L2
	   do i = 1, n
	      f(i,K,0) = f(i,K,0) - c(i,0) * f(i,K,1)
	   enddo
	   enddo
	endif

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



