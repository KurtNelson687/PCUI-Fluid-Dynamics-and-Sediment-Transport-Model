cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine bc_prep( level, ii, jj, kk, i0, k0, id, kd, p, r, b,
     <		g11, g12, g13, g21, g22, g23, g31, g32, g33, gcc )

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "para.inc"

	integer level, ii, jj, kk
	integer i0, k0, id, kd

	double precision, dimension(0:ii+1,0:jj+1,0:kk+1) :: p, r, b
	
	double precision, dimension(-1:ii+2,-1:jj+2,-1:kk+2) ::
     <		g11, g12, g13, g21, g22, g23, g31, g32, g33, gcc

	logical dowest, doeast, doback, dofrnt
	integer i, j, k, is, ie, ks, ke

	if ( id .eq. 1 ) then
	   dowest = .true.
	   doeast = .true.
	   is = 1
	   ie = ii
	else
	   if ( i0 .eq. 1 ) then
	      dowest = .true.
	      doeast = .false.
	      is  = 1
	      ie  = ii-1
	   else
	      dowest = .false.
	      doeast = .true.
	      is  = 2
	      ie  = ii
	   endif
	endif
	
	if ( kd .eq. 1 ) then
	   doback = .true.
	   dofrnt = .true.
	   ks = 1
	   ke = kk
	else
	   if ( k0 .eq. 1 ) then
	      doback = .true.
	      dofrnt = .false.
	      ks  = 1
	      ke  = kk-1
	   else
	      doback = .false.
	      dofrnt = .true.
	      ks  = 2
	      ke  = kk
	   endif
	endif
	
	if ( dowest .and. n_west .eq. MPI_PROC_NULL ) then
	   i = 0
	   do k = ks, ke, kd
	   do j = 1, jj
	      r(i,j,k) = - b(i,j,k) +
     <		( g12(i,  j,k) * ( p(i,  j+1,k) - p(i,  j-1,k)
     <		                 + p(i+1,j+1,k) - p(i+1,j-1,k) )
     <		+ g13(i,  j,k) * ( p(i,  j,k+1) - p(i,  j,k-1)
     <		                 + p(i+1,j,k+1) - p(i+1,j,k-1) ) )
 	   enddo
	   enddo
	endif

        if ( doeast .and. n_east .eq. MPI_PROC_NULL ) then
	   i = ii + 1
	   do k = ks, ke, kd
	   do j = 1, jj
	      r(i,j,k) = - b(i,j,k) +
     <		( g12(i-1,j,k) * ( p(i,  j+1,k) - p(i,  j-1,k)
     <		                 + p(i-1,j+1,k) - p(i-1,j-1,k) )
     <		+ g13(i-1,j,k) * ( p(i,  j,k+1) - p(i,  j,k-1)
     <		                 + p(i-1,j,k+1) - p(i-1,j,k-1) ) )
 	   enddo
	   enddo
	endif

        if ( doback .and. n_back .eq. MPI_PROC_NULL ) then
           k = 0
           do j = 1, jj
           do i = is, ie, id
              r(i,j,k) = - b(i,j,k) +
     <		( g31(i,j,k) * ( p(i+1,j,k  ) - p(i-1,j,k  )
     <		               + p(i+1,j,k+1) - p(i-1,j,k+1) )
     <		+ g32(i,j,k) * ( p(i,j+1,k  ) - p(i,j-1,k  )
     <		               + p(i,j+1,k+1) - p(i,j-1,k+1) ) )
	   enddo
	   enddo
        endif

        if ( dofrnt .and. n_frnt .eq. MPI_PROC_NULL ) then
           k = kk + 1
           do j = 1, jj
           do i = is, ie, id
	      r(i,j,k) = - b(i,j,k) +
     <		( g31(i,j,k-1) * ( p(i+1,j,k  ) - p(i-1,j,k  )
     <		                 + p(i+1,j,k-1) - p(i-1,j,k-1) )
     <		+ g32(i,j,k-1) * ( p(i,j+1,k  ) - p(i,j-1,k  )
     <		                 + p(i,j+1,k-1) - p(i,j-1,k-1) ) )
	   enddo
	   enddo
        endif

	return
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine bc_updt( level, ii, jj, kk, i0, k0, id, kd, p, r, b,
     <		g11, g12, g13, g21, g22, g23, g31, g32, g33, gcc )

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "para.inc"

	integer level, ii, jj, kk
	integer i0, k0, id, kd
	
	double precision, dimension(0:ii+1,0:jj+1,0:kk+1) :: p, r, b
	
	double precision, dimension(-1:ii+2,-1:jj+2,-1:kk+2) ::
     <		g11, g12, g13, g21, g22, g23, g31, g32, g33, gcc

	double precision, dimension(ii/id) ::
     <		sEmsg, rEmsg, sWmsg, rWmsg

	double precision, dimension(kk/kd) ::
     <		sBmsg, rBmsg, sFmsg, rFmsg

	double precision, dimension(0:jj+1) ::
     <		sSmsg, rSmsg, sNmsg, rNmsg

	integer status(MPI_STATUS_SIZE,4), req(4), nreq, len

	logical dowest, doeast, doback, dofrnt
	integer i, j, k, is, ie, ks, ke, it, kt

	if ( id .eq. 1 ) then
	   dowest = .true.
	   doeast = .true.
	   is = 1
	   ie = ii
	else
	   if ( i0 .eq. 1 ) then
	      dowest = .true.
	      doeast = .false.
	      is  = 1
	      ie  = ii-1
	   else
	      dowest = .false.
	      doeast = .true.
	      is  = 2
	      ie  = ii
	   endif
	endif
	
	if ( kd .eq. 1 ) then
	   doback = .true.
	   dofrnt = .true.
	   ks = 1
	   ke = kk
	else
	   if ( k0 .eq. 1 ) then
	      doback = .true.
	      dofrnt = .false.
	      ks  = 1
	      ke  = kk-1
	   else
	      doback = .false.
	      dofrnt = .true.
	      ks  = 2
	      ke  = kk
	   endif
	endif
	
	if ( dowest .and. n_west .eq. MPI_PROC_NULL ) then
	   i = 0
	   do k = ks, ke, kd
	   do j = 1, jj
	      p(i,j,k) = p(i+1,j,k) + r(i,j,k) / g11(i,  j,k)
	   enddo
	   enddo
	endif

        if ( doeast .and. n_east .eq. MPI_PROC_NULL ) then
	   i = ii + 1
	   do k = ks, ke, kd
	   do j = 1, jj
	      p(i,j,k) = p(i-1,j,k) - r(i,j,k) / g11(i-1,j,k)
	   enddo
	   enddo
	endif

        if ( doback .and. n_back .eq. MPI_PROC_NULL ) then
           k = 0
           do j = 1, jj
           do i = is, ie, id
	      p(i,j,k) = p(i,j,k+1) + r(i,j,k) / g33(i,j,k  )
	   enddo
	   enddo
        endif

        if ( dofrnt .and. n_frnt .eq. MPI_PROC_NULL ) then
           k = kk + 1
           do j = 1, jj
           do i = is, ie, id
	      p(i,j,k) = p(i,j,k-1) - r(i,j,k) / g33(i,j,k-1)
	   enddo
	   enddo
        endif
        
CCC	call p_exchange( ii, jj, kk, p ) 
	
CCC	IF ( 0 .EQ. 1 ) THEN
       
	if ( dowest .and. n_west .eq. MPI_PROC_NULL ) then
	   i = 0
	   nreq = 0
	   len = kk / kd
	   if ( n_suth .ne. MPI_PROC_NULL ) then
	      nreq = nreq + 1
	      call MPI_IRECV( rBmsg(1), len, MPI_DOUBLE_PRECISION,  
     <	                      n_suth, 0, comm3d, req(nreq), ierr )
     	      kt = 0
	      do k = ks, ke, kd
	         kt = kt + 1
	         sBmsg(kt) = p(i,1,k)
	      enddo
	      nreq = nreq + 1
	      call MPI_ISEND( sBmsg(1), len, MPI_DOUBLE_PRECISION,  
     <	                      n_suth, 1, comm3d, req(nreq), ierr )
	   endif
	   if ( n_nrth .ne. MPI_PROC_NULL ) then
	      nreq = nreq + 1
	      call MPI_IRECV( rFmsg(1), len, MPI_DOUBLE_PRECISION,  
     <	                      n_nrth, 1, comm3d, req(nreq), ierr )
     	      kt = 0
	      do k = ks, ke, kd
	         kt = kt + 1
	         sFmsg(kt) = p(i,jj,k)
	      enddo
	      nreq = nreq + 1
	      call MPI_ISEND( sFmsg(1), len, MPI_DOUBLE_PRECISION,  
     <	                      n_nrth, 0, comm3d, req(nreq), ierr )
	   endif
	   call MPI_WAITALL( nreq, req, status, ierr )
	   if ( n_suth .ne. MPI_PROC_NULL ) then
     	      kt = 0
	      do k = ks, ke, kd
	         kt = kt + 1
	         p(i,0,k) = rBmsg(kt)
	      enddo
	   endif
	   if ( n_nrth .ne. MPI_PROC_NULL ) then
     	      kt = 0
	      do k = ks, ke, kd
	         kt = kt + 1
	         p(i,jj+1,k) = rFmsg(kt)
	      enddo
	   endif
	   nreq = 0
	   len = jj + 2
	   if ( doback .and. n_back .ne. MPI_PROC_NULL ) then
	      do j = 0, jj+1
	         sSmsg(j) = p(i,j,1)
	      enddo
	      nreq = nreq + 1
	      call MPI_ISEND( sSmsg(0), len, MPI_DOUBLE_PRECISION,  
     <	                      n_back, 1, comm3d, req(nreq), ierr )
	   endif
	   if ( doback .and. n_frnt .ne. MPI_PROC_NULL ) then
	      nreq = nreq + 1
	      call MPI_IRECV( rNmsg(0), len, MPI_DOUBLE_PRECISION,  
     <	                      n_frnt, 1, comm3d, req(nreq), ierr )
	   endif
	   call MPI_WAITALL( nreq, req, status, ierr )
	   if ( doback .and. n_frnt .ne. MPI_PROC_NULL ) then
	      do j = 0, jj+1
	         p(i,j,kk+1) = rNmsg(j)
	      enddo
	   endif
	   nreq = 0
	   len = jj + 2
	   if ( dofrnt .and. n_back .ne. MPI_PROC_NULL ) then
	      nreq = nreq + 1
	      call MPI_IRECV( rSmsg(0), len, MPI_DOUBLE_PRECISION,  
     <	                      n_back, 0, comm3d, req(nreq), ierr )
	   endif
	   if ( dofrnt .and. n_frnt .ne. MPI_PROC_NULL ) then
	      do j = 0, jj+1
	         sNmsg(j) = p(i,j,kk)
	      enddo
	      nreq = nreq + 1
	      call MPI_ISEND( sNmsg(0), len, MPI_DOUBLE_PRECISION,  
     <	                      n_frnt, 0, comm3d, req(nreq), ierr )
	   endif
	   call MPI_WAITALL( nreq, req, status, ierr )
	   if ( dofrnt .and. n_back .ne. MPI_PROC_NULL ) then
	      do j = 0, jj+1
	         p(i,j,0) = rSmsg(j)
	      enddo
	   endif
	endif

	if ( doeast .and. n_east .eq. MPI_PROC_NULL ) then
	   i = ii + 1
	   nreq = 0
	   len = kk / kd
	   if ( n_suth .ne. MPI_PROC_NULL ) then
	      nreq = nreq + 1
	      call MPI_IRECV( rBmsg(1), len, MPI_DOUBLE_PRECISION,  
     <	                      n_suth, 0, comm3d, req(nreq), ierr )
     	      kt = 0
	      do k = ks, ke, kd
	         kt = kt + 1
	         sBmsg(kt) = p(i,1,k)
	      enddo
	      nreq = nreq + 1
	      call MPI_ISEND( sBmsg(1), len, MPI_DOUBLE_PRECISION,  
     <	                      n_suth, 1, comm3d, req(nreq), ierr )
	   endif
	   if ( n_nrth .ne. MPI_PROC_NULL ) then
	      nreq = nreq + 1
	      call MPI_IRECV( rFmsg(1), len, MPI_DOUBLE_PRECISION,  
     <	                      n_nrth, 1, comm3d, req(nreq), ierr )
     	      kt = 0
	      do k = ks, ke, kd
	         kt = kt + 1
	         sFmsg(kt) = p(i,jj,k)
	      enddo
	      nreq = nreq + 1
	      call MPI_ISEND( sFmsg(1), len, MPI_DOUBLE_PRECISION,  
     <	                      n_nrth, 0, comm3d, req(nreq), ierr )
	   endif
	   call MPI_WAITALL( nreq, req, status, ierr )
	   if ( n_suth .ne. MPI_PROC_NULL ) then
     	      kt = 0
	      do k = ks, ke, kd
	         kt = kt + 1
	         p(i,0,k) = rBmsg(kt)
	      enddo
	   endif
	   if ( n_nrth .ne. MPI_PROC_NULL ) then
     	      kt = 0
	      do k = ks, ke, kd
	         kt = kt + 1
	         p(i,jj+1,k) = rFmsg(kt)
	      enddo
	   endif
	   nreq = 0
	   len = jj + 2
	   if ( doback .and. n_back .ne. MPI_PROC_NULL ) then
	      do j = 0, jj+1
	         sSmsg(j) = p(i,j,1)
	      enddo
	      nreq = nreq + 1
	      call MPI_ISEND( sSmsg(0), len, MPI_DOUBLE_PRECISION,  
     <	                      n_back, 1, comm3d, req(nreq), ierr )
	   endif
	   if ( doback .and. n_frnt .ne. MPI_PROC_NULL ) then
	      nreq = nreq + 1
	      call MPI_IRECV( rNmsg(0), len, MPI_DOUBLE_PRECISION,  
     <	                      n_frnt, 1, comm3d, req(nreq), ierr )
	   endif
	   call MPI_WAITALL( nreq, req, status, ierr )
	   if ( doback .and. n_frnt .ne. MPI_PROC_NULL ) then
	      do j = 0, jj+1
	         p(i,j,kk+1) = rNmsg(j)
	      enddo
	   endif
	   nreq = 0
	   len = jj + 2
	   if ( dofrnt .and. n_back .ne. MPI_PROC_NULL ) then
	      nreq = nreq + 1
	      call MPI_IRECV( rSmsg(0), len, MPI_DOUBLE_PRECISION,  
     <	                      n_back, 0, comm3d, req(nreq), ierr )
	   endif
	   if ( dofrnt .and. n_frnt .ne. MPI_PROC_NULL ) then
	      do j = 0, jj+1
	         sNmsg(j) = p(i,j,kk)
	      enddo
	      nreq = nreq + 1
	      call MPI_ISEND( sNmsg(0), len, MPI_DOUBLE_PRECISION,  
     <	                      n_frnt, 0, comm3d, req(nreq), ierr )
	   endif
	   call MPI_WAITALL( nreq, req, status, ierr )
	   if ( dofrnt .and. n_back .ne. MPI_PROC_NULL ) then
	      do j = 0, jj+1
	         p(i,j,0) = rSmsg(j)
	      enddo
	   endif
	endif

	if ( doback .and. n_back .eq. MPI_PROC_NULL ) then
	   k = 0
	   nreq = 0
	   len = ii / id
	   if ( n_suth .ne. MPI_PROC_NULL ) then
	      nreq = nreq + 1
	      call MPI_IRECV( rWmsg(1), len, MPI_DOUBLE_PRECISION,  
     <	                      n_suth, 0, comm3d, req(nreq), ierr )
     	      it = 0
	      do i = is, ie, id
	         it = it + 1
	         sWmsg(it) = p(i,1,k)
	      enddo
	      nreq = nreq + 1
	      call MPI_ISEND( sWmsg(1), len, MPI_DOUBLE_PRECISION,  
     <	                      n_suth, 1, comm3d, req(nreq), ierr )
	   endif
	   if ( n_nrth .ne. MPI_PROC_NULL ) then
	      nreq = nreq + 1
	      call MPI_IRECV( rEmsg(1), len, MPI_DOUBLE_PRECISION,  
     <	                      n_nrth, 1, comm3d, req(nreq), ierr )
     	      it = 0
	      do i = is, ie, id
	         it = it + 1
	         sEmsg(it) = p(i,jj,k)
	      enddo
	      nreq = nreq + 1
	      call MPI_ISEND( sEmsg(1), len, MPI_DOUBLE_PRECISION,  
     <	                      n_nrth, 0, comm3d, req(nreq), ierr )
	   endif
	   call MPI_WAITALL( nreq, req, status, ierr )
	   if ( n_suth .ne. MPI_PROC_NULL ) then
     	      it = 0
	      do i = is, ie, id
	         it = it + 1
	         p(i,0,k) = rWmsg(it)
	      enddo
	   endif
	   if ( n_nrth .ne. MPI_PROC_NULL ) then
     	      it = 0
	      do i = is, ie, id
	         it = it + 1
	         p(i,jj+1,k) = rEmsg(it)
	      enddo
	   endif
	   nreq = 0
	   len = jj + 2
	   if ( dowest .and. n_west .ne. MPI_PROC_NULL ) then
	      do j = 0, jj+1
	         sSmsg(j) = p(1,j,k)
	      enddo
	      nreq = nreq + 1
	      call MPI_ISEND( sSmsg(0), len, MPI_DOUBLE_PRECISION,  
     <	                      n_west, 1, comm3d, req(nreq), ierr )
	   endif
	   if ( dowest .and. n_east .ne. MPI_PROC_NULL ) then
	      nreq = nreq + 1
	      call MPI_IRECV( rNmsg(0), len, MPI_DOUBLE_PRECISION,  
     <	                      n_east, 1, comm3d, req(nreq), ierr )
	   endif
	   call MPI_WAITALL( nreq, req, status, ierr )
	   if ( dowest .and. n_east .ne. MPI_PROC_NULL ) then
	      do j = 0, jj+1
	         p(ii+1,j,k) = rNmsg(j)
	      enddo
	   endif
	   nreq = 0
	   len = jj + 2
	   if ( doeast .and. n_west .ne. MPI_PROC_NULL ) then
	      nreq = nreq + 1
	      call MPI_IRECV( rSmsg(0), len, MPI_DOUBLE_PRECISION,  
     <	                      n_west, 0, comm3d, req(nreq), ierr )
	   endif
	   if ( doeast .and. n_east .ne. MPI_PROC_NULL ) then
	      do j = 0, jj+1
	         sNmsg(j) = p(ii,j,k)
	      enddo
	      nreq = nreq + 1
	      call MPI_ISEND( sNmsg(0), len, MPI_DOUBLE_PRECISION,  
     <	                      n_east, 0, comm3d, req(nreq), ierr )
	   endif
	   call MPI_WAITALL( nreq, req, status, ierr )
	   if ( doeast .and. n_west .ne. MPI_PROC_NULL ) then
	      do j = 0, jj+1
	         p(0,j,k) = rSmsg(j)
	      enddo
	   endif
	endif

	if ( dofrnt .and. n_frnt .eq. MPI_PROC_NULL ) then
	   k = kk + 1
	   nreq = 0
	   len = ii / id
	   if ( n_suth .ne. MPI_PROC_NULL ) then
	      nreq = nreq + 1
	      call MPI_IRECV( rWmsg(1), len, MPI_DOUBLE_PRECISION,  
     <	                      n_suth, 0, comm3d, req(nreq), ierr )
     	      it = 0
	      do i = is, ie, id
	         it = it + 1
	         sWmsg(it) = p(i,1,k)
	      enddo
	      nreq = nreq + 1
	      call MPI_ISEND( sWmsg(1), len, MPI_DOUBLE_PRECISION,  
     <	                      n_suth, 1, comm3d, req(nreq), ierr )
	   endif
	   if ( n_nrth .ne. MPI_PROC_NULL ) then
	      nreq = nreq + 1
	      call MPI_IRECV( rEmsg(1), len, MPI_DOUBLE_PRECISION,  
     <	                      n_nrth, 1, comm3d, req(nreq), ierr )
     	      it = 0
	      do i = is, ie, id
	         it = it + 1
	         sEmsg(it) = p(i,jj,k)
	      enddo
	      nreq = nreq + 1
	      call MPI_ISEND( sEmsg(1), len, MPI_DOUBLE_PRECISION,  
     <	                      n_nrth, 0, comm3d, req(nreq), ierr )
	   endif
	   call MPI_WAITALL( nreq, req, status, ierr )
	   if ( n_suth .ne. MPI_PROC_NULL ) then
     	      it = 0
	      do i = is, ie, id
	         it = it + 1
	         p(i,0,k) = rWmsg(it)
	      enddo
	   endif
	   if ( n_nrth .ne. MPI_PROC_NULL ) then
     	      it = 0
	      do i = is, ie, id
	         it = it + 1
	         p(i,jj+1,k) = rEmsg(it)
	      enddo
	   endif
	   nreq = 0
	   len = jj + 2
	   if ( dowest .and. n_west .ne. MPI_PROC_NULL ) then
	      do j = 0, jj+1
	         sSmsg(j) = p(1,j,k)
	      enddo
	      nreq = nreq + 1
	      call MPI_ISEND( sSmsg(0), len, MPI_DOUBLE_PRECISION,  
     <	                      n_west, 1, comm3d, req(nreq), ierr )
	   endif
	   if ( dowest .and. n_east .ne. MPI_PROC_NULL ) then
	      nreq = nreq + 1
	      call MPI_IRECV( rNmsg(0), len, MPI_DOUBLE_PRECISION,  
     <	                      n_east, 1, comm3d, req(nreq), ierr )
	   endif
	   call MPI_WAITALL( nreq, req, status, ierr )
	   if ( dowest .and. n_east .ne. MPI_PROC_NULL ) then
	      do j = 0, jj+1
	         p(ii+1,j,k) = rNmsg(j)
	      enddo
	   endif
	   nreq = 0
	   len = jj + 2
	   if ( doeast .and. n_west .ne. MPI_PROC_NULL ) then
	      nreq = nreq + 1
	      call MPI_IRECV( rSmsg(0), len, MPI_DOUBLE_PRECISION,  
     <	                      n_west, 0, comm3d, req(nreq), ierr )
	   endif
	   if ( doeast .and. n_east .ne. MPI_PROC_NULL ) then
	      do j = 0, jj+1
	         sNmsg(j) = p(ii,j,k)
	      enddo
	      nreq = nreq + 1
	      call MPI_ISEND( sNmsg(0), len, MPI_DOUBLE_PRECISION,  
     <	                      n_east, 0, comm3d, req(nreq), ierr )
	   endif
	   call MPI_WAITALL( nreq, req, status, ierr )
	   if ( doeast .and. n_west .ne. MPI_PROC_NULL ) then
	      do j = 0, jj+1
	         p(0,j,k) = rSmsg(j)
	      enddo
	   endif
	endif

CCC	ENDIF

	if ( n_west .eq. MPI_PROC_NULL ) then
	   if ( n_suth .eq. MPI_PROC_NULL ) then
	      do k = 0, kk+1
C******	0 - 0 - *
 	p(0,0,k) = 0.5D0 * ( p(1,0,k) + p(0,1,k) )
	      enddo
	   endif
	   if ( n_nrth .eq. MPI_PROC_NULL ) then
	      do k = 0, kk+1
C******	0 - Y - *
 	p(0,jj+1,k) = 0.5D0 * ( p(1,jj+1,k) + p(0,jj,k) )
	      enddo
	   endif
	   if ( n_back .eq. MPI_PROC_NULL ) then
	      do j = 0, jj+1
C******	0 - * - 0
	p(0,j,0) = 0.5D0 * ( p(1,j,0) + p(0,j,1) )
	      enddo
	   endif
	   if ( n_frnt .eq. MPI_PROC_NULL ) then
	      do j = 0, jj+1
C******	0 - * - Z
	p(0,j,kk+1) = 0.5D0 * ( p(1,j,kk+1) + p(0,j,kk) )
	      enddo
	   endif
	endif	      

	if ( n_east .eq. MPI_PROC_NULL ) then
	   if ( n_suth .eq. MPI_PROC_NULL ) then
	      do k = 0, kk+1
C******	X - 0 - *
CC	p(ii+1,0,k) = p(ii,0,k) + p(ii+1,1,k) - p(ii,1,k)
	p(ii+1,0,k) = 0.5D0 * ( p(ii,0,k) + p(ii+1,1,k) )
	      enddo
	   endif
	   if ( n_nrth .eq. MPI_PROC_NULL ) then
	      do k = 0, kk+1
C******	X - Y - *
	p(ii+1,jj+1,k) = 0.5D0 * ( p(ii,jj+1,k) + p(ii+1,jj,k) )
	      enddo
	   endif
	   if ( n_back .eq. MPI_PROC_NULL ) then
	      do j = 0, jj+1
C******	X - * - 0
 	p(ii+1,j,0) = 0.5D0 * ( p(ii,j,0) + p(ii+1,j,1) )
	      enddo
	   endif
	   if ( n_frnt .eq. MPI_PROC_NULL ) then
	      do j = 0, jj+1
C******	X - * - Z
 	p(ii+1,j,kk+1) = 0.5D0 * ( p(ii,j,kk+1) + p(ii+1,j,kk) )
	      enddo
	   endif
	endif	
	
	if ( n_suth .eq. MPI_PROC_NULL ) then
	   if ( n_back .eq. MPI_PROC_NULL ) then
	      do i = 0, ii+1
C******	* - 0 - 0
 	p(i,0,0) = 0.5D0 * ( p(i,1,0) + p(i,0,1) )
	      enddo
	   endif
	   if ( n_frnt .eq. MPI_PROC_NULL ) then
	      do i = 0, ii+1
C******	* - 0 - Z
 	p(i,0,kk+1) = 0.5D0 * ( p(i,1,kk+1) + p(i,0,kk) )
	      enddo
	   endif
	endif      

	if ( n_nrth .eq. MPI_PROC_NULL ) then
	   if ( n_back .eq. MPI_PROC_NULL ) then
	      do i = 0, ii+1
C******	* - Y - 0
 	p(i,jj+1,0) = 0.5D0 * ( p(i,jj,0) + p(i,jj+1,1) )
	      enddo
	   endif
	   if ( n_frnt .eq. MPI_PROC_NULL ) then
	      do i = 0, ii+1
C******	* - Y - Z
 	p(i,jj+1,kk+1) = 0.5D0 * ( p(i,jj,kk+1) + p(i,jj+1,kk) )
	      enddo
	   endif
	endif  

	return
	end

