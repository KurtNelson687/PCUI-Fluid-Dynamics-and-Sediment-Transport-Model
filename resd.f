cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine residual( resid, bbsum, ermin, ermax, resbc, sumbc,
     <		level, ii, jj, kk, p, r, b, jac,
     <		g11, g12, g13, g21, g22, g23, g31, g32, g33, gcc,
     <		t11, t12, t13, t21, t22, t23, t31, t32, t33, tcc )

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "para.inc"

	integer level, ii, jj, kk
	double precision, dimension(0:ii+1,0:jj+1,0:kk+1) :: p, r, b
	
	double precision, dimension(-1:ii+2,-1:jj+2,-1:kk+2) :: jac,
     <		g11, g12, g13, g21, g22, g23, g31, g32, g33, gcc,
     <		t11, t12, t13, t21, t22, t23, t31, t32, t33, tcc

	integer i, j, k

	double precision resid, res, bbsum, rhs
	double precision ermin, ermax, ern, erx
	double precision resbc, rbc, sumbc, sbc

	resid = 0.D0
	res = 0.D0
	bbsum = 0.D0
	rhs = 0.D0
	ermin = 0.D0
	ermax = 0.D0
	ern = 0.D0
	erx = 0.D0
	resbc = 0.D0
	rbc = 0.D0
	sumbc = 0.D0
	sbc = 0.D0

	do k = 0, kk+1
	do j = 0, jj+1
	do i = 0, ii+1
	   r(i,j,k) = 0.D0
	enddo
	enddo
	enddo

	do k = 1, kk
	do j = 1, jj
	do i = 1, ii

	   r(i,j,k) = b(i,j,k) -
     <		( t11(i,  j,  k  ) * ( p(i+1,j,  k  ) - p(i,j,k) )
     <		+ t11(i-1,j,  k  ) * ( p(i-1,j,  k  ) - p(i,j,k) )
     <		+ t22(i,  j,  k  ) * ( p(i,  j+1,k  ) - p(i,j,k) )
     <		+ t22(i,  j-1,k  ) * ( p(i,  j-1,k  ) - p(i,j,k) )
     <		+ t33(i,  j,  k  ) * ( p(i,  j,  k+1) - p(i,j,k) )
     <		+ t33(i,  j,  k-1) * ( p(i,  j,  k-1) - p(i,j,k) ) )

	   r(i,j,k) = r(i,j,k) -
     <		( t12(i,  j,k) * ( p(i,  j+1,k) - p(i,  j-1,k) 
     <		                 + p(i+1,j+1,k) - p(i+1,j-1,k) )
     <		- t12(i-1,j,k) * ( p(i,  j+1,k) - p(i,  j-1,k)
     <		                 + p(i-1,j+1,k) - p(i-1,j-1,k) )
     <		+ t13(i,  j,k) * ( p(i,  j,k+1) - p(i,  j,k-1)
     <		                 + p(i+1,j,k+1) - p(i+1,j,k-1) )
     <		- t13(i-1,j,k) * ( p(i,  j,k+1) - p(i,  j,k-1)
     <		                 + p(i-1,j,k+1) - p(i-1,j,k-1) ) )

	   r(i,j,k) = r(i,j,k) -
     <		( t23(i,j,  k) * ( p(i,j,  k+1) - p(i,j,  k-1)
     <		                 + p(i,j+1,k+1) - p(i,j+1,k-1) )
     <		- t23(i,j-1,k) * ( p(i,j,  k+1) - p(i,j,  k-1)
     <		                 + p(i,j-1,k+1) - p(i,j-1,k-1) ) 
     <		+ t21(i,j,  k) * ( p(i+1,j,  k) - p(i-1,j,  k)
     <		                 + p(i+1,j+1,k) - p(i-1,j+1,k) )
     <		- t21(i,j-1,k) * ( p(i+1,j,  k) - p(i-1,j,  k)
     <		                 + p(i+1,j-1,k) - p(i-1,j-1,k) ) )

	   r(i,j,k) = r(i,j,k) -
     <		( t31(i,j,k  ) * ( p(i+1,j,k  ) - p(i-1,j,k  )
     <		                 + p(i+1,j,k+1) - p(i-1,j,k+1) )
     <		- t31(i,j,k-1) * ( p(i+1,j,k  ) - p(i-1,j,k  )
     <		                 + p(i+1,j,k-1) - p(i-1,j,k-1) )
     <		+ t32(i,j,k  ) * ( p(i,j+1,k  ) - p(i,j-1,k  )
     <		                 + p(i,j+1,k+1) - p(i,j-1,k+1) )
     <		- t32(i,j,k-1) * ( p(i,j+1,k  ) - p(i,j-1,k  )
     <		                 + p(i,j+1,k-1) - p(i,j-1,k-1) ) )

	   res = res + ( r(i,j,k) )**2
	   rhs = rhs + ( b(i,j,k) )**2
	   ern = min(ern, r(i,j,k))
	   erx = max(erx, r(i,j,k))

	enddo
	enddo
	enddo

	if ( n_west .eq. MPI_PROC_NULL ) then
	   i = 0
	   do k = 1, kk
	   do j = 1, jj
	      r(i,j,k) = b(i,j,k) - 
     <		( g11(i,  j,k) * ( p(i+1,j,  k) - p(i,  j,  k) )
     <		+ g12(i,  j,k) * ( p(i,  j+1,k) - p(i,  j-1,k) 
     <		                 + p(i+1,j+1,k) - p(i+1,j-1,k) )
     <		+ g13(i,  j,k) * ( p(i,  j,k+1) - p(i,  j,k-1)
     <		                 + p(i+1,j,k+1) - p(i+1,j,k-1) ) )
CC	      rbc = rbc + ( r(i,j,k) )**2
CC	      sbc = sbc + ( b(i,j,k) )**2
	   enddo
	   enddo
	endif

        if ( n_east .eq. MPI_PROC_NULL ) then
	   i = ii + 1
	   do k = 1, kk
	   do j = 1, jj
	      r(i,j,k) = b(i,j,k) -
     <		( g11(i-1,j,k) * ( p(i,  j,  k) - p(i-1,j,  k) )
     <		+ g12(i-1,j,k) * ( p(i,  j+1,k) - p(i,  j-1,k)
     <		                 + p(i-1,j+1,k) - p(i-1,j-1,k) )
     <		+ g13(i-1,j,k) * ( p(i,  j,k+1) - p(i,  j,k-1)
     <		                 + p(i-1,j,k+1) - p(i-1,j,k-1) ) )
CC	      rbc = rbc + ( r(i,j,k) )**2
CC	      sbc = sbc + ( b(i,j,k) )**2
	   enddo
	   enddo
        endif

        if ( n_suth .eq. MPI_PROC_NULL ) then
           j = 0
           do k = 1, kk
           do i = 1, ii
	      r(i,j,k) = b(i,j,k) -
     <		( g22(i,j,  k) * ( p(i,j+1,k  ) - p(i,j,  k  ) )
     <		+ g23(i,j,  k) * ( p(i,j,  k+1) - p(i,j,  k-1)
     <		                 + p(i,j+1,k+1) - p(i,j+1,k-1) )
     <		+ g21(i,j,  k) * ( p(i+1,j,  k) - p(i-1,j,  k)
     <		                 + p(i+1,j+1,k) - p(i-1,j+1,k) ) )
CC	      rbc = rbc + ( r(i,j,k) )**2
CC	      sbc = sbc + ( b(i,j,k) )**2
	   enddo
	   enddo
        endif

        if ( n_nrth .eq. MPI_PROC_NULL ) then
           j = jj + 1
           do k = 1, kk
           do i = 1, ii
	      r(i,j,k) = b(i,j,k) - 
     <		( g22(i,j-1,k) * ( p(i,j,  k  ) - p(i,j-1,k  ) )
     <		+ g23(i,j-1,k) * ( p(i,j,  k+1) - p(i,j,  k-1)
     <		                 + p(i,j-1,k+1) - p(i,j-1,k-1) ) 
     <		+ g21(i,j-1,k) * ( p(i+1,j,  k) - p(i-1,j,  k)
     <		                 + p(i+1,j-1,k) - p(i-1,j-1,k) ) )
CC	      rbc = rbc + ( r(i,j,k) )**2
CC	      sbc = sbc + ( b(i,j,k) )**2
	   enddo
	   enddo
        endif

        if ( n_back .eq. MPI_PROC_NULL ) then
           k = 0
           do j = 1, jj
           do i = 1, ii
	      r(i,j,k) = b(i,j,k) -
     <		( g33(i,j,k  ) * ( p(i,  j,k+1) - p(i,  j,k  ) )
     <		+ g31(i,j,k  ) * ( p(i+1,j,k  ) - p(i-1,j,k  )
     <		                 + p(i+1,j,k+1) - p(i-1,j,k+1) )
     <		+ g32(i,j,k  ) * ( p(i,j+1,k  ) - p(i,j-1,k  )
     <		                 + p(i,j+1,k+1) - p(i,j-1,k+1) ) )
CC	      rbc = rbc + ( r(i,j,k) )**2
CC	      sbc = sbc + ( b(i,j,k) )**2
	   enddo
	   enddo
        endif

        if ( n_frnt .eq. MPI_PROC_NULL ) then
           k = kk + 1
           do j = 1, jj
           do i = 1, ii
	      r(i,j,k) = b(i,j,k) -
     <		( g33(i,j,k-1) * ( p(i,  j,k  ) - p(i,  j,k-1) )
     <		+ g31(i,j,k-1) * ( p(i+1,j,k  ) - p(i-1,j,k  )
     <		                 + p(i+1,j,k-1) - p(i-1,j,k-1) )
     <		+ g32(i,j,k-1) * ( p(i,j+1,k  ) - p(i,j-1,k  )
     <		                 + p(i,j+1,k-1) - p(i,j-1,k-1) ) )
CC	      rbc = rbc + ( r(i,j,k) )**2
CC	      sbc = sbc + ( b(i,j,k) )**2
	   enddo
	   enddo
        endif
        
	call MPI_REDUCE(res,resid,1,MPI_DOUBLE_PRECISION,
     <	                MPI_SUM,0,comm3d,ierr)
	call MPI_REDUCE(rhs,bbsum,1,MPI_DOUBLE_PRECISION,
     <	                MPI_SUM,0,comm3d,ierr)
	call MPI_REDUCE(ern,ermin,1,MPI_DOUBLE_PRECISION,
     <	                MPI_MIN,0,comm3d,ierr)
	call MPI_REDUCE(erx,ermax,1,MPI_DOUBLE_PRECISION,
     <	                MPI_MAX,0,comm3d,ierr)
CC	call MPI_REDUCE(rbc,resbc,1,MPI_DOUBLE_PRECISION,
CC     <	                MPI_SUM,0,comm3d,ierr)
CC	call MPI_REDUCE(sbc,sumbc,1,MPI_DOUBLE_PRECISION,
CC     <	                MPI_SUM,0,comm3d,ierr)

	if ( MYID .eq. 0 ) then
	   resid = resid / dble(ii*px*jj*py*kk*pz)
	   resid = dsqrt(resid)
	   bbsum = bbsum / dble(ii*px*jj*py*kk*pz) 	   
	   bbsum = dsqrt(bbsum)
CC	   if ( periodic .eq. 1 ) then
CC	      resbc = resbc / dble(2*kk*pz*ii*px+2*jj*py*kk*pz)
CC	      sumbc = sumbc / dble(2*kk*pz*ii*px+2*jj*py*kk*pz)
CC	   else
CC	      resbc = resbc 
CC     <	            / dble(2*kk*pz*ii*px+2*jj*py*kk*pz+2*ii*px*jj*py)
CC	      sumbc = sumbc 
CC     <	            / dble(2*kk*pz*ii*px+2*jj*py*kk*pz+2*ii*px*jj*py)
CC	   endif
CC	   resbc = dsqrt(resbc)
CC	   sumbc = dsqrt(sumbc)
	endif

	call MPI_BCAST(resid,1,MPI_DOUBLE_PRECISION,0,comm3d,ierr)
	call MPI_BCAST(bbsum,1,MPI_DOUBLE_PRECISION,0,comm3d,ierr)
	call MPI_BCAST(ermin,1,MPI_DOUBLE_PRECISION,0,comm3d,ierr)
	call MPI_BCAST(ermax,1,MPI_DOUBLE_PRECISION,0,comm3d,ierr)
CC	call MPI_BCAST(resbc,1,MPI_DOUBLE_PRECISION,0,comm3d,ierr)
CC	call MPI_BCAST(sumbc,1,MPI_DOUBLE_PRECISION,0,comm3d,ierr)

	return
	end
