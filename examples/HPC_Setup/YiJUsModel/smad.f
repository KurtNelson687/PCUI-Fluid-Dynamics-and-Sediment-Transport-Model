cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine smooth( resid, bbsum, ermin, ermax, resbc, sumbc,
     <		level, ii, jj, kk, q, r, b, jac,
     <		g11, g12, g13, g21, g22, g23, g31, g32, g33, gcc )

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "para.inc"
	include "sedi.inc"

	integer level, ii, jj, kk 
	double precision q(0:ii+1,0:jj+1,0:kk+1), 
     &                   r(0:ii+1,0:jj+1,0:kk+1), 
     &                   b(0:ii+1,0:jj+1,0:kk+1)
	
	double precision jac(-1:ii+2,-1:jj+2,-1:kk+2),
     &                   g11(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   g12(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   g13(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   g21(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   g22(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   g23(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   g31(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   g32(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   g33(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   gcc(-1:ii+2,-1:jj+2,-1:kk+2),
     &                   t11(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   t12(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   t13(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   t21(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   t22(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   t23(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   t31(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   t32(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   t33(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   tcc(-1:ii+2,-1:jj+2,-1:kk+2)

	integer m, n

	double precision resid, bbsum
	double precision ermin, ermax, erabs
	double precision resbc, sumbc
	double precision res_old, res_new
	     
	integer i, j, k

	do k = -1, kk+2
	do j = -1, jj+2
	do i = -1, ii+2
	   t11(i,j,k) = g11(i,j,k)
	   t12(i,j,k) = g12(i,j,k)
	   t13(i,j,k) = g13(i,j,k)
	   t21(i,j,k) = g21(i,j,k)
	   t22(i,j,k) = g22(i,j,k)
	   t23(i,j,k) = g23(i,j,k)
	   t31(i,j,k) = g31(i,j,k)
	   t32(i,j,k) = g32(i,j,k)
	   t33(i,j,k) = g33(i,j,k)
	   tcc(i,j,k) = 0.D0
	enddo
	enddo
	enddo

	if ( n_west .eq. MPI_PROC_NULL ) then
	   i = 0
	   do k = -1, kk+2
	   do j = -1, jj+2
	      t11(i,j,k) = 0.D0
	      t12(i,j,k) = 0.D0
	      t13(i,j,k) = 0.D0
	   enddo
	   enddo
	endif

	if ( n_east .eq. MPI_PROC_NULL ) then
	   i = ii
	   do k = -1, kk+2
	   do j = -1, jj+2
	      t11(i,j,k) = 0.D0
	      t12(i,j,k) = 0.D0
	      t13(i,j,k) = 0.D0
	   enddo
	   enddo
	endif

	if ( n_suth .eq. MPI_PROC_NULL ) then
	   j = 0
	   do k = -1, kk+2
	   do i = -1, ii+2
	      t21(i,j,k) = 0.D0
	      t22(i,j,k) = 0.D0
	      t23(i,j,k) = 0.D0
	   enddo
	   enddo
	endif

	if ( n_nrth .eq. MPI_PROC_NULL ) then
	   j = jj
	   do k = -1, kk+2
	   do i = -1, ii+2
	      t21(i,j,k) = 0.D0
	      t22(i,j,k) = 0.D0
	      t23(i,j,k) = 0.D0
	   enddo
	   enddo
	endif

	if ( n_back .eq. MPI_PROC_NULL ) then
	   k = 0
	   do j = -1, jj+2
	   do i = -1, ii+2
	      t31(i,j,k) = 0.D0
	      t32(i,j,k) = 0.D0
	      t33(i,j,k) = 0.D0
	   enddo
	   enddo
	endif

	if ( n_frnt .eq. MPI_PROC_NULL ) then
	   k = kk
	   do j = -1, jj+2
	   do i = -1, ii+2
	      t31(i,j,k) = 0.D0
	      t32(i,j,k) = 0.D0
	      t33(i,j,k) = 0.D0
	   enddo
	   enddo
	endif

        do k = 1, kk
        do j = 1, jj
        do i = 1, ii
           tcc(i,j,k) = t11(i, j, k)*pcorr(i,j,k,1)
     <                + t11(i-1, j,   k  )*pcorr(i-1,j,k,1)
     <                + t22(i, j, k)*pcorr(i,j,k,2)
     <                + t22(i,   j-1, k  )*pcorr(i,j-1,k,2)
     <                + t33(i, j, k)*pcorr(i,j,k,3)
     <                + t33(i,   j,   k-1)*pcorr(i,j,k-1,3)
        enddo 
        enddo
        enddo

	resid = 0.D0
	res_old = 0.D0
	res_new = 0.D0
	
	call residual_NB( resid, bbsum, ermin, ermax, resbc, sumbc,
     <		level, ii, jj, kk, q, r, b, jac,
     <		g11, g12, g13, g21, g22, g23, g31, g32, g33, gcc,
     <		t11, t12, t13, t21, t22, t23, t31, t32, t33, tcc )
	erabs = max(dabs(ermin), dabs(ermax))
	res_old = resid
     	IF (  MOD(ISTEP,NSAVE) .EQ. 0 .and. MYID .EQ. 0 ) 
     <	WRITE(*,1) LEVEL, RESID, BBSUM, ERMIN, ERMAX
    1	FORMAT(' L', I1,           ' START WITH    ', 2E9.3,
     <	       ' E ', 2E9.3 )
CC     <	WRITE(*,1) LEVEL, RESID, BBSUM, ERMIN, ERMAX, RESBC, SUMBC
CC    1	FORMAT(' L', I1,           ' START WITH    ', 2E9.3, 
CC     <	       ' E', 2E9.3, ' BC', 2E9.3)

	do n = 1, maxiter(level)

 	do m = 1, iterchk(level)

	call smoothj( level, ii, jj, kk, q, r, b,
     <		g11, g12, g13, g21, g22, g23, g31, g32, g33, gcc,
     <		t11, t12, t13, t21, t22, t23, t31, t32, t33, tcc )

	call smoothk( level, ii, jj, kk, q, r, b,
     <		g11, g12, g13, g21, g22, g23, g31, g32, g33, gcc,
     <		t11, t12, t13, t21, t22, t23, t31, t32, t33, tcc )

	call smoothi( level, ii, jj, kk, q, r, b,
     <		g11, g12, g13, g21, g22, g23, g31, g32, g33, gcc,
     <		t11, t12, t13, t21, t22, t23, t31, t32, t33, tcc )

 	enddo

	call residual_NB( resid, bbsum, ermin, ermax, resbc, sumbc,
     <		level, ii, jj, kk, q, r, b, jac,
     <		g11, g12, g13, g21, g22, g23, g31, g32, g33, gcc,
     <		t11, t12, t13, t21, t22, t23, t31, t32, t33, tcc )
 	erabs = max(dabs(ermin), dabs(ermax))
 	res_new = resid

CC	if ( res_new .lt. tol(level) .and. erabs .lt. ter(level) ) then
CC	IF ( MOD(ISTEP,NSAVE) .EQ. 0 .and. MYID .EQ. 0 ) 
CC  <	WRITE(*,2) LEVEL, N, RESID, RESID/BBSUM, ERMIN, ERMAX, 
CC  <	                     RESBC, RESBC/SUMBC
CC 2	FORMAT(' L', I1, ' #', I2, ' CONVERGED ', 2E9.3, 
CC  <	       ' E', 2E9.3, ' BC', 2E9.3)
CC	return
CC	endif
     	
 	if ( res_new / res_old .gt. slowiter(level) ) then
 	IF ( MOD(ISTEP,NSAVE) .EQ. 0 .and. MYID .EQ. 0 ) 
     <	WRITE(*,3) LEVEL, N, RESID, RESID/BBSUM, ERMIN, ERMAX
    3	FORMAT(' L', I1, ' #', I2, ' TOO SLOW! ', 2E9.3, ' E ', 2E9.3)
CC     <	WRITE(*,3) LEVEL, N, RESID, RESID/BBSUM, ERMIN, ERMAX, 
CC     <	                     RESBC, RESBC/SUMBC
CC    3	FORMAT(' L', I1, ' #', I2, ' TOO SLOW! ', 2E9.3, 
CC     <	       ' E', 2E9.3, ' BC', 2E9.3)
	   return
	else
	   res_old = res_new
	endif
     
	enddo
	
     	IF ( MOD(ISTEP,NSAVE) .EQ. 0 .and. MYID .EQ. 0 ) 
     <	WRITE(*,4) LEVEL, N-1, RESID, RESID/BBSUM, ERMIN, ERMAX 
    4	FORMAT(' L', I1, ' #', I2, ' MAX STEP! ', 2E9.3, ' E ', 2E9.3)
CC     <	WRITE(*,4) LEVEL, N-1, RESID, RESID/BBSUM, ERMIN, ERMAX, 
CC     <	                     RESBC, RESBC/SUMBC
CC    4	FORMAT(' L', I1, ' #', I2, ' MAX STEP! ', 2E9.3, 
CC     <	       ' E', 2E9.3, ' BC', 2E9.3)
    
	return
	end	

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine smoothi( level, ii, jj, kk, q, r, b,
     <		g11, g12, g13, g21, g22, g23, g31, g32, g33, gcc,
     <		t11, t12, t13, t21, t22, t23, t31, t32, t33, tcc )

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "para.inc"
	include "sedi.inc"

	integer level, ii, jj, kk
	double precision q(0:ii+1,0:jj+1,0:kk+1), 
     &                   r(0:ii+1,0:jj+1,0:kk+1), 
     &                   b(0:ii+1,0:jj+1,0:kk+1)
	
	double precision g11(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   g12(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   g13(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   g21(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   g22(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   g23(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   g31(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   g32(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   g33(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   gcc(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   t11(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   t12(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   t13(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   t21(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   t22(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   t23(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   t31(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   t32(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   t33(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   tcc(-1:ii+2,-1:jj+2,-1:kk+2)

	double precision aji(jj,0:ii+1), bji(jj,0:ii+1), 
     &                   cji(jj,0:ii+1), fji(jj,0:ii+1)

	integer i, j, k

	do k = 1, kk
	do j = 1, jj
	do i = 1, ii

	   r(i,j,k) = - b(i,j,k) + 
     <		( t22(i,  j,  k  ) * q(i,  j+1,k  )*pcorr(i,j,k,2)
     <		+ t22(i,  j-1,k  ) * q(i,  j-1,k  )*pcorr(i,j-1,k,2)
     <		+ t33(i,  j,  k  ) * q(i,  j,  k+1)*pcorr(i,j,k,3)
     <		+ t33(i,  j,  k-1) * q(i,  j,  k-1)*pcorr(i,j,k-1,3) )

	   r(i,j,k) = r(i,j,k) +
     <		( t12(i,  j,k) * ( q(i,  j+1,k) - q(i,  j-1,k) 
     <		                 + q(i+1,j+1,k) - q(i+1,j-1,k) )
     <          *0.25D0*(pcorr(i,j,k,2) + pcorr(i,j-1,k,2)
     <                  +pcorr(i+1,j,k,2) + pcorr(i+1,j-1,k,2) )        
     <		- t12(i-1,j,k) * ( q(i,  j+1,k) - q(i,  j-1,k)
     <		                 + q(i-1,j+1,k) - q(i-1,j-1,k) )
     <          *0.25D0*(pcorr(i,j,k,2) + pcorr(i,j-1,k,2)
     <                  +pcorr(i-1,j,k,2) + pcorr(i-1,j-1,k,2) )
     <		+ t13(i,  j,k) * ( q(i,  j,k+1) - q(i,  j,k-1)
     <		                 + q(i+1,j,k+1) - q(i+1,j,k-1) )
     <          *0.25D0*(pcorr(i,j,k,3) + pcorr(i,j,k-1,3)
     <                  +pcorr(i+1,j,k,3) + pcorr(i+1,j,k-1,3) )
     <		- t13(i-1,j,k) * ( q(i,  j,k+1) - q(i,  j,k-1)
     <		                 + q(i-1,j,k+1) - q(i-1,j,k-1) )
     <          *0.25D0*(pcorr(i,j,k,3) + pcorr(i,j,k-1,3)
     <                  +pcorr(i-1,j,k,3) + pcorr(i-1,j,k-1,3) ) )

	   r(i,j,k) = r(i,j,k) +
     <		( t23(i,j,  k) * ( q(i,j,  k+1) - q(i,j,  k-1)
     <		                 + q(i,j+1,k+1) - q(i,j+1,k-1) )
     <          *0.25D0*(pcorr(i,j,k,3) + pcorr(i,j,k-1,3)
     <                  +pcorr(i,j+1,k,3) + pcorr(i,j+1,k-1,3) )
     <		- t23(i,j-1,k) * ( q(i,j,  k+1) - q(i,j,  k-1)
     <		                 + q(i,j-1,k+1) - q(i,j-1,k-1) )
     <          *0.25D0*(pcorr(i,j,k,3) + pcorr(i,j,k-1,3)
     <                  +pcorr(i,j-1,k,3) + pcorr(i,j-1,k-1,3) )
     <		+ t21(i,j,  k) * ( q(i+1,j+1,k) - q(i-1,j+1,k) )
     <          *0.5D0*(pcorr(i,j+1,k,1) + pcorr(i-1,j+1,k,1))
     <		- t21(i,j-1,k) * ( q(i+1,j-1,k) - q(i-1,j-1,k) ) 
     <          *0.5D0*(pcorr(i,j-1,k,1) + pcorr(i-1,j-1,k,1)  ) )


	   r(i,j,k) = r(i,j,k) +
     <		( t31(i,j,k  ) * ( q(i+1,j,k+1) - q(i-1,j,k+1) )
     <          *0.5D0*(pcorr(i,j,k+1,1) + pcorr(i-1,j,k+1,1) )
     <		- t31(i,j,k-1) * ( q(i+1,j,k-1) - q(i-1,j,k-1) )
     <          *0.5D0*(pcorr(i,j,k-1,1) + pcorr(i-1,j,k-1,1) )
     <		+ t32(i,j,k  ) * ( q(i,j+1,k  ) - q(i,j-1,k  )
     <		                 + q(i,j+1,k+1) - q(i,j-1,k+1) )
     <          *0.25D0*(pcorr(i,j,k,2) + pcorr(i,j-1,k,2)
     <                  +pcorr(i,j,k+1,2) + pcorr(i,j-1,k+1,2) )
     <		- t32(i,j,k-1) * ( q(i,j+1,k  ) - q(i,j-1,k  )
     <		                 + q(i,j+1,k-1) - q(i,j-1,k-1) )
     <          *0.25D0*(pcorr(i,j,k,2) + pcorr(i,j-1,k,2)
     <                  +pcorr(i,j,k-1,2) + pcorr(i,j-1,k-1,2) ) )
 
	enddo
	enddo
	enddo
	
	if ( n_west .eq. MPI_PROC_NULL ) then
	   i = 0
	   do k = 1, kk
	   do j = 1, jj
	      r(i,j,k) = - b(i,j,k) + 
     <		( g12(i,  j,k) * ( q(i,  j+1,k) - q(i,  j-1,k) 
     <		                 + q(i+1,j+1,k) - q(i+1,j-1,k) )
     <	        *0.25D0*(pcorr(i,j,k,2) + pcorr(i,j-1,k,2)
     <                  +pcorr(i+1,j,k,2) + pcorr(i+1,j-1,k,2) )        
     <		+ g13(i,  j,k) * ( q(i,  j,k+1) - q(i,  j,k-1)
     <		                 + q(i+1,j,k+1) - q(i+1,j,k-1) )
     <          *0.25D0*(pcorr(i,j,k,3) + pcorr(i,j,k-1,3)
     <                  +pcorr(i+1,j,k,3) + pcorr(i+1,j,k-1,3) ))
	   enddo
	   enddo
	endif

        if ( n_east .eq. MPI_PROC_NULL ) then
	   i = ii + 1
	   do k = 1, kk
	   do j = 1, jj
	      r(i,j,k) = - b(i,j,k) +
     <		( g12(i-1,j,k) * ( q(i,  j+1,k) - q(i,  j-1,k)
     <		                 + q(i-1,j+1,k) - q(i-1,j-1,k) )
     <          *0.25D0*(pcorr(i,j,k,2) + pcorr(i,j-1,k,2)
     <                  +pcorr(i-1,j,k,2) + pcorr(i-1,j-1,k,2) )        
     <		+ g13(i-1,j,k) * ( q(i,  j,k+1) - q(i,  j,k-1)
     <		                 + q(i-1,j,k+1) - q(i-1,j,k-1) )
     <          *0.25D0*(pcorr(i,j,k,3) + pcorr(i,j,k-1,3)
     <                  +pcorr(i-1,j,k,3) + pcorr(i-1,j,k-1,3) ))
	   enddo
	   enddo
        endif

        if ( n_suth .eq. MPI_PROC_NULL ) then
           j = 0
           do k = 1, kk
           do i = 1, ii
	      r(i,j,k) = - b(i,j,k) +
     <		( g23(i,j,  k) * ( q(i,j,  k+1) - q(i,j,  k-1)
     <		                 + q(i,j+1,k+1) - q(i,j+1,k-1) )
     <          *0.25D0*(pcorr(i,j,k,3) + pcorr(i,j,k-1,3)
     <                  +pcorr(i,j+1,k,3) + pcorr(i,j+1,k-1,3) )
     <		+ g21(i,j,  k) * ( q(i+1,j,  k) - q(i-1,j,  k)
     <		                 + q(i+1,j+1,k) - q(i-1,j+1,k) ) 
     <          *0.25D0*(pcorr(i,j,k,1) + pcorr(i-1,j,k,1)
     <                  +pcorr(i,j+1,k,1) + pcorr(i-1,j+1,k,1) ) )
	   enddo
	   enddo
        endif

        if ( n_nrth .eq. MPI_PROC_NULL ) then
           j = jj + 1
           do k = 1, kk
           do i = 1, ii
	      r(i,j,k) = - b(i,j,k) + 
     <		( g23(i,j-1,k) * ( q(i,j,  k+1) - q(i,j,  k-1)
     <		                 + q(i,j-1,k+1) - q(i,j-1,k-1) )
     <          *0.25D0*(pcorr(i,j,k,3) + pcorr(i,j,k-1,3)
     <                  +pcorr(i,j-1,k,3) + pcorr(i,j-1,k-1,3) ) 
     <		+ g21(i,j-1,k) * ( q(i+1,j,  k) - q(i-1,j,  k)
     <		                 + q(i+1,j-1,k) - q(i-1,j-1,k) )
     <          *0.25D0*(pcorr(i,j,k,1) + pcorr(i-1,j,k,1)
     <                  +pcorr(i,j-1,k,1) + pcorr(i-1,j-1,k,1) ) )
	   enddo
	   enddo
        endif

        if ( n_back .eq. MPI_PROC_NULL ) then
           k = 0
           do j = 1, jj
           do i = 1, ii
	      r(i,j,k) = - b(i,j,k) +
     <		( g31(i,j,k  ) * ( q(i+1,j,k  ) - q(i-1,j,k  )
     <		                 + q(i+1,j,k+1) - q(i-1,j,k+1) )
     <          *0.25D0*(pcorr(i,j,k,1) + pcorr(i-1,j,k,1)
     <                  +pcorr(i,j,k+1,1) + pcorr(i-1,j,k+1,1) ) 
     <		+ g32(i,j,k  ) * ( q(i,j+1,k  ) - q(i,j-1,k  )
     <		                 + q(i,j+1,k+1) - q(i,j-1,k+1) ) 
     <          *0.25D0*(pcorr(i,j,k,2) + pcorr(i,j-1,k,2)
     <                  +pcorr(i,j,k+1,2) + pcorr(i,j-1,k+1,2) ) )
	   enddo
	   enddo
        endif

        if ( n_frnt .eq. MPI_PROC_NULL ) then
           k = kk + 1
           do j = 1, jj
           do i = 1, ii
	      r(i,j,k) = - b(i,j,k) +
     <		( g31(i,j,k-1) * ( q(i+1,j,k  ) - q(i-1,j,k  )
     <		                 + q(i+1,j,k-1) - q(i-1,j,k-1) )
     <          *0.25D0*(pcorr(i,j,k,1) + pcorr(i-1,j,k,1)
     <                  +pcorr(i,j,k-1,1) + pcorr(i-1,j,k-1,1) ) 
     <		+ g32(i,j,k-1) * ( q(i,j+1,k  ) - q(i,j-1,k  )
     <		                 + q(i,j+1,k-1) - q(i,j-1,k-1) ) 
     <          *0.25D0*(pcorr(i,j,k,2) + pcorr(i,j-1,k,2)
     <                  +pcorr(i,j,k-1,2) + pcorr(i,j-1,k-1,2) ) )
	   enddo
	   enddo
        endif
        
	do k = 1, kk

	do j = 1, jj
	   do i = 1, ii
	      aji(j,i) = - t11(i-1,j,k)*pcorr(i-1,j,k,1) 
     <             + t21(i,j  ,k)*0.5*(pcorr(i,j+1,k,1)+ pcorr(i,j,k,1)) 
     <             - t21(i,j-1,k)*0.5*(pcorr(i,j-1,k,1)+ pcorr(i,j,k,1)) 
     <	           + t31(i,j,k  )*0.5*(pcorr(i,j,k,1)+ pcorr(i,j,k+1,1)) 
     <             - t31(i,j,k-1)*0.5*(pcorr(i,j,k,1)+ pcorr(i,j,k-1,1)) 



	      cji(j,i) = - t11(i,  j,k)*pcorr(i ,j,k,1)
     <             - t21(i  ,j,k)*0.5*(pcorr(i,j+1,k,1)+ pcorr(i,j,k,1))
     <             + t21(i,j-1,k)*0.5*(pcorr(i,j-1,k,1)+ pcorr(i,j,k,1)) 
     <	           - t31(i,j,k  )*0.5*(pcorr(i,j,k,1)+ pcorr(i,j,k+1,1)) 
     <             + t31(i,j,k-1)*0.5*(pcorr(i,j,k,1)+ pcorr(i,j,k-1,1))
 
	      bji(j,i) =   tcc(i,j,  k)
	      fji(j,i) =     r(i,j,  k)
	   enddo
	enddo

	if ( n_west .eq. MPI_PROC_NULL ) then
	   i = 0
	   do j = 1, jj
	      aji(j,i) =   0.D0
	      bji(j,i) =   g11(i,j,k)*pcorr(i,j,k,1)
	      cji(j,i) = - g11(i,j,k)*pcorr(i,j,k,1)
	      fji(j,i) =     r(i,j,k) 
	   enddo
	endif
	
	if ( n_east .eq. MPI_PROC_NULL ) then
           i = ii + 1
	   do j = 1, jj
	      aji(j,i) =   g11(i-1,j,k)*pcorr(i-1,j,k,1)
	      bji(j,i) = - g11(i-1,j,k)*pcorr(i-1,j,k,1)
	      cji(j,i) = 0.D0
	      fji(j,i) =     r(i,  j,k)
	   enddo
	endif
	if ( periodic .eq. 1 ) then
	call trip( aji, bji, cji, fji, jj, 1, ii, n_west, n_east )
	else
	call trid( aji, bji, cji, fji, jj, 1, ii, n_west, n_east )
	endif

c	call trid( aji, bji, cji, fji, jj, 1, ii, n_west, n_east )

	do j =  1, jj
	   do i = 0, ii+1
	      r(i,j,k) = fji(j,i)
	   enddo
	enddo

	enddo

        if ( n_suth .eq. MPI_PROC_NULL ) then
           j = 0
           do k = 1, kk
           do i = 1, ii
	      r(i,j,k) = r(i,j+1,k) + r(i,j,k)
     <                              /(g22(i,j,k)*pcorr(i,j,k,2))
	   enddo
	   enddo
        endif

        if ( n_nrth .eq. MPI_PROC_NULL ) then
           j = jj + 1
           do k = 1, kk
           do i = 1, ii
	      r(i,j,k) = r(i,j-1,k) - r(i,j,k)
     <                              /(g22(i,j-1,k)*pcorr(i,j-1,k,2))
	   enddo
	   enddo
        endif

        if ( n_back .eq. MPI_PROC_NULL ) then
           k = 0
           do j = 1, jj
           do i = 1, ii
      	      r(i,j,k) = r(i,j,k+1) + r(i,j,k)
     <                              /(g33(i,j,k)*pcorr(i,j,k,3) )
	   enddo
	   enddo
        endif

        if ( n_frnt .eq. MPI_PROC_NULL ) then
           k = kk + 1
           do j = 1, jj
           do i = 1, ii
      	      r(i,j,k) = r(i,j,k-1) - r(i,j,k)
     <                              /(g33(i,j,k-1)*pcorr(i,j,k-1,3))
	   enddo
	   enddo
        endif
        
	if ( n_west .eq. MPI_PROC_NULL ) then
	   if ( n_suth .eq. MPI_PROC_NULL ) then
	      do k = 1, kk
C******	0 - 0 - *
cc	r(0,0,k) = r(1,0,k) + r(0,1,k) - r(1,1,k)
 	r(0,0,k) = 0.5D0 * ( r(1,0,k) + r(0,1,k) )
	      enddo
	   endif
	   if ( n_nrth .eq. MPI_PROC_NULL ) then
	      do k = 1, kk
C******	0 - Y - *
cc	r(0,jj+1,k) = r(1,jj+1,k) + r(0,jj,k) - r(1,jj,k)
 	r(0,jj+1,k) = 0.5D0 * ( r(1,jj+1,k) + r(0,jj,k) )
	      enddo
	   endif
	   if ( n_back .eq. MPI_PROC_NULL ) then
	      do j = 1, jj
C******	0 - * - 0
cc	r(0,j,0) = r(1,j,0) + r(0,j,1) - r(1,j,1)
 	r(0,j,0) = 0.5D0 * ( r(1,j,0) + r(0,j,1) )
	      enddo
	   endif
	   if ( n_frnt .eq. MPI_PROC_NULL ) then
	      do j = 1, jj
C******	0 - * - Z
cc	r(0,j,kk+1) = r(1,j,kk+1) + r(0,j,kk) - r(1,j,kk)
 	r(0,j,kk+1) = 0.5D0 * ( r(1,j,kk+1) + r(0,j,kk) )
	      enddo
	   endif
	endif	      

	if ( n_east .eq. MPI_PROC_NULL ) then
	   if ( n_suth .eq. MPI_PROC_NULL ) then
	      do k = 1, kk
C******	X - 0 - *
cc	r(ii+1,0,k) = r(ii,0,k) + r(ii+1,1,k) - r(ii,1,k)
 	r(ii+1,0,k) = 0.5D0 * ( r(ii,0,k) + r(ii+1,1,k) )
	      enddo
	   endif
	   if ( n_nrth .eq. MPI_PROC_NULL ) then
	      do k = 1, kk
C******	X - Y - *
cc	r(ii+1,jj+1,k) = r(ii,jj+1,k) + r(ii+1,jj,k) - r(ii,jj,k)
 	r(ii+1,jj+1,k) = 0.5D0 * ( r(ii,jj+1,k) + r(ii+1,jj,k) )
	      enddo
	   endif
	   if ( n_back .eq. MPI_PROC_NULL ) then
	      do j = 1, jj
C******	X - * - 0
cc	r(ii+1,j,0) = r(ii,j,0) + r(ii+1,j,1) - r(ii,j,1)
 	r(ii+1,j,0) = 0.5D0 * ( r(ii,j,0) + r(ii+1,j,1) )
	      enddo
	   endif
	   if ( n_frnt .eq. MPI_PROC_NULL ) then
	      do j = 1, jj
C******	X - * - Z
cc	r(ii+1,j,kk+1) = r(ii,j,kk+1) + r(ii+1,j,kk) - r(ii,j,kk)
 	r(ii+1,j,kk+1) = 0.5D0 * ( r(ii,j,kk+1) + r(ii+1,j,kk) )
	      enddo
	   endif
	endif	
	
	if ( n_suth .eq. MPI_PROC_NULL ) then
	   if ( n_back .eq. MPI_PROC_NULL ) then
	      do i = 1, ii
C******	* - 0 - 0
cc	r(i,0,0) = r(i,1,0) + r(i,0,1) - r(i,1,1)
 	r(i,0,0) = 0.5D0 * ( r(i,1,0) + r(i,0,1) )
	      enddo
	   endif
	   if ( n_frnt .eq. MPI_PROC_NULL ) then
	      do i = 1, ii
C******	* - 0 - Z
cc	r(i,0,kk+1) = r(i,1,kk+1) + r(i,0,kk) - r(i,1,kk)
 	r(i,0,kk+1) = 0.5D0 * ( r(i,1,kk+1) + r(i,0,kk) )
	      enddo
	   endif
	endif      

	if ( n_nrth .eq. MPI_PROC_NULL ) then
	   if ( n_back .eq. MPI_PROC_NULL ) then
	      do i = 1, ii
C******	* - Y - 0
cc	r(i,jj+1,0) = r(i,jj,0) + r(i,jj+1,1) - r(i,jj,1)
 	r(i,jj+1,0) = 0.5D0 * ( r(i,jj,0) + r(i,jj+1,1) )
	      enddo
	   endif
	   if ( n_frnt .eq. MPI_PROC_NULL ) then
	      do i = 1, ii
C******	* - Y - Z
cc	r(i,jj+1,kk+1) = r(i,jj,kk+1) + r(i,jj+1,kk) - r(i,jj,kk)
 	r(i,jj+1,kk+1) = 0.5D0 * ( r(i,jj,kk+1) + r(i,jj+1,kk) )
	      enddo
	   endif
	endif  
	
	call p_exchange( ii, jj, kk, r ) 

	do k = 0, kk+1
	do j = 0, jj+1
	do i = 0, ii+1
	   q(i,j,k) = r(i,j,k)
	enddo
	enddo
	enddo

	return
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine smoothj( level, ii, jj, kk, q, r, b,
     <		g11, g12, g13, g21, g22, g23, g31, g32, g33, gcc,
     <		t11, t12, t13, t21, t22, t23, t31, t32, t33, tcc )

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "para.inc"
	include "sedi.inc"

	integer level, ii, jj, kk
	double precision q(0:ii+1,0:jj+1,0:kk+1), 
     &                   r(0:ii+1,0:jj+1,0:kk+1), 
     &                   b(0:ii+1,0:jj+1,0:kk+1)
	
	double precision g11(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   g12(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   g13(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   g21(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   g22(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   g23(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   g31(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   g32(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   g33(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   gcc(-1:ii+2,-1:jj+2,-1:kk+2),
     &                   t11(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   t12(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   t13(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   t21(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   t22(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   t23(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   t31(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   t32(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   t33(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   tcc(-1:ii+2,-1:jj+2,-1:kk+2)

	double precision aij(ii,0:jj+1), bij(ii,0:jj+1), 
     &                   cij(ii,0:jj+1), fij(ii,0:jj+1)

	integer i, j, k

	do k = 1, kk
	do j = 1, jj
	do i = 1, ii

	   r(i,j,k) = - b(i,j,k) + 
     <		( t11(i,  j,  k  ) * q(i+1,j,  k  )*pcorr(i,j,k,1)
     <		+ t11(i-1,j,  k  ) * q(i-1,j,  k  )*pcorr(i-1,j,k,1)
     <		+ t33(i,  j,  k  ) * q(i,  j,  k+1)*pcorr(i,j,k,3)
     <		+ t33(i,  j,  k-1) * q(i,  j,  k-1)*pcorr(i,j,k-1,3) )

	   r(i,j,k) = r(i,j,k) +
     <		( t12(i,  j,k) * ( q(i+1,j+1,k) - q(i+1,j-1,k) )
     <          *0.5D0*(pcorr(i+1,j,k,2)+pcorr(i+1,j-1,k,2))
     <		- t12(i-1,j,k) * ( q(i-1,j+1,k) - q(i-1,j-1,k) )
     <          *0.5D0*(pcorr(i-1,j,k,2)+pcorr(i-1,j-1,k,2))
     <		+ t13(i,  j,k) * ( q(i,  j,k+1) - q(i,  j,k-1)
     <		                 + q(i+1,j,k+1) - q(i+1,j,k-1) )
     <          *0.25D0*(pcorr(i,j,k,3)+pcorr(i,j,k-1,3)
     <                  +pcorr(i+1,j,k,3)+pcorr(i+1,j,k-1,3))
     <		- t13(i-1,j,k) * ( q(i,  j,k+1) - q(i,  j,k-1)
     <		                 + q(i-1,j,k+1) - q(i-1,j,k-1) ) 
     <          *0.25D0*(pcorr(i,j,k,3)+pcorr(i,j,k-1,3)
     <                  +pcorr(i-1,j,k,3)+pcorr(i-1,j,k-1,3))   )

	   r(i,j,k) = r(i,j,k) +
     <		( t23(i,j,  k) * ( q(i,j,  k+1) - q(i,j,  k-1)
     <		                 + q(i,j+1,k+1) - q(i,j+1,k-1) )
     <          *0.25D0*(pcorr(i,j,k,3)+pcorr(i,j,k-1,3)
     <                  +pcorr(i,j+1,k,3)+pcorr(i,j+1,k-1,3))
     <		- t23(i,j-1,k) * ( q(i,j,  k+1) - q(i,j,  k-1)
     <		                 + q(i,j-1,k+1) - q(i,j-1,k-1) )
     <          *0.25D0*(pcorr(i,j,k,3)+pcorr(i,j,k-1,3)
     <                  +pcorr(i,j-1,k,3)+pcorr(i,j-1,k-1,3)) 
     <		+ t21(i,j,  k) * ( q(i+1,j,  k) - q(i-1,j,  k)
     <		                 + q(i+1,j+1,k) - q(i-1,j+1,k) )
     <          *0.25D0*(pcorr(i,j,k,1)+pcorr(i-1,j,k,1)
     <                  +pcorr(i,j+1,k,1)+pcorr(i-1,j+1,k,1)) 
     <		- t21(i,j-1,k) * ( q(i+1,j,  k) - q(i-1,j,  k)
     <		                 + q(i+1,j-1,k) - q(i-1,j-1,k) )
     <          *0.25D0*(pcorr(i,j,k,1)+pcorr(i-1,j,k,1)
     <                  +pcorr(i,j-1,k,1)+pcorr(i-1,j-1,k,1)) )

	   r(i,j,k) = r(i,j,k) +
     <		( t31(i,j,k  ) * ( q(i+1,j,k  ) - q(i-1,j,k  )
     <		                 + q(i+1,j,k+1) - q(i-1,j,k+1) )
     <          *0.25D0*(pcorr(i,j,k,1)+pcorr(i-1,j,k,1)
     <                  +pcorr(i,j,k+1,1)+pcorr(i-1,j,k+1,1))
     <		- t31(i,j,k-1) * ( q(i+1,j,k  ) - q(i-1,j,k  )
     <		                 + q(i+1,j,k-1) - q(i-1,j,k-1) )
     <          *0.25D0*(pcorr(i,j,k,1)+pcorr(i-1,j,k,1)
     <                  +pcorr(i,j,k-1,1)+pcorr(i-1,j,k-1,1))
     <		+ t32(i,j,k  ) * ( q(i,j+1,k+1) - q(i,j-1,k+1) )
     <          *0.5D0*(pcorr(i,j,k+1,2)+pcorr(i,j-1,k+1,2))
     <		- t32(i,j,k-1) * ( q(i,j+1,k-1) - q(i,j-1,k-1) )
     <          *0.5D0*(pcorr(i,j,k-1,2)+pcorr(i,j-1,k-1,2)))
	enddo
	enddo
	enddo
	
	if ( n_west .eq. MPI_PROC_NULL ) then
	   i = 0
	   do k = 1, kk
	   do j = 1, jj
	      r(i,j,k) = - b(i,j,k) + 
     <		( g12(i,  j,k) * ( q(i,  j+1,k) - q(i,  j-1,k) 
     <		                 + q(i+1,j+1,k) - q(i+1,j-1,k) )
     <            *0.25D0*(pcorr(i,j,k,2)+pcorr(i,j-1,k,2)
     <                    +pcorr(i+1,j,k,2)+pcorr(i+1,j-1,k,2))
     <		+ g13(i,  j,k) * ( q(i,  j,k+1) - q(i,  j,k-1)
     <		                 + q(i+1,j,k+1) - q(i+1,j,k-1) )
     <            *0.25D0*(pcorr(i,j,k,3)+pcorr(i,j,k-1,3)
     <                    +pcorr(i+1,j,k,3)+pcorr(i+1,j,k-1,3)) )

	   enddo
	   enddo
	endif

        if ( n_east .eq. MPI_PROC_NULL ) then
	   i = ii + 1
	   do k = 1, kk
	   do j = 1, jj
	      r(i,j,k) = - b(i,j,k) +
     <		( g12(i-1,j,k) * ( q(i,  j+1,k) - q(i,  j-1,k)
     <		                 + q(i-1,j+1,k) - q(i-1,j-1,k) )
     <            *0.25D0*(pcorr(i,j,k,2)+pcorr(i,j-1,k,2)
     <                    +pcorr(i-1,j,k,2)+pcorr(i-1,j-1,k,2))
     <		+ g13(i-1,j,k) * ( q(i,  j,k+1) - q(i,  j,k-1)
     <		                 + q(i-1,j,k+1) - q(i-1,j,k-1) )
     <            *0.25D0*(pcorr(i,j,k,3)+pcorr(i,j,k-1,3)
     <                    +pcorr(i-1,j,k,3)+pcorr(i-1,j,k-1,3)) )
	   enddo
	   enddo
        endif

        if ( n_suth .eq. MPI_PROC_NULL ) then
           j = 0
           do k = 1, kk
           do i = 1, ii
	      r(i,j,k) = - b(i,j,k) +
     <		( g23(i,j,  k) * ( q(i,j,  k+1) - q(i,j,  k-1)
     <		                 + q(i,j+1,k+1) - q(i,j+1,k-1) )
     <            *0.25D0*(pcorr(i,j,k,3)+pcorr(i,j,k-1,3)
     <                    +pcorr(i,j+1,k,3)+pcorr(i,j+1,k-1,3))
     <		+ g21(i,j,  k) * ( q(i+1,j,  k) - q(i-1,j,  k)
     <		                 + q(i+1,j+1,k) - q(i-1,j+1,k) )
     <            *0.25D0*(pcorr(i,j,k,1)+pcorr(i-1,j,k,1)
     <                    +pcorr(i,j+1,k,1)+pcorr(i-1,j+1,k,1)))
	   enddo
	   enddo
        endif

        if ( n_nrth .eq. MPI_PROC_NULL ) then
           j = jj + 1
           do k = 1, kk
           do i = 1, ii
	      r(i,j,k) = - b(i,j,k) + 
     <		( g23(i,j-1,k) * ( q(i,j,  k+1) - q(i,j,  k-1)
     <		                 + q(i,j-1,k+1) - q(i,j-1,k-1) )
     <            *0.25D0*(pcorr(i,j,k,3)+pcorr(i,j,k-1,3)
     <                    +pcorr(i,j-1,k,3)+pcorr(i,j-1,k-1,3)) 
     <		+ g21(i,j-1,k) * ( q(i+1,j,  k) - q(i-1,j,  k)
     <		                 + q(i+1,j-1,k) - q(i-1,j-1,k) ) 
     <            *0.25D0*(pcorr(i,j,k,1)+pcorr(i-1,j,k,1)
     <                    +pcorr(i,j-1,k,1)+pcorr(i-1,j-1,k,1)))
	   enddo
	   enddo
        endif

        if ( n_back .eq. MPI_PROC_NULL ) then
           k = 0
           do j = 1, jj
           do i = 1, ii
	      r(i,j,k) = - b(i,j,k) +
     <		( g31(i,j,k  ) * ( q(i+1,j,k  ) - q(i-1,j,k  )
     <		                 + q(i+1,j,k+1) - q(i-1,j,k+1) )
     <            *0.25D0*(pcorr(i,j,k,1)+pcorr(i-1,j,k,1)
     <                    +pcorr(i,j,k+1,1)+pcorr(i-1,j,k+1,1)) 
     <		+ g32(i,j,k  ) * ( q(i,j+1,k  ) - q(i,j-1,k  )
     <		                 + q(i,j+1,k+1) - q(i,j-1,k+1) )
     <            *0.25D0*(pcorr(i,j,k,2)+pcorr(i,j-1,k,2)
     <                    +pcorr(i,j,k+1,2)+pcorr(i,j-1,k+1,2)))
	   enddo
	   enddo
        endif

        if ( n_frnt .eq. MPI_PROC_NULL ) then
           k = kk + 1
           do j = 1, jj
           do i = 1, ii
	      r(i,j,k) = - b(i,j,k) +
     <		( g31(i,j,k-1) * ( q(i+1,j,k  ) - q(i-1,j,k  )
     <		                 + q(i+1,j,k-1) - q(i-1,j,k-1) )
     <            *0.25D0*(pcorr(i,j,k,1)+pcorr(i-1,j,k,1)
     <                    +pcorr(i,j,k-1,1)+pcorr(i-1,j,k-1,1)) 
     <		+ g32(i,j,k-1) * ( q(i,j+1,k  ) - q(i,j-1,k  )
     <		                 + q(i,j+1,k-1) - q(i,j-1,k-1) )
     <            *0.25D0*(pcorr(i,j,k,2)+pcorr(i,j-1,k,2)
     <                    +pcorr(i,j,k-1,2)+pcorr(i,j-1,k-1,2)) )
	   enddo
	   enddo
        endif
        
	do k = 1, kk

	do j = 1, jj
	   do i = 1, ii

	      aij(i,j) = - t22(i,j-1,k)*pcorr(i,j-1,k,2)
     <             + t12(i  ,j,k)*0.5*(pcorr(i+1,j,k,2)+ pcorr(i,j,k,2))
     <             - t12(i-1,j,k)*0.5*(pcorr(i-1,j,k,2)+ pcorr(i,j,k,2))
     <	           + t32(i,j,k  )*0.5*(pcorr(i,j,k+1,2)+ pcorr(i,j,k,2))
     <             - t32(i,j,k-1)*0.5*(pcorr(i,j,k-1,2)+ pcorr(i,j,k,2))

	      cij(i,j) = - t22(i,j,  k)*pcorr(i,j,k,2)
     <             - t12(i  ,j,k)*0.5*(pcorr(i+1,j,k,2)+ pcorr(i,j,k,2))
     <             + t12(i-1,j,k)*0.5*(pcorr(i-1,j,k,2)+ pcorr(i,j,k,2))
     <	           - t32(i,j,k  )*0.5*(pcorr(i,j,k+1,2)+ pcorr(i,j,k,2))
     <             + t32(i,j,k-1)*0.5*(pcorr(i,j,k-1,2)+ pcorr(i,j,k,2))

	      bij(i,j) =   tcc(i,j,  k)
	      fij(i,j) =     r(i,j,  k)
	   enddo
	enddo

	if ( n_suth .eq. MPI_PROC_NULL ) then
	   j = 0
	   do i = 1, ii
	      aij(i,j) =   0.D0
	      bij(i,j) =   g22(i,j,k)*pcorr(i,j,k,2)
	      cij(i,j) = - g22(i,j,k)*pcorr(i,j,k,2)
	      fij(i,j) =   r(i,j,k)
	      
	   enddo
	endif
	
	if ( n_nrth .eq. MPI_PROC_NULL ) then
           j = jj + 1
	   do i = 1, ii
	      aij(i,j) =   g22(i,j-1,k)*pcorr(i,j-1,k,2)
	      bij(i,j) = - g22(i,j-1,k)*pcorr(i,j-1,k,2)
	      cij(i,j) = 0.D0
	      fij(i,j) =     r(i,j  ,k)
	   enddo
	endif
	
	call trid( aij, bij, cij, fij, ii, 1, jj, n_suth, n_nrth )

	do j =  0, jj+1
	   do i = 1, ii
	      r(i,j,k) = fij(i,j)
	   enddo
	enddo

	enddo

	if ( n_west .eq. MPI_PROC_NULL ) then
	   i = 0
	   do k = 1, kk
	   do j = 1, jj
     	      r(i,j,k) = r(i+1,j,k) + r(i,j,k)
     <                 /(g11(i,j,k)*pcorr(i,j,k,1))
	   enddo
	   enddo
	endif

        if ( n_east .eq. MPI_PROC_NULL ) then
	   i = ii + 1
	   do k = 1, kk
	   do j = 1, jj
     	      r(i,j,k) = r(i-1,j,k) - r(i,j,k)
     <                 /(g11(i-1,j,k)*pcorr(i-1,j,k,1))
	   enddo
	   enddo
        endif

        if ( n_back .eq. MPI_PROC_NULL ) then
           k = 0
           do j = 1, jj
           do i = 1, ii
      	      r(i,j,k) = r(i,j,k+1) + r(i,j,k)
     <                 /(g33(i,j,k)*pcorr(i,j,k,3))
	   enddo
	   enddo
        endif

        if ( n_frnt .eq. MPI_PROC_NULL ) then
           k = kk + 1
           do j = 1, jj
           do i = 1, ii
      	      r(i,j,k) = r(i,j,k-1) - r(i,j,k)
     <                 /(g33(i,j,k-1)*pcorr(i,j,k-1,3))
	   enddo
	   enddo
        endif
        
	if ( n_west .eq. MPI_PROC_NULL ) then
	   if ( n_suth .eq. MPI_PROC_NULL ) then
	      do k = 1, kk
C******	0 - 0 - *
cc	r(0,0,k) = r(1,0,k) + r(0,1,k) - r(1,1,k)
 	r(0,0,k) = 0.5D0 * ( r(1,0,k) + r(0,1,k) )
	      enddo
	   endif
	   if ( n_nrth .eq. MPI_PROC_NULL ) then
	      do k = 1, kk
C******	0 - Y - *
cc	r(0,jj+1,k) = r(1,jj+1,k) + r(0,jj,k) - r(1,jj,k)
 	r(0,jj+1,k) = 0.5D0 * ( r(1,jj+1,k) + r(0,jj,k) )
	      enddo
	   endif
	   if ( n_back .eq. MPI_PROC_NULL ) then
	      do j = 1, jj
C******	0 - * - 0
cc	r(0,j,0) = r(1,j,0) + r(0,j,1) - r(1,j,1)
 	r(0,j,0) = 0.5D0 * ( r(1,j,0) + r(0,j,1) )
	      enddo
	   endif
	   if ( n_frnt .eq. MPI_PROC_NULL ) then
	      do j = 1, jj
C******	0 - * - Z
cc	r(0,j,kk+1) = r(1,j,kk+1) + r(0,j,kk) - r(1,j,kk)
 	r(0,j,kk+1) = 0.5D0 * ( r(1,j,kk+1) + r(0,j,kk) )
	      enddo
	   endif
	endif	      

	if ( n_east .eq. MPI_PROC_NULL ) then
	   if ( n_suth .eq. MPI_PROC_NULL ) then
	      do k = 1, kk
C******	X - 0 - *
cc	r(ii+1,0,k) = r(ii,0,k) + r(ii+1,1,k) - r(ii,1,k)
 	r(ii+1,0,k) = 0.5D0 * ( r(ii,0,k) + r(ii+1,1,k) )
	      enddo
	   endif
	   if ( n_nrth .eq. MPI_PROC_NULL ) then
	      do k = 1, kk
C******	X - Y - *
cc	r(ii+1,jj+1,k) = r(ii,jj+1,k) + r(ii+1,jj,k) - r(ii,jj,k)
 	r(ii+1,jj+1,k) = 0.5D0 * ( r(ii,jj+1,k) + r(ii+1,jj,k) )
	      enddo
	   endif
	   if ( n_back .eq. MPI_PROC_NULL ) then
	      do j = 1, jj
C******	X - * - 0
cc	r(ii+1,j,0) = r(ii,j,0) + r(ii+1,j,1) - r(ii,j,1)
 	r(ii+1,j,0) = 0.5D0 * ( r(ii,j,0) + r(ii+1,j,1) )
	      enddo
	   endif
	   if ( n_frnt .eq. MPI_PROC_NULL ) then
	      do j = 1, jj
C******	X - * - Z
cc	r(ii+1,j,kk+1) = r(ii,j,kk+1) + r(ii+1,j,kk) - r(ii,j,kk)
 	r(ii+1,j,kk+1) = 0.5D0 * ( r(ii,j,kk+1) + r(ii+1,j,kk) )
	      enddo
	   endif
	endif	
	
	if ( n_suth .eq. MPI_PROC_NULL ) then
	   if ( n_back .eq. MPI_PROC_NULL ) then
	      do i = 1, ii
C******	* - 0 - 0
cc	r(i,0,0) = r(i,1,0) + r(i,0,1) - r(i,1,1)
 	r(i,0,0) = 0.5D0 * ( r(i,1,0) + r(i,0,1) )
	      enddo
	   endif
	   if ( n_frnt .eq. MPI_PROC_NULL ) then
	      do i = 1, ii
C******	* - 0 - Z
cc	r(i,0,kk+1) = r(i,1,kk+1) + r(i,0,kk) - r(i,1,kk)
 	r(i,0,kk+1) = 0.5D0 * ( r(i,1,kk+1) + r(i,0,kk) )
	      enddo
	   endif
	endif      

	if ( n_nrth .eq. MPI_PROC_NULL ) then
	   if ( n_back .eq. MPI_PROC_NULL ) then
	      do i = 1, ii
C******	* - Y - 0
cc	r(i,jj+1,0) = r(i,jj,0) + r(i,jj+1,1) - r(i,jj,1)
 	r(i,jj+1,0) = 0.5D0 * ( r(i,jj,0) + r(i,jj+1,1) )
	      enddo
	   endif
	   if ( n_frnt .eq. MPI_PROC_NULL ) then
	      do i = 1, ii
C******	* - Y - Z
cc	r(i,jj+1,kk+1) = r(i,jj,kk+1) + r(i,jj+1,kk) - r(i,jj,kk)
 	r(i,jj+1,kk+1) = 0.5D0 * ( r(i,jj,kk+1) + r(i,jj+1,kk) )
	      enddo
	   endif
	endif  
	
	call p_exchange( ii, jj, kk, r ) 

	do k = 0, kk+1
	do j = 0, jj+1
	do i = 0, ii+1
	   q(i,j,k) = r(i,j,k)
	enddo
	enddo
	enddo

	return
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine smoothk( level, ii, jj, kk, q, r, b,
     <		g11, g12, g13, g21, g22, g23, g31, g32, g33, gcc,
     <		t11, t12, t13, t21, t22, t23, t31, t32, t33, tcc )

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "para.inc"
	include "sedi.inc"

	integer level, ii, jj, kk
	double precision q(0:ii+1,0:jj+1,0:kk+1), 
     &                   r(0:ii+1,0:jj+1,0:kk+1), 
     &                   b(0:ii+1,0:jj+1,0:kk+1)
	
	double precision g11(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   g12(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   g13(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   g21(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   g22(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   g23(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   g31(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   g32(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   g33(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   gcc(-1:ii+2,-1:jj+2,-1:kk+2),
     &                   t11(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   t12(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   t13(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   t21(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   t22(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   t23(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   t31(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   t32(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   t33(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   tcc(-1:ii+2,-1:jj+2,-1:kk+2)

	double precision aik(ii,0:kk+1), bik(ii,0:kk+1), 
     &                   cik(ii,0:kk+1), fik(ii,0:kk+1)

	integer i, j, k

	do k = 1, kk
	do j = 1, jj
	do i = 1, ii

	   r(i,j,k) = - b(i,j,k) + 
     <		( t11(i,  j,  k  ) * q(i+1,j,  k  )*pcorr(i,j,k,1)
     <		+ t11(i-1,j,  k  ) * q(i-1,j,  k  )*pcorr(i-1,j,k,1)
     <		+ t22(i,  j,  k  ) * q(i,  j+1,k  )*pcorr(i,j,k,2)
     <		+ t22(i,  j-1,k  ) * q(i,  j-1,k  )*pcorr(i,j-1,k,2) )
     
	   r(i,j,k) = r(i,j,k) +
     <		( t12(i,  j,k) * ( q(i,  j+1,k) - q(i,  j-1,k) 
     <		                 + q(i+1,j+1,k) - q(i+1,j-1,k) )
     <            *0.25D0*(pcorr(i,j,k,2)+pcorr(i,j,k-1,2)
     <                    +pcorr(i+1,j,k-1,2)+pcorr(i+1,j,k-1,2)) 
     <		- t12(i-1,j,k) * ( q(i,  j+1,k) - q(i,  j-1,k)
     <		                 + q(i-1,j+1,k) - q(i-1,j-1,k) )
     <            *0.25D0*(pcorr(i,j,k,2)+pcorr(i,j,k-1,2)
     <                    +pcorr(i-1,j,k-1,2)+pcorr(i-1,j,k-1,2)) 
     <		+ t13(i,  j,k) * ( q(i+1,j,k+1) - q(i+1,j,k-1) )
     <            *0.5D0*(pcorr(i+1,j,k,3)+pcorr(i+1,j,k-1,3))
     <		- t13(i-1,j,k) * ( q(i-1,j,k+1) - q(i-1,j,k-1) )
     <            *0.5D0*(pcorr(i-1,j,k,3)+pcorr(i-1,j,k-1,3)))

	   r(i,j,k) = r(i,j,k) +
     <		( t23(i,j,  k) * ( q(i,j+1,k+1) - q(i,j+1,k-1) )
     <            *0.5D0*(pcorr(i,j+1,k,3)+pcorr(i,j+1,k-1,3))
     <		- t23(i,j-1,k) * ( q(i,j-1,k+1) - q(i,j-1,k-1) )
     <            *0.5D0*(pcorr(i,j-1,k,3)+pcorr(i,j-1,k-1,3)) 
     <		+ t21(i,j,  k) * ( q(i+1,j,  k) - q(i-1,j,  k)
     <		                 + q(i+1,j+1,k) - q(i-1,j+1,k) )
     <            *0.25D0*(pcorr(i,j,k,1)+pcorr(i-1,j,k,1)
     <                    +pcorr(i,j+1,k,1)+pcorr(i-1,j+1,k,1)) 
     <		- t21(i,j-1,k) * ( q(i+1,j,  k) - q(i-1,j,  k)
     <		                 + q(i+1,j-1,k) - q(i-1,j-1,k) )
     <            *0.25D0*(pcorr(i,j,k,1)+pcorr(i-1,j,k,1)
     <                    +pcorr(i,j-1,k,1)+pcorr(i-1,j-1,k,1)) )

	   r(i,j,k) = r(i,j,k) +
     <		( t31(i,j,k  ) * ( q(i+1,j,k  ) - q(i-1,j,k  )
     <		                 + q(i+1,j,k+1) - q(i-1,j,k+1) )
     <            *0.25D0*(pcorr(i,j,k,1)+pcorr(i-1,j,k,1)
     <                    +pcorr(i,j,k+1,1)+pcorr(i-1,j,k+1,1))
     <		- t31(i,j,k-1) * ( q(i+1,j,k  ) - q(i-1,j,k  )
     <		                 + q(i+1,j,k-1) - q(i-1,j,k-1) )
     <            *0.25D0*(pcorr(i,j,k,1)+pcorr(i-1,j,k,1)
     <                    +pcorr(i,j,k-1,1)+pcorr(i-1,j,k-1,1))
     <		+ t32(i,j,k  ) * ( q(i,j+1,k  ) - q(i,j-1,k  )
     <		                 + q(i,j+1,k+1) - q(i,j-1,k+1) )
     <            *0.25D0*(pcorr(i,j,k,2)+pcorr(i,j-1,k,2)
     <                    +pcorr(i,j,k+1,2)+pcorr(i,j-1,k+1,2))
     <		- t32(i,j,k-1) * ( q(i,j+1,k  ) - q(i,j-1,k  )
     <		                 + q(i,j+1,k-1) - q(i,j-1,k-1) )
     <            *0.25D0*(pcorr(i,j,k,2)+pcorr(i,j-1,k,2)
     <                    +pcorr(i,j,k-1,2)+pcorr(i,j-1,k-1,2)))
 
	enddo
	enddo
	enddo
	
	if ( n_west .eq. MPI_PROC_NULL ) then
	   i = 0
	   do k = 1, kk
	   do j = 1, jj
	      r(i,j,k) = - b(i,j,k) + 
     <		( g12(i,  j,k) * ( q(i,  j+1,k) - q(i,  j-1,k) 
     <		                 + q(i+1,j+1,k) - q(i+1,j-1,k) )
     <            *0.25D0*(pcorr(i,j,k,2)+pcorr(i,j-1,k,2)
     <                    +pcorr(i+1,j,k,2)+pcorr(i+1,j-1,k,2))
     <		+ g13(i,  j,k) * ( q(i,  j,k+1) - q(i,  j,k-1)
     <		                 + q(i+1,j,k+1) - q(i+1,j,k-1) )
     <            *0.25D0*(pcorr(i,j,k,3)+pcorr(i,j,k-1,3)
     <                    +pcorr(i+1,j,k,3)+pcorr(i+1,j,k-1,3)))
	   enddo
	   enddo
	endif

        if ( n_east .eq. MPI_PROC_NULL ) then
	   i = ii + 1
	   do k = 1, kk
	   do j = 1, jj
	      r(i,j,k) = - b(i,j,k) +
     <		( g12(i-1,j,k) * ( q(i,  j+1,k) - q(i,  j-1,k)
     <		                 + q(i-1,j+1,k) - q(i-1,j-1,k) )
     <            *0.25D0*(pcorr(i,j,k,2)+pcorr(i,j-1,k,2)
     <                    +pcorr(i-1,j,k,2)+pcorr(i-1,j-1,k,2))
     <		+ g13(i-1,j,k) * ( q(i,  j,k+1) - q(i,  j,k-1)
     <		                 + q(i-1,j,k+1) - q(i-1,j,k-1) )
     <            *0.25D0*(pcorr(i,j,k,3)+pcorr(i,j,k-1,3)
     <                    +pcorr(i-1,j,k,3)+pcorr(i-1,j,k-1,3)))
	   enddo
	   enddo
        endif

        if ( n_suth .eq. MPI_PROC_NULL ) then
           j = 0
           do k = 1, kk
           do i = 1, ii
	      r(i,j,k) = - b(i,j,k) +
     <		( g23(i,j,  k) * ( q(i,j,  k+1) - q(i,j,  k-1)
     <		                 + q(i,j+1,k+1) - q(i,j+1,k-1) )
     <            *0.25D0*(pcorr(i,j,k,3)+pcorr(i,j,k-1,3)
     <                    +pcorr(i,j+1,k,3)+pcorr(i,j+1,k-1,3))
     <		+ g21(i,j,  k) * ( q(i+1,j,  k) - q(i-1,j,  k)
     <		                 + q(i+1,j+1,k) - q(i-1,j+1,k) )
     <            *0.25D0*(pcorr(i,j,k,1)+pcorr(i-1,j,k,1)
     <                    +pcorr(i,j+1,k,1)+pcorr(i-1,j+1,k,1)))
	   enddo
	   enddo
        endif

        if ( n_nrth .eq. MPI_PROC_NULL ) then
           j = jj + 1
           do k = 1, kk
           do i = 1, ii
	      r(i,j,k) = - b(i,j,k) + 
     <		( g23(i,j-1,k) * ( q(i,j,  k+1) - q(i,j,  k-1)
     <		                 + q(i,j-1,k+1) - q(i,j-1,k-1) )
     <            *0.25D0*(pcorr(i,j,k,3)+pcorr(i,j,k-1,3)
     <                    +pcorr(i,j-1,k,3)+pcorr(i,j-1,k-1,3)) 
     <		+ g21(i,j-1,k) * ( q(i+1,j,  k) - q(i-1,j,  k)
     <		                 + q(i+1,j-1,k) - q(i-1,j-1,k) )
     <            *0.25D0*(pcorr(i,j,k,1)+pcorr(i-1,j,k,1)
     <                    +pcorr(i,j-1,k,1)+pcorr(i-1,j-1,k,1)))
	   enddo
	   enddo
        endif

        if ( n_back .eq. MPI_PROC_NULL ) then
           k = 0
           do j = 1, jj
           do i = 1, ii
	      r(i,j,k) = - b(i,j,k) +
     <		( g31(i,j,k  ) * ( q(i+1,j,k  ) - q(i-1,j,k  )
     <		                 + q(i+1,j,k+1) - q(i-1,j,k+1) )
     <            *0.25D0*(pcorr(i,j,k,1)+pcorr(i-1,j,k,1)
     <                    +pcorr(i,j,k+1,1)+pcorr(i-1,j,k+1,1))
     <		+ g32(i,j,k  ) * ( q(i,j+1,k  ) - q(i,j-1,k  )
     <		                 + q(i,j+1,k+1) - q(i,j-1,k+1) )
     <            *0.25D0*(pcorr(i,j,k,2)+pcorr(i,j-1,k,2)
     <                    +pcorr(i,j,k+1,2)+pcorr(i,j-1,k+1,2)))
	   enddo
	   enddo
        endif

        if ( n_frnt .eq. MPI_PROC_NULL ) then
           k = kk + 1
           do j = 1, jj
           do i = 1, ii
	      r(i,j,k) = - b(i,j,k) +
     <		( g31(i,j,k-1) * ( q(i+1,j,k  ) - q(i-1,j,k  )
     <		                 + q(i+1,j,k-1) - q(i-1,j,k-1) )
     <            *0.25D0*(pcorr(i,j,k,1)+pcorr(i-1,j,k,1)
     <                    +pcorr(i,j,k-1,1)+pcorr(i-1,j,k-1,1))
     <		+ g32(i,j,k-1) * ( q(i,j+1,k  ) - q(i,j-1,k  )
     <		                 + q(i,j+1,k-1) - q(i,j-1,k-1) )
     <            *0.25D0*(pcorr(i,j,k,2)+pcorr(i,j-1,k,2)
     <                    +pcorr(i,j,k-1,2)+pcorr(i,j-1,k-1,2)))
	   enddo
	   enddo
        endif
        
	do j = 1, jj

	do k = 1, kk
	   do i = 1, ii

	      aik(i,k) = - t33(i,j,k-1)*pcorr(i,j,k-1,3)
     <             + t13(i  ,j,k)*0.5*(pcorr(i+1,j,k,3)+ pcorr(i,j,k,3))
     <             - t13(i-1,j,k)*0.5*(pcorr(i-1,j,k,3)+ pcorr(i,j,k,3))
     <	           + t23(i,  j,k)*0.5*(pcorr(i,j+1,k,3)+ pcorr(i,j,k,3))
     <             - t23(i,j-1,k)*0.5*(pcorr(i,j-1,k,3)+ pcorr(i,j,k,3))

	      cik(i,k) = - t33(i,j,k  )*pcorr(i,j,k,3)
     <             - t13(i  ,j,k)*0.5*(pcorr(i+1,j,k,3)+ pcorr(i,j,k,3))
     <             + t13(i-1,j,k)*0.5*(pcorr(i-1,j,k,3)+ pcorr(i,j,k,3))
     <	           - t23(i,  j,k)*0.5*(pcorr(i,j+1,k,3)+ pcorr(i,j,k,3))
     <             + t23(i,j-1,k)*0.5*(pcorr(i,j-1,k,3)+ pcorr(i,j,k,3))

	      bik(i,k) =   tcc(i,j,  k)
	      fik(i,k) =     r(i,j,  k)
	   enddo
	enddo

	if ( n_back .eq. MPI_PROC_NULL ) then
	   k = 0
	   do i = 1, ii
	      aik(i,k) =   0.D0
	      bik(i,k) =   g33(i,j,k)*pcorr(i,j,k,3)
	      cik(i,k) = - g33(i,j,k)*pcorr(i,j,k,3)
	      fik(i,k) =     r(i,j,k) 
	   enddo
	endif
	
	if ( n_frnt .eq. MPI_PROC_NULL ) then
           k = kk + 1
	   do i = 1, ii
	      aik(i,k) =   g33(i,j,k-1)*pcorr(i,j,k-1,3)
	      bik(i,k) = - g33(i,j,k-1)*pcorr(i,j,k-1,3)
	      cik(i,k) = 0.D0
	      fik(i,k) =     r(i,j,k  )
	   enddo
	endif
	
	if ( periodicZ .eq. 1 ) then
	call trip( aik, bik, cik, fik, ii, 1, kk, n_back, n_frnt )
	else
	call trid( aik, bik, cik, fik, ii, 1, kk, n_back, n_frnt )
	endif

	do k = 0, kk+1
	   do i = 1, ii
	      r(i,j,k) = fik(i,k)
	   enddo
	enddo

	enddo

	if ( n_west .eq. MPI_PROC_NULL ) then
	   i = 0
	   do k = 1, kk
	   do j = 1, jj
     	      r(i,j,k) = r(i+1,j,k) + r(i,j,k)
     <                 / (g11(i,j,k)*pcorr(i,j,k,1))
	   enddo
	   enddo
	endif

        if ( n_east .eq. MPI_PROC_NULL ) then
	   i = ii + 1
	   do k = 1, kk
	   do j = 1, jj
     	      r(i,j,k) = r(i-1,j,k) - r(i,j,k)
     <                 / (g11(i-1,j,k)*pcorr(i-1,j,k,1))
	   enddo
	   enddo
        endif

        if ( n_suth .eq. MPI_PROC_NULL ) then
           j = 0
           do k = 1, kk
           do i = 1, ii
	      r(i,j,k) = r(i,j+1,k) + r(i,j,k)
     <                 / (g22(i,j,k)*pcorr(i,j,k,2))
	   enddo
	   enddo
        endif

        if ( n_nrth .eq. MPI_PROC_NULL ) then
           j = jj + 1
           do k = 1, kk
           do i = 1, ii
	      r(i,j,k) = r(i,j-1,k) - r(i,j,k)
     <               / (g22(i,j-1,k)*pcorr(i,j-1,k,2))
	   enddo
	   enddo
        endif

	if ( n_west .eq. MPI_PROC_NULL ) then
	   if ( n_suth .eq. MPI_PROC_NULL ) then
	      do k = 1, kk
C******	0 - 0 - *
cc	r(0,0,k) = r(1,0,k) + r(0,1,k) - r(1,1,k)
 	r(0,0,k) = 0.5D0 * ( r(1,0,k) + r(0,1,k) )
	      enddo
	   endif
	   if ( n_nrth .eq. MPI_PROC_NULL ) then
	      do k = 1, kk
C******	0 - Y - *
cc	r(0,jj+1,k) = r(1,jj+1,k) + r(0,jj,k) - r(1,jj,k)
 	r(0,jj+1,k) = 0.5D0 * ( r(1,jj+1,k) + r(0,jj,k) )
	      enddo
	   endif
	   if ( n_back .eq. MPI_PROC_NULL ) then
	      do j = 1, jj
C******	0 - * - 0
cc	r(0,j,0) = r(1,j,0) + r(0,j,1) - r(1,j,1)
 	r(0,j,0) = 0.5D0 * ( r(1,j,0) + r(0,j,1) )
	      enddo
	   endif
	   if ( n_frnt .eq. MPI_PROC_NULL ) then
	      do j = 1, jj
C******	0 - * - Z
cc	r(0,j,kk+1) = r(1,j,kk+1) + r(0,j,kk) - r(1,j,kk)
 	r(0,j,kk+1) = 0.5D0 * ( r(1,j,kk+1) + r(0,j,kk) )
	      enddo
	   endif
	endif	      

	if ( n_east .eq. MPI_PROC_NULL ) then
	   if ( n_suth .eq. MPI_PROC_NULL ) then
	      do k = 1, kk
C******	X - 0 - *
cc	r(ii+1,0,k) = r(ii,0,k) + r(ii+1,1,k) - r(ii,1,k)
 	r(ii+1,0,k) = 0.5D0 * ( r(ii,0,k) + r(ii+1,1,k) )
	      enddo
	   endif
	   if ( n_nrth .eq. MPI_PROC_NULL ) then
	      do k = 1, kk
C******	X - Y - *
cc	r(ii+1,jj+1,k) = r(ii,jj+1,k) + r(ii+1,jj,k) - r(ii,jj,k)
 	r(ii+1,jj+1,k) = 0.5D0 * ( r(ii,jj+1,k) + r(ii+1,jj,k) )
	      enddo
	   endif
	   if ( n_back .eq. MPI_PROC_NULL ) then
	      do j = 1, jj
C******	X - * - 0
cc	r(ii+1,j,0) = r(ii,j,0) + r(ii+1,j,1) - r(ii,j,1)
 	r(ii+1,j,0) = 0.5D0 * ( r(ii,j,0) + r(ii+1,j,1) )
	      enddo
	   endif
	   if ( n_frnt .eq. MPI_PROC_NULL ) then
	      do j = 1, jj
C******	X - * - Z
cc	r(ii+1,j,kk+1) = r(ii,j,kk+1) + r(ii+1,j,kk) - r(ii,j,kk)
 	r(ii+1,j,kk+1) = 0.5D0 * ( r(ii,j,kk+1) + r(ii+1,j,kk) )
	      enddo
	   endif
	endif	
	
	if ( n_suth .eq. MPI_PROC_NULL ) then
	   if ( n_back .eq. MPI_PROC_NULL ) then
	      do i = 1, ii
C******	* - 0 - 0
cc	r(i,0,0) = r(i,1,0) + r(i,0,1) - r(i,1,1)
 	r(i,0,0) = 0.5D0 * ( r(i,1,0) + r(i,0,1) )
	      enddo
	   endif
	   if ( n_frnt .eq. MPI_PROC_NULL ) then
	      do i = 1, ii
C******	* - 0 - Z
cc	r(i,0,kk+1) = r(i,1,kk+1) + r(i,0,kk) - r(i,1,kk)
 	r(i,0,kk+1) = 0.5D0 * ( r(i,1,kk+1) + r(i,0,kk) )
	      enddo
	   endif
	endif      

	if ( n_nrth .eq. MPI_PROC_NULL ) then
	   if ( n_back .eq. MPI_PROC_NULL ) then
	      do i = 1, ii
C******	* - Y - 0
cc	r(i,jj+1,0) = r(i,jj,0) + r(i,jj+1,1) - r(i,jj,1)
 	r(i,jj+1,0) = 0.5D0 * ( r(i,jj,0) + r(i,jj+1,1) )
	      enddo
	   endif
	   if ( n_frnt .eq. MPI_PROC_NULL ) then
	      do i = 1, ii
C******	* - Y - Z
cc	r(i,jj+1,kk+1) = r(i,jj,kk+1) + r(i,jj+1,kk) - r(i,jj,kk)
 	r(i,jj+1,kk+1) = 0.5D0 * ( r(i,jj,kk+1) + r(i,jj+1,kk) )
	      enddo
	   endif
	endif  
	
	call p_exchange( ii, jj, kk, r ) 

	do k = 0, kk+1
	do j = 0, jj+1
	do i = 0, ii+1
	   q(i,j,k) = r(i,j,k)
	enddo
	enddo
	enddo

	return
	end
