cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine smooth( resid, bbsum, ermin, ermax, resbc, sumbc,
     <		level, ii, jj, kk, p, r, b, jac,
     <		g11, g12, g13, g21, g22, g23, g31, g32, g33, gcc )

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "para.inc"

	integer level, ii, jj, kk 
	double precision, dimension(0:ii+1,0:jj+1,0:kk+1) :: p, r, b
	
	double precision, dimension(-1:ii+2,-1:jj+2,-1:kk+2) :: jac,
     <		g11, g12, g13, g21, g22, g23, g31, g32, g33, gcc,
     <		t11, t12, t13, t21, t22, t23, t31, t32, t33, tcc

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
           tcc(i,j,k) = t11(i, j, k) + t11(i-1, j,   k  )
     <                + t22(i, j, k) + t22(i,   j-1, k  )
     <                + t33(i, j, k) + t33(i,   j,   k-1)
        enddo 
        enddo
        enddo

	resid = 0.D0
	res_old = 0.D0
	res_new = 0.D0
	
	call residual( resid, bbsum, ermin, ermax, resbc, sumbc,
     <		level, ii, jj, kk, p, r, b, jac,
     <		g11, g12, g13, g21, g22, g23, g31, g32, g33, gcc,
     <		t11, t12, t13, t21, t22, t23, t31, t32, t33, tcc )
	erabs = max(dabs(ermin), dabs(ermax))
	res_old = resid

     	IF ( MOD(ISTEP,NSAVE) .EQ. 0 .and. MYID .EQ. 0 ) 
     <	WRITE(*,1) LEVEL, RESID, BBSUM, ERMIN, ERMAX, RESBC, SUMBC
    1	FORMAT(' L', I1,           ' START WITH    ', 2E9.3, 
     <	       ' E', 2E9.3, ' BC', 2E9.3)

	do n = 1, maxiter(level)

	do m = 1, iterchk(level)

	call bc_prep( level, ii, jj, kk, 1, 1, 2, 2, p, r, b,
     <		g11, g12, g13, g21, g22, g23, g31, g32, g33, gcc )
	call smoothj( level, ii, jj, kk, 1, 1, p, r, b,
     <		g11, g12, g13, g21, g22, g23, g31, g32, g33, gcc,
     <		t11, t12, t13, t21, t22, t23, t31, t32, t33, tcc )
	call bc_updt( level, ii, jj, kk, 1, 1, 2, 2, p, r, b,
     <		g11, g12, g13, g21, g22, g23, g31, g32, g33, gcc )

	call bc_prep( level, ii, jj, kk, 1, 2, 2, 2, p, r, b,
     <		g11, g12, g13, g21, g22, g23, g31, g32, g33, gcc )
	call smoothj( level, ii, jj, kk, 1, 2, p, r, b,
     <		g11, g12, g13, g21, g22, g23, g31, g32, g33, gcc,
     <		t11, t12, t13, t21, t22, t23, t31, t32, t33, tcc )
	call bc_updt( level, ii, jj, kk, 1, 2, 2, 2, p, r, b,
     <		g11, g12, g13, g21, g22, g23, g31, g32, g33, gcc )

	call bc_prep( level, ii, jj, kk, 2, 1, 2, 2, p, r, b,
     <		g11, g12, g13, g21, g22, g23, g31, g32, g33, gcc )
	call smoothj( level, ii, jj, kk, 2, 1, p, r, b,
     <		g11, g12, g13, g21, g22, g23, g31, g32, g33, gcc,
     <		t11, t12, t13, t21, t22, t23, t31, t32, t33, tcc )
	call bc_updt( level, ii, jj, kk, 2, 1, 2, 2, p, r, b,
     <		g11, g12, g13, g21, g22, g23, g31, g32, g33, gcc )

	call bc_prep( level, ii, jj, kk, 2, 2, 2, 2, p, r, b,
     <		g11, g12, g13, g21, g22, g23, g31, g32, g33, gcc )
	call smoothj( level, ii, jj, kk, 2, 2, p, r, b,
     <		g11, g12, g13, g21, g22, g23, g31, g32, g33, gcc,
     <		t11, t12, t13, t21, t22, t23, t31, t32, t33, tcc )
	call bc_updt( level, ii, jj, kk, 2, 2, 2, 2, p, r, b,
     <		g11, g12, g13, g21, g22, g23, g31, g32, g33, gcc )

	enddo

	call residual( resid, bbsum, ermin, ermax, resbc, sumbc,
     <		level, ii, jj, kk, p, r, b, jac,
     <		g11, g12, g13, g21, g22, g23, g31, g32, g33, gcc,
     <		t11, t12, t13, t21, t22, t23, t31, t32, t33, tcc )
	erabs = max(dabs(ermin), dabs(ermax))
	res_new = resid

C     	if ( res_new .lt. tol(level) .and. erabs .lt. ter(level) ) then
C     	IF ( MOD(ISTEP,NSAVE) .EQ. 0 .and. MYID .EQ. 0 ) 
C    <	WRITE(*,2) LEVEL, N, RESID, RESID/BBSUM, ERMIN, ERMAX, 
C    <	                     RESBC, RESBC/SUMBC
C   2	FORMAT(' L', I1, ' #', I2, ' CONVERGED ', 2E9.3, 
C    <	       ' E', 2E9.3, ' BC', 2E9.3)
C     	   return
C     	endif
     	
     	if ( res_new / res_old .gt. slowiter(level) ) then
     	IF ( MOD(ISTEP,NSAVE) .EQ. 0 .and. MYID .EQ. 0 ) 
     <	WRITE(*,3) LEVEL, N, RESID, RESID/BBSUM, ERMIN, ERMAX, 
     <	                     RESBC, RESBC/SUMBC
    3	FORMAT(' L', I1, ' #', I2, ' TOO SLOW! ', 2E9.3, 
     <	       ' E', 2E9.3, ' BC', 2E9.3)
	   return
     	else
     	   res_old = res_new
     	endif
     
	enddo
	
     	IF ( MOD(ISTEP,NSAVE) .EQ. 0 .and. MYID .EQ. 0 ) 
     <	WRITE(*,4) LEVEL, N-1, RESID, RESID/BBSUM, ERMIN, ERMAX, 
     <	                     RESBC, RESBC/SUMBC
    4	FORMAT(' L', I1, ' #', I2, ' MAX STEP! ', 2E9.3, 
     <	       ' E', 2E9.3, ' BC', 2E9.3)
    
	return
	end	

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine smoothj( level, ii, jj, kk, i0, k0, p, r, b,
     <		g11, g12, g13, g21, g22, g23, g31, g32, g33, gcc,
     <		t11, t12, t13, t21, t22, t23, t31, t32, t33, tcc )

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "para.inc"

	integer level, ii, jj, kk, i0, k0
	double precision, dimension(0:ii+1,0:jj+1,0:kk+1) :: p, r, b
	
	double precision, dimension(-1:ii+2,-1:jj+2,-1:kk+2) ::
     <		g11, g12, g13, g21, g22, g23, g31, g32, g33, gcc,
     <		t11, t12, t13, t21, t22, t23, t31, t32, t33, tcc

	double precision, dimension(ii/2,0:jj+1) ::
     <		aij, bij, cij, fij
	double precision, dimension(kk/2,0:jj+1) ::
     <		akj, bkj, ckj, fkj

	integer i, j, k, id, kd, is, ie, ks, ke
	integer iis, iie, kks, kke

	if ( i0 .eq. 1 ) then
	   is  = 1
	   ie  = ii-1
	   iis = 1
	   iie = ii+1
	else
	   is  = 2
	   ie  = ii
	   iis = 0
	   iie = ii
	endif
	
	if ( k0 .eq. 1 ) then
	   ks  = 1
	   ke  = kk-1
	   kks = 1
	   kke = kk+1
	else
	   ks  = 2
	   ke  = kk
	   kks = 0
	   kke = kk
	endif
	
	do k = ks, ke, 2
	do j =  1, jj
	do i = is, ie, 2

	   r(i,j,k) = - b(i,j,k) + 
     <		( t11(i,  j,  k  ) * p(i+1,j,  k  )
     <		+ t11(i-1,j,  k  ) * p(i-1,j,  k  )
     <		+ t33(i,  j,  k  ) * p(i,  j,  k+1)
     <		+ t33(i,  j,  k-1) * p(i,  j,  k-1) )

	   r(i,j,k) = r(i,j,k) +
     <		( t12(i,  j,k) * ( p(i+1,j+1,k) - p(i+1,j-1,k) )
     <		- t12(i-1,j,k) * ( p(i-1,j+1,k) - p(i-1,j-1,k) )
     <		+ t13(i,  j,k) * ( p(i,  j,k+1) - p(i,  j,k-1)
     <		                 + p(i+1,j,k+1) - p(i+1,j,k-1) )
     <		- t13(i-1,j,k) * ( p(i,  j,k+1) - p(i,  j,k-1)
     <		                 + p(i-1,j,k+1) - p(i-1,j,k-1) ) )

	   r(i,j,k) = r(i,j,k) +
     <		( t23(i,j,  k) * ( p(i,j,  k+1) - p(i,j,  k-1)
     <		                 + p(i,j+1,k+1) - p(i,j+1,k-1) )
     <		- t23(i,j-1,k) * ( p(i,j,  k+1) - p(i,j,  k-1)
     <		                 + p(i,j-1,k+1) - p(i,j-1,k-1) ) 
     <		+ t21(i,j,  k) * ( p(i+1,j,  k) - p(i-1,j,  k)
     <		                 + p(i+1,j+1,k) - p(i-1,j+1,k) )
     <		- t21(i,j-1,k) * ( p(i+1,j,  k) - p(i-1,j,  k)
     <		                 + p(i+1,j-1,k) - p(i-1,j-1,k) ) )

	   r(i,j,k) = r(i,j,k) +
     <		( t31(i,j,k  ) * ( p(i+1,j,k  ) - p(i-1,j,k  )
     <		                 + p(i+1,j,k+1) - p(i-1,j,k+1) )
     <		- t31(i,j,k-1) * ( p(i+1,j,k  ) - p(i-1,j,k  )
     <		                 + p(i+1,j,k-1) - p(i-1,j,k-1) )
     <		+ t32(i,j,k  ) * ( p(i,j+1,k+1) - p(i,j-1,k+1) )
     <		- t32(i,j,k-1) * ( p(i,j+1,k-1) - p(i,j-1,k-1) ) )
 
	enddo
	enddo
	enddo
	
	do k = ks, ke, 2

	do j =  1, jj
	   id = 0
	   do i = is, ie, 2
	      id = id + 1
	      aij(id,j) = - t22(i,j-1,k) + t12(i,j,k) - t12(i-1,j,k)
     <	                                 + t32(i,j,k) - t32(i,j,k-1)
	      cij(id,j) = - t22(i,j,  k) - t12(i,j,k) + t12(i-1,j,k)
     <	                                 - t32(i,j,k) + t32(i,j,k-1)
	      bij(id,j) =   tcc(i,j,  k)
	      fij(id,j) =     r(i,j,  k)
	   enddo
	enddo

	if ( n_suth .eq. MPI_PROC_NULL ) then
	   j = 0
	   id = 0
	   do i = is, ie, 2
	      id = id + 1
	      aij(id,j) =   0.D0
	      bij(id,j) =   g22(i,j,k)
	      cij(id,j) = - g22(i,j,k)
	      fij(id,j) = -   b(i,j,k) + 
     <		( g23(i,j,  k) * ( p(i,j,  k+1) - p(i,j,  k-1)
     <		                 + p(i,j+1,k+1) - p(i,j+1,k-1) )
     <		+ g21(i,j,  k) * ( p(i+1,j,  k) - p(i-1,j,  k)
     <		                 + p(i+1,j+1,k) - p(i-1,j+1,k) ) )
	   enddo
	endif
	
	if ( n_nrth .eq. MPI_PROC_NULL ) then
           j = jj + 1
	   id = 0
	   do i = is, ie, 2
	      id = id + 1
	      aij(id,j) =   g22(i,j-1,k)
	      bij(id,j) = - g22(i,j-1,k)
	      cij(id,j) = 0.D0
	      fij(id,j) = -   b(i,j  ,k) +
     <		( g23(i,j-1,k) * ( p(i,j,  k+1) - p(i,j,  k-1)
     <		                 + p(i,j-1,k+1) - p(i,j-1,k-1) ) 
     <		+ g21(i,j-1,k) * ( p(i+1,j,  k) - p(i-1,j,  k)
     <		                 + p(i+1,j-1,k) - p(i-1,j-1,k) ) )
	   enddo
	endif
	
	call trid( aij, bij, cij, fij, ii/2, 1, jj, n_suth, n_nrth )

	do j =  0, jj+1
	   id = 0
	   do i = is, ie, 2
	      id = id + 1
	      p(i,j,k) = fij(id,j)
	   enddo
	enddo

	enddo

CCC	call p_exchange( ii, jj, kk, p ) 

CCC	IF ( 0 .EQ. 1 ) THEN

	call pj_exchange( ii, jj, kk, p, is, ie, 2, ks, ke, 2,
     <		jj,    0, n_nrth, n_suth )
	call pj_exchange( ii, jj, kk, p, is, ie, 2, ks, ke, 2,
     <		 1, jj+1, n_suth, n_nrth )
	
	if ( i0 .eq. 1 ) then
	call pi_exchange( ii, jj, kk, p, 0, jj+1, 1, kks, kke, 2,
     <		 1, ii+1, n_west, n_east )
     	else
	call pi_exchange( ii, jj, kk, p, 0, jj+1, 1, kks, kke, 2,
     <		ii,    0, n_east, n_west )
     	endif

	if ( k0 .eq. 1 ) then
	call pk_exchange( ii, jj, kk, p, iis, iie, 2, 0, jj+1, 1,
     <		 1, kk+1, n_back, n_frnt )
     	else
	call pk_exchange( ii, jj, kk, p, iis, iie, 2, 0, jj+1, 1,
     <		kk,    0, n_frnt, n_back )
     	endif

CCC	ENDIF

	return
	end
