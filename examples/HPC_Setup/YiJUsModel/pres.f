ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine pressure

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "ns.inc"
	include "metric.inc"
	include "para.inc"
	include "sedi.inc"

	double precision r(0:nni +1,0:nnj +1,0:nnk +1), 
     &                   b(0:nni +1,0:nnj +1,0:nnk +1)
	double precision p1(0:nni1+1,0:nnj1+1,0:nnk1+1), 
     &                   r1(0:nni1+1,0:nnj1+1,0:nnk1+1), 
     &                   b1(0:nni1+1,0:nnj1+1,0:nnk1+1) 
	double precision p2(0:nni2+1,0:nnj2+1,0:nnk2+1), 
     &                   r2(0:nni2+1,0:nnj2+1,0:nnk2+1), 
     &                   b2(0:nni2+1,0:nnj2+1,0:nnk2+1)
	double precision p3(0:nni3+1,0:nnj3+1,0:nnk3+1), 
     &                   r3(0:nni3+1,0:nnj3+1,0:nnk3+1), 
     &                   b3(0:nni3+1,0:nnj3+1,0:nnk3+1)
	double precision p4(0:nni4+1,0:nnj4+1,0:nnk4+1), 
     &                   r4(0:nni4+1,0:nnj4+1,0:nnk4+1), 
     &                   b4(0:nni4+1,0:nnj4+1,0:nnk4+1)
	double precision p5(0:nni5+1,0:nnj5+1,0:nnk5+1), 
     &                   r5(0:nni5+1,0:nnj5+1,0:nnk5+1), 
     &                   b5(0:nni5+1,0:nnj5+1,0:nnk5+1)

	integer i, j, k, L, n

	double precision resid, bbsum, ermin, ermax, resbc, sumbc
	double precision temp

	resid = 0.D0

	do k = 0, nnk+1
	do j = 0, nnj+1
	do i = 0, nni+1
	   p(i,j,k) = 0.D0
	   r(i,j,k) = 0.D0
	   b(i,j,k) = 0.D0
	enddo
	enddo
	enddo

	do k = 0, nnk1+1
	do j = 0, nnj1+1
	do i = 0, nni1+1
	   p1(i,j,k) = 0.D0
	   r1(i,j,k) = 0.D0
	   b1(i,j,k) = 0.D0
	enddo
	enddo
	enddo

	do k = 0, nnk2+1
	do j = 0, nnj2+1
	do i = 0, nni2+1
	   p2(i,j,k) = 0.D0
	   r2(i,j,k) = 0.D0
	   b2(i,j,k) = 0.D0
	enddo
	enddo
	enddo

	do k = 0, nnk3+1
	do j = 0, nnj3+1
	do i = 0, nni3+1
	   p3(i,j,k) = 0.D0
	   r3(i,j,k) = 0.D0
	   b3(i,j,k) = 0.D0
	enddo
	enddo
	enddo

	do k = 0, nnk4+1
	do j = 0, nnj4+1
	do i = 0, nni4+1
	   p4(i,j,k) = 0.D0
	   r4(i,j,k) = 0.D0
	   b4(i,j,k) = 0.D0
	enddo
	enddo
	enddo

	do k = 0, nnk5+1
	do j = 0, nnj5+1
	do i = 0, nni5+1
	   p5(i,j,k) = 0.D0
	   r5(i,j,k) = 0.D0
	   b5(i,j,k) = 0.D0
	enddo
	enddo
	enddo

	
        call quick
	call totalU

CBCBCBC	BCBCBCBCBCBCBCBCBCBC

	temp = 0.125D0 / dtime

	if ( n_west .eq. MPI_PROC_NULL ) then
	do k = 1, nnk
	do j = 1, nnj
	   b(0,j,k) = temp *
     <	      ( xix(0,j,k) * ( 15.D0 * upf(1,j,k,1) 
     <	                     - 10.D0 * upf(2,j,k,1) 
     <                       +  3.D0 * upf(3,j,k,1) )
     <	      + xiy(0,j,k) * ( 15.D0 * upf(1,j,k,2)
     <	                     - 10.D0 * upf(2,j,k,2) 
     <                       +  3.D0 * upf(3,j,k,2) )
     <	      + xiz(0,j,k) * ( 15.D0 * upf(1,j,k,3) 
     <	                     - 10.D0 * upf(2,j,k,3) 
     <                       +  3.D0 * upf(3,j,k,3) ) )
	   uxi(0,j,k) = 0.D0
	enddo
	enddo
	endif
	if ( n_east .eq. MPI_PROC_NULL ) then
	do k = 1, nnk
	do j = 1, nnj
	   b(nni+1,j,k) = temp *
     <	      ( xix(nni,j,k) * ( 15.D0 * upf(nni,  j,k,1) 
     <	                       - 10.D0 * upf(nni-1,j,k,1) 
     <                         +  3.D0 * upf(nni-2,j,k,1) )
     <	      + xiy(nni,j,k) * ( 15.D0 * upf(nni,  j,k,2)
     <	                       - 10.D0 * upf(nni-1,j,k,2) 
     <                         +  3.D0 * upf(nni-2,j,k,2) )
     <	      + xiz(nni,j,k) * ( 15.D0 * upf(nni,  j,k,3) 
     <	                       - 10.D0 * upf(nni-1,j,k,3) 
     <                         +  3.D0 * upf(nni-2,j,k,3) ) )
	   uxi(nni,j,k) = 0.D0
	enddo
	enddo
	endif

	if ( n_suth .eq. MPI_PROC_NULL ) then
	do k = 1, nnk
	do i = 1, nni
	   
	   b(i,0,k) = temp *
     <	      ( etx(i,0,k) * ( 15.D0 * upf(i,1,k,1) 
     <	                     - 10.D0 * upf(i,2,k,1) 
     <                       +  3.D0 * upf(i,3,k,1) )	     
     <	      + ety(i,0,k) * ( 15.D0 * upf(i,1,k,2)
     <	                     - 10.D0 * upf(i,2,k,2) 
     <                       +  3.D0 * upf(i,3,k,2) )
     <	      + etz(i,0,k) * ( 15.D0 * upf(i,1,k,3) 
     <	                     - 10.D0 * upf(i,2,k,3) 
     <                       +  3.D0 * upf(i,3,k,3) ) )
	   if (movingGrid .ne. 1) then
  	      uej(i,0,k) = 0.D0
	   else
	      if (mod(istep, ncheck) .ne. 0.D0) then
		 uej(i,0,k) = 0.D0
	      endif
	   endif
	enddo
	enddo
	endif
	if ( n_nrth .eq. MPI_PROC_NULL ) then
	do k = 1, nnk
	do i = 1, nni
	   if (open_top .ne. 0) then
	   b(i,nnj+1,k) = temp *
     <	      ( etx(i,nnj,k) * ( 15.D0 * upf(i,nnj,  k,1) 
     <	                       - 10.D0 * upf(i,nnj-1,k,1) 
     <                         +  3.D0 * upf(i,nnj-2,k,1) )
     <	      + ety(i,nnj,k) * ( 15.D0 * upf(i,nnj,  k,2)
     <	                       - 10.D0 * upf(i,nnj-1,k,2) 
     <                         +  3.D0 * upf(i,nnj-2,k,2) )
     <	      + etz(i,nnj,k) * ( 15.D0 * upf(i,nnj,  k,3) 
     <	                       - 10.D0 * u(i,nnj-1,k,3) 
     <                         +  3.D0 * u(i,nnj-2,k,3) ) )
c	   uej(i,nnj,k) = 0.D0
	   else
	   b(i,nnj+1,k) = temp *
     <	      ( etx(i,nnj,k) * ( 15.D0 * upf(i,nnj,  k,1) 
     <	                       - 10.D0 * upf(i,nnj-1,k,1) 
     <                         +  3.D0 * upf(i,nnj-2,k,1) )
     <	      + ety(i,nnj,k) * ( 15.D0 * upf(i,nnj,  k,2)
     <	                       - 10.D0 * upf(i,nnj-1,k,2) 
     <                         +  3.D0 * upf(i,nnj-2,k,2) )
     <	      + etz(i,nnj,k) * ( 15.D0 * upf(i,nnj,  k,3) 
     <	                       - 10.D0 * upf(i,nnj-1,k,3) 
     <                         +  3.D0 * upf(i,nnj-2,k,3) ) )
	   
	   uej(i,nnj,k) = 0.D0
	   endif

	enddo
	enddo
	endif

	if ( n_back .eq. MPI_PROC_NULL ) then
	do j = 1, nnj
	do i = 1, nni
	   b(i,j,0) = temp *
     <	      ( ztx(i,j,0) * ( 15.D0 * upf(i,j,1,1) 
     <	                     - 10.D0 * upf(i,j,2,1) 
     <                       +  3.D0 * upf(i,j,3,1) )
     <	      + zty(i,j,0) * ( 15.D0 * upf(i,j,1,2)
     <	                     - 10.D0 * upf(i,j,2,2) 
     <                       +  3.D0 * upf(i,j,3,2) )
     <	      + ztz(i,j,0) * ( 15.D0 * upf(i,j,1,3) 
     <	                     - 10.D0 * upf(i,j,2,3) 
     <                       +  3.D0 * upf(i,j,3,3) ) )
	   uzk(i,j,0) = 0.D0
	enddo
	enddo
	endif
	if ( n_frnt .eq. MPI_PROC_NULL ) then
	do j = 1, nnj
	do i = 1, nni
	   b(i,j,nnk+1) = temp *
     <	      ( ztx(i,j,nnk) * ( 15.D0 * upf(i,j,nnk,  1) 
     <	                       - 10.D0 * upf(i,j,nnk-1,1) 
     <                         +  3.D0 * upf(i,j,nnk-2,1) )
     <	      + zty(i,j,nnk) * ( 15.D0 * upf(i,j,nnk,  2)
     <	                       - 10.D0 * upf(i,j,nnk-1,2) 
     <                         +  3.D0 * upf(i,j,nnk-2,2) )
     <	      + ztz(i,j,nnk) * ( 15.D0 * upf(i,j,nnk,  3) 
     <	                       - 10.D0 * upf(i,j,nnk-1,3) 
     <                         +  3.D0 * upf(i,j,nnk-2,3) ) )
	   uzk(i,j,nnk) = 0.D0
	enddo
	enddo
	endif

CBCBCBC	BCBCBCBCBCBCBCBCBCBC

	temp = 1.D0 / dtime

	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni
c	   if (open_top .eq. 1) then
c	   if (istep .eq. 1) then
c           if (j .eq. 1) then
c	      uxi(i,j,k) = 0.D0
c	   endif
c	   if (j .eq. nnj) then
c	      uxi(i,nnj-1,k) = 0.D0
c	   endif
c	   endif
c	   endif
	   b(i,j,k) = temp * ( ( uxi(i,j,k) - uxi(i-1,j,k) )
     <	                     + ( uej(i,j,k) - uej(i,j-1,k) )
     <	                     + ( uzk(i,j,k) - uzk(i,j,k-1) )
     <                       + ( wej(i,j,k)*sedf(i,j,k,2)
     <                         - wej(i,j-1,k)*sedf(i,j-1,k,2) )	)   
c	   if (j .eq. 1) then
c	   write(*,*) uej(i,j,k) - uej(i,j-1,k) 
c	   endif
	   
	enddo
	enddo
	enddo

	L = 1

	call smooth( resid, bbsum, ermin, ermax, resbc, sumbc,
     <		L, nni, nnj, nnk, p, r, b, jac,
     <		g11, g12, g13, g21, g22, g23, g31, g32, g33, gcc )

	do n = 1, maxstep

	if ( mod(istep,nsave) .eq. 0 .and. MYID .eq. 0 ) then
	   write(*,*) ' V-cycle # ', n
	endif

	if ( mg_level .ge. 2 ) then

	call thrd_f2c( nni, nnj, nnk, nni1, nnj1, nnk1, r, b1 )

	L = 2

	do k = 0, nnk1+1
	do j = 0, nnj1+1
	do i = 0, nni1+1
	   p1(i,j,k) = 0.D0
	enddo
	enddo
	enddo

	call smooth( resid, bbsum, ermin, ermax, resbc, sumbc,
     <		L, nni1, nnj1, nnk1, p1, r1, b1, jaf,
     <		f11, f12, f13, f21, f22, f23, f31, f32, f33, fcc )

	if ( mg_level .ge. 3 ) then

	call thrd_f2c( nni1, nnj1, nnk1, nni2, nnj2, nnk2, r1, b2 )

	L = 3

	do k = 0, nnk2+1
	do j = 0, nnj2+1
	do i = 0, nni2+1
	   p2(i,j,k) = 0.D0
	enddo
	enddo
	enddo

	call smooth( resid, bbsum, ermin, ermax, resbc, sumbc,
     <		L, nni2, nnj2, nnk2, p2, r2, b2, jah,
     <		h11, h12, h13, h21, h22, h23, h31, h32, h33, hcc )

	if ( mg_level .ge. 4 ) then

	call thrd_f2c( nni2, nnj2, nnk2, nni3, nnj3, nnk3, r2, b3 )

	L = 4

	do k = 0, nnk3+1
	do j = 0, nnj3+1
	do i = 0, nni3+1
	   p3(i,j,k) = 0.D0
	enddo
	enddo
	enddo

	call smooth( resid, bbsum, ermin, ermax, resbc, sumbc,
     <		L, nni3, nnj3, nnk3, p3, r3, b3, jap,
     <		p11, p12, p13, p21, p22, p23, p31, p32, p33, pcc )

	if ( mg_level .ge. 5 ) then

	call thrd_f2c( nni3, nnj3, nnk3, nni4, nnj4, nnk4, r3, b4 )

	L = 5

	do k = 0, nnk4+1
	do j = 0, nnj4+1
	do i = 0, nni4+1
	   p4(i,j,k) = 0.D0
	enddo
	enddo
	enddo

	call smooth( resid, bbsum, ermin, ermax, resbc, sumbc,
     <		L, nni4, nnj4, nnk4, p4, r4, b4, jar,
     <		r11, r12, r13, r21, r22, r23, r31, r32, r33, rcc )

	if ( mg_level .ge. 6 ) then

	call thrd_f2c( nni4, nnj4, nnk4, nni5, nnj5, nnk5, r4, b5 )

	L = 6

	do k = 0, nnk5+1
	do j = 0, nnj5+1
	do i = 0, nni5+1
	   p5(i,j,k) = 0.D0
	enddo
	enddo
	enddo

	call smooth( resid, bbsum, ermin, ermax, resbc, sumbc,
     <		L, nni5, nnj5, nnk5, p5, r5, b5, jae,
     <		e11, e12, e13, e21, e22, e23, e31, e32, e33, ecc )

	call thrd_c2f( nni5, nnj5, nnk5, nni4, nnj4, nnk4, p5, p4 )


	L = 5

	call smooth( resid, bbsum, ermin, ermax, resbc, sumbc,
     <		L, nni4, nnj4, nnk4, p4, r4, b4, jar,
     <		r11, r12, r13, r21, r22, r23, r31, r32, r33, rcc )

C	6
	endif

	call thrd_c2f( nni4, nnj4, nnk4, nni3, nnj3, nnk3, p4, p3 )

	L = 4

	call smooth( resid, bbsum, ermin, ermax, resbc, sumbc,
     <		L, nni3, nnj3, nnk3, p3, r3, b3, jap,
     <		p11, p12, p13, p21, p22, p23, p31, p32, p33, pcc )

C	5
	endif

	call thrd_c2f( nni3, nnj3, nnk3, nni2, nnj2, nnk2, p3, p2 )

	L = 3

	call smooth( resid, bbsum, ermin, ermax, resbc, sumbc,
     <		L, nni2, nnj2, nnk2, p2, r2, b2, jah,
     <		h11, h12, h13, h21, h22, h23, h31, h32, h33, hcc )

C	4
	endif

	call thrd_c2f( nni2, nnj2, nnk2, nni1, nnj1, nnk1, p2, p1 )

	L = 2

	call smooth( resid, bbsum, ermin, ermax, resbc, sumbc,
     <		L, nni1, nnj1, nnk1, p1, r1, b1, jaf,
     <		f11, f12, f13, f21, f22, f23, f31, f32, f33, fcc )

C	3
	endif

	call thrd_c2f( nni1, nnj1, nnk1, nni, nnj, nnk, p1, p )

C	2
	endif

	L = 1

	call smooth( resid, bbsum, ermin, ermax, resbc, sumbc,
     <		L, nni, nnj, nnk, p, r, b, jac,
     <		g11, g12, g13, g21, g22, g23, g31, g32, g33, gcc )
	
	if ( resid .lt. tol(L) .and. resbc .lt. tol(L) .and.
     <	     max(dabs(ermin), dabs(ermax)) .lt. ter(L) .and. 
     <	     resid/bbsum .lt. factor ) then	     
	   if ( MYID .EQ. 0 ) 
     <	      write(*,*) ' Total V-cycle # ', n, resid
	   return
	endif
c        write(*,*) ' Total V-cycle # ', n, resid 
	enddo

	return
	end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine openP_bc_nrth
	
	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "ns.inc"
	include "metric.inc"
	include "para.inc"
	
	integer i,j,k

	if (n_nrth .eq. MPI_PROC_NULL) then
	   j = nnj+1
	   do k = 0, nnk+1
	   do j = 0, nnj+1
	      p(i,j,k) = 2.D0*p(i,j-1,k) - p(i,j-2,k)
	   enddo
	   enddo
	endif

	return 
	end
