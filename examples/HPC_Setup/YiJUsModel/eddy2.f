cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine eddy2

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "ns.inc"
	include "metric.inc"
	include "para.inc"
	include "eddy.inc"
	include "sedi.inc"

	double precision ss(0:nni+1,0:nnj+1,0:nnk+1,9), 
     &                   tt(0:nni+1,0:nnj+1,0:nnk+1,9), 
     &                   zz(0:nni+1,0:nnj+1,0:nnk+1,9)
     
	integer i, j, k, m, jj

	double precision tiny, alpha2, weit, fact, temp
	double precision fmin, fmax, tmpzz
	double precision hc, C, absU, az, dz, temp, sum

	tiny = 1.D-20

	alpha2 = 4.D0

        do m = 1, 9
        do k = 0, nnk+1
        do j = 0, nnj+1
        do i = 0, nni+1
	   ss(i,j,k,m) = 0.D0
	   tt(i,j,k,m) = 0.D0
	   zz(i,j,k,m) = 0.D0
	enddo
	enddo
	enddo 
	enddo 

        do m = 1, 6
        do k = 0, nnk+1
        do j = 0, nnj+1
        do i = 0, nni+1
	   rr(i,j,k,m) = 0.D0
	enddo
	enddo
	enddo 
	enddo 

C	T:
C	xix (1) etx (2) ztx (3)
C	xiy (4) ety (5) zty (6)
C	xiz (7) etz (8) ztz (9)

        do k = 0, nnk+1
        do j = 0, nnj+1
        do i = 0, nni+1
c Multiplying jac by 0.25 is equal to dividing by 2 here and 
c again in the next step.... 
	   temp = 0.25D0 * jac(i,j,k)
	   tt(i,j,k,1) = temp * ( xix(i-1,j,k) + xix(i,j,k) )
	   tt(i,j,k,2) = temp * ( etx(i,j-1,k) + etx(i,j,k) )
	   tt(i,j,k,3) = temp * ( ztx(i,j,k-1) + ztx(i,j,k) )
	   tt(i,j,k,4) = temp * ( xiy(i-1,j,k) + xiy(i,j,k) )
	   tt(i,j,k,5) = temp * ( ety(i,j-1,k) + ety(i,j,k) )
	   tt(i,j,k,6) = temp * ( zty(i,j,k-1) + zty(i,j,k) )
	   tt(i,j,k,7) = temp * ( xiz(i-1,j,k) + xiz(i,j,k) )
	   tt(i,j,k,8) = temp * ( etz(i,j-1,k) + etz(i,j,k) )
	   tt(i,j,k,9) = temp * ( ztz(i,j,k-1) + ztz(i,j,k) )
	enddo
	enddo
	enddo 

C	Z:
C	uxi (1) uet (2) uzt (3)
C	vxi (4) vet (5) vzt (6)
C	wxi (7) wet (8) wzt (9)

        do k = 0, nnk+1
        do j = 0, nnj+1
        do i = 0, nni+1
	   zz(i,j,k,1) = u(i+1,j,k,1) - u(i-1,j,k,1)
	   zz(i,j,k,2) = u(i,j+1,k,1) - u(i,j-1,k,1)
	   zz(i,j,k,3) = u(i,j,k+1,1) - u(i,j,k-1,1)
	   zz(i,j,k,4) = u(i+1,j,k,2) - u(i-1,j,k,2)
	   zz(i,j,k,5) = u(i,j+1,k,2) - u(i,j-1,k,2)
	   zz(i,j,k,6) = u(i,j,k+1,2) - u(i,j,k-1,2)
	   zz(i,j,k,7) = u(i+1,j,k,3) - u(i-1,j,k,3)
	   zz(i,j,k,8) = u(i,j+1,k,3) - u(i,j-1,k,3)
	   zz(i,j,k,9) = u(i,j,k+1,3) - u(i,j,k-1,3)
	enddo
	enddo
	enddo 

C	R:
C	rxi (1) ret (2) rzt (3)

        do k = 0, nnk+1
        do j = 0, nnj+1
        do i = 0, nni+1
c	   rr(i,j,k,1) = phi(i+1,j,k) - phi(i-1,j,k)
c	   rr(i,j,k,2) = phi(i,j+1,k) - phi(i,j-1,k)
c	   rr(i,j,k,3) = phi(i,j,k+1) - phi(i,j,k-1)
           rr(i,j,k,1) = sedc(i+1,j,k) - sedc(i-1,j,k)
	   rr(i,j,k,2) = sedc(i,j+1,k) - sedc(i,j-1,k)
	   rr(i,j,k,3) = sedc(i,j,k+1) - sedc(i,j,k-1)
c           write(*,*) i,j,k,sedc(i,j+1,k), sedc(i,j-1,k), rr(i,j,k,2)
	enddo
	enddo
	enddo 

c- Now obtain the strain rate tensor for Smargorinsky term
 
C	S:
C	S11 (1) S22 (2) S33 (3)
C	S12 (4) S23 (5) S31 (6)
C	rsx (7) rsy (8) rsz (9)

        do k = 0, nnk+1
        do j = 0, nnj+1
        do i = 0, nni+1
	   ss(i,j,k,1) = tt(i,j,k,1) * zz(i,j,k,1)
     <	               + tt(i,j,k,2) * zz(i,j,k,2)
     <	               + tt(i,j,k,3) * zz(i,j,k,3)
	   ss(i,j,k,2) = tt(i,j,k,4) * zz(i,j,k,4)
     <	               + tt(i,j,k,5) * zz(i,j,k,5)
     <	               + tt(i,j,k,6) * zz(i,j,k,6)
	   ss(i,j,k,3) = tt(i,j,k,7) * zz(i,j,k,7)
     <	               + tt(i,j,k,8) * zz(i,j,k,8)
     <	               + tt(i,j,k,9) * zz(i,j,k,9)
	   ss(i,j,k,4) = tt(i,j,k,4) * zz(i,j,k,1)
     <	               + tt(i,j,k,5) * zz(i,j,k,2)
     <	               + tt(i,j,k,6) * zz(i,j,k,3)
     <	               + tt(i,j,k,1) * zz(i,j,k,4)
     <	               + tt(i,j,k,2) * zz(i,j,k,5)
     <	               + tt(i,j,k,3) * zz(i,j,k,6)
	   ss(i,j,k,5) = tt(i,j,k,7) * zz(i,j,k,4)
     <	               + tt(i,j,k,8) * zz(i,j,k,5)
     <	               + tt(i,j,k,9) * zz(i,j,k,6)
     <	               + tt(i,j,k,4) * zz(i,j,k,7)
     <	               + tt(i,j,k,5) * zz(i,j,k,8)
     <	               + tt(i,j,k,6) * zz(i,j,k,9)
	   ss(i,j,k,6) = tt(i,j,k,1) * zz(i,j,k,7)
     <	               + tt(i,j,k,2) * zz(i,j,k,8)
     <	               + tt(i,j,k,3) * zz(i,j,k,9)
     <	               + tt(i,j,k,7) * zz(i,j,k,1)
     <	               + tt(i,j,k,8) * zz(i,j,k,2)
     <	               + tt(i,j,k,9) * zz(i,j,k,3)
	   ss(i,j,k,4) = 0.5D0 * ss(i,j,k,4)
	   ss(i,j,k,5) = 0.5D0 * ss(i,j,k,5)
	   ss(i,j,k,6) = 0.5D0 * ss(i,j,k,6)
	   ss(i,j,k,7) = tt(i,j,k,1) * rr(i,j,k,1)
     <	               + tt(i,j,k,2) * rr(i,j,k,2)
     <	               + tt(i,j,k,3) * rr(i,j,k,3)
	   ss(i,j,k,8) = tt(i,j,k,4) * rr(i,j,k,1)
     <	               + tt(i,j,k,5) * rr(i,j,k,2)
     <	               + tt(i,j,k,6) * rr(i,j,k,3)
	   ss(i,j,k,9) = tt(i,j,k,7) * rr(i,j,k,1)
     <	               + tt(i,j,k,8) * rr(i,j,k,2)
     <	               + tt(i,j,k,9) * rr(i,j,k,3)
	enddo
	enddo
	enddo 

C	|S|

        do k = 0, nnk+1
        do j = 0, nnj+1
        do i = 0, nni+1
	   sab(i,j,k) = 
     <	              ss(i,j,k,1)**2 + ss(i,j,k,2)**2 + ss(i,j,k,3)**2
     <	   + 2.D0 * ( ss(i,j,k,4)**2 + ss(i,j,k,5)**2 + ss(i,j,k,6)**2 )
	   sab(i,j,k) = dsqrt( 2.D0 * sab(i,j,k) )
	enddo
	enddo
	enddo 

C	Z:
C	|S|S11 (1) |S|S22 (2) |S|S33 (3)
C	|S|S12 (4) |S|S23 (5) |S|S31 (6)
C	|S|rsx (7) |S|rsy (8) |S|rsz (9)

        do k = 0, nnk+1
        do j = 0, nnj+1
        do i = 0, nni+1
	   zz(i,j,k,1) = sab(i,j,k) * ss(i,j,k,1)
	   zz(i,j,k,2) = sab(i,j,k) * ss(i,j,k,2)
	   zz(i,j,k,3) = sab(i,j,k) * ss(i,j,k,3)
	   zz(i,j,k,4) = sab(i,j,k) * ss(i,j,k,4)
	   zz(i,j,k,5) = sab(i,j,k) * ss(i,j,k,5)
	   zz(i,j,k,6) = sab(i,j,k) * ss(i,j,k,6)
	   zz(i,j,k,7) = sab(i,j,k) * ss(i,j,k,7)
	   zz(i,j,k,8) = sab(i,j,k) * ss(i,j,k,8)
	   zz(i,j,k,9) = sab(i,j,k) * ss(i,j,k,9)
	enddo
	enddo
	enddo 

C	Z:
C	<|S|S11> (1) <|S|S22> (2) <|S|S33> (3)
C	<|S|S12> (4) <|S|S23> (5) <|S|S31> (6)
C	<|S|rsx> (7) <|S|rsy> (8) <|S|rsz> (9)

	weit = 2.D0
	fact = 1.D0 / 64.D0 
	call filter(nni,nnj,nnk, zz(0,0,0,1),rr(0,0,0,5), weit,fact)
	call filter(nni,nnj,nnk, zz(0,0,0,2),rr(0,0,0,5), weit,fact)
	call filter(nni,nnj,nnk, zz(0,0,0,3),rr(0,0,0,5), weit,fact)
	call filter(nni,nnj,nnk, zz(0,0,0,4),rr(0,0,0,5), weit,fact)
	call filter(nni,nnj,nnk, zz(0,0,0,5),rr(0,0,0,5), weit,fact)
	call filter(nni,nnj,nnk, zz(0,0,0,6),rr(0,0,0,5), weit,fact)
	call filter(nni,nnj,nnk, zz(0,0,0,7),rr(0,0,0,5), weit,fact)
	call filter(nni,nnj,nnk, zz(0,0,0,8),rr(0,0,0,5), weit,fact)
	call filter(nni,nnj,nnk, zz(0,0,0,9),rr(0,0,0,5), weit,fact)

C	R:
C	|S| (4)

        do k = 0, nnk+1
        do j = 0, nnj+1
        do i = 0, nni+1
	   rr(i,j,k,4) = sab(i,j,k)
	enddo
	enddo
	enddo 

C	R:
C	<|S|> (4)

	weit = 2.D0
	fact = 1.D0 / 64.D0 
	call filter(nni,nnj,nnk, rr(0,0,0,4),rr(0,0,0,5), weit,fact)

C	S:
C	<S11> (1) <S22> (2) <S33> (3)
C	<S12> (1) <S23> (2) <S31> (3)
C	<rsx> (1) <rsy> (2) <rsz> (3)

	weit = 2.D0
	fact = 1.D0 / 64.D0 
	call filter(nni,nnj,nnk, ss(0,0,0,1),rr(0,0,0,5), weit,fact)
	call filter(nni,nnj,nnk, ss(0,0,0,2),rr(0,0,0,5), weit,fact)
	call filter(nni,nnj,nnk, ss(0,0,0,3),rr(0,0,0,5), weit,fact)
	call filter(nni,nnj,nnk, ss(0,0,0,4),rr(0,0,0,5), weit,fact)
	call filter(nni,nnj,nnk, ss(0,0,0,5),rr(0,0,0,5), weit,fact)
	call filter(nni,nnj,nnk, ss(0,0,0,6),rr(0,0,0,5), weit,fact)
	call filter(nni,nnj,nnk, ss(0,0,0,7),rr(0,0,0,5), weit,fact)
	call filter(nni,nnj,nnk, ss(0,0,0,8),rr(0,0,0,5), weit,fact)
	call filter(nni,nnj,nnk, ss(0,0,0,9),rr(0,0,0,5), weit,fact)


C	Z:
C	M11 (1) M22 (2) M33 (3)
C	M12 (4) M23 (5) M31 (6)
C	N1  (7) N2  (8) N3  (9)

        do k = 1, nnk
        do j = 1, nnj
        do i = 1, nni
	   zz(i,j,k,1) = alpha2 * rr(i,j,k,4) * ss(i,j,k,1) - zz(i,j,k,1) 
	   zz(i,j,k,2) = alpha2 * rr(i,j,k,4) * ss(i,j,k,2) - zz(i,j,k,2) 
	   zz(i,j,k,3) = alpha2 * rr(i,j,k,4) * ss(i,j,k,3) - zz(i,j,k,3) 
	   zz(i,j,k,4) = alpha2 * rr(i,j,k,4) * ss(i,j,k,4) - zz(i,j,k,4) 
	   zz(i,j,k,5) = alpha2 * rr(i,j,k,4) * ss(i,j,k,5) - zz(i,j,k,5) 
	   zz(i,j,k,6) = alpha2 * rr(i,j,k,4) * ss(i,j,k,6) - zz(i,j,k,6) 
	   zz(i,j,k,7) = alpha2 * rr(i,j,k,4) * ss(i,j,k,7) - zz(i,j,k,7) 
	   zz(i,j,k,8) = alpha2 * rr(i,j,k,4) * ss(i,j,k,8) - zz(i,j,k,8)  
	   zz(i,j,k,9) = alpha2 * rr(i,j,k,4) * ss(i,j,k,9) - zz(i,j,k,9)
c           write(*,*) zz(i,j,k,7), zz(i,j,k,8), zz(i,j,k,9)  
	enddo
	enddo
	enddo 

C	T:
C	uu (1) vv (2) ww (3)
C	uv (4) vw (5) wu (6)
C	ru (7) rv (8) rw (9)

        do k = 0, nnk+1
        do j = 0, nnj+1
        do i = 0, nni+1
	   tt(i,j,k,1) = u(i,j,k,1)**2
	   tt(i,j,k,2) = bsum(i,j,k,2)
	   tt(i,j,k,3) = u(i,j,k,3)**2
	   tt(i,j,k,4) = bsum(i,j,k,1)
	   tt(i,j,k,5) = bsum(i,j,k,3)
	   tt(i,j,k,6) = u(i,j,k,3) * u(i,j,k,1)
c          tt(i,j,k,7) = phi(i,j,k) * u(i,j,k,1)
c	   tt(i,j,k,8) = phi(i,j,k) * u(i,j,k,2)
c	   tt(i,j,k,9) = phi(i,j,k) * u(i,j,k,3)
	   tt(i,j,k,7) = sedc(i,j,k) * u(i,j,k,1)
	   tt(i,j,k,8) = sedc(i,j,k) * u(i,j,k,2)
	   tt(i,j,k,9) = sedc(i,j,k) * u(i,j,k,3)
	enddo
	enddo
	enddo 

C	T:
C	<uu> (1) <vv> (2) <ww> (3)
C	<uv> (4) <vw> (5) <wu> (6)
C	<ru> (7) <rv> (8) <rw> (9)

	weit = 2.D0
	fact = 1.D0 / 64.D0 
	call filter(nni,nnj,nnk, tt(0,0,0,1),rr(0,0,0,5), weit,fact)
	call filter(nni,nnj,nnk, tt(0,0,0,2),rr(0,0,0,5), weit,fact)
	call filter(nni,nnj,nnk, tt(0,0,0,3),rr(0,0,0,5), weit,fact)
	call filter(nni,nnj,nnk, tt(0,0,0,4),rr(0,0,0,5), weit,fact)
	call filter(nni,nnj,nnk, tt(0,0,0,5),rr(0,0,0,5), weit,fact)
	call filter(nni,nnj,nnk, tt(0,0,0,6),rr(0,0,0,5), weit,fact)
	call filter(nni,nnj,nnk, tt(0,0,0,7),rr(0,0,0,5), weit,fact)
	call filter(nni,nnj,nnk, tt(0,0,0,8),rr(0,0,0,5), weit,fact)
	call filter(nni,nnj,nnk, tt(0,0,0,9),rr(0,0,0,5), weit,fact)

C	R:
C	u (1) v (2) w (3) r (4)

        do k = 0, nnk+1
        do j = 0, nnj+1
        do i = 0, nni+1
	   rr(i,j,k,1) = u(i,j,k,1)
	   rr(i,j,k,2) = u(i,j,k,2)
	   rr(i,j,k,3) = u(i,j,k,3)
c	   rr(i,j,k,4) = phi(i,j,k)
           rr(i,j,k,4) = sedc(i,j,k)
	enddo
	enddo
	enddo 

C	R:
C	<u> (1) <v> (2) <w> (3) <r> (4)

	weit = 2.D0
	fact = 1.D0 / 64.D0 
	call filter(nni,nnj,nnk, rr(0,0,0,1),rr(0,0,0,5), weit,fact)
	call filter(nni,nnj,nnk, rr(0,0,0,2),rr(0,0,0,5), weit,fact)
	call filter(nni,nnj,nnk, rr(0,0,0,3),rr(0,0,0,5), weit,fact)
	call filter(nni,nnj,nnk, rr(0,0,0,4),rr(0,0,0,5), weit,fact)

	hc = 2.D0*(xp(2,1,1,1)-xp(1,1,1,1))
        C  = 0.8D0
   	do k = 0, nnk+1
        do j = 0, nnj+1
        do i = 0, nni+1
	   dz = xp(i,j,k,2)+depth(i,k)
	   if (dz .le. hc) then
	      absU = dsqrt( rr(i,j,k,1)**2.D0+rr(i,j,k,2)**2.D0
     <                     +rr(i,j,k,3)**2.D0 )
	      az = (dcos(pi*dz/(2.D0*hc)))**2.D0
	      btau(i,j,k,1) = -C*az*absU*rr(i,j,k,1)
	      btau(i,j,k,2) = -C*az*absU*rr(i,j,k,2)
              btau(i,j,k,3) = -C*az*absU*rr(i,j,k,3)
	enddo
	enddo
	enddo

	do m = 0, 3
	do k = 0, nnk+1
        do j = 0, nnj+1
        do i = 0, nni+1
	   if (btau(i,j,k,m) .ne. 0.D0) then
	      sum = 0.D0
	      do jj = 1, j
		 sum = sum + btau(i,j,k,m)*0.5D0*jac(i,jj,k)
     <                 *(xp(i,j+1,k,2)-xp(i,j-1,k,2))	   
              enddo
	      bsum(i,j,k,m) = sum
	   endif
	enddo
	enddo
	enddo
	enddo

C	T:
C	L11 (1) L22 (2) L33 (3)
C	L12 (4) L23 (5) L31 (6)
C	K1  (7) K2  (8) K3  (9)

        do k = 1, nnk
        do j = 1, nnj
        do i = 1, nni
	   tt(i,j,k,1) = tt(i,j,k,1) - rr(i,j,k,1)**2
	   tt(i,j,k,2) = tt(i,j,k,2) - bsum(i,j,k,2)
	   tt(i,j,k,3) = tt(i,j,k,3) - rr(i,j,k,3)**2
	   tt(i,j,k,4) = tt(i,j,k,4) - bsum(i,j,k,1)
	   tt(i,j,k,5) = tt(i,j,k,5) - bsum(i,j,k,3)
	   tt(i,j,k,6) = tt(i,j,k,6) - rr(i,j,k,3) * rr(i,j,k,1)
	   tt(i,j,k,7) = tt(i,j,k,7) - rr(i,j,k,4) * rr(i,j,k,1)
	   tt(i,j,k,8) = tt(i,j,k,8) - rr(i,j,k,4) * rr(i,j,k,2)
	   tt(i,j,k,9) = tt(i,j,k,9) - rr(i,j,k,4) * rr(i,j,k,3)
	enddo
	enddo
	enddo 

	call ssbc(rr(0,0,0,1))
	call ssbc(rr(0,0,0,2))
	call ssbc(rr(0,0,0,3))
	call ssbc(rr(0,0,0,4))
	call p_exchange( nni, nnj, nnk, rr(0,0,0,1) )   
	call p_exchange( nni, nnj, nnk, rr(0,0,0,2) )   
	call p_exchange( nni, nnj, nnk, rr(0,0,0,3) )   
	call p_exchange( nni, nnj, nnk, rr(0,0,0,4) )   

C	MijMij => R1

	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni
	   rr(i,j,k,1) = 
     <	              zz(i,j,k,1)**2 + zz(i,j,k,2)**2 + zz(i,j,k,3)**2
     <	   + 2.D0 * ( zz(i,j,k,4)**2 + zz(i,j,k,5)**2 + zz(i,j,k,6)**2 )
	enddo
	enddo
	enddo 

C	MijLij => R2

	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni
	   rr(i,j,k,2) = zz(i,j,k,1) * tt(i,j,k,1) 
     <	               + zz(i,j,k,2) * tt(i,j,k,2) 
     <	               + zz(i,j,k,3) * tt(i,j,k,3) 
     <	      + 2.D0 * ( zz(i,j,k,4) * tt(i,j,k,4) 
     <	               + zz(i,j,k,5) * tt(i,j,k,5) 
     <	               + zz(i,j,k,6) * tt(i,j,k,6) )
	enddo
	enddo
	enddo 

C	MijMij => Z1
C	MijLij => Z2
C	MijHij => Z3

	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni
	   zz(i,j,k,1) = rr(i,j,k,1)
	   zz(i,j,k,2) = rr(i,j,k,2)
c	   zz(i,j,k,3) = rr(i,j,k,3)
	enddo
	enddo
	enddo


C	               Z2       Z3
C	           - MijLij + MijHij
C	2 C D**2 = -------------------
C	                 MijMij
C	                   Z1

C	                  Z5     Z6
C	               - NiKi + NiJi
C	C D**2 / Prt = ---------------
C	                    NiNi
C	                     Z4

	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni
	   zz(i,j,k,1) = max(zz(i,j,k,1), tiny)
           zz(i,j,k,4) = max(zz(i,j,k,4), tiny)
	enddo
	enddo
	enddo

C	C D**2 => R1

	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni
	   rr(i,j,k,1) = zz(i,j,k,2)
	   rr(i,j,k,1) = 0.5D0 * rr(i,j,k,1) / zz(i,j,k,1)
	enddo
	enddo
	enddo

C	C/Prt D**2 => R2

	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni
	   rr(i,j,k,2) = - zz(i,j,k,5) + zz(i,j,k,6)
	   rr(i,j,k,2) = rr(i,j,k,2) / zz(i,j,k,4)
c           write(*,*) rr(i,j,k,2)
	enddo
	enddo
	enddo
	
C	1 / D**2 => R3

	temp = 2.D0 / 3.D0
	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni
	   rr(i,j,k,3) = jac(i,j,k) ** temp
	enddo
	enddo
	enddo

C	C => R1
C	C/PrT => R2

	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni
	   rr(i,j,k,1) = rr(i,j,k,1) * rr(i,j,k,3)
	   rr(i,j,k,2) = rr(i,j,k,2) * rr(i,j,k,3)
	enddo
	enddo
	enddo

	call ssbc(rr(0,0,0,1))
	call ssbc(rr(0,0,0,2))
	call p_exchange( nni, nnj, nnk, rr(0,0,0,1) )   
	call p_exchange( nni, nnj, nnk, rr(0,0,0,2) )   

C	{C} => R1
C	{C/PrT} => R2

	weit = 1.D0
	fact = 1.D0 / 27.D0
	call filter(nni,nnj,nnk, rr(0,0,0,1),rr(0,0,0,5), weit,fact)
	call filter(nni,nnj,nnk, rr(0,0,0,2),rr(0,0,0,5), weit,fact)

	if ( mod(istep, nsave) .eq. 0 ) then
	   call minmax(nni, nnj, nnk, rr(0,0,0,1), fmin, fmax)
	   if ( MYID .eq. 0 ) then
	      write(*,*) ' C    = ', fmin, fmax
	   endif
	   call minmax(nni, nnj, nnk, rr(0,0,0,2), fmin, fmax)
	   if ( MYID .eq. 0 ) then
	      write(*,*) ' C/Pr = ', fmin, fmax
	   endif   
	endif
	
C	D**2 => R3

	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni
	   rr(i,j,k,3) = 1.D0 / rr(i,j,k,3)
	enddo
	enddo
	enddo

C	muT = C D**2 |S|
C	kappaT = C/PrT D**2 |S|

	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni
	    vst(i,j,k) = rr(i,j,k,1) * rr(i,j,k,3) * sab(i,j,k)
	   akst(i,j,k) = rr(i,j,k,2) * rr(i,j,k,3) * sab(i,j,k)
	   akeq(i,j,k) = akst(i,j,k)
	enddo
	enddo
	enddo

	if ( mod(istep, nsave) .eq. 0 ) then
	   call minmax(nni, nnj, nnk, vst, fmin, fmax)
	   if ( MYID .eq. 0 ) then
	      write(*,*) ' vst = ', fmin, fmax
	   endif
	   call minmax(nni, nnj, nnk, akst, fmin, fmax)
	   if ( MYID .eq. 0 ) then
	      write(*,*) ' kst = ', fmin, fmax
	   endif
           call minmax(nni, nnj, nnk, aksd, fmin, fmax)
	   if ( MYID .eq. 0 ) then
	      write(*,*) ' kst_sedi = ', fmin, fmax
	   endif
	   
	endif

        do k = 1, nnk
        do j = 1, nnj
        do i = 1, nni
c	    vst(i,j,k) = max( vst(i,j,k), -vis)
c	   akst(i,j,k) = max(akst(i,j,k), -ak )
c	   aksd(i,j,k) = max(aksd(i,j,k), -ak )
            vst(i,j,k) = max( vst(i,j,k),  0.D0)
	   akst(i,j,k) = max(akst(i,j,k),  0.D0)
	   aksd(i,j,k) = max(aksd(i,j,k),  0.D0)
	enddo
	enddo
	enddo 

	call eddybc
	call p_exchange( nni, nnj, nnk,  vst )   
	call p_exchange( nni, nnj, nnk, akst )
	call p_exchange( nni, nnj, nnk, aksd )

	

	return
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine ssbc(s)

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "eddy.inc"

	double precision s(0:nni+1,0:nnj+1,0:nnk+1)
	integer i, j, k

	if ( n_west .eq. MPI_PROC_NULL ) then
	   do k = 0, nnk+1
	   do j = 0, nnj+1
	      s(0,j,k) = 2.D0 * s(1,j,k) - s(2,j,k)
	   enddo
	   enddo
	endif

	if ( n_east .eq. MPI_PROC_NULL ) then
	   do k = 0, nnk+1
	   do j = 0, nnj+1
	      s(nni+1,j,k) = 2.D0 * s(nni,j,k) - s(nni-1,j,k)
	   enddo
	   enddo
	endif

	if ( n_suth .eq. MPI_PROC_NULL ) then
	   do k = 0, nnk+1
	   do i = 0, nni+1
	      s(i,0,k) = 2.D0 * s(i,1,k) - s(i,2,k)
	   enddo
	   enddo
	endif

	if ( n_nrth .eq. MPI_PROC_NULL ) then
	   do k = 0, nnk+1
	   do i = 0, nni+1
	      s(i,nnj+1,k) = 2.D0 * s(i,nnj,k) - s(i,nnj-1,k)
	   enddo
	   enddo
	endif

	if ( n_back .eq. MPI_PROC_NULL ) then
	   do j = 0, nnj+1
	   do i = 0, nni+1
	      s(i,j,0) = 2.D0 * s(i,j,1) - s(i,j,2)
	   enddo
	   enddo
	endif

	if ( n_frnt .eq. MPI_PROC_NULL ) then
	   do j = 0, nnj+1
	   do i = 0, nni+1
	      s(i,j,nnk+1) = 2.D0 * s(i,j,nnk) - s(i,j,nnk-1)
	   enddo
	   enddo
	endif

	return
	end
	
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine eddybc

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "eddy.inc"

	integer i, j, k

	if ( n_west .eq. MPI_PROC_NULL ) then
	   do k = 0, nnk+1
	   do j = 0, nnj+1
	       vst(0,j,k) =   vst(1,j,k)
	      akst(0,j,k) =  akst(1,j,k)
	      aksd(0,j,k) =  aksd(1,j,k)
	   enddo
	   enddo
	endif

	if ( n_east .eq. MPI_PROC_NULL ) then
	   do k = 0, nnk+1
	   do j = 0, nnj+1
	       vst(nni+1,j,k) =   vst(nni,j,k)
	      akst(nni+1,j,k) =  akst(nni,j,k)
	      aksd(nni+1,j,k) =  aksd(nni,j,k)
	   enddo
	   enddo
	endif

	if ( n_suth .eq. MPI_PROC_NULL ) then
	   do k = 0, nnk+1
	   do i = 0, nni+1
	       vst(i,0,k) =   vst(i,1,k)
	      akst(i,0,k) =  akst(i,1,k)
	      aksd(i,0,k) =  aksd(i,1,k)
	   enddo
	   enddo
	endif

	if ( n_nrth .eq. MPI_PROC_NULL ) then
	   do k = 0, nnk+1
	   do i = 0, nni+1
c	       vst(i,nnj+1,k) = - vst(i,nnj,k)
c	      akst(i,nnj+1,k) = - akst(i,nnj,k)
c	      aksd(i,nnj+1,k) = - aksd(i,nnj,k)
               vst(i,nnj+1,k) = vst(i,nnj,k)
	      akst(i,nnj+1,k) = akst(i,nnj,k)
	      aksd(i,nnj+1,k) = aksd(i,nnj,k)
	   enddo
	   enddo
	endif

	if ( n_back .eq. MPI_PROC_NULL ) then
	   do j = 0, nnj+1
	   do i = 0, nni+1
c	       vst(i,j,0) = -  vst(i,j,1)
c	      akst(i,j,0) = - akst(i,j,1)
c	      aksd(i,j,0) = - aksd(i,j,1)
               vst(i,j,0) = vst(i,j,1)
	      akst(i,j,0) = akst(i,j,1)
	      aksd(i,j,0) = aksd(i,j,1)
	   enddo
	   enddo
	endif

	if ( n_frnt .eq. MPI_PROC_NULL ) then
	   do j = 0, nnj+1
	   do i = 0, nni+1
c	       vst(i,j,nnk+1) = -  vst(i,j,nnk)
c	      akst(i,j,nnk+1) = - akst(i,j,nnk)
c	      aksd(i,j,nnk+1) = - aksd(i,j,nnk)
               vst(i,j,nnk+1) =   vst(i,j,nnk)
	      akst(i,j,nnk+1) =  akst(i,j,nnk)
	      aksd(i,j,nnk+1) =  aksd(i,j,nnk)

	   enddo
	   enddo
	endif

	return
	end
	
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine filter(ii, jj, kk, f, r, weit, fact)

	implicit none

	integer ii, jj, kk

	double precision f(0:ii+1,0:jj+1,0:kk+1),
     &                   r(0:ii+1,0:jj+1,0:kk+1)

	double precision weit, fact

	integer i, j, k

        do k = 0, kk+1
        do j = 0, jj+1
        do i = 1, ii
	   r(i,j,k) = f(i-1,j,k) + weit * f(i,j,k) + f(i+1,j,k)
	enddo
	enddo
	enddo 

        do k = 0, kk+1
        do j = 1, jj
        do i = 1, ii
	   f(i,j,k) = r(i,j-1,k) + weit * r(i,j,k) + r(i,j+1,k)
	enddo
	enddo
	enddo 

        do k = 1, kk
        do j = 1, jj
        do i = 1, ii
	   r(i,j,k) = f(i,j,k-1) + weit * f(i,j,k) + f(i,j,k+1)
	enddo
	enddo
	enddo 

        do k = 1, kk
        do j = 1, jj
        do i = 1, ii
	   f(i,j,k) = fact * r(i,j,k)
	enddo
	enddo
	enddo 

	return
	end
	
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine minmax(ii, jj, kk, f, fmin, fmax)

	implicit none
	include "mpif.h"
	include "mpi.inc"

	integer ii, jj, kk

	double precision f(0:ii+1,0:jj+1,0:kk+1)

	double precision fmin, fmax, ffmin, ffmax

	integer i, j, k

C f77 does not recognize HUGE, so it is set to 10000.D0
c	fmin =    HUGE(1.D0)
c	fmax =   -HUGE(1.D0)
c	ffmin =   HUGE(1.D0)
c	ffmax =  -HUGE(1.D0)
	fmin =    10000.D0
	fmax =   -10000.D0
	ffmin =   10000.D0
	ffmax =  -10000.D0

	do k = 1, kk
	do j = 1, jj
	do i = 1, ii
	   ffmin = min(ffmin, f(i,j,k))
	   ffmax = max(ffmax, f(i,j,k))
	enddo
	enddo
	enddo

	call MPI_REDUCE(ffmin,fmin,1,MPI_DOUBLE_PRECISION,
     <	                MPI_MIN,0,comm3d,ierr)
	call MPI_REDUCE(ffmax,fmax,1,MPI_DOUBLE_PRECISION,
     <	                MPI_MAX,0,comm3d,ierr)

	return
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine parabolicVis

	include "size.inc"
	include "ns.inc"
	include "metric.inc"
	include "para.inc"
	include "eddy.inc"
	include "sedi.inc"

	integer i, j, k
        double precision n0, kappa, yy

        n0 = diam
        kappa = 0.41D0
        
        do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni

           yy = depth(i,k) + xp(i,j,k,2)
           aksd(i,j,k) = kappa*dsqrt(Pgrd*depth(i,k))*(yy+n0)
     <                  * (1.D0-yy/(depth(i,k)))
c           if (sedi .eq. 1) then
c              aksd(i,j,k) = vst(i,j,k)     
c           endif
        enddo
	enddo
	enddo
         
        return
        end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine NewDiffusivity
       
        include "size.inc"
        include "ns.inc"
        include "metric.inc"
        include "para.inc"
        include "eddy.inc"
	include "sedi.inc"  
          
        integer i, j, k
	double precision h, temp
          

	do k = 1, nnk
        do j = 1, nnj
        do i = 1, nni
        temp = xp(i,j+1,k,2) - xp(i,j-1,k,2)
        temp = 1.D0/temp
        temp = (u(i,j+1,k,1) - u(i,j-1,k,1))*temp
        temp = 1.D0/temp
        aksd(i,j,k) = Pgrd*(-xp(i,j,k,2)+diam)*temp
        
        enddo
        enddo
        enddo	


        return
        end
