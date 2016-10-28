cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine eddy

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "ns.inc"
	include "metric.inc"
	include "para.inc"
	include "eddy.inc"
	include "sedi.inc"

	double precision ss(0:nni+1,0:nnj+1,0:nnk+1,12), 
     &                   tt(0:nni+1,0:nnj+1,0:nnk+1,12), 
     &                   zz(0:nni+1,0:nnj+1,0:nnk+1,12)
     
	double precision kk(0:nni+1,1:cut,0:nnk+1)

	integer i, j, k, m

	double precision tiny, alpha2, weit, fact, temp
	double precision fmin, fmax, tmpzz

	tiny = 1.D-20

	alpha2 = 4.D0

        do m = 1, 12
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

c        do m = 1, 6
c        do k = 0, nnk+1
c        do j = 0, nnj+1
c        do i = 0, nni+1
c	   rr(i,j,k,m) = 0.D0
c	enddo
c	enddo
c	enddo 
c	enddo 

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
           rr(i,j,k,1) = sedc(i+1,j,k) - sedc(i-1,j,k)
	   rr(i,j,k,2) = sedc(i,j+1,k) - sedc(i,j-1,k)
	   rr(i,j,k,3) = sedc(i,j,k+1) - sedc(i,j,k-1)
	   rr(i,j,k,4) = phi(i+1,j,k) - phi(i-1,j,k)
	   rr(i,j,k,5) = phi(i,j+1,k) - phi(i,j-1,k)
	   rr(i,j,k,6) = phi(i,j,k+1) - phi(i,j,k-1)
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
	   ss(i,j,k,10) = tt(i,j,k,1) * rr(i,j,k,4)
     <	                + tt(i,j,k,2) * rr(i,j,k,5)
     <	                + tt(i,j,k,3) * rr(i,j,k,6)
	   ss(i,j,k,11) = tt(i,j,k,4) * rr(i,j,k,4)
     <	                + tt(i,j,k,5) * rr(i,j,k,5)
     <	                + tt(i,j,k,6) * rr(i,j,k,6)
	   ss(i,j,k,12) = tt(i,j,k,7) * rr(i,j,k,4)
     <	                + tt(i,j,k,8) * rr(i,j,k,5)
     <	                + tt(i,j,k,9) * rr(i,j,k,6)
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
	   zz(i,j,k,10) = sab(i,j,k) * ss(i,j,k,10)
	   zz(i,j,k,11) = sab(i,j,k) * ss(i,j,k,11)
	   zz(i,j,k,12) = sab(i,j,k) * ss(i,j,k,12)
	enddo
	enddo
	enddo 

C	Z:
C	<|S|S11> (1) <|S|S22> (2) <|S|S33> (3)
C	<|S|S12> (4) <|S|S23> (5) <|S|S31> (6)
C	<|S|rsx> (7) <|S|rsy> (8) <|S|rsz> (9)

	weit = 2.D0
	fact = 1.D0 / 64.D0 
	call filter(nni,nnj,nnk, zz(0,0,0,1),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, zz(0,0,0,2),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, zz(0,0,0,3),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, zz(0,0,0,4),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, zz(0,0,0,5),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, zz(0,0,0,6),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, zz(0,0,0,7),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, zz(0,0,0,8),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, zz(0,0,0,9),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, zz(0,0,0,10),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, zz(0,0,0,11),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, zz(0,0,0,12),rr(0,0,0,8), weit,fact)

C	R:
C	|S| (4)

        do k = 0, nnk+1
        do j = 0, nnj+1
        do i = 0, nni+1
	   rr(i,j,k,7) = sab(i,j,k)
	enddo
	enddo
	enddo 

C	R:
C	<|S|> (4)

	weit = 2.D0
	fact = 1.D0 / 64.D0 
	call filter(nni,nnj,nnk, rr(0,0,0,7),rr(0,0,0,8), weit,fact)

C	S:
C	<S11> (1) <S22> (2) <S33> (3)
C	<S12> (1) <S23> (2) <S32> (3)
C	<rsx> (1) <rsy> (2) <rsz> (3)

	weit = 2.D0
	fact = 1.D0 / 64.D0 
	call filter(nni,nnj,nnk, ss(0,0,0,1),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, ss(0,0,0,2),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, ss(0,0,0,3),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, ss(0,0,0,4),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, ss(0,0,0,5),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, ss(0,0,0,6),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, ss(0,0,0,7),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, ss(0,0,0,8),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, ss(0,0,0,9),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, ss(0,0,0,10),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, ss(0,0,0,11),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, ss(0,0,0,12),rr(0,0,0,8), weit,fact)


C	Z:
C	M11 (1) M22 (2) M33 (3)
C	M12 (4) M23 (5) M31 (6)
C	N1  (7) N2  (8) N3  (9)

        do k = 1, nnk
        do j = 1, nnj
        do i = 1, nni
	   zz(i,j,k,1) = alpha2 * rr(i,j,k,7) * ss(i,j,k,1) - zz(i,j,k,1) 
	   zz(i,j,k,2) = alpha2 * rr(i,j,k,7) * ss(i,j,k,2) - zz(i,j,k,2) 
	   zz(i,j,k,3) = alpha2 * rr(i,j,k,7) * ss(i,j,k,3) - zz(i,j,k,3) 
	   zz(i,j,k,4) = alpha2 * rr(i,j,k,7) * ss(i,j,k,4) - zz(i,j,k,4) 
	   zz(i,j,k,5) = alpha2 * rr(i,j,k,7) * ss(i,j,k,5) - zz(i,j,k,5) 
	   zz(i,j,k,6) = alpha2 * rr(i,j,k,7) * ss(i,j,k,6) - zz(i,j,k,6) 
	   zz(i,j,k,7) = alpha2 * rr(i,j,k,7) * ss(i,j,k,7) - zz(i,j,k,7) 
	   zz(i,j,k,8) = alpha2 * rr(i,j,k,7) * ss(i,j,k,8) - zz(i,j,k,8)  
	   zz(i,j,k,9) = alpha2 * rr(i,j,k,7) * ss(i,j,k,9) - zz(i,j,k,9)
	   zz(i,j,k,10) = alpha2 * rr(i,j,k,7)*ss(i,j,k,10)-zz(i,j,k,10) 
	   zz(i,j,k,11) = alpha2 * rr(i,j,k,7)*ss(i,j,k,11)-zz(i,j,k,11)  
	   zz(i,j,k,12) = alpha2 * rr(i,j,k,7)*ss(i,j,k,12)-zz(i,j,k,12)
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
	   tt(i,j,k,2) = u(i,j,k,2)**2
	   tt(i,j,k,3) = u(i,j,k,3)**2
	   tt(i,j,k,4) = u(i,j,k,1) * u(i,j,k,2)
	   tt(i,j,k,5) = u(i,j,k,2) * u(i,j,k,3)
	   tt(i,j,k,6) = u(i,j,k,3) * u(i,j,k,1)
	   tt(i,j,k,7) = sedc(i,j,k) * u(i,j,k,1)
	   tt(i,j,k,8) = sedc(i,j,k) * u(i,j,k,2)
	   tt(i,j,k,9) = sedc(i,j,k) * u(i,j,k,3)
	   tt(i,j,k,10) = phi(i,j,k) * u(i,j,k,1)
	   tt(i,j,k,11) = phi(i,j,k) * u(i,j,k,2)
	   tt(i,j,k,12) = phi(i,j,k) * u(i,j,k,3)
	enddo
	enddo
	enddo 

C	T:
C	<uu> (1) <vv> (2) <ww> (3)
C	<uv> (4) <vw> (5) <wu> (6)
C	<ru> (7) <rv> (8) <rw> (9)

	weit = 2.D0
	fact = 1.D0 / 64.D0 
	call filter(nni,nnj,nnk, tt(0,0,0,1),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, tt(0,0,0,2),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, tt(0,0,0,3),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, tt(0,0,0,4),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, tt(0,0,0,5),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, tt(0,0,0,6),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, tt(0,0,0,7),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, tt(0,0,0,8),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, tt(0,0,0,9),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, tt(0,0,0,10),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, tt(0,0,0,11),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, tt(0,0,0,12),rr(0,0,0,8), weit,fact)

C	R:
C	u (1) v (2) w (3) r (4)

        do k = 0, nnk+1
        do j = 0, nnj+1
        do i = 0, nni+1
	   rr(i,j,k,1) = u(i,j,k,1)
	   rr(i,j,k,2) = u(i,j,k,2)
	   rr(i,j,k,3) = u(i,j,k,3)
           rr(i,j,k,4) = sedc(i,j,k)
	   rr(i,j,k,5) = phi(i,j,k)
	enddo
	enddo
	enddo 

C	R:
C	<u> (1) <v> (2) <w> (3) <r> (4)

	weit = 2.D0
	fact = 1.D0 / 64.D0 
	call filter(nni,nnj,nnk, rr(0,0,0,1),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, rr(0,0,0,2),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, rr(0,0,0,3),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, rr(0,0,0,4),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, rr(0,0,0,5),rr(0,0,0,8), weit,fact)

C	T:
C	L11 (1) L22 (2) L33 (3)
C	L12 (4) L23 (5) L31 (6)
C	K1  (7) K2  (8) K3  (9)

        do k = 1, nnk
        do j = 1, nnj
        do i = 1, nni
	   tt(i,j,k,1) = tt(i,j,k,1) - rr(i,j,k,1)**2
	   tt(i,j,k,2) = tt(i,j,k,2) - rr(i,j,k,2)**2
	   tt(i,j,k,3) = tt(i,j,k,3) - rr(i,j,k,3)**2
	   tt(i,j,k,4) = tt(i,j,k,4) - rr(i,j,k,1) * rr(i,j,k,2)
	   tt(i,j,k,5) = tt(i,j,k,5) - rr(i,j,k,2) * rr(i,j,k,3)
	   tt(i,j,k,6) = tt(i,j,k,6) - rr(i,j,k,3) * rr(i,j,k,1)
	   tt(i,j,k,7) = tt(i,j,k,7) - rr(i,j,k,4) * rr(i,j,k,1)
	   tt(i,j,k,8) = tt(i,j,k,8) - rr(i,j,k,4) * rr(i,j,k,2)
	   tt(i,j,k,9) = tt(i,j,k,9) - rr(i,j,k,4) * rr(i,j,k,3)
	   tt(i,j,k,10) = tt(i,j,k,10) - rr(i,j,k,5) * rr(i,j,k,1)
	   tt(i,j,k,11) = tt(i,j,k,11) - rr(i,j,k,5) * rr(i,j,k,2)
	   tt(i,j,k,12) = tt(i,j,k,12) - rr(i,j,k,5) * rr(i,j,k,3)
	enddo
	enddo
	enddo 

	if (cut .gt. 0) then
c       MijMij => R1

	   do k = 1, nnk
	   do j = 1, cut
	   do i = 1, nni
	      rr(i,j,k,1) = 
     <	              zz(i,j,k,1)**2 + zz(i,j,k,2)**2 + zz(i,j,k,3)**2
     <	    + 2.D0 * ( zz(i,j,k,4)**2 + zz(i,j,k,5)**2 + zz(i,j,k,6)**2)
	      rr(i,j,k,1) = max(rr(i,j,k,1), tiny)
	   enddo
	   enddo
	   enddo

C	LijMij => R2

	   do k = 1, nnk
	   do j = 1, cut
	   do i = 1, nni
	      rr(i,j,k,2) = tt(i,j,k,1) * zz(i,j,k,1) 
     <	               + tt(i,j,k,2) * zz(i,j,k,2) 
     <	               + tt(i,j,k,3) * zz(i,j,k,3) 
     <	      + 2.D0 * ( tt(i,j,k,4) * zz(i,j,k,4) 
     <	               + tt(i,j,k,5) * zz(i,j,k,5) 
     <	               + tt(i,j,k,6) * zz(i,j,k,6) )
	   enddo
	   enddo
	   enddo 

C	NiNi => R3

	   do k = 1, nnk
	   do j = 1, cut
	   do i = 1, nni
	      rr(i,j,k,3) = zz(i,j,k,7)**2
     <         	   + zz(i,j,k,8)**2
     <                 + zz(i,j,k,9)**2
           rr(i,j,k,3) = max(rr(i,j,k,3), tiny)
	   enddo
	   enddo
	   enddo 

C	KiNi => R4

	   do k = 1, nnk
	   do j = 1, cut
	   do i = 1, nni
              rr(i,j,k,4) = tt(i,j,k,7) * zz(i,j,k,7)
     <                 + tt(i,j,k,8) * zz(i,j,k,8)
     <                 + tt(i,j,k,9) * zz(i,j,k,9)
           enddo
	   enddo
	   enddo
	


C	NiNi => R5 for scalar

	   do k = 1, nnk
	   do j = 1, cut
	   do i = 1, nni
	      rr(i,j,k,5) = zz(i,j,k,10)**2
     <                 + zz(i,j,k,11)**2
     <                 + zz(i,j,k,12)**2
              rr(i,j,k,5) = max(rr(i,j,k,5), tiny)
           enddo
	   enddo
	   enddo 

C	KiNi => R6 for scalar

	   do k = 1, nnk
	   do j = 1, cut
	   do i = 1, nni
	      rr(i,j,k,6) = tt(i,j,k,10) * zz(i,j,k,10)
     <                 + tt(i,j,k,11) * zz(i,j,k,11)
     <                 + tt(i,j,k,12) * zz(i,j,k,12)
           enddo
	   enddo
	   enddo 

cccccccccccccccccccccccccccccccccccccccccccccccccccccc

c       CD**2 = - 0.5 LijMij / MijMij
c         R2          R2        R1    

	do k = 1, nnk
	do j = 1, cut
	do i = 1, nni
           rr(i,j,k,2) = - 0.5D0 * rr(i,j,k,2) / rr(i,j,k,1)
	enddo
	enddo
	enddo 

c       C/PrT D**2 = -KiNi / NiNi
c            R4        R4     R3

	do k = 1, nnk
	do j = 1, cut
	do i = 1, nni
           rr(i,j,k,4) = - rr(i,j,k,4) / rr(i,j,k,3)
	enddo
	enddo
	enddo 

c       C/PrT D**2 = -KiNi / NiNi
c            R4        R4     R3

	do k = 1, nnk
	do j = 1, cut
	do i = 1, nni
           rr(i,j,k,4) = - rr(i,j,k,4) / rr(i,j,k,3)
	enddo
	enddo
	enddo 

c       C/PrT D**2 = -KiNi / NiNi     for scalar
c            R6        R6     R5 

	do k = 1, nnk
	do j = 1, cut
	do i = 1, nni
           rr(i,j,k,6) = - rr(i,j,k,6) / rr(i,j,k,5)
	enddo
	enddo
	enddo 

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc


c       1 / D**2 => R1

	temp = 2.D0 / 3.D0
	do k = 1, nnk
	do j = 1, cut
	do i = 1, nni
           rr(i,j,k,1) = jac(i,j,k) ** temp
	enddo
	enddo
	enddo 

c       C/PrT => R4

	do k = 1, nnk
	do j = 1, cut
	do i = 1, nni
           kk(i,j,k) = rr(i,j,k,4) * rr(i,j,k,1)
	enddo
	enddo
	enddo 

	endif


C	R:
C	u (1) v (2) w (3) r (4)

        do k = 0, nnk+1
        do j = 0, nnj+1
        do i = 0, nni+1
	   rr(i,j,k,1) = u(i,j,k,1)
	   rr(i,j,k,2) = u(i,j,k,2)
	   rr(i,j,k,3) = u(i,j,k,3)
	   rr(i,j,k,4) = sedc(i,j,k)
	   rr(i,j,k,5) = phi(i,j,k)
	enddo
	enddo
	enddo 

C	R:
C	[u] (1) [v] (2) [w] (3) [r] (4)

	weit = 6.D0
	fact = 1.D0 / 512.D0 
	call filter(nni,nnj,nnk, rr(0,0,0,1),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, rr(0,0,0,2),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, rr(0,0,0,3),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, rr(0,0,0,4),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, rr(0,0,0,5),rr(0,0,0,8), weit,fact)
	call ssbc(rr(0,0,0,1))
	call ssbc(rr(0,0,0,2))
	call ssbc(rr(0,0,0,3))
	call ssbc(rr(0,0,0,4))
	call ssbc(rr(0,0,0,5))
	call p_exchange( nni, nnj, nnk, rr(0,0,0,1) )   
	call p_exchange( nni, nnj, nnk, rr(0,0,0,2) )   
	call p_exchange( nni, nnj, nnk, rr(0,0,0,3) )   
	call p_exchange( nni, nnj, nnk, rr(0,0,0,4) )   
	call p_exchange( nni, nnj, nnk, rr(0,0,0,5) )   

C	S:
C	[u][u] (1) [v][v] (2) [w][w] (3)
C	[u][v] (4) [v][w] (5) [w][u] (6)
C	[r][u] (7) [r][v] (8) [r][w] (9)

        do k = 0, nnk+1
        do j = 0, nnj+1
        do i = 0, nni+1
	   ss(i,j,k,1) = rr(i,j,k,1)**2
	   ss(i,j,k,2) = rr(i,j,k,2)**2
	   ss(i,j,k,3) = rr(i,j,k,3)**2
	   ss(i,j,k,4) = rr(i,j,k,1) * rr(i,j,k,2)
	   ss(i,j,k,5) = rr(i,j,k,2) * rr(i,j,k,3)
	   ss(i,j,k,6) = rr(i,j,k,3) * rr(i,j,k,1)
	   ss(i,j,k,7) = rr(i,j,k,4) * rr(i,j,k,1)
	   ss(i,j,k,8) = rr(i,j,k,4) * rr(i,j,k,2)
	   ss(i,j,k,9) = rr(i,j,k,4) * rr(i,j,k,3)
	   ss(i,j,k,10) = rr(i,j,k,5) * rr(i,j,k,1)
	   ss(i,j,k,11) = rr(i,j,k,5) * rr(i,j,k,2)
	   ss(i,j,k,12) = rr(i,j,k,5) * rr(i,j,k,3)
	enddo
	enddo
	enddo 

C	S:
C	<[u][u]> (1) <[v][v]> (2) <[w][w]> (3)
C	<[u][v]> (4) <[v][w]> (5) <[w][u]> (6)
C	<[r][u]> (7) <[r][v]> (8) <[r][w]> (9)

	weit = 2.D0
	fact = 1.D0 / 64.D0 
	call filter(nni,nnj,nnk, ss(0,0,0,1),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, ss(0,0,0,2),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, ss(0,0,0,3),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, ss(0,0,0,4),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, ss(0,0,0,5),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, ss(0,0,0,6),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, ss(0,0,0,7),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, ss(0,0,0,8),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, ss(0,0,0,9),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, ss(0,0,0,10),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, ss(0,0,0,11),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, ss(0,0,0,12),rr(0,0,0,8), weit,fact)


C	R:
C	<[u]> (1) <[v]> (2) <[w]> (3) <[r]> (4)

	weit = 2.D0
	fact = 1.D0 / 64.D0 
	call filter(nni,nnj,nnk, rr(0,0,0,1),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, rr(0,0,0,2),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, rr(0,0,0,3),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, rr(0,0,0,4),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, rr(0,0,0,5),rr(0,0,0,8), weit,fact)

C	S:
C	H11 (1) H22 (2) H33 (3)
C	H12 (4) H23 (5) H31 (6)
C	J1  (7) J2  (8) J3  (9)

        do k = 1, nnk
        do j = 1, nnj
        do i = 1, nni
	   ss(i,j,k,1) = ss(i,j,k,1) - rr(i,j,k,1)**2
	   ss(i,j,k,2) = ss(i,j,k,2) - rr(i,j,k,2)**2
	   ss(i,j,k,3) = ss(i,j,k,3) - rr(i,j,k,3)**2
	   ss(i,j,k,4) = ss(i,j,k,4) - rr(i,j,k,1) * rr(i,j,k,2)
	   ss(i,j,k,5) = ss(i,j,k,5) - rr(i,j,k,2) * rr(i,j,k,3)
	   ss(i,j,k,6) = ss(i,j,k,6) - rr(i,j,k,3) * rr(i,j,k,1)
	   ss(i,j,k,7) = ss(i,j,k,7) - rr(i,j,k,4) * rr(i,j,k,1)
	   ss(i,j,k,8) = ss(i,j,k,8) - rr(i,j,k,4) * rr(i,j,k,2)
	   ss(i,j,k,9) = ss(i,j,k,9) - rr(i,j,k,4) * rr(i,j,k,3)
	   ss(i,j,k,10) = ss(i,j,k,10) - rr(i,j,k,5) * rr(i,j,k,1)
	   ss(i,j,k,11) = ss(i,j,k,11) - rr(i,j,k,5) * rr(i,j,k,2)
	   ss(i,j,k,12) = ss(i,j,k,12) - rr(i,j,k,5) * rr(i,j,k,3)
	enddo
	enddo
	enddo 

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

C	MijHij => R3

	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni
	   rr(i,j,k,3) = zz(i,j,k,1) * ss(i,j,k,1) 
     <	               + zz(i,j,k,2) * ss(i,j,k,2) 
     <	               + zz(i,j,k,3) * ss(i,j,k,3) 
     <	      + 2.D0 * ( zz(i,j,k,4) * ss(i,j,k,4) 
     <	               + zz(i,j,k,5) * ss(i,j,k,5) 
     <	               + zz(i,j,k,6) * ss(i,j,k,6) )
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
	   zz(i,j,k,3) = rr(i,j,k,3)
	enddo
	enddo
	enddo

C	NiNi => R1

	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni
	   rr(i,j,k,1) = zz(i,j,k,7)**2
     <	               + zz(i,j,k,8)**2
     <	               + zz(i,j,k,9)**2
	enddo
	enddo
	enddo 

C	NiKi => R2

	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni
	   rr(i,j,k,2) = zz(i,j,k,7) * tt(i,j,k,7) 
     <	               + zz(i,j,k,8) * tt(i,j,k,8) 
     <	               + zz(i,j,k,9) * tt(i,j,k,9) 
	enddo
	enddo
	enddo 

C	NiJi => R3

	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni
	   rr(i,j,k,3) = zz(i,j,k,7) * ss(i,j,k,7) 
     <	               + zz(i,j,k,8) * ss(i,j,k,8) 
     <	               + zz(i,j,k,9) * ss(i,j,k,9) 
	enddo
	enddo
	enddo 

C	NiNi => Z4
C	NiKi => Z5
C	NiJi => Z6

	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni
	   zz(i,j,k,4) = rr(i,j,k,1)
	   zz(i,j,k,5) = rr(i,j,k,2)
	   zz(i,j,k,6) = rr(i,j,k,3)
	enddo
	enddo
	enddo 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccc Now, for scalar cccccccccccccccccccccccccccccccc

C	NiNi => R1

	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni
	   rr(i,j,k,1) = zz(i,j,k,10)**2
     <	               + zz(i,j,k,11)**2
     <	               + zz(i,j,k,12)**2
	enddo
	enddo
	enddo 

C	NiKi => R2

	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni
	   rr(i,j,k,2) = zz(i,j,k,10) * tt(i,j,k,10) 
     <	               + zz(i,j,k,11) * tt(i,j,k,11) 
     <	               + zz(i,j,k,12) * tt(i,j,k,12) 
	enddo
	enddo
	enddo 

C	NiJi => R3

	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni
	   rr(i,j,k,3) = zz(i,j,k,10) * ss(i,j,k,10) 
     <	               + zz(i,j,k,11) * ss(i,j,k,11) 
     <	               + zz(i,j,k,12) * ss(i,j,k,12) 
	enddo
	enddo
	enddo 

C	NiNi => Z7
C	NiKi => Z8
C	NiJi => Z9

	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni
	   zz(i,j,k,7) = rr(i,j,k,1)
	   zz(i,j,k,8) = rr(i,j,k,2)
	   zz(i,j,k,9) = rr(i,j,k,3)
	enddo
	enddo
	enddo 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc



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
           zz(i,j,k,7) = max(zz(i,j,k,7), tiny)
	enddo
	enddo
	enddo

C	C D**2 => R1

	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni
	   rr(i,j,k,1) = - zz(i,j,k,2) + zz(i,j,k,3)
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
	enddo
	enddo
	enddo

C	C/Prt D**2 => R3 for scalar

	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni
	   rr(i,j,k,3) = - zz(i,j,k,8) + zz(i,j,k,9)
	   rr(i,j,k,3) = rr(i,j,k,3) / zz(i,j,k,7)
	enddo
	enddo
	enddo

	
C	1 / D**2 => R4

	temp = 2.D0 / 3.D0
	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni
	   rr(i,j,k,4) = jac(i,j,k) ** temp
	enddo
	enddo
	enddo




C	C => R1
C	C/PrT => R2

	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni
	   rr(i,j,k,1) = rr(i,j,k,1) * rr(i,j,k,4)
c	   rr(i,j,k,2) = rr(i,j,k,2) * rr(i,j,k,4)
	   if (j .le. cut) then
	      rr(i,j,k,2) = kk(i,j,k)
	   else
	      rr(i,j,k,2) = rr(i,j,k,2) * rr(i,j,k,4)
	   endif
	   rr(i,j,k,3) = rr(i,j,k,3) * rr(i,j,k,4)
	enddo
	enddo
	enddo

	call ssbc(rr(0,0,0,1))
	call ssbc(rr(0,0,0,2))
	call ssbc(rr(0,0,0,3))
	call p_exchange( nni, nnj, nnk, rr(0,0,0,1) )   
	call p_exchange( nni, nnj, nnk, rr(0,0,0,2) )   
	call p_exchange( nni, nnj, nnk, rr(0,0,0,3) )   

C	{C} => R1
C	{C/PrT} => R2
C	{C/PrT} => R3  for scalar

	weit = 1.D0
	fact = 1.D0 / 27.D0
	call filter(nni,nnj,nnk, rr(0,0,0,1),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, rr(0,0,0,2),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, rr(0,0,0,3),rr(0,0,0,8), weit,fact)

	if ( mod(istep, nsave) .eq. 0 ) then
	   call minmax(nni, nnj, nnk, rr(0,0,0,1), fmin, fmax)
	   if ( MYID .eq. 0 ) then
	      write(*,*) ' C    = ', fmin, fmax
	   endif

	   call minmax(nni, nnj, nnk, rr(0,0,0,2), fmin, fmax)
	   if ( MYID .eq. 0 ) then
	      write(*,*) ' C/Pr (sedi)', fmin, fmax
	   endif   

	   call minmax(nni, nnj, nnk, rr(0,0,0,3), fmin, fmax)
	   if ( MYID .eq. 0 ) then
	      write(*,*) ' C/Pr (phi)', fmin, fmax
	   endif   
	endif
	
C	D**2 => R3

	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni
	   rr(i,j,k,4) = 1.D0 / rr(i,j,k,4)
	enddo
	enddo
	enddo

C	muT = C D**2 |S|
C	kappaT = C/PrT D**2 |S|

	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni
	   vst(i,j,k)  = rr(i,j,k,1) * rr(i,j,k,4) * sab(i,j,k)
	   aksd(i,j,k) = rr(i,j,k,2) * rr(i,j,k,4) * sab(i,j,k) 
	   akst(i,j,k) = rr(i,j,k,3) * rr(i,j,k,4) * sab(i,j,k)

	enddo
	enddo
	enddo

	if ( mod(istep, nsave) .eq. 0 ) then
	   call minmax(nni, nnj, nnk, vst, fmin, fmax)
	   if ( MYID .eq. 0 ) then
	      write(*,*) ' vst = ', fmin, fmax
	   endif
	   call minmax(nni, nnj, nnk, aksd, fmin, fmax)
	   if ( MYID .eq. 0 ) then
	      write(*,*) ' sed diff = ', fmin, fmax
	   endif
	   call minmax(nni, nnj, nnk, akst, fmin, fmax)
	   if ( MYID .eq. 0 ) then
	      write(*,*) ' phi diff = ', fmin, fmax
	   endif
           	   
	endif

	call minmax(nni, nnj, nnk, vst, fmin, fmax)
	if ( MYID .eq. 0 ) then
	   write(*,*) ' vst = ', fmin, fmax
	endif
	call minmax(nni, nnj, nnk, aksd, fmin, fmax)
	if ( MYID .eq. 0 ) then
	   write(*,*) ' sed diff = ', fmin, fmax
	endif
	call minmax(nni, nnj, nnk, akst, fmin, fmax)
	if ( MYID .eq. 0 ) then
	   write(*,*) ' phi diff = ', fmin, fmax
	endif

        do k = 1, nnk
        do j = 1, nnj
        do i = 1, nni
	    vst(i,j,k) = max( vst(i,j,k), -vis)
	   akst(i,j,k) = max(akst(i,j,k), -ak )
	   aksd(i,j,k) = max(aksd(i,j,k), -ak )          
	enddo
	enddo
	enddo 

	call eddybc
	call p_exchange( nni, nnj, nnk,  vst )   
	call p_exchange( nni, nnj, nnk, akst )
	call p_exchange( nni, nnj, nnk, aksd )

	
C	T:
C	uu (1) vv (2) ww (3)
C	uv (4) vw (5) wu (6)
C	ru (7) rv (8) rw (9)

        do k = 0, nnk+1
        do j = 0, nnj+1
        do i = 0, nni+1
	   tt(i,j,k,1) = u(i,j,k,1)**2
	   tt(i,j,k,2) = u(i,j,k,2)**2
	   tt(i,j,k,3) = u(i,j,k,3)**2
	   tt(i,j,k,4) = u(i,j,k,1) * u(i,j,k,2)
	   tt(i,j,k,5) = u(i,j,k,2) * u(i,j,k,3)
	   tt(i,j,k,6) = u(i,j,k,3) * u(i,j,k,1)
           tt(i,j,k,7) = sedc(i,j,k) * u(i,j,k,1)
	   tt(i,j,k,8) = sedc(i,j,k) * u(i,j,k,2)
	   tt(i,j,k,9) = sedc(i,j,k) * u(i,j,k,3)
	   tt(i,j,k,10) = phi(i,j,k) * u(i,j,k,1)
	   tt(i,j,k,11) = phi(i,j,k) * u(i,j,k,2)
	   tt(i,j,k,12) = phi(i,j,k) * u(i,j,k,3)

	enddo
	enddo
	enddo 
        

C	T:
C	[uu] (1) [vv] (2) [ww] (3)
C	[uv] (4) [vw] (5) [wu] (6)
C	[ru] (7) [rv] (8) [rw] (9)

	weit = 6.D0
	fact = 1.D0 / 512.D0 
	call filter(nni,nnj,nnk, tt(0,0,0,1),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, tt(0,0,0,2),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, tt(0,0,0,3),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, tt(0,0,0,4),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, tt(0,0,0,5),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, tt(0,0,0,6),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, tt(0,0,0,7),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, tt(0,0,0,8),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, tt(0,0,0,9),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, tt(0,0,0,10),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, tt(0,0,0,11),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, tt(0,0,0,12),rr(0,0,0,8), weit,fact)

C	R:
C	u (1) v (2) w (3) r (4)

        do k = 0, nnk+1
        do j = 0, nnj+1
        do i = 0, nni+1
	   rr(i,j,k,1) = u(i,j,k,1)
	   rr(i,j,k,2) = u(i,j,k,2)
	   rr(i,j,k,3) = u(i,j,k,3)
           rr(i,j,k,4) = sedc(i,j,k)
	   rr(i,j,k,5) = phi(i,j,k)
	enddo
	enddo
	enddo 

C	R:
C	[u] (1) [v] (2) [w] (3) [r] (4)

	weit = 6.D0
	fact = 1.D0 / 512.D0 
	call filter(nni,nnj,nnk, rr(0,0,0,1),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, rr(0,0,0,2),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, rr(0,0,0,3),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, rr(0,0,0,4),rr(0,0,0,8), weit,fact)
	call filter(nni,nnj,nnk, rr(0,0,0,5),rr(0,0,0,8), weit,fact)

C	T:
C	L11 (1) L22 (2) L33 (3)
C	L12 (4) L23 (5) L31 (6)
C	P1  (7) P2  (8) P3  (9)

	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni
	   tt(i,j,k,1) = tt(i,j,k,1) - rr(i,j,k,1)**2
	   tt(i,j,k,2) = tt(i,j,k,2) - rr(i,j,k,2)**2
	   tt(i,j,k,3) = tt(i,j,k,3) - rr(i,j,k,3)**2
	   tt(i,j,k,4) = tt(i,j,k,4) - rr(i,j,k,1) * rr(i,j,k,2)
	   tt(i,j,k,5) = tt(i,j,k,5) - rr(i,j,k,2) * rr(i,j,k,3)
	   tt(i,j,k,6) = tt(i,j,k,6) - rr(i,j,k,3) * rr(i,j,k,1)
	   tt(i,j,k,7) = tt(i,j,k,7) - rr(i,j,k,4) * rr(i,j,k,1)
	   tt(i,j,k,8) = tt(i,j,k,8) - rr(i,j,k,4) * rr(i,j,k,2)
	   tt(i,j,k,9) = tt(i,j,k,9) - rr(i,j,k,4) * rr(i,j,k,3)
	   tt(i,j,k,10) = tt(i,j,k,10) - rr(i,j,k,5) * rr(i,j,k,1)
	   tt(i,j,k,11) = tt(i,j,k,11) - rr(i,j,k,5) * rr(i,j,k,2)
	   tt(i,j,k,12) = tt(i,j,k,12) - rr(i,j,k,5) * rr(i,j,k,3)
	enddo
	enddo
	enddo 

	call ssbc(tt(0,0,0,1))
	call ssbc(tt(0,0,0,2))
	call ssbc(tt(0,0,0,3))
	call ssbc(tt(0,0,0,4))
	call ssbc(tt(0,0,0,5))
	call ssbc(tt(0,0,0,6))
	call ssbc(tt(0,0,0,7))
	call ssbc(tt(0,0,0,8))
	call ssbc(tt(0,0,0,9))
	call ssbc(tt(0,0,0,10))
	call ssbc(tt(0,0,0,11))
	call ssbc(tt(0,0,0,12))
	call p_exchange( nni, nnj, nnk, tt(0,0,0,1) )   
	call p_exchange( nni, nnj, nnk, tt(0,0,0,2) )   
	call p_exchange( nni, nnj, nnk, tt(0,0,0,3) )   
	call p_exchange( nni, nnj, nnk, tt(0,0,0,4) )   
	call p_exchange( nni, nnj, nnk, tt(0,0,0,5) )   
	call p_exchange( nni, nnj, nnk, tt(0,0,0,6) )   
	call p_exchange( nni, nnj, nnk, tt(0,0,0,7) )   
	call p_exchange( nni, nnj, nnk, tt(0,0,0,8) )   
	call p_exchange( nni, nnj, nnk, tt(0,0,0,9) )
	call p_exchange( nni, nnj, nnk, tt(0,0,0,10) )   
	call p_exchange( nni, nnj, nnk, tt(0,0,0,11) )   
	call p_exchange( nni, nnj, nnk, tt(0,0,0,12) )
	
        do m = 1, 9
	do k = 0, nnk+1
	do j = 0, nnj+1
	do i = 0, nni+1
	   rr(i,j,k,m) = 0.D0
	enddo
	enddo
	enddo 
	enddo 

C	R:
C	K = 1 (5)
C	H = 1 (6)

	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni
	   rr(i,j,k,8) = 1.D0
	   rr(i,j,k,9) = 1.D0
	enddo
	enddo
	enddo


c        do k = 0, nnk
c	do j = 0, nnj
c	do i = 0, nni
c	   chi(i,j,k,5) =  tt(i,j,k,8)
c	   if (n_suth .eq. MPI_PROC_NULL) then
c	     if (j .eq. 0) then
c	        chi(i,j,k,5) = pick(i,k)
c	     endif
c	   endif
c	enddo
c	enddo 
c	enddo 
	
        do k = 1, nnk
	do j = 1, nnj
    	do i = 1, nni
	
	   rr(i,j,k,1) = rr(i,j,k,1)
     <	               + xix(i,  j,k) * ( tt(i,j,k,1) + tt(i+1,j,k,1) )
     <	               - xix(i-1,j,k) * ( tt(i,j,k,1) + tt(i-1,j,k,1) )
     <	               + xiy(i,  j,k) * ( tt(i,j,k,4) + tt(i+1,j,k,4) )
     <	               - xiy(i-1,j,k) * ( tt(i,j,k,4) + tt(i-1,j,k,4) )
     <	               + xiz(i,  j,k) * ( tt(i,j,k,6) + tt(i+1,j,k,6) )
     <	               - xiz(i-1,j,k) * ( tt(i,j,k,6) + tt(i-1,j,k,6) )
	   rr(i,j,k,1) = rr(i,j,k,1)
     <	               + etx(i,j,  k) * ( tt(i,j,k,1) + tt(i,j+1,k,1) )
     <	               - etx(i,j-1,k) * ( tt(i,j,k,1) + tt(i,j-1,k,1) )
     <	               + ety(i,j,  k) * ( tt(i,j,k,4) + tt(i,j+1,k,4) )
     <	               - ety(i,j-1,k) * ( tt(i,j,k,4) + tt(i,j-1,k,4) )
     <	               + etz(i,j,  k) * ( tt(i,j,k,6) + tt(i,j+1,k,6) )
     <	               - etz(i,j-1,k) * ( tt(i,j,k,6) + tt(i,j-1,k,6) )
	   rr(i,j,k,1) = rr(i,j,k,1)
     <	               + ztx(i,j,k  ) * ( tt(i,j,k,1) + tt(i,j,k+1,1) )
     <	               - ztx(i,j,k-1) * ( tt(i,j,k,1) + tt(i,j,k-1,1) )
     <	               + zty(i,j,k  ) * ( tt(i,j,k,4) + tt(i,j,k+1,4) )
     <	               - zty(i,j,k-1) * ( tt(i,j,k,4) + tt(i,j,k-1,4) )
     <	               + ztz(i,j,k  ) * ( tt(i,j,k,6) + tt(i,j,k+1,6) )
     <	               - ztz(i,j,k-1) * ( tt(i,j,k,6) + tt(i,j,k-1,6) )

	   rr(i,j,k,2) = rr(i,j,k,2)
     <	               + xix(i,  j,k) * ( tt(i,j,k,4) + tt(i+1,j,k,4) )
     <	               - xix(i-1,j,k) * ( tt(i,j,k,4) + tt(i-1,j,k,4) )
     <	               + xiy(i,  j,k) * ( tt(i,j,k,2) + tt(i+1,j,k,2) )
     <	               - xiy(i-1,j,k) * ( tt(i,j,k,2) + tt(i-1,j,k,2) )
     <	               + xiz(i,  j,k) * ( tt(i,j,k,5) + tt(i+1,j,k,5) )
     <	               - xiz(i-1,j,k) * ( tt(i,j,k,5) + tt(i-1,j,k,5) )
	   rr(i,j,k,2) = rr(i,j,k,2)
     <	               + etx(i,j,  k) * ( tt(i,j,k,4) + tt(i,j+1,k,4) )
     <	               - etx(i,j-1,k) * ( tt(i,j,k,4) + tt(i,j-1,k,4) )
     <	               + ety(i,j,  k) * ( tt(i,j,k,2) + tt(i,j+1,k,2) )
     <	               - ety(i,j-1,k) * ( tt(i,j,k,2) + tt(i,j-1,k,2) )
     <	               + etz(i,j,  k) * ( tt(i,j,k,5) + tt(i,j+1,k,5) )
     <	               - etz(i,j-1,k) * ( tt(i,j,k,5) + tt(i,j-1,k,5) )
	   rr(i,j,k,2) = rr(i,j,k,2)
     <	               + ztx(i,j,k  ) * ( tt(i,j,k,4) + tt(i,j,k+1,4) )
     <	               - ztx(i,j,k-1) * ( tt(i,j,k,4) + tt(i,j,k-1,4) )
     <	               + zty(i,j,k  ) * ( tt(i,j,k,2) + tt(i,j,k+1,2) )
     <	               - zty(i,j,k-1) * ( tt(i,j,k,2) + tt(i,j,k-1,2) )
     <	               + ztz(i,j,k  ) * ( tt(i,j,k,5) + tt(i,j,k+1,5) )
     <	               - ztz(i,j,k-1) * ( tt(i,j,k,5) + tt(i,j,k-1,5) )

	   rr(i,j,k,3) = rr(i,j,k,3)
     <	               + xix(i,  j,k) * ( tt(i,j,k,6) + tt(i+1,j,k,6) )
     <	               - xix(i-1,j,k) * ( tt(i,j,k,6) + tt(i-1,j,k,6) )
     <	               + xiy(i,  j,k) * ( tt(i,j,k,5) + tt(i+1,j,k,5) )
     <	               - xiy(i-1,j,k) * ( tt(i,j,k,5) + tt(i-1,j,k,5) )
     <	               + xiz(i,  j,k) * ( tt(i,j,k,3) + tt(i+1,j,k,3) )
     <	               - xiz(i-1,j,k) * ( tt(i,j,k,3) + tt(i-1,j,k,3) )
	   rr(i,j,k,3) = rr(i,j,k,3)
     <	               + etx(i,j,  k) * ( tt(i,j,k,6) + tt(i,j+1,k,6) )
     <	               - etx(i,j-1,k) * ( tt(i,j,k,6) + tt(i,j-1,k,6) )
     <	               + ety(i,j,  k) * ( tt(i,j,k,5) + tt(i,j+1,k,5) )
     <	               - ety(i,j-1,k) * ( tt(i,j,k,5) + tt(i,j-1,k,5) )
     <	               + etz(i,j,  k) * ( tt(i,j,k,3) + tt(i,j+1,k,3) )
     <	               - etz(i,j-1,k) * ( tt(i,j,k,3) + tt(i,j-1,k,3) )
	   rr(i,j,k,3) = rr(i,j,k,3)
     <	               + ztx(i,j,k  ) * ( tt(i,j,k,6) + tt(i,j,k+1,6) )
     <	               - ztx(i,j,k-1) * ( tt(i,j,k,6) + tt(i,j,k-1,6) )
     <	               + zty(i,j,k  ) * ( tt(i,j,k,5) + tt(i,j,k+1,5) )
     <	               - zty(i,j,k-1) * ( tt(i,j,k,5) + tt(i,j,k-1,5) )
     <	               + ztz(i,j,k  ) * ( tt(i,j,k,3) + tt(i,j,k+1,3) )
     <	               - ztz(i,j,k-1) * ( tt(i,j,k,3) + tt(i,j,k-1,3) )

	   sai(i,j,k) = tt(i,j,k,8)
	   rr(i,j,k,4) = rr(i,j,k,4)
     <	               + xix(i,  j,k) * ( tt(i,j,k,7) + tt(i+1,j,k,7) )
     <	               - xix(i-1,j,k) * ( tt(i,j,k,7) + tt(i-1,j,k,7) )
     <	               + xiy(i,  j,k) * ( tt(i,j,k,8) + tt(i+1,j,k,8) )
     <	               - xiy(i-1,j,k) * ( tt(i,j,k,8) + tt(i-1,j,k,8) )
     <	               + xiz(i,  j,k) * ( tt(i,j,k,9) + tt(i+1,j,k,9) )
     <	               - xiz(i-1,j,k) * ( tt(i,j,k,9) + tt(i-1,j,k,9) )
           rr(i,j,k,4) = rr(i,j,k,4)
     <	               + etx(i,j,  k) * ( tt(i,j,k,7) + tt(i,j+1,k,7) )
     <	               - etx(i,j-1,k) * ( tt(i,j,k,7) + tt(i,j-1,k,7) )
     <	               + ety(i,j,  k) * ( tt(i,j,k,8) + tt(i,j+1,k,8) )
     <	               - ety(i,j-1,k) * ( tt(i,j,k,8) + tt(i,j-1,k,8) )
     <	               + etz(i,j,  k) * ( tt(i,j,k,9) + tt(i,j+1,k,9) )
     <	               - etz(i,j-1,k) * ( tt(i,j,k,9) + tt(i,j-1,k,9) )
	   rr(i,j,k,4) = rr(i,j,k,4)
     <	               + ztx(i,j,k  ) * ( tt(i,j,k,7) + tt(i,j,k+1,7) )
     <	               - ztx(i,j,k-1) * ( tt(i,j,k,7) + tt(i,j,k-1,7) )
     <	               + zty(i,j,k  ) * ( tt(i,j,k,8) + tt(i,j,k+1,8) )
     <	               - zty(i,j,k-1) * ( tt(i,j,k,8) + tt(i,j,k-1,8) )
     <	               + ztz(i,j,k  ) * ( tt(i,j,k,9) + tt(i,j,k+1,9) )
     <	               - ztz(i,j,k-1) * ( tt(i,j,k,9) + tt(i,j,k-1,9) )
	   rr(i,j,k,5) = rr(i,j,k,5)
     <	               + xix(i,  j,k) * ( tt(i,j,k,10) + tt(i+1,j,k,10) )
     <	               - xix(i-1,j,k) * ( tt(i,j,k,10) + tt(i-1,j,k,10) )
     <	               + xiy(i,  j,k) * ( tt(i,j,k,11) + tt(i+1,j,k,11) )
     <	               - xiy(i-1,j,k) * ( tt(i,j,k,11) + tt(i-1,j,k,11) )
     <	               + xiz(i,  j,k) * ( tt(i,j,k,12) + tt(i+1,j,k,12) )
     <	               - xiz(i-1,j,k) * ( tt(i,j,k,12) + tt(i-1,j,k,12) )
           rr(i,j,k,5) = rr(i,j,k,5)
     <	               + etx(i,j,  k) * ( tt(i,j,k,10) + tt(i,j+1,k,10) )
     <	               - etx(i,j-1,k) * ( tt(i,j,k,10) + tt(i,j-1,k,10) )
     <	               + ety(i,j,  k) * ( tt(i,j,k,11) + tt(i,j+1,k,11) )
     <	               - ety(i,j-1,k) * ( tt(i,j,k,11) + tt(i,j-1,k,11) )
     <	               + etz(i,j,  k) * ( tt(i,j,k,12) + tt(i,j+1,k,12) )
     <	               - etz(i,j-1,k) * ( tt(i,j,k,12) + tt(i,j-1,k,12) )
	   rr(i,j,k,5) = rr(i,j,k,5)
     <	               + ztx(i,j,k  ) * ( tt(i,j,k,10) + tt(i,j,k+1,10) )
     <	               - ztx(i,j,k-1) * ( tt(i,j,k,10) + tt(i,j,k-1,10) )
     <	               + zty(i,j,k  ) * ( tt(i,j,k,11) + tt(i,j,k+1,11) )
     <	               - zty(i,j,k-1) * ( tt(i,j,k,11) + tt(i,j,k-1,11) )
     <	               + ztz(i,j,k  ) * ( tt(i,j,k,12) + tt(i,j,k+1,12) )
     <	               - ztz(i,j,k-1) * ( tt(i,j,k,12) + tt(i,j,k-1,12) )

	enddo
	enddo
	enddo 

	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni
	   rr(i,j,k,1) = 0.5D0 * rr(i,j,k,1)
	   rr(i,j,k,2) = 0.5D0 * rr(i,j,k,2)
	   rr(i,j,k,3) = 0.5D0 * rr(i,j,k,3)
	   if (n_suth .eq. MPI_PROC_NULL) then
c	      if (j .le. cut) then
	         if (j .eq. 1) then 
		    rr(i,j,k,4) = - ety(i,j-1,k)*pick(i,k)
		 else
		    rr(i,j,k,4) = 0.0D0
		 endif
c	      endif
           endif
	enddo
	enddo
	enddo

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
c	      s(i,nnj+1,k) = -s(i,nnj,k)
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
	       vst(0,j,k) =   -vst(1,j,k)
	      akst(0,j,k) =  -akst(1,j,k)
	      aksd(0,j,k) =  -aksd(1,j,k)
	   enddo
	   enddo
	endif

	if ( n_east .eq. MPI_PROC_NULL ) then
	   do k = 0, nnk+1
	   do j = 0, nnj+1
	       vst(nni+1,j,k) =   -vst(nni,j,k)
	      akst(nni+1,j,k) =  -akst(nni,j,k)
	      aksd(nni+1,j,k) =  -aksd(nni,j,k)
	   enddo
	   enddo
	endif

	if ( n_suth .eq. MPI_PROC_NULL ) then
	   do k = 0, nnk+1
	   do i = 0, nni+1
	       vst(i,0,k) =   -vst(i,1,k)
	      akst(i,0,k) = - akst(i,1,k)
	      aksd(i,0,k) = - aksd(i,1,k)
	   enddo
	   enddo
	endif

	if ( n_nrth .eq. MPI_PROC_NULL ) then
	   do k = 0, nnk+1
	   do i = 0, nni+1
	       vst(i,nnj+1,k) = - vst(i,nnj,k)
	      akst(i,nnj+1,k) = - akst(i,nnj,k)
	      aksd(i,nnj+1,k) = - aksd(i,nnj,k)
c               vst(i,nnj+1,k) = vst(i,nnj,k)
c	      akst(i,nnj+1,k) = akst(i,nnj,k)
c	      aksd(i,nnj+1,k) = aksd(i,nnj,k)
	   enddo
	   enddo
	endif

	if ( n_back .eq. MPI_PROC_NULL ) then
	   do j = 0, nnj+1
	   do i = 0, nni+1
	       vst(i,j,0) = -  vst(i,j,1)
	      akst(i,j,0) = - akst(i,j,1)
	      aksd(i,j,0) = - aksd(i,j,1)
c               vst(i,j,0) = vst(i,j,1)
c	      akst(i,j,0) = akst(i,j,1)
c	      aksd(i,j,0) = aksd(i,j,1)
	   enddo
	   enddo
	endif

	if ( n_frnt .eq. MPI_PROC_NULL ) then
	   do j = 0, nnj+1
	   do i = 0, nni+1
	       vst(i,j,nnk+1) = -  vst(i,j,nnk)
	      akst(i,j,nnk+1) = - akst(i,j,nnk)
	      aksd(i,j,nnk+1) = - aksd(i,j,nnk)
c               vst(i,j,nnk+1) =   vst(i,j,nnk)
c	      akst(i,j,nnk+1) =  akst(i,j,nnk)
c	      aksd(i,j,nnk+1) =  aksd(i,j,nnk)

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
