cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine eddy

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "ns.inc"
	include "metric.inc"
	include "para.inc"
	include "eddy.inc"

	double precision, dimension(0:nni+1,0:nnj+1,0:nnk+1,9) :: 
     <		ss, tt, zz

	integer i, j, k, m

	double precision tiny, alpha2, weit, fact, temp
	double precision fmin, fmax

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
	   rr(i,j,k,1) = phi(i+1,j,k) - phi(i-1,j,k)
	   rr(i,j,k,2) = phi(i,j+1,k) - phi(i,j-1,k)
	   rr(i,j,k,3) = phi(i,j,k+1) - phi(i,j,k-1)
	enddo
	enddo
	enddo 

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
C	<S12> (1) <S23> (2) <S32> (3)
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
	   tt(i,j,k,7) = phi(i,j,k) * u(i,j,k,1)
	   tt(i,j,k,8) = phi(i,j,k) * u(i,j,k,2)
	   tt(i,j,k,9) = phi(i,j,k) * u(i,j,k,3)
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
	   rr(i,j,k,4) = phi(i,j,k)
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
	enddo
	enddo
	enddo 

C	R:
C	u (1) v (2) w (3) r (4)

        do k = 0, nnk+1
        do j = 0, nnj+1
        do i = 0, nni+1
	   rr(i,j,k,1) = u(i,j,k,1)
	   rr(i,j,k,2) = u(i,j,k,2)
	   rr(i,j,k,3) = u(i,j,k,3)
	   rr(i,j,k,4) = phi(i,j,k)
	enddo
	enddo
	enddo 

C	R:
C	[u] (1) [v] (2) [w] (3) [r] (4)

	weit = 6.D0
	fact = 1.D0 / 512.D0 
	call filter(nni,nnj,nnk, rr(0,0,0,1),rr(0,0,0,5), weit,fact)
	call filter(nni,nnj,nnk, rr(0,0,0,2),rr(0,0,0,5), weit,fact)
	call filter(nni,nnj,nnk, rr(0,0,0,3),rr(0,0,0,5), weit,fact)
	call filter(nni,nnj,nnk, rr(0,0,0,4),rr(0,0,0,5), weit,fact)

	call ssbc(rr(0,0,0,1))
	call ssbc(rr(0,0,0,2))
	call ssbc(rr(0,0,0,3))
	call ssbc(rr(0,0,0,4))
	call p_exchange( nni, nnj, nnk, rr(0,0,0,1) )   
	call p_exchange( nni, nnj, nnk, rr(0,0,0,2) )   
	call p_exchange( nni, nnj, nnk, rr(0,0,0,3) )   
	call p_exchange( nni, nnj, nnk, rr(0,0,0,4) )   

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
	enddo
	enddo
	enddo 

C	S:
C	<[u][u]> (1) <[v][v]> (2) <[w][w]> (3)
C	<[u][v]> (4) <[v][w]> (5) <[w][u]> (6)
C	<[r][u]> (7) <[r][v]> (8) <[r][w]> (9)

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

C	R:
C	<[u]> (1) <[v]> (2) <[w]> (3) <[r]> (4)

	weit = 2.D0
	fact = 1.D0 / 64.D0 
	call filter(nni,nnj,nnk, rr(0,0,0,1),rr(0,0,0,5), weit,fact)
	call filter(nni,nnj,nnk, rr(0,0,0,2),rr(0,0,0,5), weit,fact)
	call filter(nni,nnj,nnk, rr(0,0,0,3),rr(0,0,0,5), weit,fact)
	call filter(nni,nnj,nnk, rr(0,0,0,4),rr(0,0,0,5), weit,fact)

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
	endif

        do k = 1, nnk
        do j = 1, nnj
        do i = 1, nni
	    vst(i,j,k) = max( vst(i,j,k), -vis)
	   akst(i,j,k) = max(akst(i,j,k), -ak )
	enddo
	enddo
	enddo 

	call eddybc
	call p_exchange( nni, nnj, nnk,  vst )   
	call p_exchange( nni, nnj, nnk, akst )   

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
	   tt(i,j,k,7) = phi(i,j,k) * u(i,j,k,1)
	   tt(i,j,k,8) = phi(i,j,k) * u(i,j,k,2)
	   tt(i,j,k,9) = phi(i,j,k) * u(i,j,k,3)
	enddo
	enddo
	enddo 

C	T:
C	[uu] (1) [vv] (2) [ww] (3)
C	[uv] (4) [vw] (5) [wu] (6)
C	[ru] (7) [rv] (8) [rw] (9)

	weit = 6.D0
	fact = 1.D0 / 512.D0 
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
	   rr(i,j,k,4) = phi(i,j,k)
	enddo
	enddo
	enddo 

C	R:
C	[u] (1) [v] (2) [w] (3) [r] (4)

	weit = 6.D0
	fact = 1.D0 / 512.D0 
	call filter(nni,nnj,nnk, rr(0,0,0,1),rr(0,0,0,5), weit,fact)
	call filter(nni,nnj,nnk, rr(0,0,0,2),rr(0,0,0,5), weit,fact)
	call filter(nni,nnj,nnk, rr(0,0,0,3),rr(0,0,0,5), weit,fact)
	call filter(nni,nnj,nnk, rr(0,0,0,4),rr(0,0,0,5), weit,fact)

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
	call p_exchange( nni, nnj, nnk, tt(0,0,0,1) )   
	call p_exchange( nni, nnj, nnk, tt(0,0,0,2) )   
	call p_exchange( nni, nnj, nnk, tt(0,0,0,3) )   
	call p_exchange( nni, nnj, nnk, tt(0,0,0,4) )   
	call p_exchange( nni, nnj, nnk, tt(0,0,0,5) )   
	call p_exchange( nni, nnj, nnk, tt(0,0,0,6) )   
	call p_exchange( nni, nnj, nnk, tt(0,0,0,7) )   
	call p_exchange( nni, nnj, nnk, tt(0,0,0,8) )   
	call p_exchange( nni, nnj, nnk, tt(0,0,0,9) )

	do m = 1, 6
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
	   rr(i,j,k,5) = 1.D0
	   rr(i,j,k,6) = 1.D0
	enddo
	enddo
	enddo

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

	enddo
	enddo
	enddo 

	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni
	   rr(i,j,k,1) = 0.5D0 * rr(i,j,k,1)
	   rr(i,j,k,2) = 0.5D0 * rr(i,j,k,2)
	   rr(i,j,k,3) = 0.5D0 * rr(i,j,k,3)
	   rr(i,j,k,4) = 0.5D0 * rr(i,j,k,4)
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

	double precision, dimension(0:nni+1,0:nnj+1,0:nnk+1) :: s
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
	       vst(0,j,k) = -  vst(1,j,k)
	      akst(0,j,k) = - akst(1,j,k)
	   enddo
	   enddo
	endif

	if ( n_east .eq. MPI_PROC_NULL ) then
	   do k = 0, nnk+1
	   do j = 0, nnj+1
	       vst(nni+1,j,k) = -  vst(nni,j,k)
	      akst(nni+1,j,k) = - akst(nni,j,k)
	   enddo
	   enddo
	endif

	if ( n_suth .eq. MPI_PROC_NULL ) then
	   do k = 0, nnk+1
	   do i = 0, nni+1
	       vst(i,0,k) = -  vst(i,1,k)
	      akst(i,0,k) = - akst(i,1,k)
	   enddo
	   enddo
	endif

	if ( n_nrth .eq. MPI_PROC_NULL ) then
	   do k = 0, nnk+1
	   do i = 0, nni+1
	       vst(i,nnj+1,k) = -  vst(i,nnj,k)
	      akst(i,nnj+1,k) = - akst(i,nnj,k)
	   enddo
	   enddo
	endif

	if ( n_back .eq. MPI_PROC_NULL ) then
	   do j = 0, nnj+1
	   do i = 0, nni+1
	       vst(i,j,0) = -  vst(i,j,1)
	      akst(i,j,0) = - akst(i,j,1)
	   enddo
	   enddo
	endif

	if ( n_frnt .eq. MPI_PROC_NULL ) then
	   do j = 0, nnj+1
	   do i = 0, nni+1
	       vst(i,j,nnk+1) = -  vst(i,j,nnk)
	      akst(i,j,nnk+1) = - akst(i,j,nnk)
	   enddo
	   enddo
	endif

	return
	end
	
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine filter(ii, jj, kk, f, r, weit, fact)

	implicit none

	integer ii, jj, kk

	double precision, dimension(0:ii+1,0:jj+1,0:kk+1) :: 
     <		f, r

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

	double precision, dimension(0:ii+1,0:jj+1,0:kk+1) :: f

	double precision fmin, fmax, ffmin, ffmax

	integer i, j, k

	fmin =    HUGE(1.D0)
	fmax =   -HUGE(1.D0)
	ffmin =   HUGE(1.D0)
	ffmax =  -HUGE(1.D0)

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
