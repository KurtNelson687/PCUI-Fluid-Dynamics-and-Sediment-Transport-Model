cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine convection

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "ns.inc"

	integer i, j, k, m
	double precision fne, fpe, fnw, fpw, fnn, fpn, 
     <	                 fns, fps, fnf, fpf, fnb, fpb,
     <	                 aeu, awu, anu, asu, afu, abu, apu, rms,
     <	                 bpw, bne, bps, bnn, bpb, bnf

	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni

	   fne = 0.5d0 * (uxi(i,  j,  k) - dabs(uxi(i,  j,  k  )))
	   fpe = 0.5d0 * (uxi(i,  j,  k) + dabs(uxi(i,  j,  k  )))
	   fnw = 0.5d0 * (uxi(i-1,j,  k) - dabs(uxi(i-1,j,  k  )))
	   fpw = 0.5d0 * (uxi(i-1,j,  k) + dabs(uxi(i-1,j,  k  )))
	   fnn = 0.5d0 * (uej(i,  j,  k) - dabs(uej(i,  j,  k  )))
	   fpn = 0.5d0 * (uej(i,  j,  k) + dabs(uej(i,  j,  k  )))
	   fns = 0.5d0 * (uej(i,  j-1,k) - dabs(uej(i,  j-1,k  )))
	   fps = 0.5d0 * (uej(i,  j-1,k) + dabs(uej(i,  j-1,k  )))
	   fnf = 0.5d0 * (uzk(i,  j,  k) - dabs(uzk(i,  j,  k  )))
	   fpf = 0.5d0 * (uzk(i,  j,  k) + dabs(uzk(i,  j,  k  )))
	   fnb = 0.5d0 * (uzk(i,  j,k-1) - dabs(uzk(i,  j,  k-1)))
	   fpb = 0.5d0 * (uzk(i,  j,k-1) + dabs(uzk(i,  j,  k-1)))

	   rms = ( uxi(i,j,k) - uxi(i-1,j,k) )
     <	       + ( uej(i,j,k) - uej(i,j-1,k) )
     <	       + ( uzk(i,j,k) - uzk(i,j,k-1) )
         
	   aeu = - 0.5  D0 * uxi(i,j,k) 
     <	         + 0.125D0 * ( fpe - 2.d0 * fne - fnw )
	   awu =   0.5  D0 * uxi(i-1,j,k) 
     <	         + 0.125D0 * ( fpe + 2.d0 * fpw - fnw )
	   anu = - 0.5  D0 * uej(i,j,k)
     <	         + 0.125D0 * ( fpn - 2.d0 * fnn - fns )
	   asu =   0.5  D0 * uej(i,j-1,k)
     <	         + 0.125D0 * ( fpn + 2.d0 * fps - fns )
	   afu = - 0.5  D0 * uzk(i,j,k)
     <	         + 0.125D0 * ( fpf - 2.d0 * fnf - fnb )
	   abu =   0.5  D0 * uzk(i,j,k-1)
     <	         + 0.125D0 * ( fpf + 2.d0 * fpb - fnb )
	   apu =   aeu + awu + anu + asu + afu + abu
     <	         + 0.125D0 * ( fne + fnn + fnf - fpw - fps - fpb )
     <	         + rms
	
	   do m = 1, 3

	      bpw = - fpw * u(i-2,j,k,m)
	      bne =   fne * u(i+2,j,k,m)
	      bps = - fps * u(i,j-2,k,m)
	      bnn =   fnn * u(i,j+2,k,m)
	      bpb = - fpb * u(i,j,k-2,m)
	      bnf =   fnf * u(i,j,k+2,m)
            
	      hb(i,j,k,m) = aeu * u(i+1,j,k,m) + awu * u(i-1,j,k,m)
     <	                  + anu * u(i,j+1,k,m) + asu * u(i,j-1,k,m)
     <	                  + afu * u(i,j,k+1,m) + abu * u(i,j,k-1,m)
     <	                  - apu * u(i,j,k,m)
     <		+ 0.125D0 * ( bpw + bne + bps + bnn + bpb + bnf )

	   enddo

	enddo
	enddo
	enddo

	return
	end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine quick

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "para.inc"
	include "ns.inc"
	include "metric.inc"

	integer i, j, k, L

C......	UXI

	do k = 1, nnk
	do j = 1, nnj
	do i = ius, iue
	if ( uxi(i,j,k) .ge. 0.D0 ) then
	   uxi(i,j,k) = xix(i,j,k) * 
     <	        ( 0.75 D0 * u(i,j,k,1) + 0.375D0 * u(i+1,j,k,1) 
     <          - 0.125D0 * u(i-1,j,k,1) )
     <	              + xiy(i,j,k) *
     <	        ( 0.75 D0 * u(i,j,k,2) + 0.375D0 * u(i+1,j,k,2) 
     <          - 0.125D0 * u(i-1,j,k,2) )
     <	              + xiz(i,j,k) *
     <	        ( 0.75 D0 * u(i,j,k,3) + 0.375D0 * u(i+1,j,k,3) 
     <          - 0.125D0 * u(i-1,j,k,3) )
     	else
	   uxi(i,j,k) = xix(i,j,k) * 
     <	        ( 0.75 D0 * u(i+1,j,k,1) + 0.375D0 * u(i,j,k,1) 
     <          - 0.125D0 * u(i+2,j,k,1) )
     <	              + xiy(i,j,k) *
     <	        ( 0.75 D0 * u(i+1,j,k,2) + 0.375D0 * u(i,j,k,2) 
     <          - 0.125D0 * u(i+2,j,k,2) )
     <	              + xiz(i,j,k) *
     <	        ( 0.75 D0 * u(i+1,j,k,3) + 0.375D0 * u(i,j,k,3) 
     <          - 0.125D0 * u(i+2,j,k,3) )
     	endif
	enddo
	enddo
	enddo

C......	UEJ

	do k = 1, nnk
	do j = jus, jue
	do i = 1, nni
	if ( uej(i,j,k) .ge. 0.D0 ) then
	   uej(i,j,k) = etx(i,j,k) * 
     <	        ( 0.75 D0 * u(i,j,k,1) + 0.375D0 * u(i,j+1,k,1) 
     <          - 0.125D0 * u(i,j-1,k,1) )
     <	              + ety(i,j,k) *
     <	        ( 0.75 D0 * u(i,j,k,2) + 0.375D0 * u(i,j+1,k,2) 
     <          - 0.125D0 * u(i,j-1,k,2) )
     <	              + etz(i,j,k) *
     <	        ( 0.75 D0 * u(i,j,k,3) + 0.375D0 * u(i,j+1,k,3) 
     <          - 0.125D0 * u(i,j-1,k,3) )
     	else
	   uej(i,j,k) = etx(i,j,k) * 
     <	        ( 0.75 D0 * u(i,j+1,k,1) + 0.375D0 * u(i,j,k,1) 
     <          - 0.125D0 * u(i,j+2,k,1) )
     <	              + ety(i,j,k) *
     <	        ( 0.75 D0 * u(i,j+1,k,2) + 0.375D0 * u(i,j,k,2) 
     <          - 0.125D0 * u(i,j+2,k,2) )
     <	              + etz(i,j,k) *
     <	        ( 0.75 D0 * u(i,j+1,k,3) + 0.375D0 * u(i,j,k,3) 
     <          - 0.125D0 * u(i,j+2,k,3) )
     	endif
	enddo
	enddo
	enddo

C......	UZK

	do k = kus, kue
	do j = 1, nnj
	do i = 1, nni
	if ( uzk(i,j,k) .ge. 0.D0 ) then
	   uzk(i,j,k) = ztx(i,j,k) * 
     <	        ( 0.75 D0 * u(i,j,k,1) + 0.375D0 * u(i,j,k+1,1) 
     <          - 0.125D0 * u(i,j,k-1,1) )
     <	              + zty(i,j,k) *
     <	        ( 0.75 D0 * u(i,j,k,2) + 0.375D0 * u(i,j,k+1,2) 
     <          - 0.125D0 * u(i,j,k-1,2) )
     <	              + ztz(i,j,k) *
     <	        ( 0.75 D0 * u(i,j,k,3) + 0.375D0 * u(i,j,k+1,3) 
     <          - 0.125D0 * u(i,j,k-1,3) )
     	else
	   uzk(i,j,k) = ztx(i,j,k) * 
     <	        ( 0.75 D0 * u(i,j,k+1,1) + 0.375D0 * u(i,j,k,1) 
     <          - 0.125D0 * u(i,j,k+2,1) )
     <	              + zty(i,j,k) *
     <	        ( 0.75 D0 * u(i,j,k+1,2) + 0.375D0 * u(i,j,k,2) 
     <          - 0.125D0 * u(i,j,k+2,2) )
     <	              + ztz(i,j,k) *
     <	        ( 0.75 D0 * u(i,j,k+1,3) + 0.375D0 * u(i,j,k,3) 
     <          - 0.125D0 * u(i,j,k+2,3) )
	endif
	enddo
	enddo
	enddo
	
	return
	end
