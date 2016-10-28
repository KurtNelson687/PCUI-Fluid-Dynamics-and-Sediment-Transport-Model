cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine updateGrid

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "metric.inc"
	include "mygrid.inc"

	integer i, j, k


        
        
	call updateMetric(1, nni, nnj, nnk,
     <          g11, g12, g13, g21, g22, g23, g31, g32, g33, gcc)
        
	call qf2c(nni, nnj, nnk, nni1, nnj1, nnk1, jac, jaf,
     <          g11, g12, g13, g21, g22, g23, g31, g32, g33, gcc,
     <          f11, f12, f13, f21, f22, f23, f31, f32, f33, fcc)

	call qf2c(nni1, nnj1, nnk1, nni2, nnj2, nnk2, jaf, jah,
     <          f11, f12, f13, f21, f22, f23, f31, f32, f33, fcc,
     <          h11, h12, h13, h21, h22, h23, h31, h32, h33, hcc)

	call qf2c(nni2, nnj2, nnk2, nni3, nnj3, nnk3, jah, jap,
     <          h11, h12, h13, h21, h22, h23, h31, h32, h33, hcc,
     <          p11, p12, p13, p21, p22, p23, p31, p32, p33, pcc)

	call qf2c(nni3, nnj3, nnk3, nni4, nnj4, nnk4, jap, jar,
     <          p11, p12, p13, p21, p22, p23, p31, p32, p33, pcc,
     <          r11, r12, r13, r21, r22, r23, r31, r32, r33, rcc)

	return
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine updateMetric(level, ii, jj, kk, 
     <          q11, q12, q13, q21, q22, q23, q31, q32, q33, qcc)

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "metric.inc"
	include "mygrid.inc"
	include "para.inc"

	integer level, ii, jj, kk, ii0, jj0, kk0

	double precision q11(-1:ii+2,-1:jj+2,-1:kk+2),
     &                   q12(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   q13(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   q21(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   q22(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   q23(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   q31(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   q32(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   q33(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   qcc(-1:ii+2,-1:jj+2,-1:kk+2) 
	double precision x(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   y(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   z(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   xi(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   et(-1:ii+2,-1:jj+2,-1:kk+2), 
     &                   zt(-1:ii+2,-1:jj+2,-1:kk+2)
     
        double precision jab
	double precision xxi, yxi, zxi, xet, yet, zet, xzt, yzt, zzt
	double precision xsx, xsy, xsz, esx, esy, esz, zsx, zsy, zsz

	integer i, j, k, ic
	double precision faci, facj, fack
        

	ii0 = npx * ii
	jj0 = npy * jj
	kk0 = npz * kk

	faci = 1.D0 / dble(ii*px)
        facj = 1.D0 / dble(jj*py)
	fack = 1.D0 / dble(kk*pz)
	
	do k = -1, kk+2
	do j = -1, jj+2
	do i = -1, ii+2
	   q11(i,j,k) = 0.D0
	   q12(i,j,k) = 0.D0
	   q13(i,j,k) = 0.D0
	   q21(i,j,k) = 0.D0
	   q22(i,j,k) = 0.D0
	   q23(i,j,k) = 0.D0
	   q31(i,j,k) = 0.D0
	   q32(i,j,k) = 0.D0
	   q33(i,j,k) = 0.D0
	   qcc(i,j,k) = 0.D0
	enddo
	enddo
	enddo

	do k = -1, kk+2
	do j = -1, jj+2
	do i = -1, ii+2
	   xi(i,j,k) = ( dble(ii0 + i) - 0.5D0 ) * faci
	   zt(i,j,k) = ( dble(kk0 + k) - 0.5D0 ) * fack
	   x(i,j,k) = 0.D0
	   z(i,j,k) = 0.D0
	enddo
	enddo
	enddo
        
	call coordxz(ii, jj, kk, x, z, xi, zt)
        
c        call returnDepth(ii, jj, kk, x, z, depth, tx1, tx2, tz1, tz2)
        
        do k = -1, kk+2
        do j = -1, jj+2
	do i = -1, ii+2
	   et(i,j,k) = -depth(i, k)+( dble(jj0 + j) - 0.5D0 )
     <                * facj * depth(i, k)
	   y(i,j,k) = 0.D0
	enddo
        enddo
	enddo
  
        call coordy(ii, jj, kk, x, y, z, xi, et, zt)

	if ( level .eq. 1 ) then
	   do k = -1, kk+2
	   do j = -1, jj+2
	   do i = -1, ii+2
	      xp(i,j,k,1) = x(i,j,k)
	      xp(i,j,k,2) = y(i,j,k)
	      xp(i,j,k,3) = z(i,j,k)
	   enddo
	   enddo
	   enddo
	endif	
	
C...... I-face
	
	do k =  0, kk+1
	do j =  0, jj+1
	do i = -1, ii+1
	   xxi = ( x(i+1,j,k) - x(i,j,k) )
	   yxi = ( y(i+1,j,k) - y(i,j,k) )
	   zxi = ( z(i+1,j,k) - z(i,j,k) )
	   xet = 0.25D0 * ( x(i,  j+1,k) - x(i,  j-1,k) 
     <                    + x(i+1,j+1,k) - x(i+1,j-1,k) )
	   yet = 0.25D0 * ( y(i,  j+1,k) - y(i,  j-1,k)
     <                    + y(i+1,j+1,k) - y(i+1,j-1,k) )
	   zet = 0.25D0 * ( z(i,  j+1,k) - z(i,  j-1,k)
     <                    + z(i+1,j+1,k) - z(i+1,j-1,k) )
	   xzt = 0.25D0 * ( x(i,  j,k+1) - x(i,  j,k-1)
     <                    + x(i+1,j,k+1) - x(i+1,j,k-1) )
	   yzt = 0.25D0 * ( y(i,  j,k+1) - y(i,  j,k-1)
     <                    + y(i+1,j,k+1) - y(i+1,j,k-1) )
	   zzt = 0.25D0 * ( z(i,  j,k+1) - z(i,  j,k-1)
     <                    + z(i+1,j,k+1) - z(i+1,j,k-1) )
	   jab = xxi * yet * zzt
     <         + xet * yzt * zxi
     <         + xzt * yxi * zet
     <         - xzt * yet * zxi
     <         - xet * yxi * zzt
     <         - xxi * yzt * zet
	   jab = 1.D0 / jab
	   xsx = ( yet * zzt - yzt * zet )
	   esx = ( yzt * zxi - yxi * zzt )
	   zsx = ( yxi * zet - yet * zxi )
	   xsy = ( zet * xzt - zzt * xet )
	   esy = ( zzt * xxi - zxi * xzt )
	   zsy = ( zxi * xet - zet * xxi )
	   xsz = ( xet * yzt - xzt * yet )
	   esz = ( xzt * yxi - xxi * yzt )
	   zsz = ( xxi * yet - xet * yxi )
	   q11(i,j,k) = jab * ( xsx * xsx + xsy * xsy + xsz * xsz )
	   q12(i,j,k) = jab * ( xsx * esx + xsy * esy + xsz * esz )
	   q13(i,j,k) = jab * ( xsx * zsx + xsy * zsy + xsz * zsz )
	   
	   if ( level .eq. 1 ) then
	      xix(i,j,k) = xsx
	      xiy(i,j,k) = xsy
	      xiz(i,j,k) = xsz
	   endif
	enddo
	enddo
	enddo
	
C...... J-face
	
	do k =  0, kk+1
	do j = -1, jj+1
	do i =  0, ii+1
	   xet = ( x(i,j+1,k) - x(i,j,k) )
	   yet = ( y(i,j+1,k) - y(i,j,k) )
	   zet = ( z(i,j+1,k) - z(i,j,k) )
	   xzt = 0.25D0 * ( x(i,j,  k+1) - x(i,j,  k-1)
     <                    + x(i,j+1,k+1) - x(i,j+1,k-1) )
	   yzt = 0.25D0 * ( y(i,j,  k+1) - y(i,j,  k-1)
     <                    + y(i,j+1,k+1) - y(i,j+1,k-1) )
	   zzt = 0.25D0 * ( z(i,j,  k+1) - z(i,j,  k-1)
     <                    + z(i,j+1,k+1) - z(i,j+1,k-1) )
	   xxi = 0.25D0 * ( x(i+1,j  ,k) - x(i-1,j  ,k) 
     <                    + x(i+1,j+1,k) - x(i-1,j+1,k) )
	   yxi = 0.25D0 * ( y(i+1,j  ,k) - y(i-1,j  ,k) 
     <                    + y(i+1,j+1,k) - y(i-1,j+1,k) )
	   zxi = 0.25D0 * ( z(i+1,j  ,k) - z(i-1,j  ,k) 
     <                    + z(i+1,j+1,k) - z(i-1,j+1,k) )
	   jab = xxi * yet * zzt
     <         + xet * yzt * zxi
     <         + xzt * yxi * zet
     <         - xzt * yet * zxi
     <         - xet * yxi * zzt
     <         - xxi * yzt * zet
	   jab = 1.D0 / jab
	   xsx = ( yet * zzt - yzt * zet )
	   esx = ( yzt * zxi - yxi * zzt )
	   zsx = ( yxi * zet - yet * zxi )
	   xsy = ( zet * xzt - zzt * xet )
	   esy = ( zzt * xxi - zxi * xzt )
	   zsy = ( zxi * xet - zet * xxi )
	   xsz = ( xet * yzt - xzt * yet )
	   esz = ( xzt * yxi - xxi * yzt )
	   zsz = ( xxi * yet - xet * yxi )
	   q21(i,j,k) = jab * ( esx * xsx + esy * xsy + esz * xsz )
	   q22(i,j,k) = jab * ( esx * esx + esy * esy + esz * esz )
	   q23(i,j,k) = jab * ( esx * zsx + esy * zsy + esz * zsz )
	   
	   if ( level .eq. 1 ) then
	      etx(i,j,k) = esx
	      ety(i,j,k) = esy
	      etz(i,j,k) = esz
	   endif
	enddo
	enddo
	enddo

C...... K-face
	
	do k = -1, kk+1
	do j =  0, jj+1
	do i =  0, ii+1
	   xzt = ( x(i,j,k+1) - x(i,j,k) )
	   yzt = ( y(i,j,k+1) - y(i,j,k) )
	   zzt = ( z(i,j,k+1) - z(i,j,k) )
	   xxi = 0.25D0 * ( x(i+1,j,k  ) - x(i-1,j,k  ) 
     <                    + x(i+1,j,k+1) - x(i-1,j,k+1) )
	   yxi = 0.25D0 * ( y(i+1,j,k  ) - y(i-1,j,k  ) 
     <                    + y(i+1,j,k+1) - y(i-1,j,k+1) )
	   zxi = 0.25D0 * ( z(i+1,j,k  ) - z(i-1,j,k  ) 
     <                    + z(i+1,j,k+1) - z(i-1,j,k+1) )
	   xet = 0.25D0 * ( x(i,j+1,k  ) - x(i,j-1,k  )
     <                    + x(i,j+1,k+1) - x(i,j-1,k+1) )
	   yet = 0.25D0 * ( y(i,j+1,k  ) - y(i,j-1,k  )
     <                    + y(i,j+1,k+1) - y(i,j-1,k+1) )
	   zet = 0.25D0 * ( z(i,j+1,k  ) - z(i,j-1,k  )
     <                    + z(i,j+1,k+1) - z(i,j-1,k+1) )
	   jab = xxi * yet * zzt
     <         + xet * yzt * zxi
     <         + xzt * yxi * zet
     <         - xzt * yet * zxi
     <         - xet * yxi * zzt
     <         - xxi * yzt * zet
	   jab = 1.D0 / jab
	   xsx = ( yet * zzt - yzt * zet )
	   esx = ( yzt * zxi - yxi * zzt )
	   zsx = ( yxi * zet - yet * zxi )
	   xsy = ( zet * xzt - zzt * xet )
	   esy = ( zzt * xxi - zxi * xzt )
	   zsy = ( zxi * xet - zet * xxi )
	   xsz = ( xet * yzt - xzt * yet )
	   esz = ( xzt * yxi - xxi * yzt )
	   zsz = ( xxi * yet - xet * yxi )
	   q31(i,j,k) = jab * ( zsx * xsx + zsy * xsy + zsz * xsz )
	   q32(i,j,k) = jab * ( zsx * esx + zsy * esy + zsz * esz )
	   q33(i,j,k) = jab * ( zsx * zsx + zsy * zsy + zsz * zsz )

           if ( level .eq. 1 ) then
	      ztx(i,j,k) = zsx
	      zty(i,j,k) = zsy
	      ztz(i,j,k) = zsz
	   endif
	enddo
	enddo
	enddo
	
C...... Center

	if ( level .eq. 1 ) then
	   do k = 0, kk+1
	   do j = 0, jj+1
	   do i = 0, ii+1
	      xxi = 0.5D0 * ( x(i+1,j,k) - x(i-1,j,k) )
	      yxi = 0.5D0 * ( y(i+1,j,k) - y(i-1,j,k) )
	      zxi = 0.5D0 * ( z(i+1,j,k) - z(i-1,j,k) )
	      xet = 0.5D0 * ( x(i,j+1,k) - x(i,j-1,k) )
	      yet = 0.5D0 * ( y(i,j+1,k) - y(i,j-1,k) )
	      zet = 0.5D0 * ( z(i,j+1,k) - z(i,j-1,k) )
	      xzt = 0.5D0 * ( x(i,j,k+1) - x(i,j,k-1) )
	      yzt = 0.5D0 * ( y(i,j,k+1) - y(i,j,k-1) )
	      zzt = 0.5D0 * ( z(i,j,k+1) - z(i,j,k-1) )
	      jab = xxi * yet * zzt
     <            + xet * yzt * zxi
     <            + xzt * yxi * zet
     <            - xzt * yet * zxi
     <            - xet * yxi * zzt
     <            - xxi * yzt * zet
	      jab = 1.D0 / jab
	      jac(i,j,k) = jab
	   enddo
	   enddo
	   enddo
	   do k = -1, kk+2
	   do j = -1, jj+2
	   do i = -1, ii+2
	      phi_init(i,j,k) = 0.D0
	   enddo
	   enddo
	   enddo
	endif

	do k = -1, kk+2
	do j = -1, jj+2
	do i = -1, ii+2
	   q12(i,j,k) = 0.25D0 * q12(i,j,k)
	   q13(i,j,k) = 0.25D0 * q13(i,j,k)
	   q23(i,j,k) = 0.25D0 * q23(i,j,k)
	   q21(i,j,k) = 0.25D0 * q21(i,j,k)
	   q31(i,j,k) = 0.25D0 * q31(i,j,k)
	   q32(i,j,k) = 0.25D0 * q32(i,j,k)
	enddo
	enddo
	enddo

	do k = 0, kk+2
	do j = 0, jj+2
	do i = 0, ii+2
	   qcc(i,j,k) = q11(i,j,k) + q11(i-1,j,  k  )
     <                + q22(i,j,k) + q22(i,  j-1,k  )
     <                + q33(i,j,k) + q33(i,  j,  k-1)
	enddo 
	enddo
	enddo

	return
	end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine oldGrid
	
	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "metric.inc"
	include "mygrid.inc"

	integer i, j, k, m

	do k = -1,nnk+2
	do j = -1,nnj+2
	do i = -1,nni+2
	   jac_old(i,j,k) = jac(i,j,k)
	   xix_old(i,j,k) = xix(i,j,k)
	   xiy_old(i,j,k) = xiy(i,j,k)
	   xiz_old(i,j,k) = xiz(i,j,k)
	   etx_old(i,j,k) = etx(i,j,k)
	   ety_old(i,j,k) = ety(i,j,k)
	   etz_old(i,j,k) = etz(i,j,k)
	   ztx_old(i,j,k) = ztx(i,j,k)
	   zty_old(i,j,k) = zty(i,j,k)
	   ztz_old(i,j,k) = ztz(i,j,k)
	enddo
	enddo
	enddo

	do m = 1, 3
	do k = -1,nnk+2
	do j = -1,nnj+2
	do i = -1,nni+2
	    xp_old(i,j,k,m) =  xp(i,j,k,m)
	enddo
	enddo
	enddo
	enddo

        

	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine getEtadot

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "metric.inc"
	include "mygrid.inc"
	include "para.inc"
	include "ns.inc"
	include "sedi.inc"

	integer i, j, k, m
	integer indx_east, indx_west, indx_nrth, indx_suth, 
     &          indx_back, indx_frnt
	double precision temp
c	double precision qxi(-1:nni+2, -1:nnj+2, -1:nnk+2),
c     &                   qet(-1:nni+2, -1:nnj+2, -1:nnk+2),
c     &                   qzt(-1:nni+2, -1:nnj+2, -1:nnk+2)	
	

	temp = 1.D0/dtime
	do k = -1, nnk+2
	do j = -1, nnj+2
	do i = -1, nni+2
	   jac_diff(i,j,k) = (1.D0/jac(i,j,k)-1.D0/jac_old(i,j,k))*temp
	enddo
	enddo
	enddo

	temp =0.5D0/(ncheck*dtime)
	do k = 0, nnk+1
c	do j = 0, nnj+1
	do i = 0, nni+1
	   kzk(i,0,k) = 0.D0
c	   kxi(i,j,k) = temp*((xix_old(i,j,k)+xix_old(i,j,k))
c     <                  *0.5D0*(    xp(i,j,k,1)+    xp(i+1,j,k,1)
c     <                         -xp_old(i,j,k,1)-xp_old(i+1,j,k,1))
c     <                       +(xiy_old(i,j,k)+xiy_old(i,j,k))
c     <                  *0.5D0*(    xp(i,j,k,2)+    xp(i+1,j,k,2)
c     <                         -xp_old(i,j,k,2)-xp_old(i+1,j,k,2))
c     <                       +(xiz_old(i,j,k)+xiz_old(i,j,k))
c     <                  *0.5D0*(    xp(i,j,k,3)+    xp(i+1,j,k,3)
c     <                         -xp_old(i,j,k,3)-xp_old(i+1,j,k,3)))

	   kej(i,0,k) = temp*((etx(i,0,k)+etx(i,0,k))
     <                  *0.5D0*(    xp(i,0,k,1)+    xp(i,  1,k,1)
     <                         -xp_old(i,0,k,1)-xp_old(i,  1,k,1))
     <                       +(ety(i,0,k)+ety(i,0,k))
     <                  *0.5D0*(    xp(i,0,k,2)+    xp(i,  1,k,2)
     <                         -xp_old(i,0,k,2)-xp_old(i,  1,k,2))
     <                       +(etz(i,0,k)+etz(i,0,k))
     <                  *0.5D0*(    xp(i,0,k,3)+    xp(i,  1,k,3)
     <                         -xp_old(i,0,k,3)-xp_old(i,  1,k,3)))
	   if (istep .ne. 1) then
	   kej(i,0,k) = 1.5D0*kej(i,0,k)-0.5D0*kej_old(i,0,k)
	   endif
           kzk(i,0,k) = 0.D0
c	   kzk(i,j,k) = temp*((ztx_old(i,j,k)+ztx_old(i,j,k))
c     <                  *0.5D0*(    xp(i,j,k,1)+    xp(i,j,k+1,1)
c     <                         -xp_old(i,j,k,1)-xp_old(i,j,k+1,1))
c     <                       +(zty_old(i,j,k)+zty_old(i,j,k))
c     <                  *0.5D0*(    xp(i,j,k,2)+    xp(i,j,k+1,2)
c     <                         -xp_old(i,j,k,2)-xp_old(i,j,k+1,2))
c     <                       +(ztz_old(i,j,k)+ztz_old(i,j,k))
c     <                  *0.5D0*(    xp(i,j,k,3)+    xp(i,j,k+1,3)
c     <                         -xp_old(i,j,k,3)-xp_old(i,j,k+1,3)))
	   
	enddo
c	enddo
	enddo
	
	temp = 1.D0/3.D0
	do k = 0, nnk+1
	do j = 1, nnj+1
	do i = 0, nni+1
	   kxi(i,j,k) = 0.D0
	   kzk(i,j,k) = 0.D0
	   if (istep .ne. 1) then
	   kej(i,j,k) = kej(i,j-1,k)+2.D0*temp*jac_diff(i,j,k)
     <                   + temp*(kej_old(i,j,k) - kej_old(i,j-1,k))
	   else
	   kej(i,j,k) = kej(i,j-1,k)+jac_diff(i,j,k)	   
	   endif

	enddo
	enddo
	enddo



	do k = 0,nnk+1
	do j = 0,nnj+1
	do i = 0,nni+1
	   kxi_old(i,j,k) = kxi(i,j,k)
	   kej_old(i,j,k) = kej(i,j,k)
	   kzk_old(i,j,k) = kzk(i,j,k)
c	   write(*,*) i,j,k,kej(i,j,k)-kej(i,j-1,k)-jac_diff(i,j,k)
	enddo
	enddo
	enddo

	
	return 
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine getVlid

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "metric.inc"
	include "mygrid.inc"
	include "para.inc"
	include "ns.inc"
	include "sedi.inc"

	integer i,j,k
	double precision dv, da, aa, lid

	dv = 0.D0
        aa = 0.D0
        do k = -1, nnk+2
	do i = -1, nni+2
	   if (i .eq. -1) then
              if (k .eq. -1) then 
	         da = 0.25D0*(xp(i+2,1,k,1) - xp(i,1,k,1))
     <		            *(xp(i,1,k+2,3) - xp(i,1,k,3))
	      elseif (k .eq. nnk+1) then
		 da = 0.25D0*(xp(i+2,1,k,1) - xp(i,1,k,1))
     <		            *(xp(i,1,k,3) - xp(i,1,k-2,3))
              else
		 da = 0.25D0*(xp(i+2,1,k,1) - xp(i,1,k,1))
     <		            *(xp(i,1,k+1,3) - xp(i,1,k-1,3))
              endif
	   elseif (i .eq. nni+1) then
              if (k .eq. -1) then 
	         da = 0.25D0*(xp(i,1,k,1) - xp(i-2,1,k,1))
     <		            *(xp(i,1,k+2,3) - xp(i,1,k,3))
	      elseif (k .eq. nnk+1) then
		 da = 0.25D0*(xp(i,1,k,1) - xp(i-2,1,k,1))
     <		            *(xp(i,1,k,3) - xp(i,1,k-2,3))
              else
		 da = 0.25D0*(xp(i,1,k,1) - xp(i-2,1,k,1))
     <		            *(xp(i,1,k+1,3) - xp(i,1,k-1,3))
              endif
	   elseif (k .eq. -1)    then
              	 da = 0.25D0*(xp(i+1,1,k,1) - xp(i-1,1,k,1))
     <                      *(xp(i,1,k+2,3) - xp(i,1,k  ,3))
	   elseif (k .eq. nnk+1) then
	         da = 0.25D0*(xp(i+1,1,k,1) - xp(i-1,1,k,1))
     <                      *(xp(i,1,k  ,3) - xp(i,1,k-2,3))
	   else
	   da = 0.25D0*(xp(i+1,1,k,1) - xp(i-1,1,k,1))
     <                *(xp(i,1,k+1,3) - xp(i,1,k-1,3))
	   endif
	   dv = dv + 0.5D0*(     xp(i,1,k,2) +     xp(i,0,k,2)
     <                      -xp_old(i,1,k,2) - xp_old(i,0,k,2))*da
	   aa = aa + da
	enddo
	enddo

        lid = dv/aa/(dtime*ncheck)

	do k = -1, nnk+1
	do i = -1, nni+1
	   v_lid(i,k) = lid
	enddo
	enddo

	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine update_wbc

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "metric.inc"
	include "mygrid.inc"
	include "para.inc"
	include "ns.inc"
	include "sedi.inc"

	integer i, k

	do i = -1, nni
	do k = -1, nnk

c	   u(i, 0, k, 2) = u(i, 0, k, 2) + w_bed(i, k)
c	   u(i,-1, k, 2) = u(i,-1, k, 2) + w_bed(i, k)
c	   u(i, 1, k, 2) = u(i, 1, k, 2) + w_bed(i, k)
	   u(i,  0, k, 2) = u_bed(i,k,2)
	   u(i, -1, k, 2) = u_bed(i,k,2)
	enddo
	enddo

	return
	end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine zero_wbc

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "metric.inc"
	include "mygrid.inc"
	include "para.inc"
	include "ns.inc"
	include "sedi.inc"

	integer i, k

	do i = -1, nni
	do k = -1, nnk

	   u(i, 0, k, 2) = 0.D0*u(i, 0, k, 2)
	   u(i,-1, k, 2) = 0.D0*u(i, -1, k, 2)
c	   u(i, 1, k, 2) = u(i, 1, k, 2)
	enddo
	enddo

	return
	end
