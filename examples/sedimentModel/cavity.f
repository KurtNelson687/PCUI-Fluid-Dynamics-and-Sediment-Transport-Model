cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine grid

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "metric.inc"
	include "cavity.inc"
	INCLUDE "ns.inc"

	integer i, j, k

C	These are setting the size of the domain. The streching is done assuming length 1 then scaled by these
	bx = 0.5D0 !note: bx and several other varibles in the start of this subroutine are common mapped
	by = 1.D0
	bz = 0.5D0

C	These are flags indicate if the grid is going to be streched in the x,y, or z direction
	stretchx = 0
	stretchy = 0
	stretchz = 0

C	These are parameters for streching
	dm = 0.D0    
	am = 3.2D0 
	bm = 0.48D0
	cm = dm
     *	   + dlog( dcosh(am*(1.D0-bm)) / dcosh(am*bm) ) / am
     *	   - dlog( dcosh(am*bm) / dcosh(am*(1.D0-bm)) ) / am
	cm = 1.D0 / cm

	do k = -1, nnk+2
	do j = -1, nnj+2
	do i = -1, nni+2
	   jac(i,j,k) = 0.D0
	   xix(i,j,k) = 0.D0
	   xiy(i,j,k) = 0.D0
	   xiz(i,j,k) = 0.D0
	   etx(i,j,k) = 0.D0
	   ety(i,j,k) = 0.D0
	   etz(i,j,k) = 0.D0
	   ztx(i,j,k) = 0.D0
	   zty(i,j,k) = 0.D0
	   ztz(i,j,k) = 0.D0
	enddo
	enddo
	enddo

	call metric(1, nni, nnj, nnk,
     <          g11, g12, g13, g21, g22, g23, g31, g32, g33, gcc)!nni is the number of points on each processor in the i direction i.e ni/px. It is and input

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

	subroutine metric(level, ii, jj, kk, 
     <          q11, q12, q13, q21, q22, q23, q31, q32, q33, qcc)

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "metric.inc"
	include "cavity.inc"
	include "para.inc"

	integer level, ii, jj, kk, ii0, jj0, kk0

	double precision, dimension(-1:ii+2,-1:jj+2,-1:kk+2) ::
     <          q11, q12, q13, q21, q22, q23, q31, q32, q33, qcc
	double precision, dimension(-1:ii+2,-1:jj+2,-1:kk+2) ::
     <          x, y, z, xi, et, zt
     
	double precision jab
	double precision xxi, yxi, zxi, xet, yet, zet, xzt, yzt, zzt
	double precision xsx, xsy, xsz, esx, esy, esz, zsx, zsy, zsz

     
	integer i, j, k
	double precision faci, facj, fack

C	ii here is the number of points on processors in the x direction and npx is equal to the first coordinate for the processor mapping
	ii0 = npx * ii
	jj0 = npy * jj
	kk0 = npz * kk

C	These give dx, dy, and dz unstreched
	faci = 1.D0 / dble(ii*px) !dble converts ii*px to double percision. Px is the number of processors in the x direction
	facj = 1.D0 / dble(jj*py) !These are 1/(number of points in the given direction)
	fack = 1.D0 / dble(kk*pz)

	do k = -1, kk+2 !ii, jj, and kk
	do j = -1, jj+2
	do i = -1, ii+2
	   q11(i,j,k) = 0.D0 !I'm not sure this is true - These are components of the mesh skewness tensor (Gmn - see 2.55 in Yangs defense
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
	   xi(i,j,k) = ( dble(ii0 + i) - 0.5D0 ) * faci !This is setting up the unstretched grid spacing. These are the cell center locations at this point
	   et(i,j,k) = ( dble(jj0 + j) - 0.5D0 ) * facj
	   zt(i,j,k) = ( dble(kk0 + k) - 0.5D0 ) * fack
	   x(i,j,k) = 0.D0
	   y(i,j,k) = 0.D0
	   z(i,j,k) = 0.D0
	enddo
	enddo
	enddo

	call coord(ii, jj, kk, x, y, z, xi, et, zt) !subroutine that stretches the domain if activated

	if ( level .eq. 1 .and. newrun .eq. 1 ) then !level is hard coded to 1 when metric is called. newrun is set in io.f
	   write(200+myid) x, y, z
	endif

	if ( level .eq. 1 ) then
	   do k = -1, kk+2
	   do j = -1, jj+2
	   do i = -1, ii+2
	      xp(i,j,k,1) = x(i,j,k)!xp is a 4-dimensional array storing cell centered coordinates
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
C.... these are computing the derivatives of the cartesian components with respect to the general curvilinear components on the i faces
	   xxi = ( x(i+1,j,k) - x(i,j,k) )
	   yxi = ( y(i+1,j,k) - y(i,j,k) )
	   zxi = ( z(i+1,j,k) - z(i,j,k) )
	   xet = 0.25D0 * ( x(i,  j+1,k) - x(i,  j-1,k) !This is just the x coordinate of the face
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
	   jab = xxi * yet * zzt !This is the Jacobian on the i faces
     <         + xet * yzt * zxi
     <         + xzt * yxi * zet
     <         - xzt * yet * zxi
     <         - xet * yxi * zzt
     <         - xxi * yzt * zet
	   jab = 1.D0 / jab !This is the inverse Jacobian or cell volume
	   xsx = ( yet * zzt - yzt * zet )
	   esx = ( yzt * zxi - yxi * zzt )
	   zsx = ( yxi * zet - yet * zxi )
	   xsy = ( zet * xzt - zzt * xet )
	   esy = ( zzt * xxi - zxi * xzt )
	   zsy = ( zxi * xet - zet * xxi )
	   xsz = ( xet * yzt - xzt * yet )
	   esz = ( xzt * yxi - xxi * yzt )
	   zsz = ( xxi * yet - xet * yxi )
	   q11(i,j,k) = jab * ( xsx * xsx + xsy * xsy + xsz * xsz )!These are the mesh skewness tensor values on the i faces
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

	
	do k =  0, kk+1
	do j =  0, jj+1
	do i =  0, ii+1
	   yetjface(i,j,k) = (y(i,j+1,k)-y(i,j,k))
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
	      jac(i,j,k) = jab!This is the jacobian at the cell center
	   enddo
	   enddo
	   enddo
	   do k = -1, kk+2
	   do j = -1, jj+2
	   do i = -1, ii+2
	      rho_init(i,j,k) = 0.D0
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

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine qf2c(iif, jjf, kkf, iic, jjc, kkc, jaq, jas,
     <          q11, q12, q13, q21, q22, q23, q31, q32, q33, qcc,
     <          s11, s12, s13, s21, s22, s23, s31, s32, s33, scc)

	implicit none	

	integer iif, jjf, kkf, iic, jjc, kkc

	double precision, dimension(-1:iif+2,-1:jjf+2,-1:kkf+2) :: jaq,
     <          q11, q12, q13, q21, q22, q23, q31, q32, q33, qcc
     
	double precision, dimension(-1:iic+2,-1:jjc+2,-1:kkc+2) :: jas,
     <          s11, s12, s13, s21, s22, s23, s31, s32, s33, scc

	integer i, j, k, i0, j0, k0

	do k = -1, kkc+2
	do j = -1, jjc+2
	do i = -1, iic+2
	   jas(i,j,k) = 0.D0
	   s11(i,j,k) = 0.D0
	   s12(i,j,k) = 0.D0
	   s13(i,j,k) = 0.D0
	   s21(i,j,k) = 0.D0
	   s22(i,j,k) = 0.D0
	   s23(i,j,k) = 0.D0
	   s31(i,j,k) = 0.D0
	   s32(i,j,k) = 0.D0
	   s33(i,j,k) = 0.D0
	enddo
	enddo
	enddo

	do k = 1, kkc
	   k0 = 2 * k - 1
	do j = 1, jjc
	   j0 = 2 * j - 1
	do i = 0, iic
	   i0 = 2 * i - 1
	   jas(i,j,k) = 1.D0 / jaq(i0,  j0,  k0  )
     <	              + 1.D0 / jaq(i0+1,j0,  k0  )
     <	              + 1.D0 / jaq(i0  ,j0+1,k0  )
     <	              + 1.D0 / jaq(i0+1,j0+1,k0  )
     <	              + 1.D0 / jaq(i0,  j0,  k0+1)
     <	              + 1.D0 / jaq(i0+1,j0,  k0+1)
     <	              + 1.D0 / jaq(i0,  j0+1,k0+1)
     <	              + 1.D0 / jaq(i0+1,j0+1,k0+1)
	   jas(i,j,k) = 1.D0 / jas(i,j,k)
	enddo
	enddo
	enddo

	do k = 1, kkc
	   k0 = 2 * k - 1
	do j = 1, jjc
	   j0 = 2 * j - 1
	do i = 0, iic
	   i0 = 2 * i
	      s11(i,j,k) = 0.5D0 * ( q11(i0,j0,  k0  ) 
     <	                           + q11(i0,j0+1,k0  )
     <	                           + q11(i0,j0  ,k0+1)
     <	                           + q11(i0,j0+1,k0+1) )
	      s12(i,j,k) = 0.5D0 * ( q12(i0,j0,  k0  ) 
     <	                           + q12(i0,j0+1,k0  )
     <	                           + q12(i0,j0  ,k0+1)
     <	                           + q12(i0,j0+1,k0+1) )
	      s13(i,j,k) = 0.5D0 * ( q13(i0,j0,  k0  ) 
     <	                           + q13(i0,j0+1,k0  )
     <	                           + q13(i0,j0  ,k0+1)
     <	                           + q13(i0,j0+1,k0+1) )
	enddo
	enddo
	enddo

	do k = 1, kkc
	   k0 = 2 * k - 1
	do j = 0, jjc
	   j0 = 2 * j
	do i = 1, iic
	   i0 = 2 * i - 1
	      s21(i,j,k) = 0.5D0 * ( q21(i0,  j0,k0  ) 
     <	                           + q21(i0+1,j0,k0  )
     <	                           + q21(i0,  j0,k0+1)
     <	                           + q21(i0+1,j0,k0+1) )
	      s22(i,j,k) = 0.5D0 * ( q22(i0,  j0,k0  ) 
     <	                           + q22(i0+1,j0,k0  )
     <	                           + q22(i0,  j0,k0+1)
     <	                           + q22(i0+1,j0,k0+1) )
	      s23(i,j,k) = 0.5D0 * ( q23(i0,  j0,k0  ) 
     <	                           + q23(i0+1,j0,k0  )
     <	                           + q23(i0,  j0,k0+1)
     <	                           + q23(i0+1,j0,k0+1) )
	enddo
	enddo
	enddo

	do k = 0, kkc
	   k0 = 2 * k
	do j = 1, jjc
	   j0 = 2 * j - 1
	do i = 1, iic
	   i0 = 2 * i - 1
	      s31(i,j,k) = 0.5D0 * ( q31(i0,  j0,  k0) 
     <	                           + q31(i0+1,j0,  k0)
     <	                           + q31(i0,  j0+1,k0)
     <	                           + q31(i0+1,j0+1,k0) )
	      s32(i,j,k) = 0.5D0 * ( q32(i0,  j0,  k0) 
     <	                           + q32(i0+1,j0,  k0)
     <	                           + q32(i0,  j0+1,k0)
     <	                           + q32(i0+1,j0+1,k0) )
	      s33(i,j,k) = 0.5D0 * ( q33(i0,  j0,  k0) 
     <	                           + q33(i0+1,j0,  k0)
     <	                           + q33(i0,  j0+1,k0)
     <	                           + q33(i0+1,j0+1,k0) )
	enddo
	enddo
	enddo

	do k = 1, kkc
	do j = 1, jjc
	do i = 1, iic
	   scc(i,j,k) = s11(i,j,k) + s11(i-1,j,  k  )
     <	              + s22(i,j,k) + s22(i,  j-1,k  )
     <	              + s33(i,j,k) + s33(i,  j,  k-1)
	enddo 
	enddo
	enddo

	return
	end
     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine mapx(ii, jj, kk, x, xi)

	include "size.inc"
	include "mpif.h"
	include "cavity.inc"
	include "mpi.inc"

	integer ii, jj, kk !number of points on each processor in the i, j, k direction

	double precision, dimension(-1:ii+2,-1:jj+2,-1:kk+2) :: x, xi

	integer i, j, k
	
	if ( stretchx .eq. 1 ) then

	do k = -1, kk+2
	do j = -1, jj+2
	do i = -1, ii+2
	   x(i,j,k) = dm * xi(i,j,k) 
     <              + dlog( dcosh( am * ( xi(i,j,k) - bm ) )
     <                    / dcosh( am *               bm   ) ) 
     <              / am
     <              - dlog( dcosh( am * ( xi(i,j,k) - 1.D0+bm ) )
     <                    / dcosh( am * (             1.D0-bm ) ) )
     <              / am
	   x(i,j,k) = cm * x(i,j,k)
	enddo
	enddo
	enddo

CBCBCBCBCBC
	IF ( N_WEST .EQ. MPI_PROC_NULL ) THEN
	DO K = -1, KK+2
	DO J = -1, JJ+2
	   X( 0,J,K) = - X(1,J,K) !I believe these are just setting ghost cells
	   X(-1,J,K) = - X(2,J,K)
	ENDDO
	ENDDO
	ENDIF
	IF ( N_EAST .EQ. MPI_PROC_NULL ) THEN
	DO K = -1, KK+2
	DO J = -1, JJ+2
	   X(II+1,J,K) = 2.D0 - X(II,  J,K)
	   X(II+2,J,K) = 2.D0 - X(II-1,J,K)
	ENDDO
	ENDDO
	ENDIF
CBCBCBCBCBC

	else

	do k = -1, kk+2
	do j = -1, jj+2
	do i = -1, ii+2
	   x(i,j,k) = xi(i,j,k)!This sets the unstreched coordinate to the final x value
	enddo
	enddo
	enddo

	endif

	return
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine mapy(ii, jj, kk, x, xi)

	include "size.inc"
	include "mpif.h"
	include "cavity.inc"
	include "mpi.inc"

	integer nd, ii, jj, kk

	double precision, dimension(-1:ii+2,-1:jj+2,-1:kk+2) :: x, xi

	integer i, j, k
	
	if ( stretchy .eq. 1 ) then

	do k = -1, kk+2
	do j = -1, jj+2
	do i = -1, ii+2
	   x(i,j,k) = dm * xi(i,j,k) 
     <              + dlog( dcosh( am * ( xi(i,j,k) - bm ) )
     <                    / dcosh( am *               bm   ) ) 
     <              / am
     <              - dlog( dcosh( am * ( xi(i,j,k) - 1.D0+bm ) )
     <                    / dcosh( am * (             1.D0-bm ) ) )
     <              / am
	   x(i,j,k) = cm * x(i,j,k)
	enddo
	enddo
	enddo

CBCBCBCBCBC
	IF ( N_SUTH .EQ. MPI_PROC_NULL ) THEN
	DO K = -1, KK+2
	DO I = -1, II+2
	   X(I, 0,K) = - X(I,1,K)
	   X(I,-1,K) = - X(I,2,K)
	ENDDO
	ENDDO
	ENDIF
	IF ( N_NRTH .EQ. MPI_PROC_NULL ) THEN
	DO K = -1, KK+2
	DO I = -1, II+2
	   X(I,JJ+1,K) = 2.D0 - X(I,JJ,  K)
	   X(I,JJ+2,K) = 2.D0 - X(I,JJ-1,K)
	ENDDO
	ENDDO
	ENDIF
CBCBCBCBCBC

	else

	do k = -1, kk+2
	do j = -1, jj+2
	do i = -1, ii+2
	   x(i,j,k) = xi(i,j,k)
	enddo
	enddo
	enddo

	endif

	return
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine mapz(ii, jj, kk, x, xi)

	include "size.inc"
	include "mpif.h"
	include "cavity.inc"
	include "mpi.inc"

	integer nd, ii, jj, kk

	double precision, dimension(-1:ii+2,-1:jj+2,-1:kk+2) :: x, xi

	integer i, j, k
	
	if ( stretchz .eq. 1 ) then

	do k = -1, kk+2
	do j = -1, jj+2
	do i = -1, ii+2
	   x(i,j,k) = dm * xi(i,j,k) 
     <              + dlog( dcosh( am * ( xi(i,j,k) - bm ) )
     <                    / dcosh( am *               bm   ) ) 
     <              / am
     <              - dlog( dcosh( am * ( xi(i,j,k) - 1.D0+bm ) )
     <                    / dcosh( am * (             1.D0-bm ) ) )
     <              / am
	   x(i,j,k) = cm * x(i,j,k)
	enddo
	enddo
	enddo

CBCBCBCBCBC
	IF ( N_BACK .EQ. MPI_PROC_NULL ) THEN
	DO J = -1, JJ+2
	DO I = -1, II+2
	   X(I,J, 0) = - X(I,J,1)
	   X(I,J,-1) = - X(I,J,2)
	ENDDO
	ENDDO
	ENDIF
	IF ( N_FRNT .EQ. MPI_PROC_NULL ) THEN
	DO J = -1, JJ+2
	DO I = -1, II+2
	   X(I,J,KK+1) = 2.D0 - X(I,J,KK  )
	   X(I,J,KK+2) = 2.D0 - X(I,J,KK-1)
	ENDDO
	ENDDO
	ENDIF
CBCBCBCBCBC

	else

	do k = -1, kk+2
	do j = -1, jj+2
	do i = -1, ii+2
	   x(i,j,k) = xi(i,j,k)
	enddo
	enddo
	enddo

	endif

	return
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine coord(ii, jj, kk, x, y, z, xi, et, zt)

	include "size.inc"
	include "cavity.inc"

	integer ii, jj, kk !These are the number of points on each processors in the i, j, and k direction

	double precision, dimension(-1:ii+2,-1:jj+2,-1:kk+2) :: 
     <          x, y, z

	double precision, dimension(-1:ii+2,-1:jj+2,-1:kk+2) :: 
     <          xi, et, zt

	integer i, j, k
	
	call mapx(ii, jj, kk, x, xi) !maps the unstretched x to stretched x if stretch = 1
	call mapy(ii, jj, kk, y, et)
	call mapz(ii, jj, kk, z, zt)

C	This is simply scaling the domain by what is indicated by bx, :by, and bz
	do k = -1, kk+2
	do j = -1, jj+2
	do i = -1, ii+2
	   x(i,j,k) = bx * x(i,j,k)
	   y(i,j,k) = by * y(i,j,k)
	   z(i,j,k) = bz * z(i,j,k)
	enddo
	enddo
	enddo

	return
	end
