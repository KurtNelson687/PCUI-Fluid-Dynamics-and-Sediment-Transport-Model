!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     
!c     Subroutine: particle_transport
!c     -------------------------
!c     
!c     Goncalo Gil 
!c     
!c     Stanford University, Environmental Fluid Mechanics Laboratory
!c
!c     Lagrangian Particle tracking code. Integration in time with a 
!c     fourth-order accurate Runge-Kutta algorithm. Second-order trilinear 
!c     interpolation is performed to obtain the particle's velocity at its 
!c     current location in the fixed grid. The particle is located via an 
!c     O(log N) kd-tree search algorithm available at:
!c     https://github.com/jmhodges/kdtree2
!c

      subroutine particle_transport

      implicit none 

      include 'para.inc'
      include "mpi.inc"

      if (istep.eq.1) call init_particles
      call getU
      if (myid.eq.0) call rk4
            
      end subroutine

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine init_particles

      include 'size.inc'
      include 'mpif.h'
      include 'metric.inc'
      include 'particles.inc'
      include 'cavity.inc'
      include "mpi.inc"
      include 'para.inc'

      integer :: i, j, k, ntPart, n, ci
     
      ipstep = 1

      xxL(1) = 0.D0                !(xp(0,0,0,:)+xp(1,1,1,:))/2.D0
      xxL(2) = -by
      xxL(3) = 0.D0
      xxR(1) = xxL(1)+bx
      xxR(2) = 0.D0
      xxR(3) = xxL(3)+bz

      call getX
      call getU

      if (myid.eq.0) then

         if (newrun.eq.1) then
           xPartB = 0
           xPartBT = 0
           xPartS = 0
           ntPart = 0

           do k = lbk, ubk
           do j = lbj, ubj
           do i = lbi, ubi
             ntPart = ntPart + 1         
             xPart(ntPart,:) = xxp(i,j,k,:)
             uPart(ntPart,:) = uu(i,j,k,:)
           end do 
           end do
           end do

           if (ntPart.ne.nPart) then
             print *, 'Number of particles mismatch', ntPart, nPart
             stop
           end if

           call output_particles

         else
           call input_continue_particles
           call output_particles
         end if

      end if

      end subroutine

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
!c      rk4 advances the particles in time using Runge-Kutta. The timestep
!c      of integration for rk4 is twice the timestep for the NS solver since
!c      rk4 uses intermediate values. 
!c
!c

      subroutine rk4

      include "size.inc"
      include "mpif.h"
      include "ns.inc"
      include "particles.inc"
      include "para.inc"
      include "mpi.inc"

      integer                        :: n,ci
              
      ! If the index of the NS solver is the same as the index of the particle
      ! solver where i_ns = 2*i_particle-1 then the particle solver is at a
      ! full timestep, otherwise it is at a half timestep.

      if (istep.eq.2*ipstep-1) then 
         if (istep.gt.1) then 
                      
            do n = 1,nPart
               do ci = 1,3
                  xPartC(n,ci) = xPart(n,ci)+k3(n,ci)
                  call boundary_adjustment(n,ci)   
               end do
            end do

            call interpolate3D(uPart,xPartC,xxp,uu,ni+2,nj+2,nk+2,
     <                         nPart,xxL,xxR,xPartB,xPartBT,xPartS)

            k4 = 2.D0*dtime*uPart

!     ---- Advance particle in time using RK4 and adjust for boundary crossings
            xPart = xPart + (1.D0/6.D0)*(k1 + 2.D0*(k2+k3)+k4)
         end if

!        Location of particle at current position in RK4         
         xPartC = xPart

         do n = 1,nPart
            do ci = 1,3
               call boundary_adjustment(n,ci)
            end do
         end do

         xPart = xPartC

         call interpolate3D(uPart,xPartC,xxp,uu,ni+2,nj+2,nk+2,nPart,
     <                      xxL,xxR,xPartB,xPartBT,xPartS)
        
         k1 = 2.D0*dtime*uPart  

!     ---- Store results
         if (mod(istep+1,nsave) .eq. 0) then
            print *,'Output particles.'
            call output_particles          
         endif
         if (mod(istep+1,ncont) .eq. 0) then
            print *,'Output continue particles.'
            call output_continue_particles          
         endif

      else

         do n = 1,nPart
            do ci = 1,3
               xPartC(n,ci) = xPart(n,ci)+0.5D0*k1(n,ci)
               call boundary_adjustment(n,ci)   
            end do
         end do

         call interpolate3D(uPart,xPartC,xxp,uu,ni+2,nj+2,nk+2,
     <                      nPart,xxL,xxR,xPartB,xPartBT,xPartS)

         k2 = 2.D0*dtime*uPart

         do n = 1,nPart
            do ci = 1,3
               xPartC(n,ci) = xPart(n,ci)+0.5D0*k2(n,ci)
               call boundary_adjustment(n,ci)
            end do
         end do

         call interpolate3D(uPart,xPartC,xxp,uu,ni+2,nj+2,
     <                      nk+2,nPart,xxL,xxR,xPartB,xPartBT,xPartS)     

         k3 = 2.D0*dtime*uPart

         ipstep = ipstep + 1
      end if

!10      format(4f15.12)
!11      format(4i15)
!100                 format('---... Boundary adjustment (n,ci,s) = ',
!    <                        '(', i2, ',',i2,',',i2,') ...---')

      end subroutine

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine boundary_adjustment(n,ci)

      include 'size.inc'
      include 'ns.inc'
      include 'particles.inc'
      include "metric.inc"
      include "para.inc"
      include 'mpi.inc'

      double precision, dimension(3) :: x,xW,bN,rN
      double precision               :: xL,xR,D
      integer                        :: n,ci

!     Calculate appropriate boundary value
      if (ci.eq.1.or.ci.eq.3) then
        xL = xxL(ci)
        xR = xxR(ci)
        D = xR-xL
      else
        call get_bottom(xPartC(n,1),xL)
        xR = xxR(ci)
      end if

!     Calculate appropriate boundary normal unit vector 
!     rN is vector from x to the wall in the direction of bN.
      bN = 0.D0
      rN = 0.D0
      if (ci.eq.1.or.ci.eq.3) then
        if (xPartC(n,ci).lt.xL) then
          bN(ci) = 1
          rN(ci) = bN(ci)*(xL-xPartC(n,ci))
!         print *, bN, rN
        elseif(xPartC(n,ci).gt.xR) then    
          bN(ci) = 1
          rN(ci) = bN(ci)*(xR-xPartC(n,ci))
!         print *, bN, rN
        else
          return
        end if
      else
        if (xPartC(n,ci).lt.xL) then
          call get_boundary_normal(xPartC(n,1),bN)
          xW(1) = xPartC(n,1)
          xW(2) = xL
          xW(3) = xPartC(n,3)
          rN = bN*dot_product(xW-xPartC(n,:),bN) 
          print *, bN, rN, xW-xPartC(n,:)
        elseif(xPartC(n,ci).gt.xR) then    
          bN(ci) = 1
          rN(ci) = bN(ci)*(xR-xPartC(n,ci))
!         print *, bn, rN
        else
          return
        end if
      end if

!     ADJUST PARTICLE LOCATION
!     If periodic apply symmetry condition
      if (periods(ci).eqv..true.) then
        print *, '---...Symmetry...---'
        xPartC(n,ci) = xPartC(n,ci) + bN(ci)*D

      else
!     If it hits the wall it bounces off
        print *, '---...',n,' Hit-Wall...---'
        xPartC(n,:) = xPartC(n,:) + 2*rN
      end if

      end subroutine boundary_adjustment

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!     Get vertical location of bottom given horizontal particle position

      subroutine get_bottom(x1,x2)

      include 'cavity.inc'

      double precision :: x1,x2,xc,zc,xend,zend,r,s
      xc = 1.675D0
      zc = 2.44D0
      xend = 2.537D0
      zend = -0.4335D0
      r = 3.D0
      s = 0.3D0

      if (x1.lt.xc) then
        x2 = -by
      elseif (x1.ge.xc.and.x1.lt.xend) then
        x2 = -sqrt(r**2-(x1-xc)**2) + zc
      else
        x2 = s*x1 + (zend - s*xend)
      end if

      end subroutine

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!     Get boundary normal unit vector given horizontal particle position

      subroutine get_boundary_normal(x1,bN)

      double precision               :: x1,xc,xend,r,s,theta
      double precision, dimension(3) :: bN

      xc = 1.675D0
      xend = 2.537D0
      r = 3.D0
      s = 0.3D0

      if (x1.lt.xc) then
        bN(2) = 1
      elseif (x1.ge.xc.and.x1.lt.xend) then
        s = (x1-xc)/sqrt(r**2-(x1-xc)**2)
        theta =  atan(s)
        bN(1) = -sin(theta)
        bN(2) =  cos(theta)
      else
        theta =  atan(s)
        bN(1) = -sin(theta)
        bN(2) =  cos(theta)
      end if

      end subroutine
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!     Obtain u sub-matrices from each process and assemble them into a single global matrix
!     in order to perform the particle tracking on process 0
!

      subroutine getU

      include 'size.inc'
      include "mpif.h"
      include 'ns.inc'
      include 'particles.inc'
      include 'mpi.inc'

      integer :: nprocs,i,j,k,n
      integer :: resizedrecvsubarray,sendsubarray,recvsubarray
      integer(kind=mpi_address_kind) :: start, extent
      integer, dimension(3) :: sizes, subsizes, starts, coords
      double precision, dimension(-1:nni+2,-1:nnj+2,-1:nnk+2) :: xpsend
      double precision, dimension(1:ni,1:nj,1:nk) :: xprecv 
      integer, dimension(:), allocatable :: counts, disps
      double precision :: dp

      call MPI_COMM_SIZE( MPI_COMM_WORLD, nprocs, ierr )
!      call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
      call MPI_CART_COORDS( comm3d, myid, 3, coords, ierr )

      ! ---- Create subarray structure for the send subarray
      sizes(1) = nni+4
      sizes(2) = nnj+4
      sizes(3) = nnk+4
      subsizes(1) = nni
      subsizes(2) = nnj
      subsizes(3) = nnk
      starts(1) = 2
      starts(2) = 2
      starts(3) = 2

      call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,
     <                              MPI_ORDER_FORTRAN,
     <                              MPI_DOUBLE_PRECISION,sendsubarray,
     <                              ierr)
      call MPI_TYPE_COMMIT(sendsubarray,ierr)

      call MPI_Barrier(MPI_COMM_WORLD, ierr)

      ! ---- Create subarray structure for the receive subarray
      if (myid.eq.0) then

         sizes(1) = ni
         sizes(2) = nj
         sizes(3) = nk
         subsizes(1) = nni
         subsizes(2) = nnj
         subsizes(3) = nnk
         starts(1) = 0
         starts(2) = 0
         starts(3) = 0

         call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,
     <                              MPI_ORDER_FORTRAN,
     <                              MPI_DOUBLE_PRECISION,recvsubarray,
     <                              ierr)
         call MPI_TYPE_COMMIT(recvsubarray,ierr)
         extent = sizeof(dp)
         start = 0
         call MPI_TYPE_CREATE_RESIZED(recvsubarray, start, extent, 
     <                             resizedrecvsubarray,ierr)                               
         call MPI_TYPE_COMMIT(resizedrecvsubarray,ierr)
 
         allocate(counts(nprocs),disps(nprocs))
         counts = 1
         do i = 1, px
            do j = 1, py
               do k = 1, pz
                  disps((i-1)*py*pz+(j-1)*pz+k) = nni*(i-1) +
     <                              px*py*(k-1)*nni*nnj*nnk +
     <                              px*(j-1)*nni*nnj 
               end do
            end do
         end do     
      end if

      ! ---- Assemble local matrix into a global matrix
      do n = 1, 3
         xpsend = u(:,:,:,n)           
         call MPI_GATHERV(xpsend,1,sendsubarray,xprecv,counts,
     <                 disps,resizedrecvsubarray,0,MPI_COMM_WORLD,ierr)
         uuu(:,:,:,n) = xprecv
      end do

      call MPI_TYPE_FREE(sendsubarray,ierr)
      if (myid.eq.0) then
         call MPI_TYPE_FREE(recvsubarray,ierr)
         call MPI_TYPE_FREE(resizedrecvsubarray,ierr)
         deallocate(disps,counts)


         ! No-slip
         uu = 0.D0
         ! Insides
         uu(2:ni+1,2:nj+1,2:nk+1,:) = uuu
         ! Rigid Lid - free slip
         uu(1:ni+2,nj+2,1:nk+2,1) = uu(1:ni+2,nj+1,1:nk+2,1)
         uu(1:ni+2,nj+2,1:nk+2,3) = uu(1:ni+2,nj+1,1:nk+2,3)
         ! Free slip - x
         uu(1,2:nj+1,2:nk+1,:)    = uu(2,2:nj+1,2:nk+1,:) 
         uu(ni+2,2:nj+1,2:nk+1,:) = uu(ni+1,2:nj+1,2:nk+1,:)
         ! Symmetry - z
         uu(2:ni+1,2:nj+1,1,:)= (uu(2:ni+1,2:nj+1,2,:) +
     <                           uu(2:ni+1,2:nj+1,nk+1,:))/2.D0
         uu(2:ni+1,2:nj+1,nk+2,:) = uu(2:ni+1,2:nj+1,1,:)


c$$$         if(ipstep.eq.1) then
c$$$            open(123, file='output_UVW_test.dat.500',form='unformatted',
c$$$     >           status='unknown')
c$$$         else
c$$$            open(123, file='output_UVW_test.dat.500',form='unformatted',
c$$$     >           status='old', position='append')
c$$$         endif
c$$$         write(123) uu
c$$$         close(unit = 123)

      end if
      
      end subroutine

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!     Obtain x sub-matrices from each process and assemble them into a single global matrix
!     in order to perform the particle tracking on process 0
!

      subroutine getX

      include 'size.inc'
      include "mpif.h"
      include 'metric.inc'
      include 'particles.inc'
      include 'mpi.inc'

      integer :: nprocs,i,j,k,n
      integer :: resizedrecvsubarray,sendsubarray,recvsubarray
      integer(kind=mpi_address_kind) :: start, extent
      integer, dimension(3) :: sizes, subsizes, starts, coords
      double precision, dimension(-1:nni+2,-1:nnj+2,-1:nnk+2) :: xpsend
      double precision, dimension(1:ni,1:nj,1:nk) :: xprecv 
      integer, dimension(:), allocatable :: counts, disps
      double precision :: dp, bot

      call MPI_COMM_SIZE( MPI_COMM_WORLD, nprocs, ierr )
!      call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
      call MPI_CART_COORDS( comm3d, myid, 3, coords, ierr )

      ! ---- Create subarray structure for the send subarray
      sizes(1) = nni+4
      sizes(2) = nnj+4
      sizes(3) = nnk+4
      subsizes(1) = nni
      subsizes(2) = nnj
      subsizes(3) = nnk
      starts(1) = 2
      starts(2) = 2
      starts(3) = 2

      call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,
     <                              MPI_ORDER_FORTRAN,
     <                              MPI_DOUBLE_PRECISION,sendsubarray,
     <                              ierr)
      call MPI_TYPE_COMMIT(sendsubarray,ierr)

      call MPI_Barrier(MPI_COMM_WORLD, ierr)

      ! ---- Create subarray structure for the receive subarray
      if (myid.eq.0) then

         sizes(1) = ni
         sizes(2) = nj
         sizes(3) = nk
         subsizes(1) = nni
         subsizes(2) = nnj
         subsizes(3) = nnk
         starts(1) = 0
         starts(2) = 0
         starts(3) = 0

         call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,
     <                              MPI_ORDER_FORTRAN,
     <                              MPI_DOUBLE_PRECISION,recvsubarray,
     <                              ierr)
         call MPI_TYPE_COMMIT(recvsubarray,ierr)
         extent = sizeof(dp)
         start = 0
         call MPI_TYPE_CREATE_RESIZED(recvsubarray, start, extent, 
     <                             resizedrecvsubarray,ierr)                               
         call MPI_TYPE_COMMIT(resizedrecvsubarray,ierr)
 
         allocate(counts(nprocs),disps(nprocs))
         counts = 1
         do i = 1, px
            do j = 1, py
               do k = 1, pz
                  disps((i-1)*py*pz+(j-1)*pz+k) = nni*(i-1) +
     <                              px*py*(k-1)*nni*nnj*nnk +
     <                              px*(j-1)*nni*nnj 
               end do
            end do
         end do
      
      end if

      ! ---- Assemble local matrix into a global matrix
      do n = 1, 3
         xpsend = xp(:,:,:,n)           
         call MPI_GATHERV(xpsend,1,sendsubarray,xprecv,counts,
     <                 disps,resizedrecvsubarray,0,MPI_COMM_WORLD,ierr)
         xxxp(:,:,:,n) = xprecv
      end do   

      call MPI_TYPE_FREE(sendsubarray,ierr)

      if (myid.eq.0) then
         call MPI_TYPE_FREE(recvsubarray,ierr)
         call MPI_TYPE_FREE(resizedrecvsubarray,ierr)
         deallocate(disps,counts)

!        Interior points
         do k = 2, nk+1
         do j = 2, nj+1
         do i = 2, ni+1
            xxp(i,j,k,1) = xxxp(i-1,j-1,k-1,1)
            xxp(i,j,k,2) = xxxp(i-1,j-1,k-1,2)
            xxp(i,j,k,3) = xxxp(i-1,j-1,k-1,3)
         end do
         end do
         end do

!        Boundary points
!        ---Left and right
         xxp(1,:,:,1)      = xxL(1)
         xxp(ni+2,:,:,1)   = xxR(1)
         xxp(1,:,:,2)      = xxp(2,:,:,2) 
         xxp(ni+2,:,:,2)   = xxp(ni+1,:,:,2)
         xxp(1,:,:,3)      = xxp(2,:,:,3) 
         xxp(ni+2,:,:,3)   = xxp(ni+1,:,:,3)

!        ---Top and bottom
         xxp(:,1,:,1)    = xxp(:,2,:,1)
         xxp(:,nj+2,:,1) = xxp(:,nj+1,:,1)
         do i = 1, ni+2
            call get_bottom(xxp(i,2,2,1),bot)
            xxp(i,1,:,2) = bot 
         end do
         xxp(:,nj+2,:,2) = xxR(2)
         xxp(:,1,:,3)    = xxp(:,2,:,3)
         xxp(:,nj+2,:,3) = xxp(:,nj+1,:,3)

!        ---Front and back
         xxp(:,:,1,1)      = xxp(:,:,2,1)
         xxp(:,:,nk+2,1)   = xxp(:,:,nk+1,1)
         xxp(:,:,1,2)      = xxp(:,:,2,2)
         xxp(:,:,nk+2,2)   = xxp(:,:,nk+1,2) 
         xxp(:,:,1,3)      = xxL(3)
         xxp(:,:,nk+2,3)   = xxR(3)

!        print *, xxp(:,1,1,2)

!        OLD         
!        xxp(1,:,:,1)      = xxL(1)
!        xxp(ni+2,:,:,1)   = xxR(1)
!        do i = 2, ni+1
!           xxp(i,:,:,1) = xxxp(i-1,1,1,1)
!        end do
!        xxp(:,1,:,2)      = xxL(2)
!        xxp(:,nj+2,:,2)   = xxR(2)
!        do j = 2, nj+1
!           xxp(:,j,:,2) = xxxp(1,j-1,1,2)
!        end do
!        xxp(:,:,1,3)      = xxL(3)
!        xxp(:,:,nk+2,3)   = xxR(3)
!        do k = 2, nk+1
!           xxp(:,:,k,3) = xxxp(1,1,k-1,3)
!        end do

         open(123, file='output_XYZ.dat',form='unformatted',
     >           status='unknown')
         write(123) xxp
         close(unit = 123)

c$$$         if(ipstep.eq.1) then
c$$$            open(123, file='output_XYZ_test.dat.600',form='unformatted',
c$$$     >           status='unknown')
c$$$         else
c$$$            open(123, file='output_XYZ_test.dat.600',form='unformatted',
c$$$     >           status='old', position='append')
c$$$         endif
c$$$         write(123) xxp
c$$$         close(unit = 123)

      end if

      end subroutine

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      subroutine output_particles
      
      include 'size.inc'
      include 'particles.inc'      

      if(ipstep.eq.1) then
         open(123, file='output_xPart.dat', form='unformatted',
     >          status='unknown')
      else
         open(123, file='output_xPart.dat', form='unformatted',
     >          status='old', position='append')
      endif
      
      write(123) xPart     

      close(unit = 123)

      if(ipstep.eq.1) then
         open(123, file='output_uPart.dat', form='unformatted',
     >          status='unknown')
      else
         open(123, file='output_uPart.dat', form='unformatted',
     >          status='old', position='append')
      endif
      
      write(123) uPart      

      close(unit = 123)

      if(ipstep.eq.1) then
         open(123, file='output_xPartS.dat', form='unformatted',
     >          status='unknown')
      else
         open(123, file='output_xPartS.dat', form='unformatted',
     >          status='old', position='append')
      endif
      
      write(123) xPartS

      close(unit = 123)

      end subroutine

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine output_continue_particles

      include 'size.inc'
      include 'particles.inc'

      open(123, file='continue_particles.dat', form='unformatted',
     <     status='unknown')
      write(123) xPart
      write(123) uPart
      close(123)

      end subroutine

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine input_continue_particles

      include 'size.inc'
      include 'particles.inc'

      open(123, file='continue_particles.dat', form='unformatted',
     <     status='unknown')
      read(123) xPart
      read(123) uPart
      close(123)

      end subroutine

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!c$$$      if (myid.eq.0) then
!c$$$         do n = 1,3
!c$$$            print *, ' Global (n = ', n,')'
!c$$$            print *, '==================='
!c$$$            write (*,'(32f5.2)') (uu(i,:,1,n),i=1,ni )
!c$$$         end do
!c$$$      end if
!c$$$
!c$$$      call MPI_Barrier(MPI_COMM_WORLD, ierr)
