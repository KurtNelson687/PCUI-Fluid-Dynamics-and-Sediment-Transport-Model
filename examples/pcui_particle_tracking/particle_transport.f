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

      integer :: i, j, k, ntPart, n, ci
     
      ipstep = 1

      xxL = 0.D0                !(xp(0,0,0,:)+xp(1,1,1,:))/2.D0
      xxR(1) = xxL(1)+bx
      xxR(2) = 0.D0
      xxR(3) = xxL(3)+bz

      call getX
      call getU

      if (myid.eq.0) then

         xPartB = 0
         xPartBT = 0
         xPartS = 0
         xPartW = 0.D0
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

      integer                        :: n, ci
              
      ! If the index of the NS solver is the same as the index of the particle
      ! solver where i_ns = 2*i_particle-1 then the particle solver is at a
      ! full timestep, otherwise it is at a half timestep.

      if (istep.eq.2*ipstep-1) then 
         if (istep.gt.1) then 
                      
            do n = 1,nPart
               do ci = 1,3
                  xPartC(n,ci) = xPart(n,ci)+k3(n,ci)
                  call boundary_adjustment(xPartC(n,:),n,ci)   
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
               call boundary_adjustment(xPartC(n,:),n,ci)
            end do
         end do

         call interpolate3D(uPart,xPartC,xxp,uu,ni+2,nj+2,nk+2,nPart,
     <                      xxL,xxR,xPartB,xPartBT,xPartS)
        
         k1 = 2.D0*dtime*uPart  

!     ---- Store results
         if (mod(istep+1,nsave) .eq. 0) then
            print *,'Output particles.'
            call output_particles          
         endif
      else

         do n = 1,nPart
            do ci = 1,3
               xPartC(n,ci) = xPart(n,ci)+0.5D0*k1(n,ci)
               call boundary_adjustment(xPartC(n,:),n,ci)   
            end do
         end do

         call interpolate3D(uPart,xPartC,xxp,uu,ni+2,nj+2,nk+2,
     <                      nPart,xxL,xxR,xPartB,xPartBT,xPartS)

         k2 = 2.D0*dtime*uPart

         do n = 1,nPart
            do ci = 1,3
               xPartC(n,ci) = xPart(n,ci)+0.5D0*k2(n,ci)
               call boundary_adjustment(xPartC(n,:),n,ci)
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

      subroutine boundary_adjustment(x,n,ci)

      include 'size.inc'
      include 'ns.inc'
      include 'particles.inc'
      include "metric.inc"
      include "para.inc"
      include 'mpi.inc'

      double precision, dimension(3) :: x,xW,bN
      double precision               :: xL,xR,D
      integer                        :: n,ci

!     Calculate appropriate boundary value
      if (ci.eq.1.or.ci.eq.3)  then
        xL = xxL(ci)
        xR = xxR(ci)
        D = xR-xL
      else
        xR = xxR(ci)
        call get_bottom(x(1),xL)
      end if

!     Calculate appropriate boundary normal unit vector 
!     bN points INTO the domain
      bN = 0.D0
      if (ci.eq.1.or.ci.eq.3)  then
        if (x(ci).lt.xL) then
          bN(n,ci) = 1
        elseif(x(ci).gt.xR) then    
          bN(n,ci) = -1
        else
          return
        end if
      else
        if (x(ci).lt.xL) then
          call get_boundary_normal(x(1),bN)
        elseif(x(ci).gt.xR) then    
          bN(n,ci) = -1
        else
          return
        end if
      end if

!     ADJUST PARTICLE LOCATION
!     If periodic apply symmetry condition
      if (periods(ci).eqv..true.) then
        print *, '---...Symmetry...---'
        xPartC(n,ci) = xPartC(n,ci) + bN(n,ci)*D(ci)

      else
!     If it hits the wall it bounces off
        print *, '---...',n,' Hit-Wall...---'
      end if

     
!     OLD!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     ---- Symmetry/boundary check
      if (x.lt.xL(ci)) then
         xPartB(n,ci) = -1
         xPartBT(n,ci) = xPartBT(n,ci) + 1               
      elseif(x.gt.xR(ci)) then    
         xPartB(n,ci) = 1
         xPartBT(n,ci) = xPartBT(n,ci) + 1
      else
         return
      end if
     
!     ---- If periodic it applies symmetry condition
      if (periods(ci).eqv..true.) then
         print *, '---...Symmetry...---'
         xPartS(n,ci) = xPartS(n,ci) + xPartB(n,ci)
         print *, '--> ',xPart(n,ci),x
         xPart(n,ci) = xPartI(n,ci) - xPartB(n,ci)*D(ci)
         print *, xPart(n,ci)
         print *, uPart(n,ci)
      else
!     ---- If it hits the wall it bounces off
         print *, '---...',n,' Hit-Wall...---'
!     ---- Unit vector in the direction of velocity vector
         uV = uPartI(n,:)
         uN = uPartI(n,:)/sqrt(uPartI(n,1)**2.D0 + 
     <        uPartI(n,2)**2.D0+uPartI(n,3)**2.D0)
         print *, Un
!     ---- Vector normal to boundary is just the boundary vector
         bN = xPartB(n,:)
         print *, bN
!     ---- Calculate distance from the particle to the boundary along the velocity vector
         if (xPartB(n,ci).eq.-1) then
            r1 = dot_product(xL - xPartI(n,:),bN)
         else if (xPartB(n,ci).eq.1) then
            r1 = dot_product(xR - xPartI(n,:),bN)
         end if
         r2 = dot_product(uV,bN)

         if (r2.eq.0.D0.and.r1.ne.0.D0) then
            print *, r1, r2
            print *, 'The line is outside the plane and is parallel 
     <                to it'
            xPart(n,ci) = xPartI(n,ci) + 20*r1
            return
         elseif (r2.eq.0.D0.and.r1.eq.0.D0) then
           print *, r1, r2
           print *, 'The line is inside the plane and is parallel to it'
           return
         end if                

         twall = r1/r2
!     ---- Bounce it off the wall                   
         xPartW(n,:) = xPartI(n,:) + twall*uV
         distV = xPart(n,:)-xPartW(n,:)

         if (ci.eq.1) then
            xPart(n,1) = xPartW(n,1)-distV(1)
            xPart(n,2) = xPartW(n,2)+distV(2)
            xPart(n,3) = xPartW(n,3)+distV(3)
         elseif (ci.eq.2) then
            xPart(n,1) = xPartW(n,1)+distV(1)
            xPart(n,2) = xPartW(n,2)-distV(2)
            xPart(n,3) = xPartW(n,3)+distV(3)
         elseif (ci.eq.3) then
            xPart(n,1) = xPartW(n,1)+distV(1)
            xPart(n,2) = xPartW(n,2)+distV(2)
            xPart(n,3) = xPartW(n,3)-distV(3)
         end if
      end if

      xPartB(n,ci) = 0

      end subroutine boundary_adjustment

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
         ! Rigid Lid
         uu(1:ni+2,nj+2,1:nk+2,1) = uu(1:ni+2,nj+1,1:nk+2,1)
         uu(1:ni+2,nj+2,1:nk+2,3) = uu(1:ni+2,nj+1,1:nk+2,3)
         ! Symmetry - x
         uu(1,2:nj+1,2:nk+1,:)= (uu(ni+1,2:nj+1,2:nk+1,:) +
     <                           uu(2,2:nj+1,2:nk+1,:))/2.D0
         uu(nj+2,2:nj+1,2:nk+1,:) = uu(1,2:nj+1,2:nk+1,:)

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

         xxp(1,:,:,1)      = xxL(1)
         xxp(ni+2,:,:,1)   = xxR(1)
         do i = 2, ni+1
            xxp(i,:,:,1) = xxxp(i-1,1,1,1)
         end do
         xxp(:,1,:,2)      = xxL(2)
         xxp(:,nj+2,:,2)   = xxR(2)
         do j = 2, nj+1
            xxp(:,j,:,2) = xxxp(1,j-1,1,2)
         end do
         xxp(:,:,1,3)      = xxL(3)
         xxp(:,:,nk+2,3)   = xxR(3)
         do k = 2, nk+1
            xxp(:,:,k,3) = xxxp(1,1,k-1,3)
         end do

         open(123, file='output_XYZ.dat',form='unformatted',
     >           status='unknown')
         write(123) xxxp
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


!c$$$      if (myid.eq.0) then
!c$$$         do n = 1,3
!c$$$            print *, ' Global (n = ', n,')'
!c$$$            print *, '==================='
!c$$$            write (*,'(32f5.2)') (uu(i,:,1,n),i=1,ni )
!c$$$         end do
!c$$$      end if
!c$$$
!c$$$      call MPI_Barrier(MPI_COMM_WORLD, ierr)
