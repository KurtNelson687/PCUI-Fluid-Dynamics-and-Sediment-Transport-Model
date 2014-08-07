!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     
!c     Subroutine: interpolate3d
!c     -------------------------
!c     
!c     Goncalo Gil 
!c     
!c     Stanford University, Environmental Fluid Mechanics Laboratory
!c
!c     Variables of interest:
!c     ----------------------
!c
!c     uparticle ----> ARRAY(3,Np) where Np = number of particles and
!c                     the array contains u,v,w velocities at the particles.
!c     xparticle ----> ARRAY(3,Np) where Np = number of particles and
!c                     the array contains x,y,z positions of the particles.
!c     xp -----------> ARRAY(il,jl,kl,3) contains grid point locations
!c     u ------------> ARRAY(il,jl,kl,3) contains velocity field at xp.
!c     il -----------> INTEGER contains number of gridpoints in x direction.
!c     jl -----------> INTEGER contains number of gridpoints in y direction.
!c     kl -----------> INTEGER contains number of gridpoints in z direction.
!C     Np -----------> INTEGER contains number of partciles.
!c

      subroutine interpolate3D(uparticle,xparticle,xp,u,il,jl,kl,Np,xxL,xxR,xPartB,xPartBT,xPartS)

      use kdtree2_module
      
      implicit none
      
      integer            :: corner(3,Np),il,jl,kl,ci,i,j,k,num,idx,Np,neighbors(8,3)
      double precision   :: xp(il,jl,kl,3),x(3),xL(3),xR(3),r,s,t, &
                            u(il,jl,kl,3),xparticle(Np,3),uparticle(Np,3),xxL(3),xxR(3)

      type(kdtree2_result)   :: results(1)
      type(kdtree2), pointer :: tree
      double precision       :: xList(3,(il)*(jl)*(kl))
      integer                :: xListInd(3,(il)*(jl)*(kl))
      integer,dimension(Np,3):: xPartS, xPartB, xPartBT

      ! Create a list of points from 3D grid to feed into kdtree locator
      num = 0
      do k = 1, kl
         do j = 1, jl
            do i = 1, il
               num = num + 1
                 do ci = 1, 3
                  xList(ci,num) = xp(i,j,k,ci)
               end do
               xListInd(1,num) = i
               xListInd(2,num) = j
               xListInd(3,num) = k
            end do
         end do
      end do

      ! ---- Create tree structure used to find particle's nearest neighbor on the grid
      tree => kdtree2_create(xList,sort=.false.,rearrange=.false.)

      ! ---- Loop through all particles, locate, find neighbors, interpolate.
      xL = abs(xxL)
      xR = abs(xxR)

      do num = 1, Np

!         if (xPartB(num,1).eq.0.and.xPartB(num,2).eq.0.and.xPartB(num,3).eq.0) then

            ! ---- Locate particle
            call kdtree2_n_nearest_brute_force(tree,xparticle(num,:),1,results)
            
            idx = results(1)%idx

            i = xListInd(1,idx)
            j = xListInd(2,idx)
            k = xListInd(3,idx)

            ! if (num.ge.Np-3) then
            !    print *, '#############################################'
            !    write(*,10) num,i,j,k,xparticle(num,1),uparticle(num,1)
            !    print *, '#############################################'
            ! end if
10          format(4i6,2f15.12)

            ! Make sure that closest grid point to particle is to the left of particle
            ! Example: if particle is located half way between two grid points the 
            !          grid point to the left is selected so that the correct grid 
            !          cell is chosen.

            x = xparticle(num,:)

            if (x(1)-xp(i,j,k,1).lt.0.D0) then
               i = i - 1
            end if
            if (x(2)-xp(i,j,k,2).lt.0.D0) then
               j = j - 1
            end if
            if (x(3)-xp(i,j,k,3).lt.0.D0) then
               k = k - 1
            end if

!            if (istep.eq.343) then
 !           end if

            corner(1,num) = i
            corner(2,num) = j
            corner(3,num) = k
            
            call findNeighbors(neighbors,corner(:,num))
            r = xparticle(num,1)-xp(neighbors(1,1),neighbors(1,2),neighbors(1,3),1)
            s = xparticle(num,2)-xp(neighbors(1,1),neighbors(1,2),neighbors(1,3),2)
            t = xparticle(num,3)-xp(neighbors(1,1),neighbors(1,2),neighbors(1,3),3)
          
!           if (num.eq.Np) then
!              print *, '*1**>',x(1),i,j,k,xp(i,j,k,1),u(i,j,k,1)
!              print *, xparticle(nP,1),uparticle(nP,1)
!              do ci = 1,8
!                 print *, u(neighbors(ci,1),neighbors(ci,2),neighbors(ci,3),1)
!              end do
!           end if
            do ci = 1,3
               call interpolate(uparticle(num,ci),neighbors,num,ci, &
                       xp,u,Np,il,jl,kl,r,s,t)          
            end do
!           if (num.eq.Np) then
!              print *, '*2**>',x(1),i,j,k,xp(i,j,k,1),u(i,j,k,1)
!              print *, xparticle(nP,1),uparticle(nP,1)
!           end if
        
 !        end if
      end do
      
      call kdtree2_destroy(tree)        

      end subroutine interpolate3D

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     
!c     Subroutine: findNeighbors
!c     -------------------------
!c     
!c     This subroutine finds the neighbors of particle given the 1st corner 
!c     grid point.
!c

      subroutine findNeighbors(neighbors,corner)
      implicit none
      
      integer          :: i,j,k,neighbors(8,3),corner(3)

      i = corner(1)
      j = corner(2)
      k = corner(3)

      neighbors(1,1) = i
      neighbors(1,2) = j
      neighbors(1,3) = k
      
      neighbors(2,1) = i+1
      neighbors(2,2) = j
      neighbors(2,3) = k

      neighbors(3,1) = i+1
      neighbors(3,2) = j
      neighbors(3,3) = k+1

      neighbors(4,1) = i
      neighbors(4,2) = j  
      neighbors(4,3) = k+1

      neighbors(5,1) = i
      neighbors(5,2) = j+1  
      neighbors(5,3) = k

      neighbors(6,1) = i+1
      neighbors(6,2) = j+1
      neighbors(6,3) = k  

      neighbors(7,1) = i+1
      neighbors(7,2) = j+1
      neighbors(7,3) = k+1  

      neighbors(8,1) = i  
      neighbors(8,2) = j+1
      neighbors(8,3) = k+1
      
      end subroutine findNeighbors

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
!c     Subroutine: interpolate
!c     -----------------------
!c
!c     This subroutine takes an nx3 array with an index list containing all
!c     the neighbors for the current particle. It then uses that array to 
!c     retrieve all the neighbor's xyz positions and respective scalar to be
!c     interpolated. n is the number of neighbors.
!c

      subroutine interpolate(u_interp,ni,n,ci,xp,u,Np,il,jl,kl,r,s,t)
      implicit none

      integer          :: ci,il,jl,kl,ni(8,3),n,Np
      double precision :: u_interp,xp(il,jl,kl,3), &
                          u(il,jl,kl,3), & 
                          u01,u001,x01,x001,u101,x101,u00,u000,x00,x000, &
                          u100,x100,u11,u011,x11,x011,u111,x111,u10, &
                          u010,x10,x010,u110,x110,u1,x1, &
                          u0,x0,r,s,t
      
!CCCCCCCCCCCCCCCCCCCCCCCC
!c     PANEL 1 x-dir    c
!CCCCCCCCCCCCCCCCCCCCCCCC

!CCCCCCCCCCCCCCCCCCCCCCCC
!c        Top           c
!CCCCCCCCCCCCCCCCCCCCCCCC
         
         x001=xp(ni(4,1),ni(4,2),ni(4,3),1)
         u001=u(ni(4,1),ni(4,2),ni(4,3),ci)
         x101=xp(ni(3,1),ni(3,2),ni(3,3),1)
         u101=u(ni(3,1),ni(3,2),ni(3,3),ci)
         call single_interp(u01,u001,u101,x001,x101,r)
                  
!CCCCCCCCCCCCCCCCCCCCCCC
!c        Down          c
!CCCCCCCCCCCCCCCCCCCCCCCC

         x000=xp(ni(1,1),ni(1,2),ni(1,3),1)
         u000=u(ni(1,1),ni(1,2),ni(1,3),ci)
         x100=xp(ni(2,1),ni(2,2),ni(2,3),1)
         u100=u(ni(2,1),ni(2,2),ni(2,3),ci)
         call single_interp(u00,u000,u100,x000,x100,r)

!CCCCCCCCCCCCCCCCCCCCCCCC
!c     PANEL 2 x-dir    c
!CCCCCCCCCCCCCCCCCCCCCCCC

!CCCCCCCCCCCCCCCCCCCCCCCC
!c        Top           c
!CCCCCCCCCCCCCCCCCCCCCCCC
         
         x011=xp(ni(8,1),ni(8,2),ni(8,3),1)
         u011=u(ni(8,1),ni(8,2),ni(8,3),ci)
         x111=xp(ni(7,1),ni(7,2),ni(7,3),1)
         u111=u(ni(7,1),ni(7,2),ni(7,3),ci)
         call single_interp(u11,u011,u111,x011,x111,r)
         
!CCCCCCCCCCCCCCCCCCCCCCCC
!c        Down          c
!CCCCCCCCCCCCCCCCCCCCCCCC

         x010=xp(ni(5,1),ni(5,2),ni(5,3),1)
         u010=u(ni(5,1),ni(5,2),ni(5,3),ci)
         x110=xp(ni(6,1),ni(6,2),ni(6,3),1)
         u110=u(ni(6,1),ni(6,2),ni(6,3),ci)
         call single_interp(u10,u000,u100,x000,x100,r)
         
!CCCCCCCCCCCCCCCCCCCCCCCC
!c     PANEL 3 y-dir    c
!CCCCCCCCCCCCCCCCCCCCCCCC

!CCCCCCCCCCCCCCCCCCCCCCCC
!c        Top           c
!CCCCCCCCCCCCCCCCCCCCCCCC
         
         x01=xp(ni(1,1),ni(1,2),ni(1,3),2)
         x11=xp(ni(5,1),ni(5,2),ni(5,3),2)
         call single_interp(u1,u01,u11,x01,x11,s)
         
!CCCCCCCCCCCCCCCCCCCCCCCC
!c        Down          c
!CCCCCCCCCCCCCCCCCCCCCCCC

         x00=xp(ni(1,1),ni(1,2),ni(1,3),2)
         x10=xp(ni(5,1),ni(5,2),ni(5,3),2)
         call single_interp(u0,u00,u10,x00,x10,s)
         
!CCCCCCCCCCCCCCCCCCCCCCCC
!c     Line z-dir       c
!CCCCCCCCCCCCCCCCCCCCCCCC

         x0=xp(ni(1,1),ni(1,2),ni(1,3),3)
         x1=xp(ni(8,1),ni(8,2),ni(8,3),3)
         call single_interp(u_interp,u0,u1,x0,x1,t)

      end subroutine interpolate

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
!c     Subroutine: single_interp
!c     -----------------------
!c
!c     This subroutine performs a single linear interpolation..
!c

      subroutine single_interp(y,y0,y1,x0,x1,d)
      double precision y,y0,y1,x,x0,x1,d
         
      if (x0.eq.x1) then
         if(y0.gt.y1) then
            y = y0
         else
            y = y1
         end if
      else
         y=y0+d*(y1-y0)/(x1-x0)
      end if
      
      end subroutine single_interp
