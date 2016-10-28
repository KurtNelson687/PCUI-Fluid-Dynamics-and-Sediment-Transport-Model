ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine randomU
        
        include "size.inc"
        include "ns.inc"
        include "mpif.h"
        include "mpi.inc"

        integer i, j, k, m
        double precision fac, randx, randy, old, seed  
        seed = 11.D0
        old = seed
        
        fac = 0.1
        do m = 1, 3
        do k = 1, nnk
        do j = 1, nnj
        do i = 1, nni
           randx = mod((57*old+1), 256.D0)
           randy = mod((57*randx+1), 256.D0)+myid
           randy = randy/256.D0
           old = randy
           if ( m .ne. 1 ) then
c              u(i, j, k, m) = (1.D0+(randy-0.5D0))*0.02D0
              u(i, j, k, m) = (randy-0.5D0)*0.1D0
           else
              u(i, j, k, m) = (1.D0+0.02D0*(randy-0.5D0+0.1D0*myid
     <                         ))*u(i, j, k, m)
           endif
        enddo
        enddo
        enddo
        enddo
      
        return
      
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine TimeAvrg
      
      include "size.inc"
      include "ns.inc"
      include "para.inc"

      integer i, j, k, m

      if (time .eq. 0.D0) then
        do m = 1, 3
        do k = 1, nnk
        do j = 1, nnj
        do i = 1, nni
c           uavg(i, j, k, m) = 0.D0
        enddo
        enddo
        enddo
        enddo     
      else
        
        do m = 1, 3
        do k = 1, nnk
        do j = 1, nnj
        do i = 1, nni
c           uavg(i, j, k, m) = (uavg(i, j, k, m)*(time-dtime)
c     <                       +    u(i, j, k, m)*dtime)/time
        enddo
        enddo
        enddo
        enddo
         
      endif
      
      return
      
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine randomSD
        
        include "size.inc"
        include "ns.inc"
        include "sedi.inc"
        include "mpif.h"
        include "mpi.inc"

        integer i, j, k, m
        double precision fac, randx, randy, old, seed  
        seed = 11.D0
        old = seed
        
        fac = 0.1
        
        do k = 1, nnk
        do j = 1, nnj
        do i = 1, nni
           randx = mod((57*old+1), 256.D0)
           randy = mod((57*randx+1), 256.D0)+myid
           randy = randy/256.D0
           old = randy
           
           sedc(i, j, k) = (1.D0+0.4D0*(randy-0.5D0))
     <                    *sedc(i, j, k)
           
        enddo
        enddo
        enddo
       
      
        return
      
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine smoothSD
        
        include "size.inc"
        include "ns.inc"
        include "sedi.inc"
        include "mpif.h"
        include "mpi.inc"

        integer i, j, k
        double precision sedc2(1:nni,1:nnj,1:nnk)
        
        do k = 1, nnk
        do j = 1, nnj
        do i = 1, nni
           sedc2(i,j,k) = sedc(i,j,k)
           
        enddo
        enddo
        enddo
       
      
        return
      
      end
