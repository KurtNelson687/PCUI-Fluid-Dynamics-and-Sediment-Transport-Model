cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
        subroutine updateDepth2

        include "size.inc"
        include "mpif.h"
        include "mpi.inc"
        include "ns.inc"
        include "metric.inc"
        include "para.inc"
        include "sedi.inc"
        include "eddy.inc"
        include "depth.inc"
        
        double precision temp
        double precision delf, phie, tot

        integer i, k        
        double precision tan_r, sx, sz
        double precision L, T, Ks
        
        call collectSedi

c       Source term due to suspended load
        do k = 1, nnk
        do i = 1, nni
           drhs(i,k) = -sedv(i,k)*1.D0/jad(i,k)          
c           drhs(i,k) = 0.D0          
        enddo
        enddo

        
        do k = -1, nnk+2
        do i = -1, nni+2
           h(i,k) = depth0(i,k)-depth(i,k)
        enddo
        enddo



c       Source term due to bed load
        
        if (bedload .eq. 1) then
        do k = 1, nnk
        do i = ius, iue
           if(qbxi(i,k) .ge. 0.D0) then
             qbxi(i,k) = xixd(i,k) * 
     <	        ( 0.75 D0 * qbi(i,k) + 0.375D0 * qbi(i+1,k) 
     <          - 0.125D0 * qbi(i-1,k) )
     <	              + xizd(i,k) *
     <	        ( 0.75 D0 * qbk(i,k) + 0.375D0 * qbk(i+1,k) 
     <          - 0.125D0 * qbk(i-1,k) )
            else
               qbxi(i,k) = xixd(i,k) * 
     <	        ( 0.75 D0 * qbi(i+1,k) + 0.375D0 * qbi(i,k) 
     <          - 0.125D0 * qbi(i+2,k) )
     <	              + xizd(i,k) *
     <	        ( 0.75 D0 * qbk(i+1,k) + 0.375D0 * qbk(i,k) 
     <          - 0.125D0 * qbk(i+2,k) )
              
            endif 
           
        enddo
        enddo        

        
	do k = kus, kue
	do i = 1, nni

           if ( qbzt(i,k) .ge. 0.D0 ) then
              qbzt(i,k) = ztxd(i,k) * 
     <	        ( 0.75 D0 * qbi(i,k) + 0.375D0 * qbi(i,k+1) 
     <          - 0.125D0 * qbi(i,k-1) )
     <	              + ztzd(i,k) *
     <	        ( 0.75 D0 * qbk(i,k) + 0.375D0 * qbk(i,k+1) 
     <          - 0.125D0 * qbk(i,k-1) )
           else
	       qbzt(i,k) = ztxd(i,k) * 
     <	        ( 0.75 D0 * qbi(i,k+1) + 0.375D0 * qbi(i,k) 
     <          - 0.125D0 * qbi(i,k+2) )
     <	              + ztzd(i,k) *
     <	        ( 0.75 D0 * qbk(i,k+1) + 0.375D0 * qbk(i,k) 
     <          - 0.125D0 * qbk(i,k+2) )
	    endif
	enddo
	enddo
        
         temp = 1.D0/(1.D0-por)

        do k = 1, nnk
        do i = 1, nni
           drhs(i, k)=  drhs(i,k) + temp*
     <                            ( -qbxi(i,k) + qbxi(i-1,k) 
     <                              -qbzt(i,k) + qbzt(i,k-1))
        enddo
        enddo
        endif

        if (d_diff .eq. 1) then
c....... obtain diffusion coefficient

           tan_r = dtan(resp*pi/180)
           L = 0.15
           T = dsqrt(g*(spwght-1)/spwght/0.03)
           Ks = L**2.D0/T


           do k = 1, nnk
           do i = 0, nni+1
              sx = -(depth(i+1,k)-depth(i-1,k))
     <             /(xp(i+1,1,k,1)-xp(i-1,1,k,1))
              if (dabs(sx) .ge. tan_r) then
                 af_x(i,k) = 0.0005*((dabs(sx)-tan_r)/tan_r)**3.D0
              else
                 af_x(i,k) = 0.000
              endif
*(dexp((sx/tan_r)**2.D0)-1.D0)
              
           enddo
           enddo

           do k = 0, nnk+1
           do i = 1, nni
              sz = -(depth(i,k+1)-depth(i,k-1))
     <           /(xp(i,1,k+1,3)-xp(i,1,k-1,3))
              if (dabs(sz) .ge. tan_r) then
                 af_z(i,k) = 0.0005*((dabs(sz)-tan_r)/tan_r)**3.D0
              else
                 af_z(i,k) = 0.000
              endif
c              af_z(i,k) = 0.0005*dabs(sz)**2.D0
*(dexp((sz/tan_r)**2.D0)-1.D0)
              
           enddo
           enddo

c...... Cross diffusion terms at step n-1 from Crank-Nicolson
        
           do k = 1, nnk
           do i = 1, nni

             drhs(i,k) = drhs(i,k)
     <          + 0.5D0*(af_z(i,k) + af_z(i+1,k)) *
     <            d13(i  ,k) * ( h(i,  k+1) - h(i  ,k-1)
     <                         + h(i+1,k+1) - h(i+1,k-1) )
     <          - 0.5D0*(af_z(i,k) + af_z(i-1,k)) *
     <            d13(i-1,k) * ( h(i  ,k+1) - h(i  ,k-1)
     <                         + h(i-1,k+1) - h(i-1,k-1) )

             drhs(i,k) = drhs(i,k)
     <         + 0.5D0*(af_x(i,k) + af_x(i,k+1))  *
     <             d31(i,k  ) * ( h(i+1,k  ) - h(i-1,k  )
     <                          + h(i+1,k+1) - h(i-1,k+1) )
     <         - 0.5D0*(af_x(i,k) + af_x(i,k-1))  *
     <             d31(i,k-1) * ( h(i+1,k  ) - h(i-1,k  )
     <                          + h(i+1,k-1) - h(i-1,k-1) )

            enddo
            enddo

            do k = 0, nnk
            do i = 0, nni

               drhs(i,k) = drhs(i,k)
     <          + 0.5D0*(af_x(i,k) + af_x(i+1,k)) *
     <            d11(i  ,k  ) * ( h(i+1,k  ) - h(i  ,k  ) )
     <          + 0.5D0*(af_x(i,k) + af_x(i-1,k)) *
     <            d11(i-1,k  ) * ( h(i-1,k  ) - h(i  ,k  ) )
     <          + 0.5D0*(af_z(i,k) + af_z(i,k+1)) *
     <            d33(i  ,k  ) * ( h(i,  k+1) - h(i  ,k  ) )
     <          + 0.5D0*(af_z(i,k) + af_z(i,k-1)) *
     <            d33(i,  k-1) * ( h(i,  k-1) - h(i  ,k  ) )    
            enddo
            enddo

            do k = 1, nnk
            do i = 1, nni
               drhs(i,k) = drhs(i,k)*jad(i,k)*dtime
            enddo      
            enddo     
           
            call depth_solve


        else
            do k = 1, nnk
            do i = 1, nni
               drhs(i,k) = drhs(i,k)*jad(i,k)*dtime
            enddo      
            enddo     
         endif
        
        do k = 1, nnk
        do i = 1, nni
           temp = depth(i,k)
           h(i,k) = h(i,k) + drhs(i,k)
           depth(i,k) = depth0(i,k) - h(i,k)
           u_bed(i, k, 2) = -(depth(i, k)-temp)/dtime
        enddo      
        enddo     
          
        call depth_bc

        call depth_exchange
                
        tot = 0.D0

        do k = 1, nnk
        do i = 1, nni
           
           tot = tot + u_bed(i,k,2)
c           u_bed(i,k,2) = 0.D0
        enddo
        enddo

        temp = 1.D0/(nnk*nni)

        do k = -1, nnk+2
        do i = -1, nni+2
c           write(*,*) tot
           v_lid(i,k) = tot*temp
           
        enddo
        enddo
              
        return
        end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine depth_solve
        
        include "size.inc"
        include "mpif.h"
        include "mpi.inc"
        include "ns.inc"
        include "metric.inc"
        include "para.inc"
        include "eddy.inc"
        include "sedi.inc"
        include "depth.inc"

        double precision ax(nnk,0:nni+1), bx(nnk,0:nni+1),
     &                   cx(nnk,0:nni+1), fx(nnk,0:nni+1)
        double precision az(nni,0:nnk+1), bz(nni,0:nnk+1),
     &                   cz(nni,0:nnk+1), fz(nni,0:nnk+1)
        
        integer i, j, k
        double precision coef, B, temp        
 
         coef = 0.5D0 * dtime

c....... solve for I-direction

        do k = 1, nnk
        do i = 1, nni
           ax(k,i) = -( 0.5D0*(af_x(i,k) + af_x(i-1,k)) )*
     <               coef * jad(i,k) * d11(i-1,k)
           cx(k,i) = -( 0.5D0*(af_x(i,k) + af_x(i+1,k)) )*
     <               coef * jad(i,k) * d11(i,  k)
           bx(k,i) = 1.D0 - ax(k,i) - cx(k,i)
           fx(k,i) = drhs(i,k)
           
        
        enddo

        if ( n_west .eq. MPI_PROC_NULL ) then
          
              drhs(0,k)= 0.D0
              ax(k,0) = 0.D0 
              bx(k,0) = 1.D0
              cx(k,0) =-1.D0
              fx(k,0) = drhs(0,k)
   
         endif 

        if ( n_east .eq. MPI_PROC_NULL ) then
              
              drhs(nni+1,k)= 0.D0
              ax(k,nni+1) = 1.D0 
              bx(k,nni+1) =-1.D0
              cx(k,nni+1) = 0.D0
              fx(k,nni+1) = drhs(nni+1,k)    
        endif
        
        if ( periodic .eq. 1 ) then
        call trip(ax, bx, cx, fx, nnk, 1, nni, n_west, n_east )
        else
        call trid(ax, bx, cx, fx, nnk, 1, nni, n_west, n_east )
        endif

        do i = 1, nni
           drhs(i,k) = fx(k,i)
        enddo
        
        enddo

        

c...... solve for K-direction

        do i = 1, nni
        do k = 1, nnk
           az(i,k) = -( 0.5D0*(af_z(i,k) + af_z(i,k-1)) )*
     <               coef * jad(i,k) * d33(i,k-1)
           cz(i,k) = -( 0.5D0*(af_z(i,k) + af_z(i,k+1)) )*
     <               coef * jad(i,k) * d33(i  ,k)
           bz(i,k) = 1.D0 - az(i,k) - cz(i,k)
           fz(i,k) = drhs(i,k)
           
        enddo

        if ( n_back .eq. MPI_PROC_NULL ) then          
              drhs(i,0) = 0.D0
              az(i,0) = 0.D0 
              bz(i,0) = 1.D0
              cz(i,0) =-1.D0
              fz(i,0) = drhs(i,0)    
        endif 

        if ( n_frnt .eq. MPI_PROC_NULL ) then           
              
              drhs(i,nnk+1) = 0.D0
              az(i,nnk+1) = 1.D0 
              bz(i,nnk+1) =-1.D0
              cz(i,nnk+1) = 0.D0
              fz(i,nnk+1) = drhs(i,nnk+1)
                
        endif
        
        if ( periodicZ .eq. 1 ) then
        call trip(az, bz, cz, fz, nni, 1, nnk, n_back, n_frnt )
        else
        call trid(az, bz, cz, fz, nni, 1, nnk, n_back, n_frnt )
        endif

        do k = 1, nnk
           drhs(i,k) = fz(i,k)
        enddo
                
        enddo

        return
        end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine depth_bc
      
      include "size.inc"
      include "mpif.h"
      include "para.inc"            
      include "mpi.inc"
      include "ns.inc"
      include "sedi.inc"
      include "eddy.inc"
      include "metric.inc"
      

      integer i, j, k
      double precision B, temp, ap, len
      

      if (btmBuffer .eq. 1) then
c         if ( n_west .eq. MPI_PROC_NULL) then
             do k = 1, nnk
c             ap = depth(10,k)-depth0(5,k)
             ap = depth(5,k)-0.2D0   
             len = xp(5,1,k,1)-xp(3,1,k,1)
             do i = 1, 5
                if (i .ge. 3) then
                   depth(i,k) = 0.2D0 + ap*
     <             dsin(pi/2.D0*(xp(i,1,k,1)-xp(3,1,k,1))/len)
                else
c                   depth(i,k) = depth0(i,k)
                   depth(i,k) = 0.2D0
                endif
             enddo
             enddo
c         endif

c         if ( n_east .eq. MPI_PROC_NULL) then
             do k = 1, nnk
             ap = 0.2D0-depth(nni-5,k)
             len = xp(nni-3,1,k,1)-xp(nni-5,1,k,1)
             do i = nni-5, nni
                if (i .le. nni-3) then
                   depth(i,k) = depth(nni-5,k) + ap*
     <             dsin(pi/2.D0*(xp(i,1,k,1)-xp(nni-5,1,k,1))/len)
                else
c                   depth(i,k) = depth0(i,k)
                   depth(i,k) = 0.2D0
                 endif
             enddo
             enddo
c         endif
         
c         if ( n_back .eq. MPI_PROC_NULL) then
             do i = 1, nni
c             ap = depth(i,10)-depth0(i,5)
             ap = depth(i,5)-0.2
             len = xp(i,1,5,3)-xp(i,1, 3,3)
             do k = 1, 5
                if (k .ge. 3) then
                   depth(i,k) = 0.2 + ap*
     <             dsin(pi/2.D0*(xp(i,1,k,3)-xp(i,1,3,3))/len)
                else
c                   depth(i,k) = depth0(i,k)
                    depth(i,k) = 0.2
                endif
             enddo
             enddo
c         endif

c         if ( n_frnt .eq. MPI_PROC_NULL) then
             do i = 1, nni
c             ap = depth(i,nnk-5)-depth0(i,nnk-10)
              ap = 0.2-depth(i,nnk-5)
              len = xp(i,1,nnk-3,3)-xp(i,1,nnk-5,3)
              do k = nnk-5, nnk
                 if (k .le. nnk-3) then
                    depth(i,k) = depth(i,nnk-5) + ap*
     <              dsin(pi/2.D0*(xp(i,1,k,3)-xp(i,1,nnk-5,3))/len)
                 else
c                   depth(i,k) = depth0(i,k)
                   depth(i,k) = 0.2
                endif
             enddo
             enddo
c         endif
      endif



c      if ( n_west .eq. MPI_PROC_NULL) then
         do k = -1, nnk+2 
            depth( 0,k) = 0.1
cdepth(1,k)
            depth(-1,k) = depth(0,k)
         enddo
c      endif

c      if ( n_east .eq. MPI_PROC_NULL) then
         do k = -1, nnk+2 
            depth(nni+1,k) = 0.1
cdepth(nni,k)
            depth(nni+2,k) = depth(nni+1,k)
         enddo
c      endif
             

      if ( n_back .eq. MPI_PROC_NULL) then
         do i = -1, nni+2
            depth(i,0) = 0.1
cdepth(i,1)
            depth(i,-1) = depth(i,0)
         enddo         
      endif

      if ( n_frnt .eq. MPI_PROC_NULL) then
         do i = -1, nni+2
            depth(i,nnk+1) = 0.1
cdepth(i,nnk)
            depth(i,nnk+2) = depth(i,nnk+1)
         enddo
      endif
         

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine depth_bc2
      
      include "size.inc"
      include "mpif.h"
      include "para.inc"            
      include "mpi.inc"
      include "ns.inc"
      include "sedi.inc"
      include "eddy.inc"
      include "metric.inc"
      

      integer i, j, k
      double precision B, temp, ap, len
      

      if (btmBuffer .eq. 1) then
c         if ( n_west .eq. MPI_PROC_NULL) then
             do k = 1, nnk
                ap = depth(5,k)-depth0(5,k)   
                len = xp(5,1,k,1)-xp(3,1,k,1)
             do i = 1, 5
                if (i .ge. 3) then
                   depth(i,k) = depth0(5,k) + ap*
     <             dsin(pi/2.D0*(xp(i,1,k,1)-xp(3,1,k,1))/len)
                else
                   depth(i,k) = depth0(i,k)
          
                endif
             enddo
             enddo
c         endif

c         if ( n_east .eq. MPI_PROC_NULL) then
             do k = 1, nnk
                ap = depth0(nni-5,k)-depth(nni-5,k)
                len = xp(nni-3,1,k,1)-xp(nni-5,1,k,1)
             do i = nni-5, nni
                if (i .le. nni-3) then
                   depth(i,k) = depth(nni-5,k) + ap*
     <             dsin(pi/2.D0*(xp(i,1,k,1)-xp(nni-5,1,k,1))/len)
                else
c                   depth(i,k) = depth0(i,k)
                   depth(i,k) = 0.2D0
                 endif
             enddo
             enddo
c         endif
         
c         if ( n_back .eq. MPI_PROC_NULL) then
             do i = 1, nni
             ap = depth(i,10)-depth0(i,5)
c             ap = depth(i,5)-0.2
             len = xp(i,1,5,3)-xp(i,1, 3,3)
             do k = 1, 5
                if (k .ge. 3) then
                   depth(i,k) = depth0(i,k) + ap*
     <             dsin(pi/2.D0*(xp(i,1,k,3)-xp(i,1,3,3))/len)
                else
                   depth(i,k) = depth0(i,k)
c                    depth(i,k) = 0.2
                endif
             enddo
             enddo
c         endif

c         if ( n_frnt .eq. MPI_PROC_NULL) then
             do i = 1, nni
c             ap = depth(i,nnk-5)-depth0(i,nnk-10)
              ap = depth0(i,k)-depth(i,nnk-5)
              len = xp(i,1,nnk-3,3)-xp(i,1,nnk-5,3)
              do k = nnk-5, nnk
                 if (k .le. nnk-3) then
                    depth(i,k) = depth(i,nnk-5) + ap*
     <              dsin(pi/2.D0*(xp(i,1,k,3)-xp(i,1,nnk-5,3))/len)
                 else
                   depth(i,k) = depth0(i,k)
c                   depth(i,k) = 0.2
                endif
             enddo
             enddo
c         endif
      endif



c      if ( n_west .eq. MPI_PROC_NULL) then
c         do k = -1, nnk+2 
c            depth( 0,k) = 0.1
cdepth(1,k)
c            depth(-1,k) = depth(0,k)
c         enddo
c      endif

c      if ( n_east .eq. MPI_PROC_NULL) then
c         do k = -1, nnk+2 
c            depth(nni+1,k) = 0.1
cdepth(nni,k)
c            depth(nni+2,k) = depth(nni+1,k)
c         enddo
c      endif
             

c      if ( n_back .eq. MPI_PROC_NULL) then
c         do i = -1, nni+2
c            depth(i,0) = 0.1
cdepth(i,1)
c            depth(i,-1) = depth(i,0)
c         enddo         
c      endif

c      if ( n_frnt .eq. MPI_PROC_NULL) then
c         do i = -1, nni+2
c            depth(i,nnk+1) = 0.1
cdepth(i,nnk)
c            depth(i,nnk+2) = depth(i,nnk+1)
c         enddo
c      endif
         

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine updateDepth3
        
        include "size.inc"
        include "metric.inc"
        include "para.inc"
        include "sedi.inc"
        include "ns.inc"
        include "mpif.h"
        include "mpi.inc"
        
        integer i,j,k
        double precision temp, tot, lx, lz, tmp, dx, dz

        call collectSedi
        tot = 0.D0
        lx = dabs(0.5*(xp(5,1,5,1)-xp(3,1,5,1)))
        lz = dabs(0.5*(xp(5,1,5,1)-xp(5,1,3,1)))

        tmp = dtime/(1.D0-por)
        do k = 1, nnk
        do i = 1, nni
           dx = xp(i,1,k,1)-xp(i-1,1,k,1)
           dz = xp(i,1,k,3)-xp(i,1,k-1,3)
           temp = depth(i,k)
           depth(i,k) = depth(i,k) + sedv(i,k)

           if (bedload .eq. 1) then
              depth(i,k) = depth(i,k)+tmp*( 
     <                     1.D0/dx*(qbxi(i,k)-qbxi(i-1,k))
c     <                   + etx(i,1,k)*(qbzt(i,k)-qbzt(i-1,k))
c     <                   + xiz(i,1,k)*(qbxi(i,k)-qbxi(i,k-1))
     <                   + 1.D0/dz*(qbzt(i,k)-qbzt(i,k-1)))         
           endif
           u_bed(i, k, 2) = -(depth(i, k)-temp)/dtime
           
        enddo
        enddo

        call depth_bc
c        call depth_exchange

        do k = 1, nnk
        do i = 1, nni
           
           tot = tot + u_bed(i,k,2)
           
        enddo
        enddo

        temp = 1.D0/(nnk*nni)

        do k = -1, nnk+2
        do i = -1, nni+2
c           write(*,*) tot
           v_lid(i,k) = tot*temp
           
        enddo
        enddo



        return
        end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine bedGeometry
        
        include "size.inc"
        include "metric.inc"
        include "para.inc"
        include "sedi.inc"
        include "ns.inc"
        include "mpif.h"
        include "mpi.inc"

        integer i,k
        double precision h_x, h_z

        do k = 0, nnk+1
        do i = 0, nni+1
           h_x = -(depth(i+1, k)-depth(i-1, k))
     <           /(xp(i+1,1,k,1)-xp(i-1,1,k,1))      
           cos_x(i,k)=1.0/dsqrt(1.0+h_x**2)
           sin_x(i,k)=h_x*cos_x(i,k)
           
           h_z = -(depth(i, k+1)-depth(i, k-1))
     <           /(xp(i,1,k+1,3)-xp(i,1,k-1,3))      
           cos_z(i,k) = 1.0/dsqrt(1.0+h_z**2)
           sin_z(i,k) = h_z*cos_z(i,k)
        enddo
        enddo

        return
        end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine collectSedi

        include "size.inc"
        include "metric.inc"
        include "para.inc"
        include "sedi.inc"
        include "ns.inc"
        include "mpif.h"
        include "mpi.inc"


        integer i,j,k
        double precision temp


        temp = 1.D0/(1.D0-por)
        do k = -1, nnk+2
        do i = -1, nni+2
           sedv(i,k) = temp*(pick(i,k)+ ws*sedc(i,1,k))
        enddo
        enddo


        return 
        end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine btm_metric

        include "size.inc"
        include "metric.inc"
        include "para.inc"
        include "sedi.inc"
        include "ns.inc"
        include "mpif.h"
        include "mpi.inc"
        include "depth.inc"


        integer i,k
        double precision jab
        double precision xxi, zxi, xzt, zzt
	double precision xsx, xsz, zsx, zsz        

        do k = -1, nnk+2
	do i = -1, nni+2
	   jad(i,k) = 0.D0
	   xixd(i,k) = 0.D0
	   xizd(i,k) = 0.D0
	   ztxd(i,k) = 0.D0
	   ztzd(i,k) = 0.D0
            d11(i,k) = 0.D0 
            d13(i,k) = 0.D0 
            d31(i,k) = 0.D0 
            d33(i,k) = 0.D0 
	enddo
	enddo

        do k =  0, nnk+1
	do i = -1, nni+1
	   xxi = ( xp(i+1,1,k,1) - xp(i,1,k,1) )
	   zxi = ( xp(i+1,1,k,3) - xp(i,1,k,3) )
           xzt = 0.25D0 * ( xp(i,  1,k+1,1) - xp(i,  1,k-1,1)
     <                   +  xp(i+1,1,k+1,1) - xp(i+1,1,k-1,1) )
	   zzt = 0.25D0 * ( xp(i,  1,k+1,3) - xp(i,  1,k-1,3)
     <                     + xp(i+1,1,k+1,3) - xp(i+1,1,k-1,3) )
	   jab = xxi * zzt - xzt * zxi
	   jab = 1.D0 / jab
	   xsx =  zzt
	   zsx = -zxi
	   xsz = -xzt
	   zsz =  xxi
	   d11(i,k) = jab * ( xsx * xsx + xsz * xsz )
	   d13(i,k) = jab * ( xsx * zsx + xsz * zsz )
         
           xixd(i,k) = xsx
           xizd(i,k) = xsz
	   
	enddo
	enddo

        	
	do k = -1, nnk+1
	do i =  0, nni+1
	   xzt = ( xp(i,1,k+1,1) - xp(i,1,k,1) )
	   zzt = ( xp(i,1,k+1,3) - xp(i,1,k,3) )
	   xxi = 0.25D0 * ( xp(i+1,1,k  ,1) - xp(i-1,1,k  ,1) 
     <                    + xp(i+1,1,k+1,1) - xp(i-1,1,k+1,1) )
	   zxi = 0.25D0 * ( xp(i+1,1,k  ,3) - xp(i-1,1,k  ,3) 
     <                    + xp(i+1,1,k+1,3) - xp(i-1,1,k+1,3) )	  
           jab = xxi * zzt - xzt * zxi
	   jab = 1.D0 / jab
	   xsx =  zzt
	   zsx = -zxi
	   xsz = -xzt
	   zsz =  xxi

	   d31(i,k) = jab * ( zsx * xsx + zsz * xsz )
	   d33(i,k) = jab * ( zsx * zsx + zsz * zsz )
        
           ztxd(i,k) = zsx
           ztzd(i,k) = zsz

	enddo
	enddo

	
        do k = 0, nnk+1
	do i = 0, nni+1
	      xxi = 0.5D0 * ( xp(i+1,1,k,1) - xp(i-1,1,k,1) )
	      zxi = 0.5D0 * ( xp(i+1,1,k,3) - xp(i-1,1,k,3) )
	      xzt = 0.5D0 * ( xp(i,1,k+1,1) - xp(i,1,k-1,1) )
	      zzt = 0.5D0 * ( xp(i,1,k+1,3) - xp(i,1,k-1,3) )
	     
              jab = xxi * zzt - xzt * zxi
	      jab = 1.D0 / jab
	      jad(i,k) = jab
              
	enddo
	enddo

        return 
        end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine q_gravity

        include "size.inc"
        include "mpif.h"
        include "mpi.inc"
        include "ns.inc"
        include "metric.inc"
        include "para.inc"
        include "sedi.inc"
        include "eddy.inc"
        include "depth.inc"
        
        integer i, k
        double precision tan_r, sin_r, slope, gp

        gp = g/spwght
        tan_r = dtan(resp*pi/180)
        sin_r = dsin(resp*pi/180)
        do k = 0, nnk+1
        do i = 0, nni+1
           if (dabs(sin_x(i,k)) .ge. sin_r) then
              if (dabs(ug(i,k,1)) .ge. 1.D-8 ) then
                 ug(i,k,1) = ug(i,k,1)-dtime*(gp*sin_x(i,k)
     <                      -ug(i,k,1)/dabs(ug(i,k,1))*gp
     <                       *cos_x(i,k)*tan_r)
              else
                 ug(i,k,1) = ug(i,k,1)-dtime*(gp*sin_x(i,k)
     <                      -gp*cos_x(i,k)*tan_r)
              endif
           else
              ug(i,k,1) = 0.D0
           endif

c           write(*,*) dabs(sin_z(i,k)), sin_r
           if (dabs(sin_z(i,k)) .ge. sin_r) then
              
              if (dabs(ug(i,k,2)) .ge. 1.D-8 ) then
                 ug(i,k,2) = ug(i,k,2)-dtime*(gp*sin_z(i,k)
     <                      -ug(i,k,2)/dabs(ug(i,k,2))*gp
     <                       *cos_z(i,k)*tan_r)
              else
                 ug(i,k,2) = ug(i,k,2)-dtime*(gp*sin_z(i,k)
     <                      -gp*cos_z(i,k)*tan_r)
              endif
           else

              ug(i,k,2) = 0.D0
           endif


c           slope = dabs((depth(i+1,k)-depth(i-1,k))
c     <         /(xp(i+1,1,k,1)-xp(i-1,1,k,1)))
c           if (slope .le. tan_resp) then
c             if (dabs(u(i,1,k,1)) .ge. 1.D0-9) then 
c                ghi(i, k) = 0.5D0*g*(datum-depth(i,k))**2.D0
c     <                    - slope/dabs(slope)
c     <                    * tan(resp*pi/180)*g
c     <                    * 0.5D0*(xp(i+1,1,k,1)-xp(i-1,1,k,1))
c             else
c                ghi(i, k) = 0.5D0*g*(datum-depth(i,k))**2.D0
c             endif
c           endif
c           slope = dabs((depth(i,k+1)-depth(i,k-1))
c     <         /(xp(i,1,k+1,3)-xp(i,1,k-1,3)))
c           if ( slope .le. tan_resp) then
c             if (dabs(u(i,1,k,3)) .ge. 1.D0-9) then 
c                ghk(i, k) = 0.5D0*g*(datum-depth(i,k))**2.D0
c     <                    - u(i,1,k,3)/dabs(u(i,1,k,3))
c     <                    * tan(resp*pi/180)*g
c     <                    * 0.5D0*(xp(i,1,k+1,3)-xp(i,1,k-1,3))              
c             else
c                ghk(i, k) = 0.5D0*g*(datum-depth(i,k))**2.D0
c             endif
c           endif
        enddo
        enddo
              
        return
        end
