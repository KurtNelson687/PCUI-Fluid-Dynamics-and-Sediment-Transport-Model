ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
        subroutine sedi_rhs

        include "size.inc"
        include "mpif.h"
        include "mpi.inc"
        include "ns.inc"
        include "metric.inc"
        include "para.inc"
        include "eddy.inc"
        include "sedi.inc"

        double precision coef, temp
        double precision delf, phie

        double precision cf(0:nni+1, 0:nnj+1, 0:nnk+1, 1:3),
     &                   u1(0:nni+1, 0:nnj+1, 0:nnk+1),
     &                   u2(0:nni+1, 0:nnj+1, 0:nnk+1),
     &                   u3(0:nni+1, 0:nnj+1, 0:nnk+1)        

        integer i, j, k
        
        double precision phimin, phimax, pmin, pmax
        double precision uep, uen, uwp, uwn,
     &                   unp, unn, usp, usn,
     &                   ufp, ufn, ubp, ubn
        
        
        coef = 1.5d0
c.......Add the settling velocity to the hydrodynamics velocity
        do k = 0, nnk+1
        do j = 0, nnj+1
        do i = 0, nni+1
           u1(i, j, k) = uxi(i, j, k) + wxi(i, j, k)
           u2(i, j, k) = uej(i, j, k) + wej(i, j, k)
           u3(i, j, k) = uzk(i, j, k) + wzk(i, j, k)
        enddo
        enddo
        enddo

c.......Take an Euler step on the first step

        if ( istep .eq. 1 ) then
           coef = 1.d0
           do k = 0, nnk+1
           do j = 0, nnj+1
           do i = 0, nni+1
              hbd(i, j, k) = 0.D0
              
           enddo
           enddo
           enddo
        endif

c.......First put in the part of Adams-Bashforth from step n-2
        
        do k = 0, nnk+1
        do j = 0, nnj+1
        do i = 0, nni+1
            sud(i, j, k) = -0.5D0 * hbd(i, j, k)
        enddo
        enddo
        enddo
        
c.......I direction

        do k = 0, nnk
        do j = 0, nnj
        do i = 0, nni
        if ( u1(i, j, k) .ge. 0.D0 ) then
           delf = sedc(i+1, j, k) - sedc(i-1, j, k)
           if( dabs(delf) .lt. 1.D-5 ) then
              cf(i, j, k, 1) = 0.125D0
           else
              phie = ( sedc(i,j,k)-sedc(i-1,j,k) ) / delf
              if ( phie .le. -1.D0 .or. phie .ge. 1.5D0 ) then
                 cf(i,j,k,1) = 0.125D0
              else
                 if ( phie .le. 0.25D0 ) then
                    if ( phie .le. 0.D0 ) then
                       cf(i,j,k,1) = 0.5D0 + 0.375D0 * phie
                    else
                       cf(i,j,k,1) = 0.5D0 - 0.625D0 * dsqrt(phie)
                    endif
                 else
                    if ( phie .le. 1.D0 ) then
                       cf(i,j,k,1) = 0.25D0 * ( 1.D0 - phie )
                    else
                       cf(i,j,k,1) =-0.25D0 * ( 1.D0 - phie )
                    endif
                endif
             endif
          endif
           sedf(i,j,k,1) = 0.5D0 * ( sedc(i,j,k) + sedc(i+1,j,k) )
     <                   -cf(i,j,k,1) * ( sedc(i-1,j,k)
     <                                  - sedc(i  ,j,k)
     <                                  - sedc(i  ,j,k)
     <                                  + sedc(i+1,j,k) )
       else
           delf = sedc(i, j, k) - sedc(i+2, j, k)
           if( dabs(delf) .lt. 1.D-5 ) then
              cf(i, j, k, 1) = 0.125D0
          else
              phie = ( sedc(i+1,j,k)-sedc(i+2,j,k) ) / delf
              if ( phie .le. -1.D0 .or. phie .ge. 1.5D0 ) then
                 cf(i,j,k,1) = 0.125D0
              else
                 if ( phie .le. 0.25D0 ) then
                    if ( phie .le. 0.D0 ) then
                       cf(i,j,k,1) = 0.5D0 + 0.375D0 * phie
                    else
                       cf(i,j,k,1) = 0.5D0 - 0.625D0 * dsqrt(phie)
                    endif
                 else
                    if ( phie .le. 1.D0 ) then
                       cf(i,j,k,1) = 0.25D0 * ( 1.D0 - phie )
                    else
                       cf(i,j,k,1) =-0.25D0 * ( 1.D0 - phie )
                    endif
                 endif
              endif
           endif
           sedf(i,j,k,1) = 0.5D0 * ( sedc(i,j,k) + sedc(i+1,j,k) )
     <                   -cf(i,j,k,1) * ( sedc(i  ,j,k)
     <                                  - sedc(i+1,j,k)
     <                                  - sedc(i+1,j,k)
     <                                  + sedc(i+2,j,k) )
        endif

        enddo
        enddo
        enddo

c.......J direction

        do k = 0, nnk
        do j = 0, nnj
        do i = 0, nni
        if ( u2(i, j, k) .ge. 0.D0 ) then
           delf = sedc(i, j+1, k) - sedc(i, j-1, k)
           if( dabs(delf) .lt. 1.D-5 ) then
              cf(i,j,k,2) = 0.125D0
           else
              phie = ( sedc(i,j,k)-sedc(i,j-1,k) ) / delf
              if ( phie .le. -1.D0 .or. phie .ge. 1.5D0 ) then
                 cf(i,j,k,2) = 0.125D0
              else
                 if ( phie .le. 0.25D0 ) then
                    if ( phie .le. 0.D0 ) then
                       cf(i,j,k,2) = 0.5D0 + 0.375D0 * phie
                    else
                       cf(i,j,k,2) = 0.5D0 - 0.625D0 * dsqrt(phie)
                    endif
                 else
                    if ( phie .le. 1.D0 ) then
                       cf(i,j,k,2) = 0.25D0 * ( 1.D0 - phie )
                    else
                       cf(i,j,k,2) =-0.25D0 * ( 1.D0 - phie )
                    endif
                 endif
              endif
           endif
           sedf(i,j,k,2) = 0.5D0 * ( sedc(i,j,k) + sedc(i,j+1,k) )
     <                   -cf(i,j,k,2) * ( sedc(i,j-1,k)
     <                                  - sedc(i,j,k)
     <                                  - sedc(i,j,k)
     <                                  + sedc(i,j+1,k) )
        else
           delf = sedc(i, j, k) - sedc(i, j+2, k)
           if( dabs(delf) .lt. 1.D-5 ) then
              cf(i,j,k,2) = 0.125D0
           else
              phie = ( sedc(i,j+1,k)-sedc(i,j+2,k) ) / delf
              if ( phie .le. -1.D0 .or. phie .ge. 1.5D0 ) then
                 cf(i,j,k,2) = 0.125D0
              else
                 if ( phie .le. 0.25D0 ) then
                    if ( phie .le. 0.D0 ) then
                      cf(i,j,k,2) = 0.5D0 + 0.375D0 * phie
                    else
                       cf(i,j,k,2) = 0.5D0 - 0.625D0 * dsqrt(phie)
                    endif
                 else
                    if ( phie .le. 1.D0 ) then
                       cf(i,j,k,2) = 0.25D0 * ( 1.D0 - phie )
                    else
                       cf(i,j,k,2) =-0.25D0 * ( 1.D0 - phie )
                    endif
                 endif
              endif
           endif
           sedf(i,j,k,2) = 0.5D0 * ( sedc(i,j,k) + sedc(i,j+1,k) )
     <                   -cf(i,j,k,2) * ( sedc(i,j  ,k)
     <                                  - sedc(i,j+1,k)
     <                                  - sedc(i,j+1,k)
     <                                  + sedc(i,j+2,k) )
        endif
           
        enddo
        enddo
        enddo


c           call get_btm_coef
c       if ( n_suth .eq. MPI_PROC_NULL ) then
c        do k = 0, nnk
c        do i = 0, nni
c           temp = (xp(i,0,k,2)+xp(i,1,k,2))/2.D0 + 2.D0*diam
c           call get_btm_coef
c           phif(i,0,k,2) = c1(i,k)*temp**2.D0+c2(i,k)*temp+c3(i,k)
c        enddo
c        enddo
c        endif
c.......K direction
        
        do k = 0, nnk
        do j = 0, nnj
        do i = 0, nni
        if ( u3(i, j, k) .ge. 0.D0 ) then
           delf = sedc(i, j, k+1) - sedc(i, j, k-1)
           if( dabs(delf) .lt. 1.D-5 ) then
              cf(i,j,k,3) = 0.125D0
           else
              phie = ( sedc(i,j,k)-sedc(i,j,k-1) ) / delf
              if ( phie .le. -1.D0 .or. phie .ge. 1.5D0 ) then
                 cf(i,j,k,3) = 0.125D0
              else
                 if ( phie .le. 0.25D0 ) then
                    if ( phie .le. 0.D0 ) then
                       cf(i,j,k,3) = 0.5D0 + 0.375D0 * phie
                    else
                       cf(i,j,k,3) = 0.5D0 - 0.625D0 * dsqrt(phie)
                    endif
                 else
                    if ( phie .le. 1.D0 ) then
                       cf(i,j,k,3) = 0.25D0 * ( 1.D0 - phie )
                    else
                       cf(i,j,k,3) =-0.25D0 * ( 1.D0 - phie )
                    endif
                 endif
              endif
           endif
           sedf(i,j,k,3) = 0.5D0 * ( sedc(i,j,k) + sedc(i,j,k+1) )
     <                   -cf(i,j,k,3) * ( sedc(i,j,k-1)
     <                                  - sedc(i,j,k  )
     <                                  - sedc(i,j,k  )
     <                                  + sedc(i,j,k+1) )
        else
           delf = sedc(i, j, k) - sedc(i, j, k+2)
           if( dabs(delf) .lt. 1.D-5 ) then
              cf(i, j, k, 3) = 0.125D0
           else
              phie = ( sedc(i,j,k+1)-sedc(i,j,k+2) ) / delf
              if ( phie .le. -1.D0 .or. phie .ge. 1.5D0 ) then
                 cf(i,j,k,3) = 0.125D0
              else
                 if ( phie .le. 0.25D0 ) then
                    if ( phie .le. 0.D0 ) then
                       cf(i,j,k,3) = 0.5D0 + 0.375D0 * phie
                    else
                       cf(i,j,k,3) = 0.5D0 - 0.625D0 * dsqrt(phie)
                    endif
                 else
                    if ( phie .le. 1.D0 ) then
                       cf(i,j,k,3) = 0.25D0 * ( 1.D0 - phie )
                    else
                       cf(i,j,k,3) =-0.25D0 * ( 1.D0 - phie )
                    endif
                 endif
              endif
           endif
           sedf(i,j,k,3) = 0.5D0 * ( sedc(i,j,k) + sedc(i,j,k+1) )
     <                   -cf(i,j,k,3) * ( sedc(i,j,k  )
     <                                  - sedc(i,j,k+1)
     <                                  - sedc(i,j,k+1)
     <                                  + sedc(i,j,k+2) )
        endif
        enddo
        enddo
        enddo

        do k = 1, nnk
        do j = 1, nnj
        do i = 1, nni
           hbd(i, j, k)=  - u1(i,  j,k) * sedf(i,  j,k,1)
     <                    + u1(i-1,j,k) * sedf(i-1,j,k,1)
     <                    - u2(i,j,  k) * sedf(i,j,  k,2)
     <                    + u2(i,j-1,k) * sedf(i,j-1,k,2)
     <                    - u3(i,j,k  ) * sedf(i,j,k,  3)
     <                    + u3(i,j,k-1) * sedf(i,j,k-1,3)

           HBD(I,J,K) = HBD(I,J,K) + SEDC(I,J,K) *
     <                ( U1(I,J,K) - U1(I-1,J,K)
     <                + U2(I,J,K) - U2(I,J-1,K)
     <                + U3(I,J,K) - U3(I,J,K-1) )
           

        enddo
        enddo
        enddo


c        do k = 0, nnk+1
c        do j = 0, nnj+1
c        do i = 0, nni+1
c           chi(i,j,k,1) = phif(i,j,k,2)*uej(i,j,k)/ety(i,j,k)
c           chi(i,j,k,2) = phif(i,j,k,2)*wej(i,j,k)/ety(i,j,k)
c           chi(i,j,k,3) = bchi(i,j,k)/ety(i,j,k)
c           chi(i,j,k,4) = (aksd(i,j,k)+ak)*(sedc(i,j+1,k)-sedc(i,j,k))/
c     <                    (xp(i,j+1,k,2)-xp(i,j,k,2))     
           
c        enddo
c        enddo
c        enddo


        do k = 1, nnk
        do j = 1, nnj
        do i = 1, nni
           hbd(i,j,k) = hbd(i,j,k) - rr(i,j,k,4)          
        enddo
        enddo
        enddo

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        if (movingGrid .eq. 1) then
           do k = 0, nnk
           do j = 0, nnj
           do i = 0, nni
              if ( -kej(i, j, k) .ge. 0.D0 ) then
                 delf = sedc(i, j+1, k) - sedc(i, j-1, k)
                 if( dabs(delf) .lt. 1.D-5 ) then
                    cf(i,j,k,2) = 0.125D0
                 else
                    phie = ( sedc(i,j,k)-sedc(i,j-1,k) ) / delf
                    if ( phie .le. -1.D0 .or. phie .ge. 1.5D0 ) then
                       cf(i,j,k,2) = 0.125D0
                    else
                       if ( phie .le. 0.25D0 ) then
                          if ( phie .le. 0.D0 ) then
                             cf(i,j,k,2) = 0.5D0 + 0.375D0 * phie
                          else
                             cf(i,j,k,2) = 0.5D0 - 0.625D0 * dsqrt(phie)
                          endif
                       else
                          if ( phie .le. 1.D0 ) then
                             cf(i,j,k,2) = 0.25D0 * ( 1.D0 - phie )
                          else
                             cf(i,j,k,2) =-0.25D0 * ( 1.D0 - phie )
                          endif
                       endif
                    endif
                 endif
                 sedf(i,j,k,2) = 0.5D0 * ( sedc(i,j,k) + sedc(i,j+1,k) )
     <                -cf(i,j,k,2) * ( sedc(i,j-1,k)
     <                - sedc(i,j,k)
     <                - sedc(i,j,k)
     <                                  + sedc(i,j+1,k) )  
              else
                 delf = sedc(i, j, k) - sedc(i, j+2, k)
                 if( dabs(delf) .lt. 1.D-5 ) then
                    cf(i,j,k,2) = 0.125D0
                 else
                    phie = ( sedc(i,j+1,k)-sedc(i,j+2,k) ) / delf
                    if ( phie .le. -1.D0 .or. phie .ge. 1.5D0 ) then
                       cf(i,j,k,2) = 0.125D0
                    else
                       if ( phie .le. 0.25D0 ) then
                          if ( phie .le. 0.D0 ) then
                             cf(i,j,k,2) = 0.5D0 + 0.375D0 * phie
                          else
                             cf(i,j,k,2) = 0.5D0 - 0.625D0 * dsqrt(phie)
                          endif
                       else
                          if ( phie .le. 1.D0 ) then
                             cf(i,j,k,2) = 0.25D0 * ( 1.D0 - phie )
                          else
                             cf(i,j,k,2) =-0.25D0 * ( 1.D0 - phie )
                          endif
                       endif
                    endif
                 endif
                 sedf(i,j,k,2) = 0.5D0 * ( sedc(i,j,k) + sedc(i,j+1,k) )
     <                -cf(i,j,k,2) * ( sedc(i,j  ,k)
     <                - sedc(i,j+1,k)
     <                - sedc(i,j+1,k)
     <                + sedc(i,j+2,k) )
              endif
           enddo
           enddo
           enddo
      
           do k = 1, nnk
           do j = 1, nnj
           do i = 1, nni
              hbd(i,j,k) = hbd(i,j,k)+kej(i,j,k)*sedf(i,j,k,2)
     <             -kej(i,j-1,k)*sedf(i,j-1,k,2)
c     <                           +sedc(i,j,k)*(kej(i,j-1,k)-kej(i,j,k))     
           enddo
           enddo
           enddo
      endif
c...... Cross viscous terms at step n-1 from Crank-Nicolson
        
        do k = 1, nnk
        do j = 1, nnj
        do i = 1, nni

           hbd(i,j,k) = hbd(i,j,k)
     <        + ( ak + 0.5D0*(aksd(i,j,k) + aksd(i+1,j,k)) ) *
     <          ( g12(i,  j,k) * ( sedc(i,  j+1,k) - sedc(i,  j-1,k)
     <                           + sedc(i+1,j+1,k) - sedc(i+1,j-1,k) )
     <          + g13(i,  j,k) * ( sedc(i,  j,k+1) - sedc(i,  j,k-1)
     <                           + sedc(i+1,j,k+1) - sedc(i+1,j,k-1) ) )
     <        - ( ak + 0.5D0*(aksd(i,j,k) + aksd(i-1,j,k)) ) *
     <          ( g12(i-1,j,k) * ( sedc(i,  j+1,k) - sedc(i,  j-1,k)
     <                           + sedc(i-1,j+1,k) - sedc(i-1,j-1,k) )
     <          + g13(i-1,j,k) * ( sedc(i,  j,k+1) - sedc(i,  j,k-1)
     <                           + sedc(i-1,j,k+1) - sedc(i-1,j,k-1) ) )

           hbd(i,j,k) = hbd(i,j,k)
     <        + ( ak + 0.5D0*(aksd(i,j,k) + aksd(i,j+1,k)) ) *
     <          ( g23(i,j  ,k) * ( sedc(i,j,  k+1) - sedc(i,j,  k-1)
     <                           + sedc(i,j+1,k+1) - sedc(i,j+1,k-1) )
     <          + g21(i,j  ,k) * ( sedc(i+1,j,  k) - sedc(i-1,j,  k)
     <                           + sedc(i+1,j+1,k) - sedc(i-1,j+1,k) ) )
     <        - ( ak + 0.5D0*(aksd(i,j,k) + aksd(i,j-1,k)) ) *
     <          ( g23(i,j-1,k) * ( sedc(i,j,  k+1) - sedc(i,j,  k-1)
     <                           + sedc(i,j-1,k+1) - sedc(i,j-1,k-1) )
     <          + g21(i,j-1,k) * ( sedc(i+1,j,  k) - sedc(i-1,j,  k)
     <                           + sedc(i+1,j-1,k) - sedc(i-1,j-1,k) ) )

           hbd(i,j,k) = hbd(i,j,k)
     <        + ( ak + 0.5D0*(aksd(i,j,k) + aksd(i,j,k+1)) ) *
     <          ( g31(i,j,k  ) * ( sedc(i+1,j,k  ) - sedc(i-1,j,k  )
     <                           + sedc(i+1,j,k+1) - sedc(i-1,j,k+1) )
     <          + g32(i,j,k  ) * ( sedc(i,j+1,k  ) - sedc(i,j-1,k  )
     <                           + sedc(i,j+1,k+1) - sedc(i,j-1,k+1) ) )
     <        - ( ak + 0.5D0*(aksd(i,j,k) + aksd(i,j,k-1)) ) *
     <          ( g31(i,j,k-1) * ( sedc(i+1,j,k  ) - sedc(i-1,j,k  )
     <                           + sedc(i+1,j,k-1) - sedc(i-1,j,k-1) )
     <          + g32(i,j,k-1) * ( sedc(i,j+1,k  ) - sedc(i,j-1,k  )
     <                           + sedc(i,j+1,k-1) - sedc(i,j-1,k-1) ) )

        enddo
        enddo
        enddo
        
        do k = 1, nnk
        do j = 1, nnj
        do i = 1, nni
           sud(i,j,k) = sud(i,j,k) + coef * hbd(i,j,k)           
        enddo
        enddo
        enddo
        
       
c...... Cross viscous terms at step n-1 from Crank-Nicolson

        do k = 0, nnk
        do j = 0, nnj
        do i = 0, nni

           sud(i,j,k) = sud(i,j,k)
     <        + ( ak + 0.5D0*(aksd(i,j,k) + aksd(i+1,j,k)) ) *
     <          g11(i,  j,  k  ) * ( sedc(i+1,j,  k  ) - sedc(i,j,k) )
     <        + ( ak + 0.5D0*(aksd(i,j,k) + aksd(i-1,j,k)) ) *
     <          g11(i-1,j,  k  ) * ( sedc(i-1,j,  k  ) - sedc(i,j,k) )
     <        + ( ak + 0.5D0*(aksd(i,j,k) + aksd(i,j+1,k)) ) *
     <          g22(i,  j,  k  ) * ( sedc(i,  j+1,k  ) - sedc(i,j,k) )
     <        + ( ak + 0.5D0*(aksd(i,j,k) + aksd(i,j-1,k)) ) *
     <          g22(i,  j-1,k  ) * ( sedc(i,  j-1,k  ) - sedc(i,j,k) )
     <        + ( ak + 0.5D0*(aksd(i,j,k) + aksd(i,j,k+1)) ) *
     <          g33(i,  j,  k  ) * ( sedc(i,  j,  k+1) - sedc(i,j,k) )
     <        + ( ak + 0.5D0*(aksd(i,j,k) + aksd(i,j,k-1)) ) *
     <          g33(i,  j,  k-1) * ( sedc(i,  j,  k-1) - sedc(i,j,k) )    
        enddo
        enddo
        enddo
        if (movingGrid .eq. 1) then
           do k = 1, nnk
           do j = 1, nnj 
           do i = 1, nni
             sud(i,j,k) = dtime * jac(i,j,k) * sud(i,j,k)
             sud(i,j,k) = sud(i,j,k)
     <                +(jac(i,j,k)/jac_old(i,j,k)-1.D0)*sedc(i,j,k)
           enddo
           enddo
           enddo
       else
           do k = 1, nnk
           do j = 1, nnj 
           do i = 1, nni
             sud(i,j,k) = dtime * jac(i,j,k) * sud(i,j,k)
           enddo
           enddo
           enddo
       endif


c        do k = 1, nnk
c        do j = 1, nnj
c        do i = 1, nni
c           sedc(i,j,k) = sedc(i,j,k) + sud(i,j,k)
c           write(*,*) sud(i,j,k)-(jac(i,j,k)/jac_old(i,j,k)-1.D0)
c     <               *sedc(i,j,k)
c     <     -(kej(i,j,k)*phif(i,j,k,2)-kej(i,j-1,k)*phif(i,j-1,k,2))      
c     <     *dtime*jac(i,j,k)  
c           write(*,*) sedc(i,j,k)
c           write(*,*) sud(i,j,k)
c        enddo
c        enddo      
c        enddo
              
        return
        end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine sedi_solve
        
        include "size.inc"
        include "mpif.h"
        include "mpi.inc"
        include "ns.inc"
        include "metric.inc"
        include "para.inc"
        include "eddy.inc"
        include "sedi.inc"

        double precision ax(nnj,0:nni+1), bx(nnj,0:nni+1),
     &                   cx(nnj,0:nni+1), fx(nnj,0:nni+1)
        double precision ay(nni,0:nnj+1), by(nni,0:nnj+1),
     &                   cy(nni,0:nnj+1), fy(nni,0:nnj+1)
        double precision az(nni,0:nnk+1), bz(nni,0:nnk+1),
     &                   cz(nni,0:nnk+1), fz(nni,0:nnk+1)
        
        integer i, j, k
        double precision coef, B, temp
        double precision dh, dzB
        
        coef = 0.5D0 * dtime

c....... solve for I-direction

        

        do k = 1, nnk

        do j = 1, nnj
        do i = 1, nni
           ax(j,i) = -( ak + 0.5D0*(aksd(i,j,k) + aksd(i-1,j,k)) )*
     <               coef * jac(i,j,k) * g11(i-1,j,k)
           cx(j,i) = -( ak + 0.5D0*(aksd(i,j,k) + aksd(i+1,j,k)) )*
     <               coef * jac(i,j,k) * g11(i,  j,k)
           bx(j,i) = 1.D0 - ax(j,i) - cx(j,i)
           fx(j,i) = sud(i,j,k)
           
        enddo
        enddo
c        do i = 1, nni
c           write(*,*) chi(i,1,k,2)
c        enddo

        if ( n_west .eq. MPI_PROC_NULL ) then
           do j = 1, nnj
              hbd(0,j,k)=( g12(0,j,k)*( sedc(0,j+1,k) - sedc(0,j-1,k)
     <                                + sedc(1,j+1,k) - sedc(1,j-1,k) ) 
     <                   + g13(0,j,k)*( sedc(0,j,k+1) - sedc(0,j,k-1)
     <                                + sedc(1,j,k+1) - sedc(1,j,k-1) ))
     <                   / g11(0,j,k)
              ax(j,0) = 0.D0 
              bx(j,0) = 1.D0
              cx(j,0) =-1.D0
              fx(j,0) = hbd(0,j,k)
           enddo    
        endif 

        if ( n_east .eq. MPI_PROC_NULL ) then
           do j = 1, nnj
              hbd(nni+1,j,k)=( g12(nni,j,k)
     <                       * ( sedc(nni,  j+1,k)- sedc(nni,  j-1,k)
     <                         + sedc(nni+1,j+1,k)- sedc(nni+1,j-1,k) ) 
     <                       + g13(nni,j,k)
     <                       * ( sedc(nni,  j,k+1)- sedc(nni,  j,k-1)  
     <                         + sedc(nni+1,j,k+1)- sedc(nni+1,j,k-1) ))
     <                       /g11(nni,j,k)                     
              ax(j,nni+1) = 1.D0 
              bx(j,nni+1) =-1.D0
              cx(j,nni+1) = 0.D0
              fx(j,nni+1) = hbd(nni+1,j,k)
           enddo    
        endif
        
        if ( periodic .eq. 1 ) then
        call trip(ax, bx, cx, fx, nnj, 1, nni, n_west, n_east )
        else
        call trid(ax, bx, cx, fx, nnj, 1, nni, n_west, n_east )
        endif

        do j = 1, nnj
        do i = 1, nni
           sud(i,j,k) = fx(j,i)
        enddo
        enddo
        
        enddo

c....... solve for J-direction

        do k = 1, nnk

        do j = 1, nnj
        do i = 1, nni
           ay(i,j) = -( ak + 0.5D0*(aksd(i,j,k) + aksd(i,j-1,k)) )*
     <               coef * jac(i,j,k) * g22(i,j-1,k)
           cy(i,j) = -( ak + 0.5D0*(aksd(i,j,k) + aksd(i,j+1,k)) )*
     <               coef * jac(i,j,k) * g22(i,  j,k)
           by(i,j) = 1.D0 - ay(i,j) - cy(i,j)
           temp = (xp(i,0,k,2)-xp(i,2,k,2))/(xp(i,1,k,2)-xp(i,2,k,2))
           if (j .eq. 1) then
c              by(i,j) = ay(i,j)*temp + by(i,j)
c              cy(i,j) = ay(i,j)*(1.D0-temp) + cy(i,j)
c              ay(i,j) = 0.D0
               by(i,j) = 1.D0
               cy(i,j) = 0.D0
               ay(i,j) = 0.D0
           endif
           if (j .eq. nnj) then
               B = (ws)*(xp(i,nnj+1,k,2)-xp(i,nnj,k,2))
     <                /(vis+aksd(i,j,k))
               by(i,j) = by(i,j)+cy(i,j)/(1.D0-B)
               cy(i,j) = 0.D0
           endif
c           write(*,*) ak
           fy(i,j) = sud(i,j,k)
        enddo
        enddo

        if ( n_suth .eq. MPI_PROC_NULL ) then
     
           do i = 1, nni
c              hbd(i,0,k)=( g23(i,0,k)*( sedc(i,0,k+1) - sedc(i,0,k-1)
c     <                                + sedc(i,1,k+1) - sedc(i,1,k-1) ) 
c     <                   + g21(i,0,k)*( sedc(i+1,0,k) - sedc(i-1,0,k)
c     <                                + sedc(1+1,1,k) - sedc(i-1,1,k) ))
c     <                   / g22(i,0,k) 
c     <             + pick(i,k)/(ak + 0.5D0*(aksd(i,0,k) + aksd(i,1,k)))
c              hbd(i,0,k) = 0.D0*pick(i,k)/(ak + aksd(i,1,k))
c               hbd(i,0,k) = pick(i,k)/(ak+aksd(i,1,k))                      
c     <                               *(xp(i,1,k,2)-xp(i,0,k,2))
c     <                   + (g23(i,0,k)*( sedc(i,0,k+1) - sedc(i,0,k-1)
c     <                                + sedc(i,1,k+1) - sedc(i,1,k-1) ) 
c     <                   + g21(i,0,k)*( sedc(i+1,0,k) - sedc(i-1,0,k)
c     <                                + sedc(1+1,1,k) - sedc(i-1,1,k) ))
c     <                   / g22(i,0,k)

               hbd(i,0,k) =0.D0*pick(i,k)*(xp(i,1,k,2)-xp(i,0,k,2))
c     <                     /(ak + 0.5*(aksd(i,1,k)+aksd(i,0,k)))
     <                   + (g23(i,0,k)*( sedc(i,0,k+1) - sedc(i,0,k-1)
     <                                + sedc(i,1,k+1) - sedc(i,1,k-1) ) 
     <                   + g21(i,0,k)*( sedc(i+1,0,k) - sedc(i-1,0,k)
     <                                + sedc(1+1,1,k) - sedc(i-1,1,k) ))
     <                   / g22(i,0,k)
              ay(i,0) = 0.D0 
              by(i,0) = 1.D0
              cy(i,0) = 0.D0
              fy(i,0) = hbd(i,0,k)
c              write(*,*) 0.5*(aksd(i,1,k)+aksd(i,0,k))
           enddo    
        endif 

        if ( n_nrth .eq. MPI_PROC_NULL ) then
           do i = 1, nni
              B = ws*(xp(i,nnj+1,k,2)-xp(i,nnj,k,2))
     <                /(100.D0*(vis+aksd(i,nnj,k)))
              hbd(i,nnj+1,k)=  0.D0 
     <                       + ( g23(i,nnj,k)
     <                       * ( sedc(i,nnj,  k+1)- sedc(i,nnj,  k-1)
     <                         + sedc(i,nnj+1,k+1)- sedc(i,nnj+1,k-1) ) 
     <                       + g21(i,nnj,k)
     <                       * ( sedc(i+1,nnj,  k)- sedc(i-1,nnj,  k)  
     <                         + sedc(i+1,nnj+1,k)- sedc(i-1,nnj+1,k) ))
     <                       /g22(i,nnj,k)
                            
              ay(i,nnj+1) =   0
              by(i,nnj+1) =   1.D0
              cy(i,nnj+1) = 0.D0
              fy(i,nnj+1) = hbd(i,nnj+1,k)
           enddo    
        endif

        call trid(ay, by, cy, fy, nni, 1, nnj, n_suth, n_nrth )
        
        do j = 1, nnj
        do i = 1, nni
           sud(i,j,k) = fy(i,j)
        enddo
        enddo
        
        enddo


c...... solve for K-direction
        do j = 1, nnj

        do k = 1, nnk
        do i = 1, nni
           az(i,k) = -( ak + 0.5D0*(aksd(i,j,k) + aksd(i,j,k-1)) )*
     <               coef * jac(i,j,k) * g33(i,j,k-1)
           cz(i,k) = -( ak + 0.5D0*(aksd(i,j,k) + aksd(i,j,k+1)) )*
     <               coef * jac(i,j,k) * g33(i,  j,k)
           bz(i,k) = 1.D0 - az(i,k) - cz(i,k)
           fz(i,k) = sud(i,j,k)
c           write(*,*) 'sedc = ', sedc(i,1,k)
        enddo
        enddo

        if ( n_back .eq. MPI_PROC_NULL ) then
           do i = 1, nni
              hbd(i,j,0)=( g31(i,j,0)*( sedc(i+1,j,0) - sedc(i-1,j,0)
     <                                + sedc(i+1,j,1) - sedc(i-1,j,1) ) 
     <                   + g32(i,j,0)*( sedc(i,j+1,0) - sedc(i,j-1,0)
     <                                + sedc(1,j+1,1) - sedc(1,j-1,1) ))
     <                   / g33(i,j,0)
              az(i,0) = 0.D0 
              bz(i,0) = 1.D0
              cz(i,0) =-1.D0
              fz(i,0) = hbd(i,j,0)
           enddo    
        endif 

        if ( n_frnt .eq. MPI_PROC_NULL ) then
           do i = 1, nni
              hbd(i,j,nnk+1)=( g31(i,j,nnk)
     <                       * ( sedc(i+1,j,nnk  )- sedc(i-1,j,nnk  )
     <                         + sedc(i+1,j,nnk+1)- sedc(i-1,j,nnk+1) ) 
     <                       + g32(i,j,nnk)
     <                       * ( sedc(i,j+1,nnk  )- sedc(i,j-1,nnk  )  
     <                         + sedc(i,j+1,nnk+1)- sedc(i,j-1,nnk+1) ))
     <                       /g33(i,j,nnk)
            
              az(i,nnk+1) = 1.D0 
              bz(i,nnk+1) =-1.D0
              cz(i,nnk+1) = 0.D0
              fz(i,nnk+1) = hbd(i,j,nnk+1)
                
           enddo    
        endif
        
        if ( periodicZ .eq. 1 ) then
        call trip(az, bz, cz, fz, nni, 1, nnk, n_back, n_frnt )
        else
        call trid(az, bz, cz, fz, nni, 1, nnk, n_back, n_frnt )
        endif

        do k = 1, nnk
        do i = 1, nni
           sud(i,j,k) = fz(i,k)
        enddo
        enddo
        
        enddo

        do k = 1, nnk
        do j = 1, nnj
        do i = 1, nni
           sedc(i,j,k) = sedc(i,j,k) + sud(i,j,k)
           if (sedc(i,j,k) .le. 0.D0) then
              sedc(i,j,k) = 0
           endif
        enddo
        enddo      
        enddo
        

        if (movingGrid .ne. 1) then
           if (j .eq. 1) then
              dzB = 2.0*(depth(i,k) + xp(i,j,k,2))
              dh = -ws*sedc(i,j,k)*dtime/(1.D0-por)
              if (dzB - dh .le. 0) then
                 write(*,*) 'Bed raises too fast-----dh = ', dh
              else
                 sedc(i,j,k) = sedc(i,j,k)*dzB/(dzB-dh)
     <                - (1.D0-por)*dh/(dzB-dh)
              endif                 
           endif
           
        endif
        
           
        


 
        call sedi_bc
        
        call sedi_exchange


        do k = -1, nnk+2
        do j = -1, nnj+2
        do i = -1, nni+2
           rho_tot(i,j,k) = sedc(i,j,k)*spwght + (1.D0-sedc(i,j,k))
        enddo
        enddo      
        enddo


        return
        end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine sedi_bc
      
      include "size.inc"
      include "mpif.h"
      include "para.inc"            
      include "mpi.inc"
      include "ns.inc"
      include "sedi.inc"
      include "eddy.inc"
      include "metric.inc"
      

      integer i, j, k
      double precision B, temp
      
      if ( n_west .eq. MPI_PROC_NULL) then
         do k = -1, nnk+2 
         do j = -1, nnj+2
            sedc( 0,j,k) = sedc(1,j,k) + hbd(0,j,k)
c             sedc( 0,j,k) = sedc(1,j,k)
c            sedc( 0,j,k) = 0.D0
         enddo
         enddo
         do k = -1, nnk+2 
         do j = -1, nnj+2
            sedc(-1,j,k) = 3.D0 * ( sedc(0,j,k) - sedc(1,j,k) )
     <                   + sedc(2,j,k)
c            sedc( -1,j,k) = sedc(0,j,k)
            
         enddo
         enddo
      endif

      if ( n_east .eq. MPI_PROC_NULL) then
         do k = -1, nnk+2 
         do j = -1, nnj+2
            sedc(nni+1,j,k) = sedc(nni,j,k) - hbd(nni+1,j,k)
c            sedc(nni+1,j,k) = sedc(nni,j,k)
c            sedc(nni+1,j,k) = 0.D0
         enddo
         enddo
         do k = -1, nnk+2 
         do j = -1, nnj+2
            sedc(nni+2,j,k) = 3.D0 * ( sedc(nni+1,j,k) - sedc(nni,j,k) )
     <                      + sedc(nni-1,j,k)
c            sedc(nni+2,j,k) = sedc(nni+1,j,k)
         enddo
         enddo
      endif
             
      if ( n_suth .eq. MPI_PROC_NULL) then
         do k = -1, nnk+2 
         do i = -1, nni+2
c             sedc(i, 0,k) = sedc(i,1,k)
c             sedc(i, 0,k) = 1.D0
              temp = (xp(i,0,k,2)+xp(i,1,k,2))/2.0
c              sedc(i,0,k) = c1(i,k)*xp(i,-1,k,2)**2.D0
c     <                    + c2(i,k)*xp(i,-1,k,2) + c3(i,k)
c              sedc(i, 0,k) = c1(i,k)*temp**2.D0 
c     <                     + c2(i,k)*temp + c3(i,k)
c              write(*,*) c1(i,k), c2(i,k), c3(i,k)
c+ sedc(i,1,k)-sedc(i,2,k)
               sedc(i,0,k) = sedc(i,1,k)
         enddo
         enddo
         do k = -1, nnk+2 
         do i = -1, nni+2
             sedc(i,-1,k) = 3.D0 * ( sedc(i,0,k) - sedc(i,1,k) )
     <                   + sedc(i,2,k)
c             sedc(i,-1,k) = sedc(i,0,k) + (sedc(i,0,k) - sedc(i,1,k))
c            sedc(i,-1,k) = 1.D0
             sedc(i, -1, k) = sedc(i,0, k)
         enddo
         enddo
      endif

      if ( n_nrth .eq. MPI_PROC_NULL) then
         do k = -1, nnk+2 
         do i = -1, nni+2
c             B = (vis+vst(i,nnj,k))/(xp(i,nnj+1,k,2)-xp(i,nnj,k,2))
              B = (ws)*(xp(i,nnj+1,k,2)-xp(i,nnj,k,2))
     <                /(30.D0*vis)
c              write(*,*) B
c            sedc(i,nnj+1,k) = sedc(i,nnj,k) - hbd(i,nnj+1,k)
c            sedc(i,nnj+1,k) = sedc(i,nnj,k) + hbd(i,nnj+1,k)
             sedc(i,nnj+1,k) = 1.D0/(1.D0-B)*sedc(i,nnj,k)
c             sedc(i,nnj+1,k) = 1.D0
         enddo
         enddo
         do k = -1, nnk+2 
         do i = -1, nni+2
            sedc(i,nnj+2,k) = 3.D0 * ( sedc(i,nnj+1,k) - sedc(i,nnj,k) )
     <                      + sedc(i,nnj-1,k)
c             sedc(i,nnj+2,k) = sedc(i, nnj+1,k)
c     <                         ( sedc(i,nnj+1,k) - sedc(i,nnj,k))
             
         enddo
         enddo
      endif

      if ( n_back .eq. MPI_PROC_NULL) then
         do j = -1, nnj+2 
         do i = -1, nni+2
            sedc(i,j, 0) = sedc(i,j,1) + hbd(i,j,0)
c             sedc(i,j, 0) = 1.D0
         enddo
         enddo
         do j = -1, nnj+2 
         do i = -1, nni+2
c            sedc(i,j,-1) = sedc(i,j,0)
            sedc(i,j,-1) = 3.D0 * ( sedc(i,j,0) - sedc(i,j,1) )
     <                   + sedc(i,j,2)
         enddo
         enddo
      endif

      if ( n_frnt .eq. MPI_PROC_NULL) then
         do j = -1, nnj+2 
         do i = -1, nni+2
            sedc(i,j,nnk+1) = sedc(i,j,nnk) - hbd(i,j,nnk+1)
c            sedc(i,j,nnk+1) = 1.D0
         enddo
         enddo
         do j = -1, nnj+2 
         do i = -1, nni+2
c            sedc(i,j,nnk+2) = sedc(i,j,nnk+1)
            sedc(i,j,nnk+2) = 3.D0 * ( sedc(i,j,nnk+1) - sedc(i,j,nnk) )
     <                      + sedc(i,j,nnk-1)
         enddo
         enddo
      endif
      
      if ( n_west .eq. MPI_PROC_NULL ) then
         
         if ( n_suth .eq. MPI_PROC_NULL ) then
            do k = 1, nnk
c            sedc(0, 0, k) = 1.D0
            sedc(0, 0, k) = 0.5D0 * ( sedc(1,0,k) + sedc(0,1,k) )
            enddo
         endif
         if ( n_nrth .eq. MPI_PROC_NULL ) then
            do k = 1, nnk
c            sedc(0,nnj+1,k) = 1.D0   
            sedc(0,nnj+1,k) = 0.5D0 * (sedc(1,nnj+1,k) + sedc(0,nnj,k) )
            enddo
         endif
         

         if ( n_back .eq. MPI_PROC_NULL ) then
            do j = 1, nnj
c            sedc(0, j, 0) = 1.D0
            sedc(0, j, 0) = 0.5D0 * ( sedc(1,j,0) + sedc(0,j,1) )
            enddo
         endif
         if ( n_frnt .eq. MPI_PROC_NULL ) then
            do j = 1, nnj
c            sedc(0,j,nnk+1) = 1.D0   
            sedc(0,j,nnk+1) = 0.5D0 * (sedc(1,j,nnk+1) + sedc(0,j,nnk) )
            enddo
         endif
         
      endif

      if ( n_east .eq. MPI_PROC_NULL ) then         
         if ( n_suth .eq. MPI_PROC_NULL ) then
            do k = 1, nnk
c            sedc(nni+1,0,k) = 1.D0
            sedc(nni+1,0,k) = 0.5D0 * (sedc(nni,0,k) + sedc(nni+1,1,k))
            enddo
         endif
         if ( n_nrth .eq. MPI_PROC_NULL ) then
            do k = 1, nnk            
c            sedc(nni+1,nnj+1,k) = 1.D0               
            sedc(nni+1,nnj+1,k) = 0.5D0 * ( sedc(nni,nnj+1,k) 
     <                                    + sedc(nni+1,nnj,k) )
            enddo
         endif
         

         if ( n_back .eq. MPI_PROC_NULL ) then
            do j = 1, nnj
c            sedc(nni+1,j,0) = 1.D0
            sedc(nni+1,j,0) = 0.5D0 * ( sedc(nni,j,0) 
     <                                 + sedc(nni+1,j,1) )
            enddo
         endif
         if ( n_frnt .eq. MPI_PROC_NULL ) then
            do j = 1, nnj
c            sedc(nni+1,j,nnk+1) = 1.D0   
            sedc(nni+1,j,nnk+1) = 0.5D0 * ( sedc(nni,j,nnk+1)
     <                                    + sedc(nni+1,j,nnk) )
            enddo
         endif        
      endif

      if ( n_suth .eq. MPI_PROC_NULL ) then
         
         if ( n_back .eq. MPI_PROC_NULL ) then
            do i = 1, nni
c            sedc(i, 0, 0) = 1.D0
            sedc(i, 0, 0) = 0.5D0 * ( sedc(i,1,0) + sedc(i,0,1) )
            enddo
         endif
         if ( n_frnt .eq. MPI_PROC_NULL ) then
            do i = 1, nnk
c            sedc(i,0,nnk+1) = 1.D0   
            sedc(i,0,nnk+1) = 0.5D0 * (sedc(i,1,nnk+1) + sedc(i,0,nnk) )
            enddo
         endif
      endif

      if ( n_nrth .eq. MPI_PROC_NULL ) then
         if ( n_back .eq. MPI_PROC_NULL ) then
            do i = 1, nni
c               sedc(i, nnj+1, 0) = 1.D0 
            sedc(i, nnj+1, 0) = 0.5D0 * ( sedc(i,nnj,0) 
     <                                  + sedc(i,nnj+1,1) )
            enddo
         endif
         if ( n_frnt .eq. MPI_PROC_NULL ) then
            do i = 1, nni
c            sedc(i,nnj+1,nnk+1) = 1.D0   
            sedc(i,nnj+1,nnk+1) = 0.5D0 * ( sedc(i,nnj,nnk+1) 
     <                                    + sedc(i,nnj+1,nnk) )
            enddo
         endif    
      endif

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine qb_bc
      
      include "size.inc"
      include "mpif.h"
      include "para.inc"            
      include "mpi.inc"
      include "ns.inc"
      include "sedi.inc"
      include "eddy.inc"
      include "metric.inc"
      

      integer i, j, k
      double precision B, temp
      
      if ( n_west .eq. MPI_PROC_NULL) then
         do k = -1, nnk+2 
            qbi( 0,k) = qbi(1,k)
            qbk( 0,k) = qbk(1,k)
            qbi( -1,k) = qbi(1,k)
            qbk( -1,k) = qbk(1,k)

         enddo
      endif

      if ( n_east .eq. MPI_PROC_NULL) then
         do k = -1, nnk+2 
            qbi(nni+1,k) = qbi(nni,k)
            qbi(nni+2,k) = qbi(nni,k)
            qbk(nni+1,k) = qbk(nni,k)
            qbk(nni+2,k) = qbk(nni,k)
         enddo         
      endif
             
 

      if ( n_back .eq. MPI_PROC_NULL) then
         do i = -1, nni+2
            qbi(i,0)  = qbi(i,1)
            qbi(i,-1) = qbi(i,1)
            qbk(i,0)  = qbk(i,1)
            qbk(i,-1) = qbk(i,1)
         enddo
      endif

      if ( n_frnt .eq. MPI_PROC_NULL) then
         do i = -1, nni+2
            qbi(i,nnk+1) = qbi(i,nnk)
            qbi(i,nnk+2) = qbi(i,nnk)
            qbk(i,nnk+1) = qbk(i,nnk)
            qbk(i,nnk+2) = qbk(i,nnk)
         enddo
      endif
    

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine bed_shear

        include "size.inc"
        include "mpif.h"
        include "mpi.inc"
        include "ns.inc"
        include "metric.inc"
        include "para.inc"
        include "eddy.inc"
        include "sedi.inc"
        
        double precision utx, utz, n2, n0, n1
        integer i,j,k

        n0 = diam
        if (n_suth .eq. MPI_PROC_NULL) then
        
        do k = 0, nnk+1   
        do i = 0, nni+1
           n2 = dabs(xp(i, 1, k, 2) + depth(i,k))
           utx = (u(i,1,k,1)-u_bed(i,k,1)) * cos_x(i,k)
           utz = u(i,1,k,3) * cos_z(i,k)
           sin_xz(i,k) = (dabs(utx*sin_x(i,k))+dabs(utz*sin_z(i,k)))
     <                   /dsqrt(utx**2.0+utz**2.0)
           cos_xz(i,k) = dsqrt(1.0-sin_xz(i,k))

           utan(i,k) = dsqrt(utx**2+utz**2)
c            utan(i,k) = abs(u(i,1,k,1))
c           ushr(i,k) = 0.41D0*utan(i,k)/dlog((n2+n0)/n0)
c           ushr(i,k) = (u(i,2,k,1) - u(i,1,k,1))
c     <               / (dlog((xp(i,2,k,2)+depth(i,k))
c     <                      /(xp(i,1,k,2)+depth(i,k))))*0.41
c           thta(i,k) = ushr(i,k)**2.D0/(dble(spwght-1.D0)*g*diam)
           Cdxi(i,k)   = 1.D0/(1.D0/0.41D0
     <                  *dlog((n2*cos_x(i,k)+n0)/n0))**2.D0
c           Cdxi(i,k)   = 1.D0/(1.D0/0.41D0
c     <                  *dlog((n2+n0)/n0))**2.D0
           Cdzt(i,k)   = 1.D0/(1.D0/0.41D0
     <                  *dlog((n2*cos_z(i,k)+n0)/n0))**2.D0
           ushr(i,k) = dsqrt(Cdxi(i,k)*utx**2
     <                  +    Cdzt(i,k)*utz**2)
           thta(i,k) = ushr(i,k)**2.D0
     <                 /(dble(spwght-1.D0)*g*diam)
           
        enddo
        enddo
        endif

        return
        end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
        subroutine pickup

        include "size.inc"
        include "metric.inc"
        include "para.inc"
        include "sedi.inc"
        include "ns.inc"
        
        integer i, k
        double precision excs, dstar, sin_r, cos_r, thta0

        sin_r = dsin(resp*pi/180.0)
        cos_r = dcos(resp*pi/180.0)

        dstar   = ((spwght-1)*g / vis**2)**(1.D0/3.D0)
     <          * diam
        thta0 = 0.24D0/(dstar) + 0.055
     <          *(1.D0 - dexp(-0.02D0*dstar))
        
        do k = 0, nnk+1          
        do i = 0, nni+1
      
           if (utan(i,k)*cos_xz(i,k) .ge. 0.0) then
              thta_cr = thta0*(sin_r*cos_xz(i,k)
     <                          +cos_r*sin_xz(i,k))/sin_r
           else
              thta_cr = thta0*(sin_r*cos_xz(i,k)
     <                          -cos_r*sin_xz(i,k))/sin_r
           endif
c           if (sin_r*cos_x(i,k)-cos_r*sin_x(i,k) .le. 0.D0) then
c              write(*,*) 'Bed slope larger than the angle of repose'
c              stop
c           endif
c           excs = thta(i,k)-thta_cr
            excs = tau_mean - 0.05
           if (excs .gt. 0.0) then
c              pick(i,k) = 0.00033D0*(excs/thta_cr)**1.5
c     <               * ((spwght-1)*g)**0.6 * diam**0.8
c     <               / vis**0.2
              pick(i,k) = 0.005/2650*dexp(6.5*dsqrt(excs))             
           else
              pick(i,k) = 0.D0
           endif

           if (bedload .eq. 1) then
              if (excs .gt. 0.0) then              
                 qbi(i,k) = 8.D0*excs**1.5D0*diam
     <                  *dsqrt((spwght-thta_cr)*diam)
     <                  *u(i,1,k,1)* dabs(cos_x(i,k))/utan(i,k)
                 qbk(i,k) = 8.D0*excs**1.5D0*diam
     <                  *dsqrt((spwght-thta_cr)*diam)
     <                  *u(i,1,k,3)* dabs(cos_z(i,k))/utan(i,k)

              else
                 qbi(i,k) = 0.D0
                 qbk(i,k) = 0.D0
              endif
           endif

        enddo
        enddo

        call qb_bc
        call qb_exchange

        return
        end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    
        subroutine bottomC

        include "size.inc"
        include "metric.inc"
        include "para.inc"
        include "sedi.inc"
        include "ns.inc"
        include "mpif.h"
        include "mpi.inc"
        integer i,j,k
        double precision temp, z1, z2

        if ( n_suth .eq. MPI_PROC_NULL) then


        do k = 0, nnk+1
        do i = 1, nni
           z1 = depth(i,k) + xp(i,1,k,2)
           z2 = depth(i,k) + xp(i,2,k,2)
           temp = (sedc(i,1,k)/sedc(i,2,k)*
     <            (((depth(i,k)-z2)/z2)**z2)/
     <             (((depth(i,k)-z1)/z1)**z1))**(1.D0/(z1-z2))
           sedc(i,0,k) = sedc(i,1,k)*
     <      (0.0001/(depth(i,k)-0.0001)*(depth(i,k)-z1)/z1)**(-z1)
           
        enddo
        enddo

        endif
        return 
        end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine updateDepth
        
        include "size.inc"
        include "metric.inc"
        include "para.inc"
        include "sedi.inc"
        include "ns.inc"
        include "mpif.h"
        include "mpi.inc"
        
        integer i,j,k
        double precision temp, checkTime, dx, dz, thta_b,
     &                   thta_c, dstar, excs, force,
     &                   cb, rp_height, ustar, sig
        
        checkTime = dtime*ncheck
        ustar = amp*omega*dsin(omega*time)
        rp_height = (hf-h0) - (hf-h0)*dexp(-grow_b*time) + h0

        if (ustar .le. 10D-8) then
           sig = 0.0
        else
           sig = dabs(ustar)/ustar
        endif

        do k = -1, nnk+2
        do i = -1, nni+2
           thta_b = 0.5*fwc*ustar**2/((spwght-1)*g*diam)
           dstar = ((spwght-1.0)*g/vis**2)**(1.0/3.0)*diam
           thta_c = 0.24/dstar + 0.055*(1-dexp(-0.02*dstar))
           excs = thta_b-thta_c
c           if (excs .le. 0.0) then
c              Qb = 0.0
c           else
c              Qb = 12*(excs)*dsqrt(thta_b)*diam
c     <             *dsqrt((spwght-1)*g*diam)*sig              
c           endif

c           cb = 0.0*Qb/(1.0-por)/(rp_height/2)
           temp = depth(i,k)
           
           depth(i,k) = 0.2D0 + rp_height/2*
     <          dcos(2.0*pi/rp_length*(xp(i,1,k,1)-cb*time))
           
           u_bed(i, k, 2) = -(depth(i, k)-temp)/checkTime
           u_bed(i, k, 1) = cb
        enddo
        enddo

        dx = xp(2,1,1,1)-xp(1,1,1,1)
        dz = xp(1,1,2,3)-xp(1,1,1,3)
        do k = 0, nnk+1
        do i = 0, nni+1

c           write(*,*) w_bed(i,k)
           u_bed(i, k,2) = u_bed(i,k,2)
c     <                   + u(i,1,k,1)*(-0.5D0*(depth(i+1,k)
c     <                                        -depth(i-1,k))/dx)
c     <                   + u(i,1,k,3)*(-0.5D0*(depth(i,k+1)
c     <                                        -depth(i,k-1))/dz)
c           write(*,*) w_bed(i,k)
        enddo
        enddo

        do k = 0, nnk+1

          
           u_bed(-1, k,2) = u_bed(-1,k,2)
c     <                   + u(-1,1,k,1)*(-(depth(0,  k)
c     <                                  -depth(-1 ,k))/dx)
c     <                   + u(-1,1,k,3)*(-0.5D0*(depth(-1,k+1)
c     <                                        -depth(-1,k-1))/dz)
           temp = depth(nni+2,k)
                     
           u_bed(nni+2, k,2) = u_bed(nni+2,k,2)
c     <                   + u(nni+2,1,k,1)*(-(depth(nni+2, k)
c     <                                      -depth(nni+1 ,k))/dx)
c     <                   + u(nni+2,1,k,3)*(-0.5D0*(depth(nni+2,k+1)
c     <                                        -depth(nni+2,k-1))/dz)
        enddo

        do i = 0, nni+1

           
          
           u_bed(i, -1,2) = u_bed(i,-1,2)
c    <                   + u(i,1,-1,1)*(-0.5D0*(depth(i+1, -1)
c    <                                  -       depth(i-1 ,-1))/dx)
c    <                   + u(i,1,-1,3)*(-(depth(i, 0)
c    <                                   -depth(i,-1))/dz)


           
           
           u_bed(i, nnk+2,2) = u_bed(i,nnk+2,2)
c     <                 + u(i,1,nnk+2,1)*(-0.5D0*(depth(i+1, nnk+2)
c     <                                          -depth(i-1, nnk+1))/dx)
c     <                   + u(i,1,nnk+2,3)*(-(depth(i,nnk+2)
c     <                                      -depth(i,nnk+1))/dz)
        enddo

           
     
           u_bed(-1, -1,2) = u_bed(-1,1,2)
c     <                   + u(-1,1,-1,1)*(-(depth(0, -1)
c     <                                    -depth(-1,-1))/dx)
c     <                   + u(-1,1,-1,3)*(-(depth(-1, 0)
c     <                                    -depth(-1,-1))/dz)

           
           
           u_bed(nni+2, nnk+2,2) = u_bed(nni+2, nnk+2,2)
c     <                 + u(nni+2,1,nnk+2,1)*(-(depth(nni+2, nnk+2)
c     <                                        -depth(nni+1, nnk+2))/dx)
c     <                 + u(nni+2,1,nnk+2,3)*(-(depth(nni+2, nnk+2)
c     <                                        -depth(nni+2, nnk+1))/dz)

           



           u_bed(-1, nnk+2,2) = u_bed(-1, nnk+2,2)
c     <                   + u(-1,1,nnk+2,1)*(-(depth(0,  nnk+2)
c     <                                       -depth(-1 ,nnk+2))/dx)
c     <                   + u(-1,1,nnk+2,3)*(-(depth(-1, nnk+2)
c     <                                       -depth(-1, nnk+1))/dz)

           



           u_bed(nni+2, -1,2) = u_bed(nni+2,-1,2)
c     <                   + u(nni+2,1,-1,1)*(-(depth(nni+2, -1)
c     <                                       -depth(nni+1 ,-1))/dx)
c     <                   + u(nni+2,1,-1,3)*(-(depth(nni+2,  0)
c     <                                       -depth(nni+2, -1))/dz)
        return
        end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        
        subroutine get_sedi_BC_coef

        include "size.inc"
        include "metric.inc"
        include "para.inc"
        include "sedi.inc"
        include "ns.inc"
        include "mpif.h"
        include "mpi.inc"

        integer i,k
        double precision temp, a, ae

        ae = 2.D0*diam
        do k = -1, nnk+2
        do i = -1, nni+2

           a = depth(i,k) + 0.5D0*(xp(i,2,k,2)+xp(i,1,k,2))
           if (a .le. ae) then
           sd_coef(i,k) = 1.D0
           else
           temp = -ws/(0.4D0*dsqrt(depth(i,k)*Pgrd))
           sd_coef(i,k) = ((depth(i,k)-a)/a
     <                     *ae/(depth(i,k)-ae))**(-temp)
           endif
c           write(*,*) sd_coef(i,k), xp(i,1,k,2)

        enddo
        enddo

        return 
        end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine get_btm_coef

        include "size.inc"
        include "metric.inc"
        include "para.inc"
        include "sedi.inc"
        include "ns.inc"
        include "mpif.h"
        include "mpi.inc"

        integer i, k
        double precision A(1:3, 1:3), b(1:3), x(1:3)
        
        do k = -1, nnk+2
        do i = -1, nni+2
           b(1) = sedc(i, 1, k)
           b(2) = sedc(i, 2, k)
           b(3) = sedc(i, 3, k)

           x(1) = 0.D0
           x(2) = 0.D0          
           x(3) = 0.D0         
           
           A(1, 1) = xp(i,1,k,2)**2.D0
           A(1, 2) = xp(i,1,k,2)
           A(1, 3) = 1.D0
           A(2, 1) = xp(i,2,k,2)**2.D0
           A(2, 2) = xp(i,2,k,2)
           A(2, 3) = 1.D0
           A(3, 1) = xp(i,3,k,2)**2.D0
           A(3, 2) = xp(i,3,k,2)
           A(3, 3) = 1.D0
        

           call lusolve(A, x, b, 3)

           c1(i, k) = x(1)
           c2(i, k) = x(2)
           c3(i, k) = x(3)
        enddo
        enddo
        return
        end
ccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine totalU

        include "size.inc"
        include "metric.inc"
        include "para.inc"
        include "sedi.inc"
        include "ns.inc"
        include "mpif.h"
        include "mpi.inc"

        integer i,j,m,k
        
        do m = 1, 3
        do k = 1, nnk+1
        do j = 1, nnj+1
        do i = 1, nni+1
           if (m .eq. 2) then
              upf(i,j,k,m) = u(i,j,k,m) + sedc(i,j,k)*ws
           else
              upf(i,j,k,m) = u(i,j,k,m)
           endif
        enddo
        enddo
        enddo
        enddo

        return
        end
