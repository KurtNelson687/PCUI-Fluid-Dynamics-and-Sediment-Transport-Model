ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      
      subroutine energy

      include "size.inc"
      include "mpif.h"
      include "mpi.inc"
      include "metric.inc"
      include "cavity.inc"
      include "para.inc"
      INCLUDE "ns.inc"

      double precision rho_local(nni*nnj*nnk)
      double precision volume_local(nni*nnj*nnk)
      double precision wksp(nni*nnj*nnk)
      double precision samples_local(px*py*pz)
      double precision samples_global(px*px*py*py*pz*pz)
      double precision pivots(px*py*pz-1)
      double precision Ep, Eb, height, area

      integer i, j, k, m, N_global, Np, indx
      integer iwksp(nni*nnj*nnk)

C..... rearrange density into 1D array on each proc
c      and calculate Ep
      m = 0
      Ep = 0
      do k = 1, nnk
      do j = 1, nnj
      do i = 1, nni
         m = m + 1
         rho_local(m) = phi(i,j,k)
         volume_local(m) = 1/jac(i,j,k)
         Ep = Ep + phi(i,j,k)*xp(i,j,k,2)*(1/jac(i,j,k))
      enddo
      enddo
      enddo

      Ep = g*Ep
      call global_sum(Ep)

C..... sort initial arrays on each proc
      call sort2(m,rho_local,volume_local,wksp,iwksp)

C..... aggregate samples of local arrays on root
      N_global = ni*nj*nk
      Np = px*py*pz
      do i = 1, Np
         indx = i*N_global/(Np*Np)
         samples_local(i) = rho_local(indx)
      enddo
      call MPI_GATHER( samples_local,  Np, MPI_DOUBLE_PRECISION,
     <                 samples_global, Np, MPI_DOUBLE_PRECISION,
     <                 0, MPI_COMM_WORLD, ierr )
      call MPI_BARRIER( MPI_COMM_WORLD, ierr)

C..... select pivots from global samples on root and broadcast them

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine global_sum(a)

      include "mpif.h"
      include "mpi.inc"

      double precision a, total

      call MPI_REDUCE( a, total, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, 
     <                 MPI_COMM_WORLD, ierr )
      call MPI_BCAST( total, 1, MPI_DOUBLE_PRECISION, 0,
     <                 MPI_COMM_WORLD, ierr )
      a = total

      return
      end
     



