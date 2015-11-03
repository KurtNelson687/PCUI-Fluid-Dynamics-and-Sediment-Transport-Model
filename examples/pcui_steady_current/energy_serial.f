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
      double precision rho_global(ni*nj*nk)
      double precision volume_global(ni*nj*nk)
      double precision Ep, Eb, height, area

      integer i, j, k, m, mmax 
      integer indx(ni*nj*nk)

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

C..... aggregate local 1D arrays into global array on proc 0
      call MPI_GATHER( rho_local,  m, MPI_DOUBLE_PRECISION,
     <                 rho_global, m, MPI_DOUBLE_PRECISION,
     <                 0, MPI_COMM_WORLD, ierr )
      call MPI_GATHER( volume_local,  m, MPI_DOUBLE_PRECISION,
     <                 volume_global, m, MPI_DOUBLE_PRECISION,
     <                 0, MPI_COMM_WORLD, ierr )
      call MPI_BARRIER( MPI_COMM_WORLD, ierr)

C..... sort global density array on proc 0 and calculate Eb
      mmax = ni*nj*nk
      if ( MYID .EQ. 0) then
        call indexx(mmax,rho_global,indx)

c      if ( MYID .EQ. 0 ) then
c        do m = 1, mmax
c          write(*,*) rho_global(indx(m)) 
c        enddo
c      endif

        area = bx*bz
        height = 0.5*volume_global(indx(mmax))/area
        Eb = rho_global(indx(mmax))*volume_global(indx(mmax))*height
        do m = mmax-1, 1, -1
           height = height + volume_global(indx(m))/area
           Eb = Eb + rho_global(indx(m))*volume_global(indx(m))*height
        enddo
        Eb = g*Eb
      endif

C..... show value on screen
      if ( MYID .EQ. 0 ) 
     <   write(*,*) 'Ep = ', Ep, " Eb = ", Eb

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
     



