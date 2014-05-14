ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine stat_init

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "para.inc"
	include "ns.inc"
	include "metric.inc"
	include "eddy.inc"
	include "stat.inc"

	integer :: i, j, k, m

	nz = 10 
	
	read(400+myid) na
	read(400+myid) time

	read(400+myid) us
	read(400+myid) ts
	
	return
	end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine stat_coll

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "para.inc"
	include "ns.inc"
	include "metric.inc"
	include "eddy.inc"
	include "stat.inc"

	integer :: i, j, k, m

	if ( mod(istep, nz) .eq. 0 ) then

	na = na + 1

	do m = 1, 3
	do k = -1, nnk+2
	do j = -1, nnj+2
	do i = -1, nni+2
	   us(i,j,k,m) = us(i,j,k,m) + u(i,j,k,m)
	enddo
	enddo
	enddo
	enddo
	
	do m = 1, 3
	do k = -1, nnk+2
	do j = -1, nnj+2
	do i = -1, nni+2
	   ts(i,j,k,m) = ts(i,j,k,m) + u(i,j,k,m)**2
	enddo
	enddo
	enddo
	enddo

	do k = -1, nnk+2
	do j = -1, nnj+2
	do i = -1, nni+2
	   ts(i,j,k,4) = ts(i,j,k,4) + u(i,j,k,1) * u(i,j,k,2)
	   ts(i,j,k,5) = ts(i,j,k,5) + u(i,j,k,2) * u(i,j,k,3)
	   ts(i,j,k,6) = ts(i,j,k,6) + u(i,j,k,3) * u(i,j,k,1)
	enddo
	enddo
	enddo

	endif
	
	return
	end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine stat_post

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "para.inc"
	include "ns.inc"
	include "metric.inc"
	include "eddy.inc"
	include "stat.inc"

	integer :: i, j, k, m

	write(500+myid) na
	write(500+myid) time

	write(500+myid) us
	write(500+myid) ts

	return
	end

