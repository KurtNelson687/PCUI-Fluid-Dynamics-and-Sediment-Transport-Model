ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine initial

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "para.inc"
	include "ns.inc"
	include "metric.inc"
	include "eddy.inc"

	integer :: i, j, k, m, L

C...... lid velocities u_lid and w_lid

	if ( case .eq. 1 .and. n_nrth .eq. MPI_PROC_NULL ) then
	   do k = -1, nnk+2
	   do i = -1, nni+2
	      u_lid(i,k) = - omg_lid 
     <                   * 0.5D0 * (xp(i,nnj,k,3)+xp(i,nnj+1,k,3))
	      w_lid(i,k) =   omg_lid  
     <                   * 0.5D0 * (xp(i,nnj,k,1)+xp(i,nnj+1,k,1))
	   enddo
	   enddo
	endif

	if ( case .eq. 0 .and. n_nrth .eq. MPI_PROC_NULL  ) then
	   do k = -1, nnk+2
	   do i = -1, nni+2
	      u_lid(i,k) = 1.D0
	      w_lid(i,k) = 0.D0
	   enddo
	   enddo
	endif

	if ( n_nrth .eq. MPI_PROC_NULL ) then
	   if ( n_west .eq. MPI_PROC_NULL ) then
	      do k = -1, nnk+2
	         u_lid( 0,k) = - u_lid(1,k)
	         u_lid(-1,k) = - u_lid(2,k)
	      enddo
	   endif
	   if ( n_east .eq. MPI_PROC_NULL ) then
	      do k = -1, nnk+2
	         u_lid(nni+1,k) = - u_lid(nni,  k)
	         u_lid(nni+2,k) = - u_lid(nni-1,k)
	      enddo
	   endif
	   if ( n_back .eq. MPI_PROC_NULL ) then
	      do i = -1, nni+2
	         u_lid(i, 0) = - u_lid(i,1)
	         u_lid(i,-1) = - u_lid(i,2)
	      enddo
	   endif
	   if ( n_frnt .eq. MPI_PROC_NULL ) then
	      do i = -1, nni+2
	         u_lid(i,nnk+1) = - u_lid(i,nnk  )
	         u_lid(i,nnk+2) = - u_lid(i,nnk-1)
	      enddo
	   endif
	endif

	do m = 1, 3
	do k = 0, nnk+1
	do j = 0, nnj+1
	do i = 0, nni+1
	   hb(i,j,k,m) = 0.D0
	enddo
	enddo
	enddo
	enddo

	do k = 0, nnk+1
	do j = 0, nnj+1
	do i = 0, nni+1
	   hbs(i,j,k) = 0.D0
	   sus(i,j,k) = 0.D0
	enddo
	enddo
	enddo

	do k = 0, nnk+1
	do j = 0, nnj+1
	do i = 0, nni+1
	   vst(i,j,k) = 0.D0
	  akst(i,j,k) = 0.D0
	   sab(i,j,k) = 0.D0
	enddo
	enddo
	enddo

	do m = 1, 6
	do k = 0, nnk+1
	do j = 0, nnj+1
	do i = 0, nni+1
	   rr(i,j,k,m) = 0.D0
	enddo
	enddo
	enddo 
	enddo 

	if ( iscalar .eq. 1 ) then
	   do k = -1, nnk+2
	   do j = -1, nnj+2
	   do i = -1, nni+2
	      phi_init(i,j,k) = 0.5D0 * ( phi1 + phi2 ) 
     <	                      + 0.5D0 * ( phi1 - phi2 ) 
     <		* dtanh( aphi * ( xp(i,j,k,2) - yphi ) )
	   enddo
	   enddo
	   enddo

	else

	   do k = -1, nnk+2
	   do j = -1, nnj+2
	   do i = -1, nni+2
	      phi_init(i,j,k) = 0.D0
	      phi     (i,j,k) = 0.D0
	   enddo
	   enddo
	   enddo

	endif

	if (  newrun .eq. 1  ) then

	kount = 1
	time = 0.D0

	do m = 1, 3
	do k = -1, nnk+2
	do j = -1, nnj+2
	do i = -1, nni+2
	   u(i,j,k,m) = 0.D0
	enddo
	enddo
	enddo
	enddo
	
	call u_bc

	do k = 0, nnk+1
	do j = 0, nnj+1
	do i = 0, nni+1
	   uxi(i,j,k) = 0.D0
	   uej(i,j,k) = 0.D0
	   uzk(i,j,k) = 0.D0
	enddo
	enddo
	enddo

	if ( iscalar .eq. 1 ) then

	   do k = -1, nnk+2
	   do j = -1, nnj+2
	   do i = -1, nni+2
	      phi(i,j,k) = phi_init(i,j,k)
	   enddo
	   enddo
	   enddo

CCC	   call phi_bc

	endif

	else

	   call input

	endif

	return
	end
