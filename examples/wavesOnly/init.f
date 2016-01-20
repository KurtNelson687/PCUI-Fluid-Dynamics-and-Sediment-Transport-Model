ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine initial

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "para.inc"
	include "ns.inc"
	include "metric.inc"
	include "eddy.inc"

	integer :: i, j, k
	logical :: iostat
	character*4 :: ID

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

	if ( case .eq. 0 .and. n_nrth .eq. MPI_PROC_NULL  ) then !This sets lid velocities to zero if case equals 0
	   do k = -1, nnk+2
	   do i = -1, nni+2
	      u_lid(i,k) = 0.D0
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

	hb = 0.D0
	
	hbs = 0.D0
	sus = 0.D0
	
	vst = 0.D0
	akst = 0.D0
	sab = 0.D0
	   
	rr = 0.D0

	write(ID, fmt='(I3)') 700+myid !initializes the density field from matlab
	if ( iscalar .eq. 1 ) then
	   inquire(file='rho_init_from_matlab.'//ID, exist=iostat)
	   if (iostat.eqv..true..and.grid_only.ne.1) then	      
	      open(700+myid, file = 'rho_init_from_matlab.'//ID,
     <                    form='unformatted',status='unknown')	 
	      read(700+myid) phi_init
	      close(700+myid)
	   end if
	end if


	if (  newrun .eq. 1  ) then

	   kount = 1
	   time = 0.D0
	   inquire(file='uvw_init_from_matlab.'//ID, exist=iostat) 
	   if (iostat.eqv..true..and.grid_only.ne.1) then
	      open(700+myid, file = 'uvw_init_from_matlab.'//ID,
     <                       form='unformatted',status='unknown')
	      read(700+myid) u
	      close(700+myid)
	   end if
	   call u_bc

	   uxi = 0.D0
	   uej = 0.D0
	   uzk = 0.D0
	
	   if ( iscalar .eq. 1 ) then
	      inquire(file='rho_full_from_matlab.'//ID, exist=iostat) 
	      if (iostat.eqv..true..and.grid_only.ne.1) then
		 open(700+myid, file = 'rho_full_from_matlab.'//ID,
     <                          form='unformatted',status='unknown')
		 read(700+myid) phi
		 close(700+myid)
	      end if
             call phi_bc
	   end if

C...... lid velocities u_lid and w_lid

	else

	   call input_continue_run

	endif

	return
	end
