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

c...... Velocity at west bc
        do k = -1, nnk+2
        do j = -1, nnj+2
           u_west(j,k) = 0.05D0/0.3D0*xp(1,j,k,2)
        enddo
        enddo

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

	write(ID, fmt='(I3)') 700+myid
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

C......	UXI

	do k = 1, nnk
	do j = 1, nnj
	do i = ius, iue
	   uxi(i,j,k) = xix(i,j,k) * 
     <	        ( 0.75 D0 * u(i,j,k,1) + 0.375D0 * u(i+1,j,k,1) 
     <          - 0.125D0 * u(i-1,j,k,1) )
     <	              + xiy(i,j,k) *
     <	        ( 0.75 D0 * u(i,j,k,2) + 0.375D0 * u(i+1,j,k,2) 
     <          - 0.125D0 * u(i-1,j,k,2) )
     <	              + xiz(i,j,k) *
     <	        ( 0.75 D0 * u(i,j,k,3) + 0.375D0 * u(i+1,j,k,3) 
     <          - 0.125D0 * u(i-1,j,k,3) )
	enddo
	enddo
	enddo
         if ( n_west .eq. MPI_PROC_NULL ) then
            do k = 1, nnk
          	do j = 1, nnj
               uxi(0,j,k) = u_west(j,k) * xix(0,j,k) 
            enddo
            enddo
         endif

         if ( n_east .eq. MPI_PROC_NULL ) then
            do k = 1, nnk
          	do j = 1, nnj
            uxi(nni,j,k) = uxi(nni-1,j,k)
c           uxi(nni,j,k) = (u(nni,j,k,1)+u(nni+1,j,k,1))/2.D0 
c    <                     * xix(nni,j,k)
            enddo
            enddo
         endif

C......	UEJ

	do k = 1, nnk
	do j = jus, jue
	do i = 1, nni
	   uej(i,j,k) = etx(i,j,k) * 
     <	        ( 0.75 D0 * u(i,j,k,1) + 0.375D0 * u(i,j+1,k,1) 
     <          - 0.125D0 * u(i,j-1,k,1) )
     <	              + ety(i,j,k) *
     <	        ( 0.75 D0 * u(i,j,k,2) + 0.375D0 * u(i,j+1,k,2) 
     <          - 0.125D0 * u(i,j-1,k,2) )
     <	              + etz(i,j,k) *
     <	        ( 0.75 D0 * u(i,j,k,3) + 0.375D0 * u(i,j+1,k,3) 
     <          - 0.125D0 * u(i,j-1,k,3) )
	enddo
	enddo
	enddo

C......	UZK

	do k = kus, kue
	do j = 1, nnj
	do i = 1, nni
	   uzk(i,j,k) = ztx(i,j,k) * 
     <	        ( 0.75 D0 * u(i,j,k,1) + 0.375D0 * u(i,j,k+1,1) 
     <          - 0.125D0 * u(i,j,k-1,1) )
     <	              + zty(i,j,k) *
     <	        ( 0.75 D0 * u(i,j,k,2) + 0.375D0 * u(i,j,k+1,2) 
     <          - 0.125D0 * u(i,j,k-1,2) )
     <	              + ztz(i,j,k) *
     <	        ( 0.75 D0 * u(i,j,k,3) + 0.375D0 * u(i,j,k+1,3) 
     <          - 0.125D0 * u(i,j,k-1,3) )
	enddo
	enddo
	enddo

	
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
