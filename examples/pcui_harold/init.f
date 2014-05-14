ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine initial

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "para.inc"
	include "ns.inc"
	include "metric.inc"
	include "eddy.inc"

	integer :: i, j, k, m, idum
        double precision :: zeta, ran0
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

	write(ID, fmt='(I4)') 1000+myid
	if ( iscalar .eq. 1 ) then
c	   inquire(file='rho_init_from_matlab.'//ID, exist=iostat)
c	   if (iostat.eqv..true..and.grid_only.ne.1) then	      
c	      open(1000+myid, file = 'rho_init_from_matlab.'//ID,
c     <                    form='unformatted',status='unknown')	 
c	      read(1000+myid) phi_init
c	      close(1000+myid)
c	   end if

           do k = -1, nnk+2
           do j = -1, nnj+2
           do i = -1, nni+2
              phi_init(i,j,k) = 1.D0
           enddo
           enddo
           enddo
	end if


	if (  newrun .eq. 1  ) then

	   kount = 1
	   time = 0.D0
c	   inquire(file='uvw_init_from_matlab.'//ID, exist=iostat) 
c	   if (iostat.eqv..true..and.grid_only.ne.1) then
c	      open(1000+myid, file = 'uvw_init_from_matlab.'//ID,
c     <                       form='unformatted',status='unknown')
c	      read(1000+myid) u
c	      close(1000+myid)
c	   end if
           
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

	   uxi = 0.D0
	   uej = 0.D0
	   uzk = 0.D0
	
	   if ( iscalar .eq. 1 ) then
c	      inquire(file='rho_full_from_matlab.'//ID, exist=iostat) 
c	      if (iostat.eqv..true..and.grid_only.ne.1) then
c		 open(1000+myid, file = 'rho_full_from_matlab.'//ID,
c     <                          form='unformatted',status='unknown')
c		 read(1000+myid) phi
c		 close(1000+myid)
c	      end if

              idum = 8654
              do k = -1, nnk+2
              do j = -1, nnj+2
              do i = -1, nni+2
                 zeta = -0.1D0*exp(-(xp(i,j,k,1)/0.7D0)**2) 
     <                  + 0.001D0*ran0(idum)
                 phi(i,j,k) = 1.D0-0.5D0*0.03D0
     <            *tanh(2D0*(xp(i,j,k,2)-zeta+0.3D0)/0.02D0*2.64665D0)
              enddo
              enddo
              enddo
              call phi_bc

c             inquire(file='phi_init_from_matlab.'//ID, exist=iostat) 
c	      if (iostat.eqv..true..and.grid_only.ne.1) then
c		 open(1000+myid, file = 'phi_init_from_matlab.'//ID,
c     <                          form='unformatted',status='unknown')
c		 read(1000+myid) phi2
c		 close(1000+myid)
c	      end if

              do k = -1, nnk+2
              do j = -1, nnj+2
              do i = -1, nni+2
                 phi2(i,j,k) = xp(i,j,k,1) 
              enddo
              enddo
              enddo
              call phi2_bc

              do k = -1, nnk+2
              do j = -1, nnj+2
              do i = -1, nni+2
                 phi3(i,j,k) = xp(i,j,k,2) 
              enddo
              enddo
              enddo
              call phi3_bc

              do k = -1, nnk+2
              do j = -1, nnj+2
              do i = -1, nni+2
                 phi4(i,j,k) = xp(i,j,k,3) 
              enddo
              enddo
              enddo
              call phi4_bc

	   end if

C...... lid velocities u_lid and w_lid

	else

	   call input_continue_run

	endif

	return
	end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        function ran0(idum)

        integer idum,IA,IM,IQ,IR,MASK
        double precision ran0,AM
        parameter(IA=16807,IM=2147483647,AM=1./IM,
     <       IQ=127773,IR=2836,MASK=123459876)
        integer k

        idum=ieor(idum,MASK)
        k=idum/IQ
        idum=IA*(idum-k*IQ)-IR*k
        if (idum.lt.0) idum=idum+IM
        ran0=AM*idum
        idum=ieor(idum,MASK)
        return
        end
