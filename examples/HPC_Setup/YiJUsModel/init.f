ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine initial

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "para.inc"
	include "ns.inc"
	include "metric.inc"
	include "eddy.inc"
	include "sedi.inc"
	include "depth.inc"

	integer i, j, k, m, L
        double precision n0
        
        Pgrd = 0.00095D0
C...... lid velocities u_lid and w_lid

	if ( case .eq. 1 .and. n_nrth .eq. MPI_PROC_NULL ) then
	   do k = -1, nnk+2
	   do i = -1, nni+2
c	      u_lid(i,k) = - omg_lid 
c     <                   * 0.5D0 * (xp(i,nnj,k,3)+xp(i,nnj+1,k,3))
c	      w_lid(i,k) =   omg_lid  
c     <                   * 0.5D0 * (xp(i,nnj,k,1)+xp(i,nnj+1,k,1))
c	      v_lid(i,k) = 0.050
	   enddo
	   enddo
	endif

	if ( case .eq. 0 .and. n_nrth .eq. MPI_PROC_NULL  ) then
	   do k = -1, nnk+2
	   do i = -1, nni+2
	      u_lid(i,k) = 0.D0      
c             u_lid(i,k) = 1.D0
	      w_lid(i,k) = 0.D0
	      v_lid(i,k) = 0.D0
	   enddo
	   enddo
	endif

	if ( case .eq. 0. and. n_west .eq. MPI_PROC_NULL ) then
	   do k = -1, nnk+2
	   do j = -1, nnj+2
	      u_flux(j,k) = 
     &           0.D0*dexp(-200.d0*(xp(1,j,k,2)-0.75D0)**2)
	   enddo
	   enddo
	end if

	if ( case .eq. 0. and. n_east .eq. MPI_PROC_NULL ) then
	   do k = -1, nnk+2
	   do j = -1, nnj+2
	      u_flux(j,k)= 
     &           0.D0*dexp(-200.d0*(xp(1,j,k,2)-0.25D0)**2)
	   enddo
	   enddo
	end if

	if ( n_nrth .eq. MPI_PROC_NULL ) then
	   if ( n_west .eq. MPI_PROC_NULL ) then
	      do k = -1, nnk+2
	         u_lid( 0,k) = - u_lid(1,k)*0.D0
	         u_lid(-1,k) = - u_lid(2,k)*0.D0
	      enddo
	   endif
	   if ( n_east .eq. MPI_PROC_NULL ) then
	      do k = -1, nnk+2
	         u_lid(nni+1,k) = - u_lid(nni,  k)*0.D0
	         u_lid(nni+2,k) = - u_lid(nni-1,k)*0.D0
	      enddo
	   endif
	   if ( n_back .eq. MPI_PROC_NULL ) then
	      do i = -1, nni+2
	         u_lid(i, 0) = - u_lid(i,1)*0.D0
	         u_lid(i,-1) = - u_lid(i,2)*0.D0
	      enddo
	   endif
	   if ( n_frnt .eq. MPI_PROC_NULL ) then
	      do i = -1, nni+2
	         u_lid(i,nnk+1) = - u_lid(i,nnk  )*0.D0
	         u_lid(i,nnk+2) = - u_lid(i,nnk-1)*0.D0
	      enddo
	   endif
	endif

	do m = 1, 3
	do k = 0, nnk+1
	do j = 0, nnj+1
	do i = 0, nni+1
	   hb(i,j,k,m) = 0.D0
	   pcorr(i,j,k,m) = 1.D0
	   sedf(i,j,k,m) = 0.0
	enddo
	enddo
	enddo
	enddo

	do k = 0, nnk+1
	do j = 0, nnj+1
	do i = 0, nni+1
	   hbs(i,j,k) = 0.D0
	   sus(i,j,k) = 0.D0
           hbd(i,j,k) = 0.D0
           sud(i,j,k) = 0.D0
	    
	enddo
	enddo
	enddo

        if ( sedi .eq. 1) then

	do k = -1, nnk+2
        do i = -1, nni+2
	   pick(i,k) = 0.D0
	   thta(i,k) = 0.D0
	   ushr(i,k) = 0.D0
	   utan(i,k) = 0.D0
	   sedv(i,k) = 0.D0
	sd_coef(i,k) = 0.D0
	    qbi(i,k) = 0.D0
	    qbk(i,k) = 0.D0
	      h(i,k) = 0.D0
c....set initial fct = 2 to make no slip condition
           fct(i,k)  = -1.D0
	    c1(i,k)  = 0.D0  
	    c2(i,k)  = 0.D0  
	    c3(i,k)  = 0.D0
	  drhs(i,k)  = 0.D0    
	  af_x(i,k)  = 0.D0
	  af_z(i,k)  = 0.D0
        enddo
	enddo

	do k = 0, nnk+1
        do i = 0, nni+1
	   qbxi(i,k) = 0.D0
	   qbzt(i,k) = 0.D0
	   ugxi(i,k) = 0.D0
	   ugzt(i,k) = 0.D0	   
        enddo
	enddo

	do m =1, 2
	do k = -1, nnk+2
        do i = -1, nni+2
	   ug(i,k,m) = 0.D0
	enddo
	enddo
	enddo
        
	endif

	do k = 0, nnk+1
	do j = 0, nnj+1
	do i = 0, nni+1
	   vst(i,j,k) = 0.D0
	  akst(i,j,k) = 0.D0
	   sab(i,j,k) = 0.D0
	enddo
	enddo
	enddo

        do m = 1, 3
	do k = 0, nnk+1
	do j = 0, nnj+1
	do i = 0, nni+1
	   btau(i,j,k,m) = 0.D0
           chi(i,j,k,m) = 0.D0
	   sut(i,j,k,m) = 0.D0
	enddo
	enddo
	enddo
	enddo

        do m = 1, 5
	do k = 0, nnk+1
	do j = 0, nnj+1
	do i = 0, nni+1
           chi(i,j,k,m) = 0.D0
	enddo
	enddo
	enddo
	enddo

	do m = 1, 3
	do k = -1, nnk+2
	do i = -1, nni+2
           u_bed(i,k,m) = 0.D0
	enddo
	enddo
	enddo

	do k = -1, nnk+2
	do j = -1, nnj+2
	do i = -1, nni+2
	   rho_tot(i,j,k) = 1.D0
	   sai(i,j,k) = 0.D0
	enddo
	enddo
	enddo
	



	if ( sedi .eq. 1 ) then
	do k = 0, nnk+1
	do j = 0, nnj+1
	do i = 0, nni+1
	   aksd(i,j,k) = 0.D0

	   bchi(i,j,k) = 0.D0
	enddo
	enddo
	enddo
	endif


	do m = 1, 9
	do k = 0, nnk+1
	do j = 0, nnj+1
	do i = 0, nni+1
	   rr(i,j,k,m) = 0.D0
	enddo
	enddo
	enddo 
	enddo 

        do m = 1, 4
	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni
	   taus(i,j,k,m) = 0.D0
	   
	enddo
	enddo
	enddo 
	enddo 

C
C This is the initial scalar profile
C
	if ( iscalar .eq. 1 ) then
	   do k = -1, nnk+2
	   do j = -1, nnj+2
	   do i = -1, nni+2        
     
	     phi_init(i,j,k) = 1.0
        
	   enddo
	   enddo
	   enddo
	else
	   do k = -1, nnk+2
	   do j = -1, nnj+2
	   do i = -1, nni+2       
	      phi_init(i,j,k) = 0.D0
	      phi(i,j,k) = 0.D0
        
	   enddo
	   enddo
	   enddo

	endif


	if (  newrun .eq. 1  ) then

	   kount = 1
	   time = 0.D0
        
	   do k = -1, nnk+2
	   do j = -1, nnj+2
           do i = -1, nni+2
	      sedc(i,1,k) = 0.D0	      
	   enddo
	   enddo
	   enddo


	
c	   if (sedi .eq. 1) then

c	      call sedi_init
c              call sedi_bc        
c	      call sedi_exchange

c	   endif


	   do m = 1, 3
	   do k = -1, nnk+2
	   do j = -1, nnj+2
	   do i = -1, nni+2
	      u(i,j,k,m) = 0.D0
	      upf(i,j,k,m) = 0.D0
	   enddo
	   enddo
	   enddo
	   enddo
        
	   n0 = diam
	   do k = 1, nnk
	   do j = 1, nnj
	   do i = 1, nni
   
	   u(i,j,k,1) = dsqrt(Pgrd*depth(i,k))/0.4D0
     <                * dlog((xp(i,j,k,2)+depth(i,k)+n0)/n0)
c	   u(i,j,k,1) = 1.D0
	   enddo
	   enddo
	   enddo
	
	   call randomU
	
	   call u_bc

	   do k = 0, nnk+1
	   do j = 0, nnj+1
	   do i = 0, nni+1
	      uxi(i,j,k) = 0.D0
	      uej(i,j,k) = 0.D0
	      uzk(i,j,k) = 0.D0
	      wxi(i,j,k) = 0.D0
	      wej(i,j,k) = 0.D0
	      wzk(i,j,k) = 0.D0
	      kxi(i,j,k) = 0.D0
	      kej(i,j,k) = 0.D0
	      kzk(i,j,k) = 0.D0
	   enddo
	   enddo
	   enddo


c	   if ( iscalar .eq. 1 ) then
c	      do k = -1, nnk+2
c              do j = -1, nnj+2
c	      do i = -1, nni+2
c		 if (xp(i,j,k,1) .ge. 50.0) then
c		    phi_init(i,j,k) = 1.D0 + 0.76D0/1000.D0*30.D0                    
c		 else
c		    phi_init(i,j,k) = 1.D0
c		 endif
c     		    phi(i,j,k) = phi_init(i,j,k)
c	      enddo
c	      enddo
c              enddo
        



c	   call phi_bc

c	   endif

	else
	   
c	   call input
	   
           if (sedi .eq. 1) then
c	      
c              call sedi_init
c              call sedi_bc        
c              call sedi_exchange
           endif
	endif
        
	return
	end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine sedi_init

        include "size.inc"
	include "metric.inc"
	include "para.inc"
	include "sedi.inc"
        include "ns.inc"

	integer i, j, k

        do k = 1, nnk
        do j = 1, nnj
        do i = 1, nni
c     set sediment BC
c	   if (xp(i,j,k,1) .le. 10.0) then
	      sedc(i,j,k) = 0.D0
c	   endif
        enddo
        enddo
        enddo
        
        
	return 
	end
