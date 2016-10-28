cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	program ns

	include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "para.inc"
	include "sedi.inc"
        include "eddy.inc"


	double precision ta, tt, t0, t1, t2
     <                   , t3, t4, t5, t6
        integer i,j,k
	
	call MPI_INIT( ierr )


	call MPI_Barrier(MPI_COMM_WORLD, ierr)
	ta = MPI_Wtime()

	t0 = 0.D0
	t1 = 0.D0
	t2 = 0.D0
	t3 = 0.D0
	t4 = 0.D0
	t5 = 0.D0
	t6 = 0.D0

	call parameter       

        call mpi_initial
	
	if (newrun .ne. 1) then
	   call input
	endif

	
	call grid
        
        call initial        
	
c	call get_sedi_BC_coef
	call oldGrid

	imove = 0
	do istep = 1, nstep

c	   call sedi_init
c	   call sedi_bc
c	   call sedi_exchange
C	   if ( mod(istep, nsave) .eq. 0 .and. MYID .EQ. 0 )
	   if ( MYID .EQ. 0 )
     <	      write(*,*) ' istep = ', istep, ' kount  = ', kount
        call MPI_Barrier(MPI_COMM_WORLD, ierr)
	if (sedi  .eq. 1) then
           call transformWs   
	endif
        call MPI_Barrier(MPI_COMM_WORLD, ierr)

        call MPI_Barrier(MPI_COMM_WORLD, ierr)
	tt = MPI_Wtime()
	   if ( ieddy .eq. 1 ) call eddy	   
	call MPI_Barrier(MPI_COMM_WORLD, ierr)
	t0 = t0 + MPI_Wtime() - tt

	call MPI_Barrier(MPI_COMM_WORLD, ierr)
	tt = MPI_Wtime()
	   if ( iscalar .eq. 1 ) call scalar_rhs
	call MPI_Barrier(MPI_COMM_WORLD, ierr)
	t5 = t5 + MPI_Wtime() - tt

	call MPI_Barrier(MPI_COMM_WORLD, ierr)	
        tt = MPI_Wtime()
	call bed_shear
	if (sedi .eq. 1) then
	   call bottom_mean_stress
	   call pickup
	endif
	call MPI_Barrier(MPI_COMM_WORLD, ierr)
	t5 = t5 + MPI_Wtime() - tt
        
	call MPI_Barrier(MPI_COMM_WORLD, ierr)
	tt = MPI_Wtime()
	   call predictor
	t1 = t1 + MPI_Wtime() - tt
	
	call MPI_Barrier(MPI_COMM_WORLD, ierr)	
	if ( sedi .eq. 1 ) then		     
             call sedi_rhs
	endif
	call MPI_Barrier(MPI_COMM_WORLD, ierr)
	
	call MPI_Barrier(MPI_COMM_WORLD, ierr)
	tt = MPI_Wtime()
	   call pres_corr
	   call pressure
	call MPI_Barrier(MPI_COMM_WORLD, ierr)
	t2 = t2 + MPI_Wtime() - tt

	call MPI_Barrier(MPI_COMM_WORLD, ierr)
	tt = MPI_Wtime()
	   call corrector
	call MPI_Barrier(MPI_COMM_WORLD, ierr)
	t3 = t3 + MPI_Wtime() - tt
	
	call MPI_Barrier(MPI_COMM_WORLD, ierr)
	tt = MPI_Wtime()
	   if ( iscalar .eq. 1 ) call scalar_solve
	call MPI_Barrier(MPI_COMM_WORLD, ierr)
	t4 = t4 + MPI_Wtime() - tt
		
        call MPI_Barrier(MPI_COMM_WORLD, ierr)	
	t5 = t5 + MPI_Wtime() - tt

	call MPI_Barrier(MPI_COMM_WORLD, ierr)
	tt = MPI_Wtime()
	  if ( sedi .eq. 1 ) call sedi_solve	
	call MPI_Barrier(MPI_COMM_WORLD, ierr)
	t4 = t4 + MPI_Wtime() - tt
	
	if(mod(istep,nsave) .eq. 0) then
	   call MPI_Barrier(MPI_COMM_WORLD, ierr)
	   tt = MPI_Wtime()
	   call output
	   call MPI_Barrier(MPI_COMM_WORLD, ierr)
	   t6 = t6 + MPI_Wtime() - tt
	end if
c        call TimeAvrg
c        call getStress
	kount = kount + 1
	time = time + dtime
	call cfl_check

	enddo
        
        call storeData
	call MPI_Barrier(MPI_COMM_WORLD, ierr)
	ta = MPI_Wtime() - ta

	if ( MYID .eq. 0 ) then
     	   WRITE(*,*) ' Total CPU time = ', ta, ta/nstep		
     	   WRITE(*,*) ' eddy:      ', t0/nstep, int(t0/ta*100.D0), '%'		
     	   WRITE(*,*) ' Predictor: ', t1/nstep, int(t1/ta*100.D0), '%'		
     	   WRITE(*,*) ' Pressure:  ', t2/nstep, int(t2/ta*100.D0), '%'		
     	   WRITE(*,*) ' Corrector: ', t3/nstep, int(t3/ta*100.D0), '%'		
     	   WRITE(*,*) ' Scalarp:   ', t5/nstep, int(t5/ta*100.D0), '%' 		
     	   WRITE(*,*) ' Scalars:   ', t4/nstep, int(t4/ta*100.D0), '%' 		
     	   WRITE(*,*) ' Output:    ', t6/nstep, int(t6/ta*100.D0), '%' 		
	endif

	call MPI_FINALIZE( ierr )

	stop
	end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine getStress
         

        include "size.inc"
	include "mpif.h"
	include "mpi.inc"
	include "para.inc"
	include "sedi.inc"
        include "eddy.inc"
        include "metric.inc"
        include "ns.inc"
 
        integer i,j,k

        do k = 1, nnk
	do j = 1, nnj
        do i = 1, nni
           
           taus(i,j,k,2) = taus(i,j,k,2) + vst(i,j,k)
     <               * ( xiy(i,j,k)*0.5D0*(u(i+1,j,k,1)-u(i-1,j,k,1))
     <                 + ety(i,j,k)*0.5D0*(u(i,j+1,k,1)-u(i,j-1,k,1))
     <                 + zty(i,j,k)*0.5D0*(u(i,j,k+1,1)-u(i,j,k-1,1)))
         
           taus(i,j,k,1) = taus(i,j,k,1) + vis
     <               * ( xiy(i,j,k)*0.5D0*(u(i+1,j,k,1)-u(i-1,j,k,1))
     <                 + ety(i,j,k)*0.5D0*(u(i,j+1,k,1)-u(i,j-1,k,1))
     <                 + zty(i,j,k)*0.5D0*(u(i,j,k+1,1)-u(i,j,k-1,1)))       
c	   write(*,*) taus(i,j,k,2), taus(i,j,k,1)
        enddo
	enddo
	enddo

        return
        end
