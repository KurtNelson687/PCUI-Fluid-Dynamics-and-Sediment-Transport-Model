cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine parameter

	include "size.inc"
        include "mpif.h"
        include "mpi.inc"
	include "para.inc"
	include "sed.inc"

	integer :: i
	character(len=100) :: buffer, label
	integer :: pos
	integer, parameter :: fh = 15
	integer :: ios = 0
	integer :: line = 0

	call MPI_COMM_RANK( MPI_COMM_WORLD, MYID, ierr )

	if ( MYID .eq. 0 ) then

	open(fh,file='pcuiRunParams.txt')
	
	do while(ios == 0)
	  read(fh, '(A)',iostat=ios) buffer
	  if (ios == 0) then
	    line = line +1
	    pos = scan(buffer,' ')
	    label = buffer(1:pos)
	    pos = scan(buffer, '=')
	    buffer = buffer(pos+1:)
	    
	    select case (label)
	    case ('dtime')
	      read(buffer, *, iostat=ios) dtime
	      write(*,*) 'dtime = ', dtime
	    case ('runCase')
	      read(buffer, *, iostat=ios) runCase
	      write(*,*) 'runCase = ', runCase
	    case ('newrun')
	      read(buffer, *, iostat=ios) newrun
	      write(*,*) 'newrun  = ', newrun
	    case ('periodic')
	      read(buffer, *, iostat=ios) periodic
	      write(*,*) 'periodic  = ', periodic
	    case ('iscalar')
	      read(buffer, *, iostat=ios) iscalar
	      write(*,*) 'iscalar  = ', iscalar
	    case ('ieddy')
	      read(buffer, *, iostat=ios) ieddy
	      write(*,*) 'ieddy  = ', ieddy
	    case ('ised')
	      read(buffer, *, iostat=ios) ised
	      write(*,*) 'ised  = ', ised
	    case ('waves')
	      read(buffer, *, iostat=ios) waves
	      write(*,*) 'waves  = ', waves
	    case ('mg_level')
	      read(buffer, *, iostat=ios) mg_level
	      write(*,*) 'mg_level  = ', mg_level
	    case ('nstep')
	      read(buffer, *, iostat=ios) nstep
	      write(*,*) 'nstep = ', nstep
	    case ('nsave')
	      read(buffer, *, iostat=ios) nsave
	      write(*,*) 'nsave  = ', nsave
	    case ('maxstep')
	      read(buffer, *, iostat=ios) maxstep
	      write(*,*) 'maxstep  = ', maxstep
	    case ('iterchk')
	      read(buffer, *, iostat=ios) iterchk
	      write(*,*) 'iterchk  = ', iterchk
	    case ('maxiter')
	      read(buffer, *, iostat=ios) maxiter
	      write(*,*) 'maxiter = ', maxiter
	    case ('tol')
	      read(buffer, *, iostat=ios) tol
	      write(*,*) 'tol  = ', tol
	    case ('ter')
	      read(buffer, *, iostat=ios) ter
	      write(*,*) 'ter  = ', ter
	    case ('slowiter')
	      read(buffer, *, iostat=ios) slowiter
	      write(*,*) 'slowiter  = ', slowiter
	    case ('vis')
	      read(buffer, *, iostat=ios) vis
	      write(*,*) 'vis = ', vis
	    case ('rhoWater')
	      read(buffer, *, iostat=ios) rhoWater
	      write(*,*) 'rhoWater  = ', rhoWater
	    case ('ak')
	      read(buffer, *, iostat=ios) ak
	      write(*,*) 'ak  = ', ak
	    case ('g')
	      read(buffer, *, iostat=ios) g
	      write(*,*) 'g  = ', g
	    case ('dpdxSteady')
	      read(buffer, *, iostat=ios) dpdxSteady
	      write(*,*) 'dpdxSteady = ', dpdxSteady
	    case ('dpdxLaminar')
	      read(buffer, *, iostat=ios) dpdxLaminar
	      write(*,*) 'dpdxLaminar = ', dpdxLaminar
	    case ('driveFac')
	      read(buffer, *, iostat=ios) driveFac
	      write(*,*) 'driveFac = ', driveFac
	    case ('waveMag')
	      read(buffer, *, iostat=ios) waveMag
	      write(*,*) 'waveMag  = ', waveMag
	    case ('Twave')
	      read(buffer, *, iostat=ios) Twave
	      write(*,*) 'Twave  = ', Twave
	    case ('omg_cyl')
	      read(buffer, *, iostat=ios) omg_cyl
	      write(*,*) 'omg_cyl  = ', omg_cyl
	    case ('omg_lid')
	      read(buffer, *, iostat=ios) omg_lid
	      write(*,*) 'omg_lid = ', omg_lid
	    case ('factor')
	      read(buffer, *, iostat=ios) factor
	      write(*,*) 'factor  = ', factor
	    case ('phi1')
	      read(buffer, *, iostat=ios) phi1
	      write(*,*) 'phi1  = ', phi1
	    case ('phi2')
	      read(buffer, *, iostat=ios) phi2
	      write(*,*) 'phi2  = ', phi2
	    case ('yphi')
	      read(buffer, *, iostat=ios) yphi
	      write(*,*) 'yphi = ', yphi
	    case ('aphi')
	      read(buffer, *, iostat=ios) aphi
	      write(*,*) 'aphi  = ', aphi
	    case ('ws')
	      read(buffer, *, iostat=ios) ws
	      write(*,*) 'ws = ', ws
	    case ('rhoSed')
	      read(buffer, *, iostat=ios) rhoSed
	      write(*,*) 'rhoSed  = ', rhoSed
	    case ('pAdjust')
	      read(buffer, *, iostat=ios) pAdjust
	      write(*,*) 'pAdjust  = ', pAdjust
	    case default
	      write(*,*) 'skipping invalid label at line', line
	    end select

	  end if
	end do

	endif

	call MPI_BCAST(runCase,     1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(newrun,      1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(periodic,    1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(pAdjust,     1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(waves,       1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(iscalar,     1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(ieddy,       1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(mg_level,    1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(nstep,       1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(nsave,       1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

	call MPI_BCAST(maxstep,     1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(iterchk(1),  5,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(maxiter(1),  5,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        call MPI_BCAST(dtime,       1,MPI_DOUBLE_PRECISION,0,
     <	                            MPI_COMM_WORLD,ierr)
	call MPI_BCAST(vis,         1,MPI_DOUBLE_PRECISION,0,
     <                              MPI_COMM_WORLD,ierr)
	call MPI_BCAST(rhoWater,         1,MPI_DOUBLE_PRECISION,0,
     <                              MPI_COMM_WORLD,ierr)
	call MPI_BCAST(ak,          1,MPI_DOUBLE_PRECISION,0,
     <                              MPI_COMM_WORLD,ierr)
	call MPI_BCAST(g,           1,MPI_DOUBLE_PRECISION,0,
     <                              MPI_COMM_WORLD,ierr)
	call MPI_BCAST(dpdxSteady,  1,MPI_DOUBLE_PRECISION,0,
     <                              MPI_COMM_WORLD,ierr)
	call MPI_BCAST(dpdxLaminar,  1,MPI_DOUBLE_PRECISION,0,
     <                              MPI_COMM_WORLD,ierr)
	call MPI_BCAST(driveFac,  1,MPI_DOUBLE_PRECISION,0,
     <                              MPI_COMM_WORLD,ierr)
	call MPI_BCAST(Twave,       1,MPI_DOUBLE_PRECISION,0,
     <                              MPI_COMM_WORLD,ierr)
	call MPI_BCAST(waveMag,    1,MPI_DOUBLE_PRECISION,0,
     <                              MPI_COMM_WORLD,ierr)
	call MPI_BCAST(omg_cyl,     1,MPI_DOUBLE_PRECISION,0,
     <                              MPI_COMM_WORLD,ierr)
	call MPI_BCAST(omg_lid,     1,MPI_DOUBLE_PRECISION,0,
     <                              MPI_COMM_WORLD,ierr)

	call MPI_BCAST(factor,      1,MPI_DOUBLE_PRECISION,0,
     <                              MPI_COMM_WORLD,ierr)
	call MPI_BCAST(tol(1),      5,MPI_DOUBLE_PRECISION,0,
     <                              MPI_COMM_WORLD,ierr)
	call MPI_BCAST(ter(1),      5,MPI_DOUBLE_PRECISION,0,
     <                              MPI_COMM_WORLD,ierr)
	call MPI_BCAST(slowiter(1), 5,MPI_DOUBLE_PRECISION,0,
     <                              MPI_COMM_WORLD,ierr)

	call MPI_BCAST(phi1,        1,MPI_DOUBLE_PRECISION,0,
     <                              MPI_COMM_WORLD,ierr)
	call MPI_BCAST(phi2,        1,MPI_DOUBLE_PRECISION,0,
     <                              MPI_COMM_WORLD,ierr)
	call MPI_BCAST(yphi,        1,MPI_DOUBLE_PRECISION,0,
     <                              MPI_COMM_WORLD,ierr)
	call MPI_BCAST(aphi,        1,MPI_DOUBLE_PRECISION,0,
     <                              MPI_COMM_WORLD,ierr)
        call MPI_BCAST(ws,          1,MPI_DOUBLE_PRECISION,0,
     <                              MPI_COMM_WORLD,ierr)
	call MPI_BCAST(rhoSed,         1,MPI_DOUBLE_PRECISION,0,
     <                              MPI_COMM_WORLD,ierr)
        call MPI_BCAST(ised,        1,MPI_DOUBLE_PRECISION,0,
     <                              MPI_COMM_WORLD,ierr)
	omg2 = omg_cyl * 2.D0
	return
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C	This subroutine is simply writing the x,y, and z coordinates from each processor
        subroutine output_xyz
        
        include "size.inc"
        include "mpif.h"
        include "mpi.inc"
	include "metric.inc"

	character*4 ID 

	if (100+myid .gt. 999) then
	write(ID, fmt='(I4)') 100+myid
	else
	write(ID, fmt='(I3)') 100+myid
	endif
	
	open(500+myid,file='xyz.'//ID,form='unformatted',status='unknown')
	write(500+myid) xp
        close(unit = 500+myid)
	
	return
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine output

	include "size.inc"
        include "mpif.h"
        include "mpi.inc"
	include "para.inc"
	include "ns.inc"
	include "metric.inc"
	include "eddy.inc"
	include "sed.inc"
	
	character*4 :: ID

	if (100+myid .gt. 999) then
	write(ID, fmt='(I4)') 100+myid
	else
	write(ID, fmt='(I3)') 100+myid
	endif

	if (myid .eq. 0) print *, 'Recording to output files ...'

!	if (kount.gt.1) then !appends to existing file - use this if you want continue run to append of existing output
	if (istep .gt.1) then !appends to existing file - use this if you want continue run to create new files
	   open(50+myid, file='output_rho.'//ID, form='unformatted',
     >          status='old',position='append')
	   write(50+myid) rho
	   close(unit = 50+myid)

	   open(200+myid, file='output_UVW.'//ID, form='unformatted',
     >          status='old',position='append')
	   write(200+myid) u
	   close(unit = 200+myid)

	   open(800+myid, file='output_p.'//ID, form='unformatted',
     >          status='old',position='append')
	   write(800+myid) p
	   close(unit = 800+myid)

	   if (ised .eq. 1) then
	       open(900+myid, file='output_Csed.'//ID, form='unformatted',
     >             status='old',position='append')
	       write(900+myid) Csed
	       close(unit = 900+myid)
	   endif

	else
	   open(50+myid, file='output_rho.'//ID, form='unformatted',
     >          status='unknown')	 
	   write(50+myid) rho
	   close(unit = 50+myid)

	   open(200+myid, file='output_UVW.'//ID, form='unformatted',
     >          status='unknown')	
	   write(200+myid) u
	   close(unit = 200+myid)

	   open(800+myid, file='output_p.'//ID, form='unformatted',
     >          status='unknown')	
	   write(800+myid) p
	   close(unit = 800+myid)

	   if (ised .eq. 1) then
	        open(900+myid, file='output_Csed.'//ID, 
     >         form='unformatted',status='unknown') 
	        write(900+myid) Csed
	        close(unit = 900+myid)
	   endif   
	end if

	return
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine input_continue_run

	include "size.inc"
        include "mpif.h"
        include "mpi.inc"
	include "para.inc"
	include "ns.inc"
	include "metric.inc"
	include "eddy.inc"
	include "sed.inc"
	include "padjust.inc"

	character*4 :: ID

	if (100+myid .gt. 999) then
	write(ID, fmt='(I4)') 100+myid
	else
	write(ID, fmt='(I3)') 100+myid
	endif
	
	open(200+myid, file='continue_run.'//ID, form='unformatted',
     >          status='old')	   

	read(200+myid) kount
	read(200+myid) vis, ak
	read(200+myid) time
	read(200+myid) u
	read(200+myid) uxi
	read(200+myid) uej
	read(200+myid) uzk
	read(200+myid) p
	read(200+myid) rho
	read(200+myid) steadyPall

	if ( iscalar .eq. 1 ) read(200+myid) phi

	if ( ieddy .gt. 0 ) then
	   read(200+myid) vst
	   read(200+myid) akst
	   read(200+myid) sab
	   read(200+myid) jac
	endif

	if ( ised .eq. 1 ) read(200+myid) Csed

	close(200+myid)

	return
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine output_continue_run

	include "size.inc"
        include "mpif.h"
        include "mpi.inc"
	include "para.inc"
	include "ns.inc"
	include "metric.inc"
	include "eddy.inc"
	include "sed.inc"
	include "padjust.inc"


	character*4 :: ID

	if (100+myid .gt. 999) then
	write(ID, fmt='(I4)') 100+myid
	else
	write(ID, fmt='(I3)') 100+myid
	endif

	write(ID, fmt='(I3)') 100+myid

	open(200+myid, file='continue_run.'//ID, form='unformatted',
     >          status='unknown')	 

	write(200+myid) kount
	write(200+myid) vis, ak
	write(200+myid) time
	write(200+myid) u
	write(200+myid) uxi
	write(200+myid) uej
	write(200+myid) uzk
	write(200+myid) p
	write(200+myid) rho
	write(200+myid) steadyPall

	if ( iscalar .eq. 1 ) write(200+myid) phi

	if ( ieddy .gt. 0 ) then
	   write(200+myid) vst
	   write(200+myid) akst
	   write(200+myid) sab
	   write(200+myid) jac
	endif
	
	if (ised .eq. 1 ) write(200+myid) Csed
	
	close(200+myid)

	return
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine output_Utheo
	include "size.inc"
        include "mpif.h"
        include "mpi.inc"
	include "padjust.inc"
	
	character*4 :: ID

	if (100+vert_id .gt. 999) then
	write(ID, fmt='(I4)') 100+vert_id
	else
	write(ID, fmt='(I3)') 100+vert_id
	endif

	if (myid .eq. 0) print *, 'writing log law to files ...'
	
	if (hor_id .eq. 0) then
	   open(50+myid, file='output_uTheo.'//ID, form='unformatted',
     >          status='unknown')	 
	   write(50+myid) uTheo 
	   close(unit = 50+myid)
	endif

	return
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine output_profiles
	include "size.inc"
	include "ns.inc"
        include "mpif.h"
        include "mpi.inc"
	include "padjust.inc"
	include "para.inc"
	include "eddy.inc"
	character*4 :: ID
	double precision, dimension(1:nnj) ::
     <     uTurb, vTurb, wTurb,uvRey, uwRey, vwRey, vstMean, rruMean,
     <     kineticMean, dissipationMean
	double precision kineticDepth, drive 
	double precision, dimension(-1:nni+2,-1:nnj+2,-1:nnk+2) :: 
     <     kinetic

C	call horizontalAverage(rr(:,:,:,1), rruMean, 1)
C	call horizontalAverage(vst, vstMean, 1)
	
	call get_turbIntensity(uTurb, vTurb, wTurb)
	call get_reynoldsStress(uvRey, uwRey, vwRey)
	call compute_totalkinetic(u,kinetic)
	call horizontalAverage(kinetic, kineticMean, 2)
	call depthAverage(kineticMean, kineticDepth)

	call compute_dissipation(dissipationMean)

	drive = steadyPall(1) 
	if (100+vert_id .gt. 999) then
	write(ID, fmt='(I4)') 100+vert_id
	else
	write(ID, fmt='(I3)') 100+vert_id
	endif


	if (myid .eq. 0) then
	   if (istep .gt.1) then !appends to existing file - use this if you want continue run to create new files
	     open(50+myid, file='outputp_time.'//ID, form='unformatted',
     >          status='old',position='append')
	     write(50+myid) time
	     close(unit = 50+myid)
	     
	open(50+myid, file='outputpVal_kineticdepth.'//ID, 
     >    form='unformatted',status='old',position='append')
	     write(50+myid) kineticDepth
	     close(unit = 50+myid)

	     open(50+myid, file='outputpVal_drive.'//ID,
     >   form='unformatted', status='old',position='append')
	     write(50+myid) drive
	     close(unit = 50+myid)
	   else

	     open(50+myid, file='outputp_time.'//ID, form='unformatted',
     >          status='unknown')
	     write(50+myid) time
	     close(unit = 50+myid)

	    open(50+myid, file='outputpVal_kineticdepth.'//ID,
     >           form='unformatted', status='unknown')
	     write(50+myid) kineticDepth
	     close(unit = 50+myid)
	     
	     open(50+myid, file='outputpVal_drive.'//ID, 
     >            form='unformatted', status='unknown')
	     write(50+myid) drive
	     close(unit = 50+myid)

	   endif
	endif

	
	if (hor_id .eq. 0) then
	   if (istep .gt.1) then !appends to existing file - use this if you want continue run to create new files
	     open(50+myid, file='outputp_dissmean.'//ID, 
     >          form='unformatted', status='old',position='append')
	     write(50+myid) dissipationMean
	     close(unit = 50+myid)

C	     open(50+myid, file='outputp_dpdx.'//ID, form='unformatted',
C     >          status='old',position='append')
C	     write(50+myid) steadyPall
C	     close(unit = 50+myid)

	     open(50+myid, file='outputp_uTurb.'//ID, form='unformatted',
     >          status='old',position='append')
	     write(50+myid) uTurb
	     close(unit = 50+myid)

	     open(50+myid, file='outputp_vTurb.'//ID, form='unformatted',
     >          status='old',position='append')
	     write(50+myid) vTurb
	     close(unit = 50+myid)

	     open(50+myid, file='outputp_wTurb.'//ID, form='unformatted',
     >          status='old',position='append')
	     write(50+myid) wTurb
	     close(unit = 50+myid)

	     open(50+myid, file='outputp_uvRey.'//ID, form='unformatted',
     >          status='old',position='append')
	     write(50+myid) uvRey
	     close(unit = 50+myid)

	     open(50+myid, file='outputp_uwRey.'//ID, form='unformatted',
     >          status='old',position='append')
	     write(50+myid) uwRey
	     close(unit = 50+myid)

	     open(50+myid, file='outputp_vwRey.'//ID, form='unformatted',
     >          status='old',position='append')
	     write(50+myid) vwRey
	     close(unit = 50+myid)

C	     open(50+myid, file='outputp_rruMean.'//ID, form='unformatted',
C     >          status='old',position='append')
C	     write(50+myid) rruMean
C	     close(unit = 50+myid)

C	     open(50+myid, file='outputp_vstMean.'//ID, form='unformatted',
C     >          status='old',position='append')
C	     write(50+myid) vstMean
C	     close(unit = 50+myid)

	   else
	     open(50+myid, file='outputp_dissmean.'//ID, 
     >          form='unformatted',status='unknown')
	     write(50+myid) dissipationMean
	     close(unit = 50+myid)

C	     open(50+myid, file='outputp_dpdx.'//ID, form='unformatted',
C     >          status='unknown')
C	     write(50+myid) steadyPall
C	     close(unit = 50+myid)

	     open(50+myid, file='outputp_uTurb.'//ID, form='unformatted',
     >          status='unknown')
	     write(50+myid) uTurb
	     close(unit = 50+myid)

	     open(50+myid, file='outputp_vTurb.'//ID, form='unformatted',
     >          status='unknown')
	     write(50+myid) vTurb
	     close(unit = 50+myid)

	     open(50+myid, file='outputp_wTurb.'//ID, form='unformatted',
     >          status='unknown')
	     write(50+myid) wTurb
	     close(unit = 50+myid)

	     open(50+myid, file='outputp_uvRey.'//ID, form='unformatted',
     >          status='unknown')
	     write(50+myid) uvRey
	     close(unit = 50+myid)

	     open(50+myid, file='outputp_uwRey.'//ID, form='unformatted',
     >          status='unknown')
	     write(50+myid) uwRey
	     close(unit = 50+myid)

	     open(50+myid, file='outputp_vwRey.'//ID, form='unformatted',
     >          status='unknown')
	     write(50+myid) vwRey
	     close(unit = 50+myid)

C	     open(50+myid, file='outputp_rruMean.'//ID, form='unformatted',
C     >          status='unknown')
C	     write(50+myid) rruMean
C	     close(unit = 50+myid)

C	     open(50+myid, file='outputp_vstMean.'//ID, form='unformatted',
C     >          status='unknown')
C	     write(50+myid) vstMean
C	     close(unit = 50+myid)

	   endif
	endif

	return
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine cfl_check

	include "size.inc"
        include "mpif.h"
	include "mpi.inc"
	include "para.inc"
	include "ns.inc"
	include "metric.inc"

	integer*4 today(3), now(3)
	integer i, j, k
	double precision cflmax, cfl, temp 

	cflmax = 0.D0
	cfl = 0.D0

	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni
           temp = dabs( uxi(i-1,j,k) + uxi(i,j,k) ) 
     <          + dabs( uej(i,j-1,k) + uej(i,j,k) )
     <          + dabs( uzk(i,j,k-1) + uzk(i,j,k) )
	   temp = temp * jac(i,j,k)
           cfl = max(cfl, temp)
	enddo
	enddo
	enddo

	call MPI_REDUCE(cfl,cflmax,1,MPI_DOUBLE_PRECISION,
     <	                MPI_MAX,0,comm3d,ierr)

	if ( myid .eq. 0 ) cflmax = 0.5D0 * dtime * cflmax

	call MPI_BCAST(cflmax,1,MPI_DOUBLE_PRECISION,0,comm3d,ierr)


	if ( MYID .eq. 0 ) then

	   if (istep .gt.1) then !appends to existing file - use this if you want continue run to create new files
	     open(50+myid, file='outputpval_cfl.100', form='unformatted',
     >          status='old',position='append')
	     write(50+myid) cflmax
	     close(unit = 50+myid)
	   else
	     open(50+myid, file='outputpval_cfl.100', form='unformatted',
     >          status='unknown')
	     write(50+myid) cflmax
	     close(unit = 50+myid)
	   endif
	endif


	if ( cflmax .gt. 0.9 .and. istep .gt. 5) then
	   if ( MYID .eq. 0 )
     <	      write(*,*) ' CFL = ', cflmax, ' > 0.9, Stop! '
	   stop
	else
	   if ( MYID .eq. 0 ) then
	      print *, ''
	      print *, '************************************************' 
     	      write(*,1) cflmax, dtime, time
	      print *, '************************************************'
	      print *, ''
	   end if
	endif
    1	format(' CFL = ', e15.9, '  dtime = ', e9.3, '  time = ', e9.3)

	if ( MYID .eq. 0) then
	   call idate(today)	! today(1)=day, (2)=month, (3)=year
	   call itime(now)	! now(1)=hour, (2)=minute, (3)=second
	   write ( *, 1000 )  today(2), today(1), today(3), now
	end if
 1000	format ( 'Date ', i2.2, '/', i2.2, '/', i4.4, '; time ',
     &	               i2.2, ':', i2.2, ':', i2.2 )

	   if (MYID .eq. 0) then
	      if(istep.eq.1) then
		 open(123, file='qoutput',status='replace')
	      else
		 open(123, file='qoutput',position='append')
	      endif

	      write(123,105) istep, nstep, time, cflmax	   
	      call idate(today)	! today(1)=day, (2)=month, (3)=year
	      call itime(now)	! now(1)=hour, (2)=minute, (3)=second
 105          format('Time step = ',i7,'/',i7', time = ',e10.3, 
     &                ', cfl = ',e10.3)
	      close(unit=123)
	   endif

	return
	end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine eqstate

	include "size.inc"
	include "para.inc"
	include "ns.inc"
	include "sed.inc"

	integer i, j, k
	

	do k = 1, nnk
	do j = 1, nnj
	do i = 1, nni
	   rho(i,j,k) = rhoWater+(1-rhoWater/rhoSed)*Csed(i,j,K)
	enddo
	enddo
	enddo

	
	return
	end

