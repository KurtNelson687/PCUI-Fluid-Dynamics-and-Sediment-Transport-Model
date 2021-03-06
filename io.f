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
	    case ('newWave')
	      read(buffer, *, iostat=ios) newWave
	      write(*,*) 'newWave  = ', newWave
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
	    case ('irho')
	      read(buffer, *, iostat=ios) irho
	      write(*,*) 'ised  = ', irho
	    case ('iTKE')
	      read(buffer, *, iostat=ios) iTKE
	      write(*,*) 'iTKE  = ', iTKE
	    case ('waves')
	      read(buffer, *, iostat=ios) waves
	      write(*,*) 'waves  = ', waves
	    case ('mg_level')
	      read(buffer, *, iostat=ios) mg_level
	      write(*,*) 'mg_level  = ', mg_level
	    case ('nstep')
	      read(buffer, *, iostat=ios) nstep
	      write(*,*) 'nstep = ', nstep
	    case ('nsavePro')
	      read(buffer, *, iostat=ios) nsavePro
	      write(*,*) 'nsavePro  = ', nsavePro
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
	    case ('dpdxWave')
	      read(buffer, *, iostat=ios) dpdxWave
	      write(*,*) 'dpdxWave  = ', dpdxWave
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
	    case ('dryBulk')
	      read(buffer, *, iostat=ios) dryBulk
	      write(*,*) 'dryBulk  = ', dryBulk
	    case ('tauCrit')
	      read(buffer, *, iostat=ios) tauCrit
	      write(*,*) 'tauCrit  = ', tauCrit
	    case ('Ased')
	      read(buffer, *, iostat=ios) Ased
	      write(*,*) 'Ased  = ', Ased
	    case ('nsed')
	      read(buffer, *, iostat=ios) nsed
	      write(*,*) 'nsed  = ', nsed
	    case ('pAdjust')
	      read(buffer, *, iostat=ios) pAdjust
	      write(*,*) 'pAdjust  = ', pAdjust
	    case ('proFileName')
	      read(buffer, *, iostat=ios) proFileName
	      write(*,*) 'proFileName = ', proFileName
	    case default
	      write(*,*) 'skipping invalid label at line', line
	    end select

	  end if
	end do

	endif

	call MPI_BCAST(runCase,     1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(newWave,      1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(newrun,      1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(periodic,    1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(pAdjust,     1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(waves,       1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(iscalar,     1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(ieddy,       1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(mg_level,    1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(nstep,       1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(nsave,       1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(nsavePro,    1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

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
	call MPI_BCAST(Twave,       1,MPI_DOUBLE_PRECISION,0,
     <                              MPI_COMM_WORLD,ierr)
	call MPI_BCAST(dpdxWave,    1,MPI_DOUBLE_PRECISION,0,
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
	call MPI_BCAST(dryBulk,         1,MPI_DOUBLE_PRECISION,0,
     <                              MPI_COMM_WORLD,ierr)
	call MPI_BCAST(tauCrit,         1,MPI_DOUBLE_PRECISION,0,
     <                              MPI_COMM_WORLD,ierr)
	call MPI_BCAST(Ased,         1,MPI_DOUBLE_PRECISION,0,
     <                              MPI_COMM_WORLD,ierr)
	call MPI_BCAST(nsed,         1,MPI_DOUBLE_PRECISION,0,
     <                              MPI_COMM_WORLD,ierr)
        call MPI_BCAST(ised,        1,MPI_DOUBLE_PRECISION,0,
     <                              MPI_COMM_WORLD,ierr)
        call MPI_BCAST(irho,        1,MPI_DOUBLE_PRECISION,0,
     <                              MPI_COMM_WORLD,ierr)
        call MPI_BCAST(iTKE,        1,MPI_DOUBLE_PRECISION,0,
     <                              MPI_COMM_WORLD,ierr)
        call MPI_BCAST(proFileName, LEN(proFileName),
     <             MPI_CHARACTER,0, MPI_COMM_WORLD,ierr)
	omg2 = omg_cyl * 2.D0
	phasePerT = int(Twave/dtime)/nsavePro
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

	subroutine output(PI)

	include "size.inc"
        include "mpif.h"
        include "mpi.inc"
	include "para.inc"
	include "ns.inc"
	include "metric.inc"
	include "eddy.inc"
	include "sed.inc"

	double precision, intent(in) :: PI
	character*4 :: ID
	double precision :: phase

	if (100+myid .gt. 999) then
	write(ID, fmt='(I4)') 100+myid
	else
	write(ID, fmt='(I3)') 100+myid
	endif

	if (myid .eq. 0) then
	   print *, 'Recording to output files ...'
	
	   if(waves .eq. 1) then
	     phase = cos(2*PI*time/Twave)
	   if (istep .gt.1) then !appends to existing file - use this if you want continue run to create new files
	     open(50+myid, file='output_wavephase.'//ID,
     >   form='unformatted', status='old',position='append')
	     write(50+myid) phase
	     close(unit = 50+myid)
	     open(50+myid, file='output_time.'//ID, form='unformatted',
     >          status='old',position='append')
	     write(50+myid) time
	     close(unit = 50+myid)
	   else
	     open(50+myid, file='output_wavephase.'//ID,
     >            form='unformatted', status='unknown')
	     write(50+myid) phase
	     close(unit = 50+myid)
	     open(50+myid, file='output_time.'//ID,
     >            form='unformatted', status='unknown')
	     write(50+myid) time
	     close(unit = 50+myid)
	   endif
	   endif
	endif

!	if (kount.gt.1) then !appends to existing file - use this if you want continue run to append of existing output
	if (istep .gt.1) then !appends to existing file - use this if you want continue run to create new files
C	   open(50+myid, file='output_rho.'//ID, form='unformatted',
C     >          status='old',position='append')
C	   write(50+myid) rho
C	   close(unit = 50+myid)

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
C	   open(50+myid, file='output_rho.'//ID, form='unformatted',
C     >          status='unknown')	 
C	   write(50+myid) rho
C	   close(unit = 50+myid)

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
	read(200+myid) numPhase

	if ( iscalar .eq. 1 ) read(200+myid) phi

	if ( ieddy .gt. 0 ) then
	   read(200+myid) vst
	   read(200+myid) akst
	   read(200+myid) sab
	   read(200+myid) jac
	endif
	
	if(newWave .ne. 1) then
	   if ( ised .eq. 1 ) read(200+myid) Csed
	endif

	close(200+myid)
	if(myid .eq. 0) write(*,*) 'time = ', time
	return
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine readMeanProfiles
C This subroutine reads in a binary file contianing a mean streamwise velocity profile and assigns
C a each processor a portion of the profile depending on its vertical position
	include "size.inc"
        include "mpif.h"
        include "mpi.inc"
	include "para.inc"
	include "ns.inc"
	include "metric.inc"
	include "padjust.inc"
	double precision, dimension(0:nj+1,1:phasePerT) ::
     <     uMeanAveAll
	integer jj0, j

	jj0 = npy * nnj !index pointing to first non-ghost cell mean velocity for processor 
	
	!Read bindary file containning mean velcoity profile for streamwise direction and top and bottom ghost cells
	open(700+myid, file = proFileName,
     <          form='unformatted',status='unknown')
	read(700+myid) uMeanAveAll
	close(700+myid)

	do j = 0, nnj+1 ! Assign mean velocity profile to processor
	   uMeanAve(j,:) = uMeanAveAll(jj0+j,:) 
	enddo
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
	write(200+myid) numPhase

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

	subroutine output_profiles(PI)
	include "size.inc"
	include "ns.inc"
        include "mpif.h"
        include "mpi.inc"
	include "padjust.inc"
	include "para.inc"
	include "eddy.inc"
	include "sed.inc"

	double precision, intent(in) :: PI
	character*4 :: ID
	double precision, dimension(1:nnj) ::
     <     uTurb, vTurb, wTurb,uvRey, uwRey, vwRey, vstMean, rruMean,
     <     kineticMean, dissipationMean, sedMean, vCsed,RiMean, 
     <     BruntNSqr, rhoMean, rhoSqrMean, rhoPrimeMean 
	double precision kineticDepth, drive, sedTotal1, sedTotal2,
     <     phase 
	double precision, dimension(-1:nni+2,-1:nnj+2,-1:nnk+2) :: 
     <     kinetic

	call horizontalAverage(rr(:,:,:,1), rruMean, 1)
	call horizontalAverage(vst, vstMean, 1)

	if(iTKE .ne.1) then
	call computeMeanAndPrimes
	endif
	call get_turbIntensity(uTurb, vTurb, wTurb)
	call get_reynoldsStress(uvRey, uwRey, vwRey)
	call compute_totalkinetic(u,kinetic)
	call horizontalAverage(kinetic, kineticMean, 2)
	call depthAverage(kineticMean, kineticDepth)
	
	if (ised .eq. 1) then
	call horizontalAverage(Csed, sedMean, 2)
	call sedMass1(Csed,sedTotal1,2)
	call SedMass2(Csed,sedTotal2,2)
	call get_sedTurbFlux(sedMean,vCsed)
	endif
	
	call get_BruntN(rho,BruntNSqr,rhoMean,rhoSqrMean,rhoPrimeMean)


	if (waves .eq. 1) then
	     phase = cos(2*PI*time/Twave)
	endif

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

	     open(50+myid, file='outputpVal_kount.'//ID,
     >   form='unformatted', status='old',position='append')
	     write(50+myid) DBLE(kount)
	     close(unit = 50+myid)
	     
	open(50+myid, file='outputpVal_kineticdepth.'//ID, 
     >    form='unformatted',status='old',position='append')
	     write(50+myid) kineticDepth
	     close(unit = 50+myid)

	     open(50+myid, file='outputpVal_drive.'//ID,
     >   form='unformatted', status='old',position='append')
	     write(50+myid) drive
	     close(unit = 50+myid)

	   if(ised .eq. 1) then
	     open(50+myid, file='outputpVal_sedTotal1.'//ID,
     >   form='unformatted', status='old',position='append')
	     write(50+myid) sedTotal1
	     close(unit = 50+myid)

	     open(50+myid, file='outputpVal_sedTotal2.'//ID,
     >   form='unformatted', status='old',position='append')
	     write(50+myid) sedTotal2
	     close(unit = 50+myid)
	  endif

	   if(waves .eq. 1) then
	     open(50+myid, file='outputpVal_wavephase.'//ID,
     >   form='unformatted', status='old',position='append')
	     write(50+myid) phase
	     close(unit = 50+myid)
	  endif

	   else

	     open(50+myid, file='outputp_time.'//ID, form='unformatted',
     >          status='unknown')
	     write(50+myid) time
	     close(unit = 50+myid)

	     open(50+myid, file='outputpVal_kount.'//ID, 
     >            form='unformatted', status='unknown')
	     write(50+myid) DBLE(kount)
	     close(unit = 50+myid)

	    open(50+myid, file='outputpVal_kineticdepth.'//ID,
     >           form='unformatted', status='unknown')
	     write(50+myid) kineticDepth
	     close(unit = 50+myid)
	     
	     open(50+myid, file='outputpVal_drive.'//ID, 
     >            form='unformatted', status='unknown')
	     write(50+myid) drive
	     close(unit = 50+myid)

	   if(ised .eq. 1) then
	     open(50+myid, file='outputpVal_sedTotal1.'//ID,
     >            form='unformatted', status='unknown')
	     write(50+myid) sedTotal1
	     close(unit = 50+myid)

	     open(50+myid, file='outputpVal_sedTotal2.'//ID,
     >            form='unformatted', status='unknown')
	     write(50+myid) sedTotal2
	     close(unit = 50+myid)
	  endif

	   if(waves .eq. 1) then
	     open(50+myid, file='outputpVal_wavephase.'//ID,
     >            form='unformatted', status='unknown')
	     write(50+myid) phase
	     close(unit = 50+myid)
	  endif

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

	     if(ieddy .eq. 1) then
	     open(50+myid, file='outputp_rruMean.'//ID, form='unformatted',
     >          status='old',position='append')
	     write(50+myid) rruMean
	     close(unit = 50+myid)

	     open(50+myid, file='outputp_vstMean.'//ID, form='unformatted',
     >          status='old',position='append')
	     write(50+myid) vstMean
	     close(unit = 50+myid)
	     endif

	     if(ised .eq. 1) then
	     open(50+myid, file='outputp_vCsed.'//ID, form='unformatted',
     >          status='old',position='append')
	     write(50+myid) vCsed
	     close(unit = 50+myid)

	     open(50+myid, file='outputp_cSed.'//ID, form='unformatted',
     >          status='old',position='append')
	     write(50+myid) sedMean
	     close(unit = 50+myid)
	     endif

	     open(50+myid, file='outputp_BruntN.'//ID, form='unformatted',
     >          status='old',position='append')
	     write(50+myid) BruntNSqr
	     close(unit = 50+myid)
	     
	     open(50+myid, file='outputp_rhoMean.'//ID, form='unformatted',
     >          status='old',position='append')
	     write(50+myid) rhoMean
	     close(unit = 50+myid)

	    
	 open(50+myid, file='outputp_rhoPrimeMean.'//ID, form='unformatted',
     >          status='old',position='append')
	     write(50+myid) rhoPrimeMean
	     close(unit = 50+myid)
 
	     open(50+myid, file='outputp_rhoSqrMean.'//ID, form='unformatted',
     >          status='old',position='append')
	     write(50+myid) rhoSqrMean
	     close(unit = 50+myid)

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

	     if(ieddy .eq. 1) then
	     open(50+myid, file='outputp_rruMean.'//ID, form='unformatted',
     >          status='unknown')
	     write(50+myid) rruMean
	     close(unit = 50+myid)

	     open(50+myid, file='outputp_vstMean.'//ID, form='unformatted',
     >          status='unknown')
	     write(50+myid) vstMean
	     close(unit = 50+myid)
	     endif

	     if(ised .eq. 1) then
	     open(50+myid, file='outputp_vCsed.'//ID, form='unformatted',
     >          status='unknown')
	     write(50+myid) vCsed
	     close(unit = 50+myid)

	     open(50+myid, file='outputp_cSed.'//ID, form='unformatted',
     >          status='unknown')
	     write(50+myid) sedMean
	     close(unit = 50+myid)
	     endif

	     open(50+myid, file='outputp_BruntN.'//ID, form='unformatted',
     >          status='unknown')
	     write(50+myid) BruntNSqr
	     close(unit = 50+myid)

	     open(50+myid, file='outputp_rhoMean.'//ID, form='unformatted',
     >          status='unknown')
	     write(50+myid) rhoMean
	     close(unit = 50+myid)

	 open(50+myid, file='outputp_rhoPrimeMean.'//ID, form='unformatted',
     >          status='unknown')
	     write(50+myid) rhoPrimeMean
	     close(unit = 50+myid)

	     open(50+myid, file='outputp_rhoSqrMean.'//ID, form='unformatted',
     >          status='unknown')
	     write(50+myid) rhoSqrMean
	     close(unit = 50+myid)
	   endif
	endif

	return
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine output_profilesTwo(PI)
	include "size.inc"
	include "ns.inc"
        include "mpif.h"
        include "mpi.inc"
	include "padjust.inc"
	include "para.inc"
	include "eddy.inc"
	include "sed.inc"

	double precision, intent(in) :: PI
	character*4 :: ID
	double precision, dimension(1:nnj) ::
     <     uTurb, vTurb, wTurb,uvRey, uwRey, vwRey,
     <     sedMean, vCsed,BruntNsqr, rhoMean, rhoSqrMean, rhoPrimeMean 
	double precision drive, phase 
	double precision, dimension(0:nnj+1) :: uMean, vMean, wMean

	if(iTKE .ne.1) then
	call getMeanVelPro(u(:,:,:,1), uMean)
	call getMeanVelPro(u(:,:,:,2), vMean)
	call getMeanVelPro(u(:,:,:,3), wMean)
	call MPI_Barrier(MPI_COMM_WORLD, ierr)
	call get_pertVel(u, uMean, velPrimes)
	endif

	call get_turbIntensity(uTurb, vTurb, wTurb)
	call get_reynoldsStress(uvRey, uwRey, vwRey)
	
	if (ised .eq. 1) then
	call horizontalAverage(Csed, sedMean, 2)
	call get_sedTurbFlux(sedMean,vCsed)
	endif
	
	call get_BruntN(rho,BruntNSqr,rhoMean,rhoSqrMean,rhoPrimeMean)


	if (waves .eq. 1) then
	     phase = cos(2*PI*time/Twave)
	endif

	drive = steadyPall(1) 
	if (100+vert_id .gt. 999) then
	write(ID, fmt='(I4)') 100+vert_id
	else
	write(ID, fmt='(I3)') 100+vert_id
	endif


	if (myid .eq. 0) then
	   if (istep .gt.2) then !appends to existing file - use this if you want continue run to create new files
	     open(50+myid, file='outputp_time2.'//ID, form='unformatted',
     >          status='old',position='append')
	     write(50+myid) time
	     close(unit = 50+myid)

	     open(50+myid, file='outputpVal_kount2.'//ID,
     >   form='unformatted', status='old',position='append')
	     write(50+myid) DBLE(kount)
	     close(unit = 50+myid)
	     
	     open(50+myid, file='outputpVal_drive2.'//ID,
     >   form='unformatted', status='old',position='append')
	     write(50+myid) drive
	     close(unit = 50+myid)

	   if(waves .eq. 1) then
	     open(50+myid, file='outputpVal_wavephase2.'//ID,
     >   form='unformatted', status='old',position='append')
	     write(50+myid) phase
	     close(unit = 50+myid)
	  endif

	   else

	     open(50+myid, file='outputp_time2.'//ID, form='unformatted',
     >          status='unknown')
	     write(50+myid) time
	     close(unit = 50+myid)

	     open(50+myid, file='outputpVal_kount2.'//ID, 
     >            form='unformatted', status='unknown')
	     write(50+myid) DBLE(kount)
	     close(unit = 50+myid)
	     
	     open(50+myid, file='outputpVal_drive2.'//ID, 
     >            form='unformatted', status='unknown')
	     write(50+myid) drive
	     close(unit = 50+myid)

	   if(waves .eq. 1) then
	     open(50+myid, file='outputpVal_wavephase2.'//ID,
     >            form='unformatted', status='unknown')
	     write(50+myid) phase
	     close(unit = 50+myid)
	  endif
	   endif
	endif

	
	if (hor_id .eq. 0) then
	   if (istep .gt.2) then !appends to existing file - use this if you want continue run to create new files
	     open(50+myid, file='outputp_umean2.'//ID, form='unformatted',
     >          status='old',position='append')
	     write(50+myid) uMean(1:nnj)
	     close(unit = 50+myid)

	     open(50+myid, file='outputp_vmean2.'//ID, form='unformatted',
     >          status='old',position='append')
	     write(50+myid) vMean(1:nnj)
	     close(unit = 50+myid)

	     open(50+myid, file='outputp_wmean2.'//ID, form='unformatted',
     >          status='old',position='append')
	     write(50+myid) wMean(1:nnj)
	     close(unit = 50+myid)

	     open(50+myid, file='outputp_uTurb2.'//ID, form='unformatted',
     >          status='old',position='append')
	     write(50+myid) uTurb
	     close(unit = 50+myid)

	     open(50+myid, file='outputp_vTurb2.'//ID, form='unformatted',
     >          status='old',position='append')
	     write(50+myid) vTurb
	     close(unit = 50+myid)

	     open(50+myid, file='outputp_wTurb2.'//ID, form='unformatted',
     >          status='old',position='append')
	     write(50+myid) wTurb
	     close(unit = 50+myid)

	     open(50+myid, file='outputp_uvRey2.'//ID, form='unformatted',
     >          status='old',position='append')
	     write(50+myid) uvRey
	     close(unit = 50+myid)

	     open(50+myid, file='outputp_uwRey2.'//ID, form='unformatted',
     >          status='old',position='append')
	     write(50+myid) uwRey
	     close(unit = 50+myid)

	     open(50+myid, file='outputp_vwRey2.'//ID, form='unformatted',
     >          status='old',position='append')
	     write(50+myid) vwRey
	     close(unit = 50+myid)

	     if(ised .eq. 1) then
	     open(50+myid, file='outputp_vCsed2.'//ID, form='unformatted',
     >          status='old',position='append')
	     write(50+myid) vCsed
	     close(unit = 50+myid)

	     open(50+myid, file='outputp_cSed2.'//ID, form='unformatted',
     >          status='old',position='append')
	     write(50+myid) sedMean
	     close(unit = 50+myid)
	     endif

	     open(50+myid, file='outputp_rhoMean2.'//ID, form='unformatted',
     >          status='old',position='append')
	     write(50+myid) rhoMean
	     close(unit = 50+myid)
	    
	 open(50+myid, file='outputp_rhoPrimeMean2.'//ID, form='unformatted',
     >          status='old',position='append')
	     write(50+myid) rhoPrimeMean
	     close(unit = 50+myid)
 
	   else
	     open(50+myid, file='outputp_umean2.'//ID, form='unformatted',
     >          status='unknown')
	     write(50+myid) uMean(1:nnj)
	     close(unit = 50+myid)

	     open(50+myid, file='outputp_vmean2.'//ID, form='unformatted',
     >          status='unknown')
	     write(50+myid) vMean(1:nnj)
	     close(unit = 50+myid)

	     open(50+myid, file='outputp_wmean2.'//ID, form='unformatted',
     >          status='unknown')
	     write(50+myid) wMean(1:nnj)
	     close(unit = 50+myid)

	     open(50+myid, file='outputp_uTurb2.'//ID, form='unformatted',
     >          status='unknown')
	     write(50+myid) uTurb
	     close(unit = 50+myid)

	     open(50+myid, file='outputp_vTurb2.'//ID, form='unformatted',
     >          status='unknown')
	     write(50+myid) vTurb
	     close(unit = 50+myid)

	     open(50+myid, file='outputp_wTurb2.'//ID, form='unformatted',
     >          status='unknown')
	     write(50+myid) wTurb
	     close(unit = 50+myid)

	     open(50+myid, file='outputp_uvRey2.'//ID, form='unformatted',
     >          status='unknown')
	     write(50+myid) uvRey
	     close(unit = 50+myid)

	     open(50+myid, file='outputp_uwRey2.'//ID, form='unformatted',
     >          status='unknown')
	     write(50+myid) uwRey
	     close(unit = 50+myid)

	     open(50+myid, file='outputp_vwRey2.'//ID, form='unformatted',
     >          status='unknown')
	     write(50+myid) vwRey
	     close(unit = 50+myid)

	     if(ised .eq. 1) then
	     open(50+myid, file='outputp_vCsed2.'//ID, form='unformatted',
     >          status='unknown')
	     write(50+myid) vCsed
	     close(unit = 50+myid)

	     open(50+myid, file='outputp_cSed2.'//ID, form='unformatted',
     >          status='unknown')
	     write(50+myid) sedMean
	     close(unit = 50+myid)
	     endif

	     open(50+myid, file='outputp_rhoMean2.'//ID, form='unformatted',
     >          status='unknown')
	     write(50+myid) rhoMean
	     close(unit = 50+myid)

	 open(50+myid, file='outputp_rhoPrimeMean2.'//ID, form='unformatted',
     >          status='unknown')
	     write(50+myid) rhoPrimeMean
	     close(unit = 50+myid)
	   endif
	endif

	return
	end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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


C	if ( MYID .eq. 0 ) then
C
C	   if (istep .gt.1) then !appends to existing file - use this if you want continue run to create new files
C	     open(50+myid, file='outputpval_cfl.100', form='unformatted',
C     >          status='old',position='append')
C	     write(50+myid) cflmax
C	     close(unit = 50+myid)
C	   else
C	     open(50+myid, file='outputpval_cfl.100', form='unformatted',
C     >          status='unknown')
C	     write(50+myid) cflmax
C	     close(unit = 50+myid)
C	   endif
C	endif


	if ( cflmax .gt. 0.9 .and. istep .gt. 5) then
	   if ( MYID .eq. 0 )
     <	      write(*,*) ' CFL = ', cflmax, ' > 0.9, Stop! '
	   stop
	else
	   if ( MYID .eq. 0 .and. mod(istep,200) .eq. 0) then
	      print *, ''
	      print *, '************************************************' 
     	      write(*,1) cflmax, dtime, time
	      print *, '************************************************'
	      print *, ''
	   end if
	endif
    1	format(' CFL = ', e15.9, '  dtime = ', e9.3, '  time = ', e9.3)

C	if ( MYID .eq. 0) then
C	   call idate(today)	! today(1)=day, (2)=month, (3)=year
C	   call itime(now)	! now(1)=hour, (2)=minute, (3)=second
C	   write ( *, 1000 )  today(2), today(1), today(3), now
C	end if
C 1000	format ( 'Date ', i2.2, '/', i2.2, '/', i4.4, '; time ',
C     &	               i2.2, ':', i2.2, ':', i2.2 )
C
C	   if (MYID .eq. 0) then
C	      if(istep.eq.1) then
C		 open(123, file='qoutput',status='replace')
C	      else
C		 open(123, file='qoutput',position='append')
C	      endif

C	      write(123,105) istep, nstep, time, cflmax	   
C	      call idate(today)	! today(1)=day, (2)=month, (3)=year
C	      call itime(now)	! now(1)=hour, (2)=minute, (3)=second
C 105          format('Time step = ',i7,'/',i7', time = ',e10.3, 
C     &                ', cfl = ',e10.3)
C	      close(unit=123)
C	   endif

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

