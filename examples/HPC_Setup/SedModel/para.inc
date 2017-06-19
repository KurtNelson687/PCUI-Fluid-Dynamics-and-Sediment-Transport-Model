	integer runCase, newrun, periodic, iscalar, ieddy, mg_level
	integer ised, iTKE, newWave, irho
	integer waves, pAdjust, nsavePro
	integer nstep, nsave, istep, kount, nunit, grid_only
	integer maxstep, iterchk(5), maxiter(5)
	common/ctl/
     <		runCase, newrun, periodic, iscalar, ieddy, mg_level,
     <		istep, nstep, nsave, kount, nunit, nsavePro,
     <		maxstep, iterchk, maxiter, grid_only,
     <		waves, ised, pAdjust, iTKE, newWave, irho

	double precision time, dtime, vis, ak, g
	double precision rhoWater
	double precision dpdxSteady, dpdxWave, Twave
	double precision omg_cyl, omg_lid, omg2
	double precision factor, tol(5), ter(5), slowiter(5)
	common/par/
     <		time, dtime, vis, ak, g,
     <          rhoWater,
     <		dpdxSteady,dpdxWave, Twave,
     <		omg_cyl, omg_lid, omg2,                     
     <		factor, tol, ter, slowiter

	double precision phi1, phi2, yphi, aphi
	common/sci/
     <		phi1, phi2, yphi, aphi
