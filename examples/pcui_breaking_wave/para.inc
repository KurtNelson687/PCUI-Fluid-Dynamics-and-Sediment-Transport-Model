	integer case, newrun, periodic, iscalar, ieddy, mg_level
	integer nstep, nsave, istep, kount, nunit, grid_only
	integer maxstep, iterchk(6), maxiter(6)
	common/ctl/
     <		case, newrun, periodic, iscalar, ieddy, mg_level,
     <		istep, nstep, nsave, kount, nunit,
     <		maxstep, iterchk, maxiter, grid_only

	double precision time, dtime, vis, ak, g
	double precision omg_cyl, omg_lid, omg2
	double precision factor, tol(6), ter(6), slowiter(6)
	common/par/
     <		time, dtime, vis, ak, g,
     <		omg_cyl, omg_lid, omg2,                     
     <		factor, tol, ter, slowiter

	double precision phi1, yphi, aphi
	common/sci/
     <		phi1, yphi, aphi
