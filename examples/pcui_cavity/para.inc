	integer case, newrun, periodic, iscalar, ieddy, mg_level
	integer nstep, nsave, istep, kount, nunit
	integer maxstep, iterchk(5), maxiter(5)
	common/ctl/
     <		case, newrun, periodic, iscalar, ieddy, mg_level,
     <		istep, nstep, nsave, kount, nunit,
     <		maxstep, iterchk, maxiter

	double precision time, dtime, vis, ak, g
	double precision omg_cyl, omg_lid, omg2
	double precision factor, tol(5), ter(5), slowiter(5)
	common/par/
     <		time, dtime, vis, ak, g,
     <		omg_cyl, omg_lid, omg2,                     
     <		factor, tol, ter, slowiter

	double precision phi1, phi2, yphi, aphi
	common/sci/
     <		phi1, phi2, yphi, aphi