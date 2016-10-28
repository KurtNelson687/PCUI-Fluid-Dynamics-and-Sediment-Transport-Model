        integer case, newrun, periodic, iscalar, ieddy, mg_level
	integer nstep, nsave, istep, kount, nunit, ncheck, imove
	integer maxstep, iterchk(6), maxiter(6), btmBuffer, d_diff
	integer sedi, logU, periodicZ, btmAug, movingGrid, bedload
	integer d_grav, open_top
	common/ctl/
     &		case, newrun, periodic, iscalar, ieddy, mg_level,
     &		istep, nstep, nsave, kount, nunit, ncheck, d_diff,
     &		maxstep, iterchk, maxiter, sedi, logU, bedload, 
     &          d_grav, periodicZ, btmAug, movingGrid, open_top,
     &          imove, btmBuffer

	double precision time, dtime, vis, ak, g, rho, ws
	double precision omg_cyl, omg_lid, omg2
	double precision factor, tol(6), ter(6), slowiter(6)
	double precision yb(2), datum
        double precision amp, omega, cflsuper
	common/par/
     &		time, dtime, vis, ak, g, datum,
     &		omg_cyl, omg_lid, omg2,                     
     &		factor, tol, ter, slowiter, 
     & 	        yb, amp, omega, rho, ws, cflsuper

	double precision phi1, phi2, yphi, aphi
	common/sci/
     &		phi1, phi2, yphi, aphi

	double precision pi
	parameter(pi=3.14159265358979)

        integer openWest, openEast
	common/BCs/
     &		openWest, openEast

	double precision h0, hf, rp_length, const_b, grow_b 
        double precision fwc, por, resp
	common/bload/ h0, hf, rp_length, const_b, grow_b,
     &                 fwc, por, resp