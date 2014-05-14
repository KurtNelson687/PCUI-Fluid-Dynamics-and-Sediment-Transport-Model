	double precision, dimension(-1:nni+2,-1:nnj+2,-1:nnk+2,1:3) ::
     <		u
	double precision, dimension(0:nni+1,0:nnj+1,0:nnk+1,1:3) ::
     <		hb
	double precision, dimension(0:nni+1,0:nnj+1,0:nnk+1) ::
     <		uxi, uej, uzk, p
	common/uvw/
     <		u, hb, uxi, uej, uzk, p

	double precision, dimension(-1:nni+2,-1:nnj+2,-1:nnk+2) ::
     <  	phi, phi2, phi3, phi4
	double precision, dimension(0:nni+1,0:nnj+1,0:nnk+1) ::
     <	        hbs, sus, hbs2, sus2, hbs3, sus3, hbs4, sus4
	common/scl/
     <		phi, phi2, phi3, phi4, 
     <          hbs, sus, hbs2, sus2, hbs3, sus3, hbs4, sus4
