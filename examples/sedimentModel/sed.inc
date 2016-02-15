	double precision ws, rhoSed
	double precision, dimension(-1:nni+2,-1:nnj+2,-1:nnk+2) ::
     <          Csed
	double precision, dimension(0:nni+1,0:nnj+1,0:nnk+1) ::
     <          hbCsed, suCsed
        common/sed/
     <          Csed, hbCsed, suCsed, ws, rhoSed
