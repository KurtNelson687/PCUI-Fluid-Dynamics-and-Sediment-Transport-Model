	double precision, dimension(0:nni+1,0:nnj+1,0:nnk+1) ::
     <		vst, akst, sab
	double precision, dimension(0:nni+1,0:nnj+1,0:nnk+1,6) :: 
     <		rr
         double precision, dimension(-1:nni+2,-1:nnj+2,-1:nnk+2) ::
     <    vst_o, akst_o, diss_sgs
	common/les/
     <		vst, akst, sab, rr, vst_o, akst_o, diss_sgs
