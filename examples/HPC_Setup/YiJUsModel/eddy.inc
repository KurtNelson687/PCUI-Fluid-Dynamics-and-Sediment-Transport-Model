
	double precision vst(0:nni+1,0:nnj+1,0:nnk+1), 
     &                   akst(0:nni+1,0:nnj+1,0:nnk+1), 
     &                  aksd(0:nni+1,0:nnj+1,0:nnk+1)
	double precision akeq(0:nni+1,0:nnj+1,0:nnk+1)
        double precision sab(0:nni+1,0:nnj+1,0:nnk+1)
	double precision rr(0:nni+1,0:nnj+1,0:nnk+1,9)
        double precision btau(0:nni+1,0:nnj+1,0:nnk+1,3)
	double precision bchi(0:nni+1,0:nnj+1,0:nnk+1)
        double precision chi(0:nni+1,0:nnj+1,0:nnk+1,5)
        double precision sut(0:nni+1,0:nnj+1,0:nnk+1,3)
	integer cut
	common/les/ vst, akst, aksd, sab, rr, btau, 
     &               bchi, chi, sut, akeq, cut
