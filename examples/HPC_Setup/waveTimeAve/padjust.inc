	double precision steadyPall(1:nnj), uTheo(1:nnj)
	double precision Ubulk
	double precision SnRun(0:1498), uRun(0:1499), Sn
	double precision, dimension(0:nni+1,0:nnj+1,0:nnk+1,1:3) ::
     <		velPrimes
	double precision, dimension(1:nni,1:nnj,1:nnk,1:3) ::
     <		ConNew, ConOld, DisNew, DisOld, DisINew
	common/pvars/
     <		steadyPall, uTheo, Ubulk, ConNew,
     <          ConOld, DisNew, DisOld, DisINew,
     <          velPrimes, Sn, SnRun, uRun
