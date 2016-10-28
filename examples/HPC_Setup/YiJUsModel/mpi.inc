	integer comm3d, ierr, myid, stat(MPI_STATUS_SIZE)
	integer n_west, n_east, n_suth, n_nrth, n_back, n_frnt
	integer npx, npy, npz
	integer ius, iue, jus, jue, kus, kue
	logical periods(3)
	common /mpi/ 
     <		comm3d, ierr, myid,
     <		n_west, n_east, n_suth, n_nrth, n_back, n_frnt,
     <		npx, npy, npz,
     <	        ius, iue, jus, jue, kus, kue,
     <		periods
