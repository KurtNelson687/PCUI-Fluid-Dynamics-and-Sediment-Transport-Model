	integer comm3d, ierr, myid
	integer vert_comm, vert_id
	integer hor_comm, hor_id, hor_nprocs
	integer n_west, n_east, n_suth, n_nrth, n_back, n_frnt
	integer npx, npy, npz
	integer ius, iue, jus, jue, kus, kue
	logical periods(3)
	common /mpi/ 
     <		comm3d, ierr, myid,
     <		n_west, n_east, n_suth, n_nrth, n_back, n_frnt,
     <		npx, npy, npz,
     <	        ius, iue, jus, jue, kus, kue,
     <		periods,
     <		hor_comm, hor_id, hor_nprocs,
     <		vert_comm, vert_id
