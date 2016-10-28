cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	parameter(mi=32, mj=32, mk=64)
	dimension x(mi,mj,mk), y(mi,mj,mk), z(mi,mj,mk)
	dimension u(mi,mj,mk,5)

	read(1) ni, nj, nk
	write(*,*)  'ni, nj, nk = ', ni, nj, nk
 	read(1) (((x(i,j,k),i=1,ni),j=1,nj),k=1,nk),
     *	        (((y(i,j,k),i=1,ni),j=1,nj),k=1,nk),
     *	        (((z(i,j,k),i=1,ni),j=1,nj),k=1,nk)

	read(2) ni, nj, nk
	read(2) a, b, c, d
	read(2) ((((u(i,j,k,m),i=1,ni),j=1,nj),k=1,nk),m=1,5)

	i = int(ni/2)
	k = int(nk/2)
	write(*,*) i, k

	do j = 1,nj
	   write(11,*) u(i,j,k,2), y(i,j,k)
	   write(12,*) x(i,j,k), u(i,j,k,3)
	enddo	

	stop

	end
