!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!C...... iPS,jPS,kPS determine the number of particles in a direction
!C......	based on ni, nj or nk. The total number of particles is:
!C...... nPart = floor((ni-1)/iPS) *
!C......         floor((nj-1)/jPS) *
!C......         floor((nk-1)/kPS)

   	double precision, parameter :: lbz = 0.D0
   	double precision, parameter :: ubz = 1.D0
   	double precision, parameter :: lbx = 0.D0
   	double precision, parameter :: ubx = 1.D0
   	double precision, parameter :: lby = 0.D0
   	double precision, parameter :: uby = 1.D0

c	integer, parameter :: lbi = (1 + (ni-1+2)*lbx)
c	integer, parameter :: ubi = (1 + (ni-1+2)*ubx)
c	integer, parameter :: lbj = (1 + (nj-1+2)*lby)
c	integer, parameter :: ubj = (1 + (nj-1+2)*uby)
c	integer, parameter :: lbk = (1 + (nk-1+2)*lbz)
c	integer, parameter :: ubk = (1 + (nk-1+2)*ubz)
 	integer, parameter :: lbi = 10
 	integer, parameter :: ubi = 50
       integer, parameter :: iin = 5
 	integer, parameter :: lbj = 2
 	integer, parameter :: ubj = nj+1
 	integer, parameter :: lbk = 8
 	integer, parameter :: ubk = 8

	integer, parameter :: nPart = ((ubi-lbi)/iin+1)*(ubj-lbj+1)*(ubk-lbk+1)
	integer            :: iPart, ipstep

	double precision,dimension(1:ni+2,1:nj+2,1:nk+2,3) :: uu,xxp
	double precision,dimension(1:ni,1:nj,1:nk,3) :: uuu,xxxp

	double precision,dimension(nPart, 3) :: xPart,uPart,k1,k2,k3,k4,
     <                                    xPartC
	double precision,dimension(3)        :: lb,rb,xxL,xxR
	integer,dimension(nPart,3)           :: xPartS,xPartB,xPartBT

	common/particles/ 
     <			  xPart, uPart, xPartS, uu, xxp, iPart, ipstep, 
     <                    lb,rb,k1,k2,k3,k4,xxL,xxR,xxxp,uuu,
     <                    xPartB,xPartBT,xPartC
