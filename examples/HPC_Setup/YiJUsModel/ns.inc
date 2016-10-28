	double precision u(-1:nni+2,-1:nnj+2,-1:nnk+2,1:3)
c	double precision ug(-1:nni+2,-1:nnj+2,-1:nnk+2,1:3)
	double precision uavg(-1:nni+2,-1:nnj+2,-1:nnk+2,1:3)
	double precision hb(0:nni+1,0:nnj+1,0:nnk+1,1:3)
	double precision uxi(0:nni+1,0:nnj+1,0:nnk+1), 
     & 	                 uej(0:nni+1,0:nnj+1,0:nnk+1), 
     &                   uzk(0:nni+1,0:nnj+1,0:nnk+1),
     &                   kxi(0:nni+1,0:nnj+1,0:nnk+1), 
     & 	                 kej(0:nni+1,0:nnj+1,0:nnk+1), 
     &                   kzk_old(0:nni+1,0:nnj+1,0:nnk+1),
     &                   kxi_old(0:nni+1,0:nnj+1,0:nnk+1), 
     & 	                 kej_old(0:nni+1,0:nnj+1,0:nnk+1), 
     &                   kzk(0:nni+1,0:nnj+1,0:nnk+1), 
     &                     p(0:nni+1,0:nnj+1,0:nnk+1),
     &                 drive(0:nni+1,0:nnj+1,0:nnk+1),
     &			 wxi(0:nni+1,0:nnj+1,0:nnk+1),  
     &			 wej(0:nni+1,0:nnj+1,0:nnk+1),  
     &			 wzk(0:nni+1,0:nnj+1,0:nnk+1),
     &                 u_bed(-1:nni+2,-1:nnk+2, 1:3),    		  
     &                 upf(-1:nni+2,-1:nnj+2,-1:nnk+2, 1:3),
     &                 sai(-1:nni+2,-1:nnj+2,-1:nnk+2)
    		  
	common/uvw/ u, uavg, hb, uxi, uej, uzk, p, wxi, wej, wzk, 
     &              kxi, kej, kzk, kxi_old, kej_old, kzk_old,
     &              u_bed, upf, sai
	double precision phi(-1:nni+2,-1:nnj+2,-1:nnk+2)
	double precision hbs(0:nni+1,0:nnj+1,0:nnk+1), 
     & 	                 sus(0:nni+1,0:nnj+1,0:nnk+1)
	common/scl/ phi, hbs, sus

 	double precision Pgrd
        common/grd/ Pgrd
         
        double precision taus(1:nni, 1:nnj,1:nnk,1:4)
        common/str/ taus