	double precision drhs(-1:nni+2, -1:nnk+2),
     &	                 af_x(-1:nni+2, -1:nnk+2),
     &	                 af_z(-1:nni+2, -1:nnk+2)
		 
        common/dphy/ drhs, af_x, af_z

        double precision jad(-1:nni+2, -1:nnk+2),
     &	                 d11(-1:nni+2, -1:nnk+2),
     &                   d13(-1:nni+2, -1:nnk+2),
     &                   d31(-1:nni+2, -1:nnk+2),
     &                   d33(-1:nni+2, -1:nnk+2),
     &		         xixd(-1:nni+2, -1:nnk+2),
     &                   ztxd(-1:nni+2, -1:nnk+2),
     &			 xizd(-1:nni+2, -1:nnk+2),
     &                   ztzd(-1:nni+2, -1:nnk+2),
     &                   h(-1:nni+2, -1:nnk+2)		 
        common/dmet/ jad,d11,d13,d31,d33,xixd,ztxd,xizd,ztzd,h
	
