	double precision sedc(-1:nni+2, -1:nnj+2, -1:nnk+2)
	double precision hbd(0:nni+1, 0:nnj+1, 0:nnk+1),
     &                   sud(0:nni+1, 0:nnj+1, 0:nnk+1),
     &                   sedf(0:nni+1, 0:nnj+1, 0:nnk+1,1:3),	
     &                   pcorr(0:nni+1, 0:nnj+1, 0:nnk+1,1:3),
     &                   rho_tot(-1:nni+2, -1:nnj+2, -1:nnk+2)		
        common/sdm/ sedc, hbd, sud, sedf, pcorr, rho_tot
	
	double precision pick(-1:nni+2, -1:nnk+2),
     &                   thta(-1:nni+2, -1:nnk+2),
     &			 ushr(-1:nni+2, -1:nnk+2),
     &			 utan(-1:nni+2, -1:nnk+2),
     &			 Cdxi(-1:nni+2, -1:nnk+2),
     &                   Cdzt(-1:nni+2, -1:nnk+2),		
     &                    fct(-1:nni+2, -1:nnk+2),
     &                   sedv(-1:nni+2, -1:nnk+2),
     & 	              sd_coef(-1:nni+2, -1:nnk+2),
     &			   c1(-1:nni+2, -1:nnk+2),
     &			   c2(-1:nni+2, -1:nnk+2),
     &			   c3(-1:nni+2, -1:nnk+2)
		
        double precision schd, tau_mean
	common/mov/ pick, thta, ushr, utan, schd, Cdxi, Cdzt,
     &              fct, sedv, sd_coef, c1, c2, c3, tau_mean			
	double precision diam, spwght, taub, thta_cr, delta  
        common/par/ diam, spwght, taub, thta_cr, delta

	double precision qbi(-1:nni+2, -1:nnk+2),
     &                   qbk(-1:nni+2, -1:nnk+2), 		
     &                   qbxi(0:nni+1, 0:nnk+1),
     &                   qbzt(0:nni+1, 0:nnk+1),
     &                   ugxi(0:nni+1, 0:nnk+1),
     &                   ugzt(0:nni+1, 0:nnk+1)
		
	common/bdl/ qbi,qbk, qbxi, qbzt, ugxi, ugzt

	double precision ug(-1:nni+2, -1:nnk+2, 1:2)

	common/ug/ ug