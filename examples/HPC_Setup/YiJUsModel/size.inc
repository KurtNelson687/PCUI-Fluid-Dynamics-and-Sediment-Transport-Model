C...... ni/px must be divisible by 16
C...... nj/py must be divisible by 16
C...... nk/pz must be divisible by 16

        implicit none

	integer px,py,pz,ni,nj,nk,ni1,nj1,nk1,ni2,nj2,nk2,
     & 	        ni3,nj3,nk3,ni4,nj4,nk4,nni,nnj,nnk,nni1,
     &          nnj1,nnk1,nni2,nnj2,nnk2,nni3,nnj3,nnk3,
     &          nni4,nnj4,nnk4,ni5,nj5,nk5,nni5,nnj5,nnk5 
	
C py is the vertical direction
	parameter(px = 8, py = 3,pz = 4,
     &            ni = 512, nj = 192, nk = 256,
     &            ni1 = ni / 2,nj1 = nj / 2,nk1 = nk / 2, 
     &            ni2 = ni / 4,nj2 = nj / 4,nk2 = nk / 4,
     &            ni3 = ni / 8,nj3 = nj / 8,nk3 = nk / 8,
     &            ni4 = ni / 16,nj4 = nj / 16,nk4 = nk / 16,
     &            ni5 = ni / 32,nj5 = nj / 32,nk5 = nk / 32,
     &            nni = ni / px,nnj = nj / py,nnk = nk / pz,
     &            nni1 = ni1 / px,nnj1 = nj1 / py,nnk1 = nk1 / pz,
     &            nni2 = ni2 / px,nnj2 = nj2 / py,nnk2 = nk2 / pz,
     &            nni3 = ni3 / px,nnj3 = nj3 / py,nnk3 = nk3 / pz,
     &            nni4 = ni4 / px,nnj4 = nj4 / py,nnk4 = nk4 / pz,
     &            nni5 = ni5 / px,nnj5 = nj5 / py,nnk5 = nk5 / pz)

	integer pcoord(px*py*pz,2)
	common/pccom/ pcoord
