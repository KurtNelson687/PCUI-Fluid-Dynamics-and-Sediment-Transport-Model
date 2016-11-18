cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	double precision, dimension(0:nni+1, 0:nnj+1, 0:nnk+1) :: 
     <		yetjface
	double precision, dimension(-1:nni+2, -1:nnj+2, -1:nnk+2) ::
     <		jac, xix, xiy, xiz, etx, ety, etz, ztx, zty, ztz
	common/met/
     <		jac, xix, xiy, xiz, etx, ety, etz, ztx, zty, ztz,
     <          yetjface
	double precision, dimension(-1:nni+2, -1:nnj+2, -1:nnk+2) ::
     <		g11, g12, g13, g21, g22, g23, g31, g32, g33, gcc
	common/ggg/
     <          g11, g12, g13, g21, g22, g23, g31, g32, g33, gcc

	double precision, dimension(-1:nni1+2, -1:nnj1+2, -1:nnk1+2) ::
     <		jaf, f11, f12, f13, f21, f22, f23, f31, f32, f33, fcc
	common/fff/
     <		jaf, f11, f12, f13, f21, f22, f23, f31, f32, f33, fcc

	double precision, dimension(-1:nni2+2, -1:nnj2+2, -1:nnk2+2) ::
     <		jah, h11, h12, h13, h21, h22, h23, h31, h32, h33, hcc
	common/hhh/
     <		jah, h11, h12, h13, h21, h22, h23, h31, h32, h33, hcc

	double precision, dimension(-1:nni3+2, -1:nnj3+2, -1:nnk3+2) ::
     <		jap, p11, p12, p13, p21, p22, p23, p31, p32, p33, pcc
	common/ppp/
     <		jap, p11, p12, p13, p21, p22, p23, p31, p32, p33, pcc

	double precision, dimension(-1:nni4+2, -1:nnj4+2, -1:nnk4+2) ::
     <		jar, r11, r12, r13, r21, r22, r23, r31, r32, r33, rcc
	common/rrr/
     <		jar, r11, r12, r13, r21, r22, r23, r31, r32, r33, rcc

	double precision, dimension(-1:nni+2, -1:nnj+2, -1:nnk+2, 3) ::
     <		xp
	double precision, dimension(-1:nni+2, -1:nnj+2, -1:nnk+2) ::
     <		rho_init
	double precision, dimension(-1:nni+2, -1:nnk+2) ::
     <		w_lid, u_lid
	common/xyz/
     <		xp, rho_init, w_lid, u_lid
