cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	double precision jac(-1:nni+2, -1:nnj+2, -1:nnk+2),
     &                   xix(-1:nni+2, -1:nnj+2, -1:nnk+2),
     &                   xiy(-1:nni+2, -1:nnj+2, -1:nnk+2),
     &                   xiz(-1:nni+2, -1:nnj+2, -1:nnk+2),
     &                   etx(-1:nni+2, -1:nnj+2, -1:nnk+2),
     &                   ety(-1:nni+2, -1:nnj+2, -1:nnk+2),
     &                   etz(-1:nni+2, -1:nnj+2, -1:nnk+2),
     &                   ztx(-1:nni+2, -1:nnj+2, -1:nnk+2),
     &                   zty(-1:nni+2, -1:nnj+2, -1:nnk+2),
     &                   ztz(-1:nni+2, -1:nnj+2, -1:nnk+2)
	common/met/ jac, xix, xiy, xiz, etx, ety, etz, ztx, zty, ztz

	double precision jac_old(-1:nni+2, -1:nnj+2, -1:nnk+2),
     &                   xix_old(-1:nni+2, -1:nnj+2, -1:nnk+2),
     &                   xiy_old(-1:nni+2, -1:nnj+2, -1:nnk+2),
     &                   xiz_old(-1:nni+2, -1:nnj+2, -1:nnk+2),
     &                   etx_old(-1:nni+2, -1:nnj+2, -1:nnk+2),
     &                   ety_old(-1:nni+2, -1:nnj+2, -1:nnk+2),
     &                   etz_old(-1:nni+2, -1:nnj+2, -1:nnk+2),
     &                   ztx_old(-1:nni+2, -1:nnj+2, -1:nnk+2),
     &                   zty_old(-1:nni+2, -1:nnj+2, -1:nnk+2),
     &                   ztz_old(-1:nni+2, -1:nnj+2, -1:nnk+2)
	common/old/ jac_old, xix_old, xiy_old, xiz_old, etx_old, 
     &              ety_old, etz_old, ztx_old, zty_old, ztz_old
	
	double precision jac_diff(-1:nni+2, -1:nnj+2, -1:nnk+2)
	common/diff/ jac_diff  

	double precision g11(-1:nni+2, -1:nnj+2, -1:nnk+2),
     &                   g12(-1:nni+2, -1:nnj+2, -1:nnk+2),
     &                   g13(-1:nni+2, -1:nnj+2, -1:nnk+2),
     &                   g21(-1:nni+2, -1:nnj+2, -1:nnk+2),
     &                   g22(-1:nni+2, -1:nnj+2, -1:nnk+2),
     &                   g23(-1:nni+2, -1:nnj+2, -1:nnk+2),
     &                   g31(-1:nni+2, -1:nnj+2, -1:nnk+2),
     &                   g32(-1:nni+2, -1:nnj+2, -1:nnk+2),
     &                   g33(-1:nni+2, -1:nnj+2, -1:nnk+2),
     &                   gcc(-1:nni+2, -1:nnj+2, -1:nnk+2)
	common/ggg/ g11, g12, g13, g21, g22, g23, g31, g32, g33, gcc

	double precision jaf(-1:nni1+2, -1:nnj1+2, -1:nnk1+2),
     &                   f11(-1:nni1+2, -1:nnj1+2, -1:nnk1+2),
     &                   f12(-1:nni1+2, -1:nnj1+2, -1:nnk1+2),
     &                   f13(-1:nni1+2, -1:nnj1+2, -1:nnk1+2),
     &                   f21(-1:nni1+2, -1:nnj1+2, -1:nnk1+2),
     &                   f22(-1:nni1+2, -1:nnj1+2, -1:nnk1+2),
     &                   f23(-1:nni1+2, -1:nnj1+2, -1:nnk1+2),
     &                   f31(-1:nni1+2, -1:nnj1+2, -1:nnk1+2),
     &                   f32(-1:nni1+2, -1:nnj1+2, -1:nnk1+2),
     &                   f33(-1:nni1+2, -1:nnj1+2, -1:nnk1+2),
     &                   fcc(-1:nni1+2, -1:nnj1+2, -1:nnk1+2)
	common/fff/ jaf, f11, f12, f13, f21, f22, f23, f31, f32, f33, fcc

	double precision jah(-1:nni2+2, -1:nnj2+2, -1:nnk2+2),
     &                   h11(-1:nni2+2, -1:nnj2+2, -1:nnk2+2),
     &                   h12(-1:nni2+2, -1:nnj2+2, -1:nnk2+2),
     &                   h13(-1:nni2+2, -1:nnj2+2, -1:nnk2+2),
     &                   h21(-1:nni2+2, -1:nnj2+2, -1:nnk2+2),
     &                   h22(-1:nni2+2, -1:nnj2+2, -1:nnk2+2),
     &                   h23(-1:nni2+2, -1:nnj2+2, -1:nnk2+2),
     &                   h31(-1:nni2+2, -1:nnj2+2, -1:nnk2+2),
     &                   h32(-1:nni2+2, -1:nnj2+2, -1:nnk2+2),
     &                   h33(-1:nni2+2, -1:nnj2+2, -1:nnk2+2),
     &                   hcc(-1:nni2+2, -1:nnj2+2, -1:nnk2+2)
	common/hhh/ jah, h11, h12, h13, h21, h22, h23, h31, h32, h33, hcc

	double precision jap(-1:nni3+2, -1:nnj3+2, -1:nnk3+2),
     &                   p11(-1:nni3+2, -1:nnj3+2, -1:nnk3+2),
     &                   p12(-1:nni3+2, -1:nnj3+2, -1:nnk3+2),
     &                   p13(-1:nni3+2, -1:nnj3+2, -1:nnk3+2),
     &                   p21(-1:nni3+2, -1:nnj3+2, -1:nnk3+2),
     &                   p22(-1:nni3+2, -1:nnj3+2, -1:nnk3+2),
     &                   p23(-1:nni3+2, -1:nnj3+2, -1:nnk3+2),
     &                   p31(-1:nni3+2, -1:nnj3+2, -1:nnk3+2),
     &                   p32(-1:nni3+2, -1:nnj3+2, -1:nnk3+2),
     &                   p33(-1:nni3+2, -1:nnj3+2, -1:nnk3+2),
     &                   pcc(-1:nni3+2, -1:nnj3+2, -1:nnk3+2)
	common/ppp/ jap, p11, p12, p13, p21, p22, p23, p31, p32, p33, pcc

	double precision jar(-1:nni4+2, -1:nnj4+2, -1:nnk4+2),
     &                   r11(-1:nni4+2, -1:nnj4+2, -1:nnk4+2),
     &                   r12(-1:nni4+2, -1:nnj4+2, -1:nnk4+2),
     &                   r13(-1:nni4+2, -1:nnj4+2, -1:nnk4+2),
     &                   r21(-1:nni4+2, -1:nnj4+2, -1:nnk4+2),
     &                   r22(-1:nni4+2, -1:nnj4+2, -1:nnk4+2),
     &                   r23(-1:nni4+2, -1:nnj4+2, -1:nnk4+2),
     &                   r31(-1:nni4+2, -1:nnj4+2, -1:nnk4+2),
     &                   r32(-1:nni4+2, -1:nnj4+2, -1:nnk4+2),
     &                   r33(-1:nni4+2, -1:nnj4+2, -1:nnk4+2),
     &                   rcc(-1:nni4+2, -1:nnj4+2, -1:nnk4+2)
	common/rrr/ jar, r11, r12, r13, r21, r22, r23, r31, r32, r33, rcc

	double precision jae(-1:nni5+2, -1:nnj5+2, -1:nnk5+2),
     &                   e11(-1:nni5+2, -1:nnj5+2, -1:nnk5+2),
     &                   e12(-1:nni5+2, -1:nnj5+2, -1:nnk5+2),
     &                   e13(-1:nni5+2, -1:nnj5+2, -1:nnk5+2),
     &                   e21(-1:nni5+2, -1:nnj5+2, -1:nnk5+2),
     &                   e22(-1:nni5+2, -1:nnj5+2, -1:nnk5+2),
     &                   e23(-1:nni5+2, -1:nnj5+2, -1:nnk5+2),
     &                   e31(-1:nni5+2, -1:nnj5+2, -1:nnk5+2),
     &                   e32(-1:nni5+2, -1:nnj5+2, -1:nnk5+2),
     &                   e33(-1:nni5+2, -1:nnj5+2, -1:nnk5+2),
     &                   ecc(-1:nni5+2, -1:nnj5+2, -1:nnk5+2)
	common/eee/ jae, e11, e12, e13, e21, e22, e23, e31, e32, e33, ecc
	

	double precision xp(-1:nni+2, -1:nnj+2, -1:nnk+2, 3),
     &                   xp_old(-1:nni+2, -1:nnj+2, -1:nnk+2, 3)
	double precision phi_init(-1:nni+2, -1:nnj+2, -1:nnk+2)
	double precision w_lid(-1:nni+2, -1:nnk+2),
     &                   u_lid(-1:nni+2, -1:nnk+2),
     &			 v_lid(-1:nni+2, -1:nnk+2),
     & 	                 u_flux(-1:nnj+2, -1:nnk+2)
	common/xyz/ xp, xp_old, phi_init, w_lid, u_lid, u_flux, v_lid

	double precision depth(-1:nni+2, -1:nnk+2),
     &                   depth0(-1:nni+2, -1:nnk+2), 
     &			 cos_x(-1:nni+2, -1:nnk+2),
     & 			 sin_x(-1:nni+2, -1:nnk+2),	
     &			 cos_z(-1:nni+2, -1:nnk+2),
     &			 sin_z(-1:nni+2, -1:nnk+2),
     &			 sin_xz(-1:nni+2, -1:nnk+2),
     &			 cos_xz(-1:nni+2, -1:nnk+2)
        common/dep/ depth,cos_x,sin_x,cos_z,sin_z,sin_xz,cos_xz,depth0