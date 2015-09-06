      subroutine allocate_matrices

      use mainvar_module

! these values will be the local dimensions of the files
	tt=endt
	nx=endx  
	ny=endy
	topo_size=endx
	topo_size_y=endy

      ALLOCATE(bl_jump(nx,ny),wvmk_loc_i(ny),
     - bl_hor_wall_topo(topo_size,topo_size_y),
     - level_topo(topo_size,topo_size_y),
     - bl_hor_wall(nx,ny),bl_x_wall(nx,ny),
     - bl_y_wall(nx,ny),level(nx,ny),
     - t(tt),x(topo_size),x_wvmk(tt),u_wvmk(tt),
     -  y(topo_size_y),v_iter(nx,ny),
     -  zeta(nx,ny,nt,num_levels+1),
     -  E(nx,ny,nt,num_levels),bl(nx,ny),
     -  zeta_iter(nx,ny),ho(topo_size,topo_size_y),
     -  u_iter(nx,ny),u(nx,ny,nt,num_levels),
     -  v(nx,ny,nt,num_levels),
     -  hxx(nx,ny),hyy(nx,ny),
     -  a1(nx,ny,num_levels),a2(nx,ny,num_levels),
     -  b1(nx,ny,num_levels),b2(nx,ny,num_levels),
     -  bet(num_levels),d1v(nx,ny,num_levels),
     -  d2v(nx,ny,num_levels),
     -  z_alp(nx,ny,nt,num_levels),
     -  WB_x(nx,ny),WB_y(nx,ny),
     -  UU(nx,ny,nt,num_levels),
     -  VV(nx,ny,nt,num_levels),
     -  F(nx,ny,nt,num_levels),
     -  F1(nx,ny,nt,num_levels),
     -  G(nx,ny,nt,num_levels),
     -  G1(nx,ny,nt,num_levels),
     -  dhzvdy(nx,ny,num_levels),dExdx(nx,ny,num_levels),
     -  dEydy(nx,ny,num_levels),dhzudx(nx,ny,num_levels),
     -  dvdxy(nx,ny,num_levels),dhvdxy(nx,ny,num_levels),
     -  dudxy(nx,ny,num_levels),dhudxy(nx,ny,num_levels),
     -  dE2xdx1(nx,ny),dE2ydy1(nx,ny),mass(1),
     -  hxy(nx,ny),
     -  dE2xdx2(nx,ny),dE2ydy2(nx,ny),dF2xdx1(nx,ny),
     -  dF2xdx2(nx,ny),dF2xdx3(nx,ny),
     -  dG2ydy1(nx,ny),dG2ydy2(nx,ny),dG2ydy3(nx,ny),
     -  dhdt(nx,ny,nt),h(nx,ny,nt),
     -  hx(nx,ny),hy(nx,ny),
     -  dhdxt(nx,ny,nt),dhdyt(nx,ny,nt),
     -  dhdtt(nx,ny,nt),dhdxtt(nx,ny,nt),dhdytt(nx,ny,nt),
     -  dz_alpdt(nx,ny,nt,num_levels),
     -  dudxx(nx,ny,num_levels),dhudxx(nx,ny,num_levels),
     -  dvdyy(nx,ny,num_levels),dhvdyy(nx,ny,num_levels),
     -  d3u(nx,ny,num_levels),
     -  d3v(nx,ny,num_levels),
     -  hp(nx,ny,nt),c_rad(nx,ny),hzu(nx,ny),hzv(nx,ny),hv(nx,ny),
     -  hu(nx,ny),E_x_grp(nx,ny),E_y_grp(nx,ny),
     -  E2_x_grp1(nx,ny),E2_x_grp2(nx,ny),E2_y_grp1(nx,ny),
     -  E2_y_grp2(nx,ny),
     -  dudx(nx,ny,num_levels),dvdy(nx,ny,num_levels),
     -  dhudx(nx,ny,num_levels),dhvdy(nx,ny,num_levels),
     -  uxpvy(nx,ny),huxphvy(nx,ny),
     -  F2_grp1(nx,ny),cutoff_mat(nx,ny),cutoff_mat_d(nx,ny),
     -  dudt(nx,ny,nt),dvdt(nx,ny,nt),dudxt(nx,ny),dvdyt(nx,ny),
     -  dhudxt(nx,ny),dhvdyt(nx,ny),
     -  x_tmp(topo_size),y_tmp(topo_size_y),
     -  ho_tmp(topo_size,topo_size_y),dzdx(nx,ny,num_levels),
     -  bl_x_tmp(nx,ny),dzdy(nx,ny,num_levels),
     -  bl_y_tmp(nx,ny),dudy(nx,ny,num_levels),dvdx(nx,ny,num_levels),
     -  mov_sh_lezs(nx,ny),
     -  mov_sh_gtzs(nx,ny),d1u(nx,ny,num_levels),
     -  d2u(nx,ny,num_levels),
     -  eddy_visc(nx,ny),dhzudy(nx,ny,num_levels),
     -  dhzvdx(nx,ny,num_levels),
     -  t_breaking(nx,ny),VG(nx,ny,2),
     -  net_sed_trans(nx,ny),roller(nx,ny),deltar(nx,ny),
     -  loc_sum(nx,ny),k_turb(nx,ny,nt),ep_turb(nx,ny,nt),
     -  array(nx,ny),array_2(nx,ny),
     -  stress_xx(4,nx,ny),stress_xy(nx,ny),stress_yy(nx,ny),
     -  f_int_src(nx,ny),f_int_src_v(nx,ny),max_vel(nx,ny),
     -  max_fs(nx,ny),bet2(num_levels),v1_iter(nx,ny),
     -  u1_iter(nx,ny),zeta1_iter(nx,ny),fs(tt),fa(tt),
     -  S_grp(nx,ny,num_levels),T_grp(nx,ny,num_levels),
     -  dSdx(nx,ny,num_levels),dSdy(nx,ny,num_levels),
     -  dTdx(nx,ny,num_levels),dTdy(nx,ny,num_levels),
     -  mean_z(nx,ny,4),cur_x(nx,ny,4),
     -  cur_y(nx,ny,4),mean_hgt(nx,ny,4),
     -  xi(num_levels),xi2(num_levels),
     -  NL_a1(nx,ny,num_levels),NL_a2(nx,ny,num_levels),
     -  NL_b1(nx,ny,num_levels),NL_b2(nx,ny,num_levels),S_tmp(tt),
     -  T_tmp(tt),h_tmp(endx,endy),wA(freq_file,dir_file),
     -  kA(freq_file,dir_file),
     -  time_bm(1000000),fse(1000000),eta_in(tt),
     -  time_cur(1000000),curf(1000000),current_in(tt),z_br(nx,ny),
     -  k_br(nx,ny),u0_b(nx,ny,2),vort(nx,ny),u_z(nx,ny),
     -  v_z(nx,ny),w_z(nx,ny),dhvdx(nx,ny,num_levels),
     -  dhudy(nx,ny,num_levels),tmp_mat(nx,ny),visc_rot_x(nx,ny),
     -  visc_rot_y(nx,ny),a_x(nx,ny),b_x(nx,ny),c_x(nx,ny),
     -  d_x(nx,ny),a_y(nx,ny),b_y(nx,ny),c_y(nx,ny),d_y(nx,ny),
     -  v_tmp(ny,nx),kt(nx,ny,4),nut(nx,ny,4),k_RHS(nx,ny,4),
     -  B_mat(nx,ny,4),B_RHS(nx,ny,4),psix(nx,ny,4),psiy(nx,ny,4) )
		

	    x0=x0+x_glob(overlap+1)
		
		do i=1,endx
			x(i)=x_glob(sx_X+i-1)
		enddo

		do j=1,endy
			y(j)=y_glob(sy_X+j-1)
		enddo

          do i=1,endx
			do j=1,endy
				ho(i,j)=ho_glob(sx_X+i-1,sy_X+j-1)
			enddo
		enddo

		DEALLOCATE(ho_glob)

          e_int_ind=endx
          end_x=end_x_t
          end_y=end_y_t

      return

      end

