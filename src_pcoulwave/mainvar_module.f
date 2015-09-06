      module mainvar_module
      implicit none
C...........................Some Variable Descriptions.......................
C      alpha = nonlinearity parameter
C      c = wave celerity
C      chan_length, chan_width = channel (domain) length (x) and width (y)
C      co = linear shallow water wave celerity
C      count1 = current # of iterations
C      dt,dy,dx = grid steps in time, y, and x
C      depth = initial water depth
C      depth_2 = constant water depth after slope
C      E = variable used in open-water zeta equations (nx by ny by nt)
C      end_t,end_y,end_x = # of time, y, and x steps in the acutal numerical domain
C      ep = convergence error - corrector iterated until all err's < ep
C      err = error at local gridpoint (max of zeta, phi or piez)
C      F = variable used in phi equations (nx by ny by nt)
C      F1 = variable used in phi equations (nx by ny by nt)
C      G = variable used in phi equations (nx by ny by nt)
C      G1 = variable used in phi equations (nx by ny by nt)      
C      h = water depth matrix (nx by my)
C      itr = max # of iterations in corrector loop (if exceeded, corrector loop will end)
C      k = wavenumber
C      L = wavelength of wave
C      maxerr = max error in current iteration
C      per = period of wave
C      slope_reg = region of constant slope
C      t = time matrix (tt by 1)
C      u = u (x velocity) matrix (nx by ny by nt)
C      u_iter = u values from previous iteration in corrector step (nx by ny by nt)
C      UU = intermediate velocity grouping = u+dudxx (the capital U in Wei & Kirby)
C      v = v (y velocity) matrix (nx by ny by nt)
C      v_iter = v values from previous iteration in corrector step (nx by ny by nt)
C      VV = intermediate velocity grouping = v+dvdyy (the capital V in Wei & Kirby)
C      wave_hgt = initial wave amplitude
C      x = x matrix (nx by 1)
C      x0 = initial location of first wave crest
C      y = y matrix (ny by 1)
C      zeta = zeta (free surface) matrix (nx by ny by nt)
C      zeta_iter = zeta values from previous iteration in corrector step (nx by ny by nt)
C      z_alp = layer where velocity is evaluated, z_alp=beta*h
C
C      nt=# time steps kept in memory (only require 4) 
C      nx, ny = max # of steps in x or y-direction (for declaration purposes
C                                    only, actual # of steps is a different variable)
C      endx,endt,endt = non-dimenional end domain values
C      a = interval at which time steps are written to file, i.e. if a=10, then free surface 
C                  values are written to file every tenth time step
C      tt = max # of time steps 
C      bc_1,bc_2 = decision variables for right and left boundary conditions
C      q,w,e,rttt = dummy variables
C      dim = decision variable for 1 or 2 dimensoins
C      nonlin_ind = decision variables for weakly or fully nonlinear equations
C      disp_prop = decision variable for nwogu's or depth-averaged dipserion properties
C      min_itr = minimum # of iterations for Corrector stage
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	integer myrank,nprocs,ierr,cur_rank,upwind_shore,step_shore,
     -pleft,pright,ptop,pbottom,comm2d,dims(2),coords(2),
     -LRstride,BTstride,xrank,yrank,xcomm,ycomm,row,col

      integer nt,nx,ny,endt,a,tt,bc_1,bc_2,itr_or,sh_mov,jj_start,
     -q,w,ee,r,ttt,dim,nonlin_ind,min_itr,disp_prop,cur_level_z,
     -shift,bc_3,bc_4,wave_type,conv,stem_wv,filt,
     -sim_opt,filtsld_int,deriv_order_ind,display,nn,i,j,n,ID,
     -slide_type,load_topo,amp_stp,y1_ind,y2_ind,e_int_ind,
     -spng_1,spng_2,spng_3,spng_4,yy,ii,m,mm,rotate,smooth,
     -ht_ind,hm_ind,hb_ind,jj,screen_output,runup,num_nodes,
     -pts_wvl,choice,slide_node,spng_inc,int_src,hor_wall_ind,
     -y_harb_ind_1,y_harb_ind_2,harbor,half_harb,mv_for,
     -mv_bc,switch_cnt,converge,check,run,bottom_fric,decomp_type,
     -wave_breaking,sed_trans,UU_separation,bw,sum3,sum4,
     -k_ep_turb,swash_filt,use_av,nn_start,n_tmp,rotationality,
     -batch,instability,instab_time,smooth_bathy,overflow,
     -itr_or_input,NaN_check,topo_size,count_min,num_freq,
     -hrs,minutes,seconds,num_theta,cc_g,g_pri_high,shift_cur,
     -g_pri_low,spec_type,n_filt,n_s,per_count,t_ind_c,nt_per,
     -num_ts,ind_ts(10000,3),cur_ts,current,vert_vel_prof,k,
     -freq_file,dir_file,pass_count,topo_size_y,h_ind,h_nodes,
     -low_bound,inc,n_mean_st,n_mean_end,overlap,total_count,match,
     -endx_glob,endy_glob,sx,ex,sy,ey,sx_X,ex_X,sy_X,ey_X,js1,je1,
     -js2,je2,js3,je3,js4,je4,js5,je5,js6,je6,js7,je7,backscatter,
     -spec_pp,numerical_scheme,overflow_catch,breaker_type,
     -offshore_cur,plus_tide,imag_leng,left_wvmk,windstress



      integer num_levels
      parameter(nt=4,pass_count=10,inc=200,overlap=6) !overlap must be >= 6
 
      integer cur_level,endx,endy,extrap,is_oreint,bf_type,elder_length,
     - int_count,total_count_s,bc_1_s,bc_2_s,limiter_on,
     - bc_3_s,bc_4_s,endx_s,endy_s,count,js,freq,theta,endx_sponge,
     - endy_sponge,startx_sponge,starty_sponge,count1,kk,n_deriv

      real err,maxerr,alpha,o,sum,spec_ts,
     -  ep,depth,wave_hgt,dt,dy,dx,end_t,end_y,c,co,x0,
     -  zc,itr,slope,slope_reg,depth_2,U_10,
     -  depth2_reg,init_reg,cpu_dt_inc,
     -  chan_width,upwind_baseline,CB,Ch_length,
     -  k_dim,per,filt_int,erru,errv,errz,
     -  alp,mod_co,x_c,y_c,beta,loc,L_r,eps

	parameter(eps=1e-6)

      real init_mass,step,aspect_ratio,inc_ang,
     -  energy,init_energy,vc,vc1,x02,slope_2,te,ts,ro,rr,ao,
     -  rad_depth,LL,TA(2),amp_top,mk,mk2,CK,CE,yt,yt2,rp1,
     -  r2,L,xso,sponge_width,sigma,cdamp1,cdamp2,delta,lambda,
     -  slot_ho,z,za,d_tmp1,d_tmp2,ramp,dom_wave,end_wave,rf,rd,
     -  uc,uc1,zc1,dom_loc,ramp_loc,cur_loc,th1,th2


      real s,t1,t2,t3,t4,t5,t6,tmp,endx_tmp,endy_tmp,
     -  y1,y2,depth_top,depth_bottom,depth_mid,
     -  slope_grid,grav,sum1,sum2,dudx_EFG,dvdy_EFG,du2dx,dv2dx,
     -  du2dy,dv2dy,dhdx,dhudx_EFG,swash_d,
     -  dhdy,dhvdy_EFG,w1x,w2x,w1y,w2y,dudyy,w1x_zeta,w1y_zeta,
     -  dvdxx,x_rel,x_l,dampfx,dampfy,x_r,corner(4,100),
     -  z_star,kappa,smooth_top,h_cutoff,end_reg,
     -  node_1,node_2,node_3,node_4,node_5,node_6,node_7,node_8,
     -  x_1,h_1,x_2,h_2,x_3,h_3,x_4,h_4,x_5,h_5,x_6,h_6,x_7,h_7,
     -  x_8,h_8,x_9,h_9,slope_1,slope_3,slope_4,slope_5,slope_6,
     -  slope_7,slope_8,writ_inc,courant,per_sld,mx_z,radius_coef,
     -  I_src_1,k_1,k_2,w_1,w_2,wave_hgt_2,L_2,per_2,hor_wall_loc,
     -  harbor_width,beta_src_1,beta_src_2,D_1,D_2,I_src_2,
     -  bas_rad,r_cur,xo,yo,r_o,AA,ww,f1_int_src,f2_int_src,
     -  dzdyy,dzdxx,param_E,param_F,param_G,conv_v,conv_u,
     -  f_BF,BF_x,BF_y,delta_breaking,dzdt_I,dzdt_star,dzdt,B,
     -  deriv_1x,deriv_2x,deriv_3x,deriv_1y,
     -  deriv_1,deriv_2,deriv_3,deriv_4,dSdt,dTdt,pres,
     -  deriv_2y,deriv_3y,t_star,dzdt_f,
     -  dyn_fric,static_fric,delta_sed,f_w,
     -  diameter,slope_x,beta_tide,D_tide,
     -  slope_y,shields_1,shields_2,p1,p2,Q1,Q2,ho_c,
     -  theta_B,theta_0,f_shape,
     -  theta_cur,del_elev,i_err,j_err,slope_e,
     -  amp,slide_wd,time_tmp2,t_o,t_o2,time_tmp,
     -  x1_o,h1_o,x2_o,d_o,h2_o,sld_wd,xm,hm,slope_ang,
     -  xc,hc,th1_o,th2_o,rate,vel,x1,x2,h1,h2,
     -  th,back_look,ho_amp,df_sl,zeta_o,tc,
     -  d,so,to,width_block,hs,end_slope,del_h,depth_dec,co_c,
     -  sum_u,sum_u_dif,sum_z,sum_z_dif,
     -  sum_v,sum_v_dif,erru2,errv2,errz2,err2,maxerr2,ep2,
     -  C_d,sigma_k,visc_water,
     -  sigma_ep,dkdx,dedx,dkdxx,dedxx,dvisctotal_kdx,dvisctotal_epdx,
     -  strain_xx,engy_xx,dkdy,dedy,dkdyy,dedyy,dvisctotal_kdy,
     -  dvisctotal_epdy,strain_xy,strain_yy,engy_xy,engy_yy,
     -  C_1ep,C_2ep,KE_x,
     -  KE_y,jmp_slp,h_swash,cutoff,mod_ang,a_abs,cpu_0,cpu_c,
     -  batch_time_limit,ep3,o_input,mass_1,
     -  cutoff_def,ep_def,ep_shift,cpu_dt,start_x,mass_tmp,k_tmp,w_tmp,
     -  L_min,L_max,per_min,per_max,shift_rn,cpu_last,cpu_per_step

      real sum_per_step,ave_per_step,end_x,end_x_t,end_y_t,
     -  tmp1,tmp2,count_cpu,start_y,frac,
     -  x_s,x_e,y_s,y_e,wave_hgt_corr,t_q_ss,gamma,Cd,Cm,Frd,
     -  U_mean(100,inc,2),W_mean(100,inc),
     -  x_ts(10000),y_ts(10000),o_c,hxx_max,hx_max,FR_cur,bf_ratio,
     -  dz,d1,d2,d3,d4,d5,d6,d7,d8,H_total,
     -  Num2,Num3,Den1,Den2,u_vec,coef_A,flux_x_tmp,
     -  coef_B,wid, Num4,Den3,max_hxx,fs_shift,amp_nl,I_src_nl,
     -  z_interp1,z_interp2,z_interp3,distance,Ch,
     -  z_interp4,z_interp5,z_interp6,z_interp7,z_interp8,z_interp9,
     -  C_1,C_2,C_3,C_4,a_1,a_3,is_coef,is_sin,c_roller,H_i,f_u,
     -  f_c,visc_coef,cos_cur,sin_cur,inc_ang_cur,x0_cur,D_src_cur,
     -  beta_src_cur,L_spec_cur,I_src_cur,amp_cur,pres_force,sum_mom,
     -  x_c2,y_c2,loc_line,y_coord, depth_tmp, Cr_c,Cr_max,max_h,
     -  mx_zeta,mx_zeta_g,mass_g,init_mass_g,u_z1,v_z1,w_z1,pi,ihvor

	character*12 filename
	character*4 ranknum

      integer, ALLOCATABLE :: bl_jump(:,:),
     - bl_hor_wall_topo(:,:),wvmk_loc_i(:),
     - level_topo(:,:),x_h(:),y_h(:),
     - bl_hor_wall(:,:),bl_x_wall(:,:),
     - bl_y_wall(:,:),level(:,:)

      real, ALLOCATABLE :: t(:),x(:),
     -  y(:),v_iter(:,:),
     -  zeta(:,:,:,:),
     -  E(:,:,:,:),bl(:,:),
     -  zeta_iter(:,:),ho(:,:),
     -  u_iter(:,:),u(:,:,:,:),
     -  v(:,:,:,:),wm_d_c(:),wm_u_c(:),
     -  hxx(:,:),hyy(:,:),
     -  a1(:,:,:),a2(:,:,:),
     -  b1(:,:,:),b2(:,:,:),
     -  bet(:),d1v(:,:,:),
     -  d2v(:,:,:),WB_x(:,:),WB_y(:,:),
     -  z_alp(:,:,:,:),
     -  UU(:,:,:,:),u_z(:,:),v_z(:,:),w_z(:,:)

	real, ALLOCATABLE ::  VV(:,:,:,:),
     -  F(:,:,:,:),
     -  F1(:,:,:,:),
     -  G(:,:,:,:),
     -  G1(:,:,:,:),
     -  dhzvdy(:,:,:),dExdx(:,:,:),
     -  dEydy(:,:,:),dhzudx(:,:,:),
     -  dvdxy(:,:,:),dhvdxy(:,:,:),
     -  dudxy(:,:,:),dhudxy(:,:,:),
     -  dE2xdx1(:,:),dE2ydy1(:,:),mass(:),
     -  hxy(:,:),
     -  dE2xdx2(:,:),dE2ydy2(:,:),dF2xdx1(:,:),
     -  dF2xdx2(:,:),dF2xdx3(:,:),
     -  dG2ydy1(:,:),dG2ydy2(:,:),dG2ydy3(:,:),
     -  dhdt(:,:,:),h(:,:,:),
     -  hx(:,:),hy(:,:),
     -  dhdxt(:,:,:),dhdyt(:,:,:),
     -  dhdtt(:,:,:),dhdxtt(:,:,:),dhdytt(:,:,:),
     -  dz_alpdt(:,:,:,:),
     -  dudxx(:,:,:),dhudxx(:,:,:),
     -  dvdyy(:,:,:),dhvdyy(:,:,:),
     -  d3u(:,:,:),
     -  d3v(:,:,:),
     -  hp(:,:,:),c_rad(:,:),hzu(:,:),hzv(:,:),hv(:,:),
     -  hu(:,:),E_x_grp(:,:),E_y_grp(:,:),
     -  E2_x_grp1(:,:),E2_x_grp2(:,:),E2_y_grp1(:,:),
     -  E2_y_grp2(:,:),
     -  dudx(:,:,:),dvdy(:,:,:),
     -  dhudx(:,:,:),dhvdy(:,:,:),
     -  uxpvy(:,:),huxphvy(:,:),
     -  F2_grp1(:,:),
     -  dudt(:,:,:),dvdt(:,:,:),dudxt(:,:),dvdyt(:,:),
     -  dhudxt(:,:),dhvdyt(:,:)

	real, ALLOCATABLE :: x_tmp(:),y_tmp(:),
     -  ho_tmp(:,:),dzdx(:,:,:),
     -  bl_x_tmp(:,:),dzdy(:,:,:),
     -  bl_y_tmp(:,:),dudy(:,:,:),dvdx(:,:,:),
     -  mov_sh_lezs(:,:),
     -  mov_sh_gtzs(:,:),d1u(:,:,:),
     -  d2u(:,:,:),
     -  eddy_visc(:,:),dhzudy(:,:,:),dhzvdx(:,:,:),
     -  t_breaking(:,:),VG(:,:,:),
     -  net_sed_trans(:,:),roller(:,:),deltar(:,:),
     -  loc_sum(:,:),k_turb(:,:,:),ep_turb(:,:,:),
     -  array(:,:),array_2(:,:),
     -  stress_xx(:,:,:),stress_xy(:,:),stress_yy(:,:),
     -  amp_spec(:,:),per_spec(:,:),
     -  L_spec(:,:),
     -  beta_src(:,:),I_src(:,:),
     -  D_src(:,:),
     -  shift_spec(:,:),f_int_src(:,:),
     -  f_int_src_v(:,:),cutoff_mat(:,:),cutoff_mat_d(:,:)

	real,  ALLOCATABLE ::  inc_ang_spec(:,:),
     -  max_fs(:,:),max_vel(:,:),bet2(:),v1_iter(:,:),
     -  u1_iter(:,:),zeta1_iter(:,:),fs(:),fa(:),
     -  S_grp(:,:,:),T_grp(:,:,:),
     -  dSdx(:,:,:),dSdy(:,:,:),
     -  dTdx(:,:,:),dTdy(:,:,:),
     -  mean_z(:,:,:),cur_x(:,:,:),
     -  cur_y(:,:,:),mean_hgt(:,:,:),
     -  xi(:),xi2(:),cosA(:,:),h_tmp(:,:),
     -  sinA(:,:),D_src_nl(:,:),
     -  NL_a1(:,:,:),NL_a2(:,:,:),
     -  NL_b1(:,:,:),NL_b2(:,:,:),S_tmp(:),T_tmp(:),
     -  wA(:,:),kA(:,:),tmpA(:,:,:,:,:),
     -  fse(:),time_bm(:),eta_in(:),
     -  curf(:),time_cur(:),current_in(:),
     -  z_br(:,:),k_br(:,:),u0_b(:,:,:),vort(:,:)

	real,  ALLOCATABLE ::  x_glob(:),
     -  y_glob(:),x_topo(:),y_topo(:),topo(:,:),
     -  topo_tmp(:,:),ho_glob(:,:)

	real,  ALLOCATABLE ::  dhvdx(:,:,:),dhudy(:,:,:),tmp_mat(:,:),
     -  visc_rot_x(:,:),visc_rot_y(:,:),a_x(:,:),b_x(:,:),c_x(:,:),
     -  d_x(:,:),a_y(:,:),b_y(:,:),c_y(:,:),d_y(:,:),v_tmp(:,:),
     -  kt(:,:,:),nut(:,:,:),k_RHS(:,:,:),B_mat(:,:,:),B_RHS(:,:,:),
     -  psix(:,:,:),psiy(:,:,:),loc_mat(:,:,:,:),x_wvmk(:),u_wvmk(:)


      real, ALLOCATABLE :: wv_signal(:,:),wm_d(:,:)

	integer,  ALLOCATABLE ::  xs(:),xe(:),ys(:),ye(:)

	integer*2 sta
	character clave*4,flou4*12
	logical isperiodic(2), reorder

	parameter(pi=3.1415927,grav=9.81)
		 
      end module mainvar_module



        
         



            





        
