      subroutine FV_allocate_matrices
      use mainvar_module, only:nx,ny,endx,endy,visc_coef,dim,overlap,
     -	bottom_fric,f_BF,limiter_on,backscatter,Ch,CB,ihvor,
     -	num_theta,num_freq
      use FV_var_module

	sx=4
	ex=endx-3

	if(dim.eq.1)then
		sy=overlap+1
		ey=endy-overlap	
	else
		sy=4
		ey=endy-3
	endif

	alpha = -0.531
	
	ilim = limiter_on  ! 1=use flux limiter scheme, with b0 and b1 coefficeints
      b0 = 2.0
      b1 = 2.0
!	ihvor = 1  ! 1=use horizontal vorticity model, speciefied now through user interface
	icor = 0  ! 0=use centered CV difference for dispersive terms with interface interpolation given by interp2, 1=use limiter, with limiter given by limiter2
!	backscatter =1 ! 0=no backscatter model, 1=backscatter model used, for 2HD only, speciefied now through user interface
!	CB=70.           ! backscatter coefficient , speciefied now through user interface
	interp2 = 3 ! interpolation to use for cell interface evaluation
!      IF(interp2.EQ.1) explicit interpolation 4th order
!      IF(interp2.EQ.2) implicit (compact) interpolation 4th order
!      IF(interp2.EQ.3) explicit interpolation 2nd order

	limiter2 = 3 ! limiter to use with icor=1 option
	limiter3 = 3 ! limiter to use for 2nd order (near shoreline) flux scheme ! not used in current model version

!      IF(limiter.EQ.1) 1st order
!      IF(limiter.EQ.2) Centred 2nd order
!      IF(limiter.EQ.3) Superbee like limiter
!      IF(limiter.EQ.4) Minimod like limiter
!      IF(limiter.EQ.5) Minimod like limiter 2

	iriemann = 2 ! type of riemann solver, 2 is generally the fastest
!	Ch = 0.07, speciefied now through user interface
	cm = visc_coef**2.      
	ks = f_BF     

      idissip	= bottom_fric         ! dissipation flag (0=no dissipation, 1=with dissipation)
      ibotfric = bottom_fric       ! bottom friction flag (0=no friction, 1=with friction)


      ALLOCATE(z_FV(nx,ny),d_FV(nx,ny),u_FV(nx,ny),
     - v_FV(nx,ny),E_FV(nx,ny),F_FV(nx,ny),F1_FV(nx,ny),G_FV(nx,ny),
     - G1_FV(nx,ny),LOWERx(nx,ny),DIAGx(nx,ny),UPPERx(nx,ny),
     - RHSx(nx,ny), LOWERy(nx,ny),DIAGy(nx,ny),UPPERy(nx,ny),
     - RHSy(nx,ny), ix(nx,ny),iy(nx,ny),ix2(nx,ny),iy2(nx,ny),
     - id(nx,ny),zx(nx,ny),zy(nx,ny),rand(nx,ny,2),ish_mat(nx,ny))
      return

      end

