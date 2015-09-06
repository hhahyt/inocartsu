      program coulwave_MPI

c      INCLUDE 'FLIB.FI'
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC PROGRAM NOTES CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      COULWAVE - Copyright 2008 by Patrick Lynett, Texas A&M University
C
C      File originally created on 8/3/99 by Pat Lynett, with regular revisions
C
C      This program is primarily a recreation of the Wei and Kirby, Fully-Nonlinear model given in
C      JFM, 1995, with some major and minor additions and revisions.  Modifcations include the ability
C	 to simulate landslide generated waves, use of the multi-layer theory, parallelization,
C	 conversion to a finite volume method, and a turbulence-vorticity closure model
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	use mainvar_module
	use FV_var_module, only: id
	include 'mpif.h'


c	! INITIALIZE MPI AND RELATED MPI SETUP 
	call mpi_init(ierr)                              ! start mpi
	call mpi_comm_size(mpi_comm_world,nprocs,ierr)   ! find number of process
	call mpi_comm_rank(mpi_comm_world,myrank,ierr)   ! find rank 

	call mpi_barrier(mpi_comm_world,ierr)  !wait for all processes to get here

	if(myrank.eq.0)then
		call user_interface  !only lead process runs interface
	endif

	call mpi_barrier(mpi_comm_world,ierr)

C............	all processes load inputs from interface
	do i=0,nprocs-1
	  if(myrank.eq.i)then
		    call load_inputs 
	  endif
	  call mpi_barrier(mpi_comm_world,ierr) ! barrier to prevent numerous simult accesses
	enddo

C............	Create files for all processes
	call create_file(myrank) 


C	SETUP DOMAIN DECOMPOSITION
	if(dim.eq.1)then
		dims(1) = nprocs   ! dimension in x direction
		dims(2) = 1   ! dimension in y direction
	else
		if(decomp_type.eq.0)then
				dims(1)=0
				dims(2)=0
			call mpi_dims_create(nprocs,2,dims,ierr )
		endif
	endif

	isperiodic(1) = .false. ! periodicity of row
	isperiodic(2) = .false. ! periodicity of column
	reorder = .true.

	call mpi_cart_create(mpi_comm_world,2,dims,
     -	isperiodic,reorder,comm2d,ierr)
	call mpi_cart_coords(comm2d,myrank,2,coords,ierr)
	call mpi_cart_shift(comm2d,0,1,pleft,pright,ierr)
	call mpi_cart_shift(comm2d,1,1,pbottom,ptop,ierr)

C............	find wave lengths, periods, and speeds for the various types of IC's
	call find_wavelengths 

C............	find global grid spacings
	call find_dxdydt

C............	determine the dimensions of the global (un-decomposed) domain
      call global_dims

C............	Load or create global depth grid
	  if(batch.eq.2)then
		call create_depth_wvmk
	  else
		call create_depth
	  endif
	call mpi_barrier(mpi_comm_world,ierr) ! wait for all processes to load depth data

C............	Here is where the global depth
C............	domain and the x and y's must be split up amongst the processes
	call decomp2d(dims,coords,endx_glob,endy_glob,sx,ex,sy,ey,
     -		bc_1,bc_2,bc_3,bc_4,overlap,sx_X,ex_X,sy_X,ey_X,ho_glob)

	endx=ex_X-sx_X+1
	endy=ey_X-sy_X+1

c	! CREATE NEW MPI VARIABLE TYPE
	call mpi_type_vector(endy,overlap,endx,
     -		mpi_real,LRstride,ierr)
	call mpi_type_vector(overlap,endx,endx,
     -		mpi_real,BTstride,ierr)
	call mpi_type_commit(LRstride,ierr)
	call mpi_type_commit(BTstride,ierr)


c	! CREATE NEW COMMUNICATORS
	if(dims(1).ne.1)then
	   row = coords(2)
	   call mpi_comm_split(comm2d,row,myrank,xcomm,ierr)
	   call mpi_comm_rank(xcomm,xrank,ierr)
	endif

	if(dims(2).ne.1)then
	   col = coords(1)
	   call mpi_comm_split(comm2d,col,myrank,ycomm,ierr)
	   call mpi_comm_rank(ycomm,yrank,ierr)
	endif

C............	allocate all the variables for each process
	call allocate_matrices

CCCCCCCCCCCCCC Display a few parameters on the screen CCCCCCCCCCCCCCCCCCCCCCCC
      if(screen_output.eq.1)then
	if(myrank.eq.0)then
      print*,'**************************************************'
      print*, 'X Start','Left Sponge Boundary',
     -            'Right Sponge Boundary','X End'
      print*,0,sponge_width*spng_1,end_x_t-sponge_width*spng_2,end_x_t
       print*,'**************************************************'
      print*, 'Y Start','Bottom Sponge Boundary',
     -            'Top Sponge Boundary','Y End'
      print*,0,sponge_width*spng_3,end_y_t-sponge_width*spng_4,end_y_t
       print*,'**************************************************'

      print*,'dx (m)=     ',dx
      print*,'dt (s)=     ',dt,'  co= ',co
	  print*,'dx/depth=  1/',ho(1,1)/dx, 
     -            ',should be <~1/3 to ensure stability'
	print*,'***************************************************'
	print*,'Processor Grid Divisions (rank, global i&j bounds)'
	endif
	
	call mpi_barrier(mpi_comm_world,ierr)
	print*,myrank,sx,ex,sy,ey
	endif

	call mpi_barrier(mpi_comm_world,ierr)

C******************************************************************************
C............sets the j loop limits, depending on whether simulation is 1D or 2D
	call set_loop_limits
	
C............create x y t vectors and h matrix and write to file
	call create_write_xyth 
	
C............Set up vertical walls surrounding numerical domain
	call bl_define 
	
C............Calculates hx and hxx values, and finds the maximum slope and curvature of the seafloor for each process
	call find_max_hxx  
	
C............Assemble dispersion properties matrices, bet, etc, calc z_alpha, zeta2 values, and calc linear tridiag coefficents
	call set_dispersion_coefs 
	
C............Assemble internal source matrices
	call set_internal_source_coefs 
	
C............Find i,j indeces for all time series
	call find_ts_indices 
	
C............ zero-out main variables, for compilers that may not do this
	call init_variables  
	
C............ if using a "current" boundary, load input file here
	if(offshore_cur.ne.0) call load_current_in
	
CCCCCCCCCCCCCCCCCCCCCCC  Initial Conditions CCCCCCCCCCCCCCCCCCC
    
      do n=1,4
	  
C............If not using internal source, calculates the hot start initial condition
		call solit_cnoidal_ic

C..................subroutine 'dhdt_calc' determines the change in
C.................bottom shape due to landslide or other phenomenon
		if(sim_opt.eq.2) call dhdt_calc(n,n)

		call move_shoreline(n)  ! adjust values of bl_hor_wall (wet-dry boolean)

cZZZZZZZZZZ
		do cur_level=1,num_levels
cZZZZZZZZZZ
C.........subroutine 'u2UU' calculates UU based on the initial
C..........  velocity profile  (see procedure for detail)
			call u2UU(nx,ny,nt,n,endx,endy,u,UU,h,hxx,b1,b2,dx,
     -                                    v,VV,hyy,dy,hx,hy,
     -                                    cur_level,num_levels,zeta,
     -                                    z_alp,nonlin_ind,dim)
C.............Enforce boundary conditions 
			call bc(n,n)
			if(left_wvmk.eq.1) call bc_wvmk(n,n)
cZZZZZZZZZZ
		enddo
cZZZZZZZZZZ

C.........subroutine 'FV_var_grp' groups the physical variables into the 
C........... finite volume based form of Kim & Lynett  
            if(n.eq.1) call FV_allocate_matrices		
		    call FV_var_grp(n,n)		    
	enddo

C.........subroutine 'mass' calculates the total mass
C............  in the numerical domain  (see procedure for detail)
C............ argument 'init_mass' is used to compare the mass at later times
      call calc_mass(init_mass,3) 

	call mpi_barrier(mpi_comm_world,ierr)

	call MPI_ALLREDUCE(init_mass, init_mass_g, 1,MPI_REAL, 
     -					MPI_SUM, mpi_comm_world,ierr) 

      if(init_mass_g.lt.1e-6)then
            init_mass_g=1.
      endif

C..........Write initial conditions to file
	call write_surfaces(3,3)
	call write_timeseries(3,3)


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCC   Begin finite difference method to determine zeta,u,v   CCCCCCCCCCCCC
      
      do nn=3,endt-1                  ! Time loop
      
49    continue
      
		n=3 ! the n, or working time, index should always be =3

		if(mod(nn,display).eq.0)then ! display info

C.......... Check to make sure that mass is being conserved
			call calc_mass(mass_1,3) 

C.......... Find the max free surface value in the domain
			call find_max_zeta(mx_zeta,3)

			call mpi_barrier(mpi_comm_world,ierr)

			call MPI_ALLREDUCE(mx_zeta, mx_zeta_g, 1,MPI_REAL, 
     -					MPI_MAX, mpi_comm_world,ierr) 

			call MPI_ALLREDUCE(mass_1, mass_g, 1,MPI_REAL, 
     -					MPI_SUM, mpi_comm_world,ierr) 
 
C..........Display and write to file some data
			if(myrank.eq.0) call write_disp_info(nn)

		endif

		call mpi_barrier(mpi_comm_world,ierr)

		if(sh_mov.eq.0) call move_shoreline(n)  ! adjust values of bl_hor_wall (wet-dry boolean)

		call exchange2d_int(endx,endy,overlap,pleft,pright,pbottom,
     -		ptop,LRstride,BTstride,comm2d,dims,bl_hor_wall(:,:),
     -		dim)  ! exhange shoreline locations

		call exchange2d_int(endx,endy,overlap,pleft,pright,pbottom,
     -		ptop,LRstride,BTstride,comm2d,dims,level(:,:),
     -		dim)  ! exhange dry cell "depth"


CCCCCCCCCCCCCCCCCCCCCCCCCCC  Begin Predictor Calculations CCCCCCCCCCCCCCCCC

C.............. all the work
	  do cur_level=1,num_levels
C.............Enforce boundary conditions 
			call bc(n,nn)
			if(left_wvmk.eq.1) call bc_wvmk(n,nn)
	  enddo

C.........subroutine 'FV_var_grp' groups the physical variables into the 
C........... finite volume based form of Kim & Lynett  		
		    call FV_var_grp(n,nn)		    

C.......exchange boundary variables - domain decomposition
cZZZZZZZZZZ
	  do cur_level=1,num_levels
cZZZZZZZZZZ
		call exchange2d(endx,endy,overlap,pleft,pright,pbottom,
     -		ptop,LRstride,BTstride,comm2d,dims,zeta(:,:,4,cur_level),
     -		dim)

		call exchange2d(endx,endy,
     -		overlap,pleft,pright,pbottom,
     -		ptop,LRstride,BTstride,comm2d,dims,u(:,:,4,cur_level),
     -		dim)
		if(dim.eq.2)then
			call exchange2d(endx,endy,
     -		overlap,pleft,pright,pbottom,
     -		ptop,LRstride,BTstride,comm2d,dims,v(:,:,4,cur_level),
     -		dim)
		endif
cZZZZZZZZZZ
	  enddo
cZZZZZZZZZZ

	  if(wave_breaking.eq.1)then
	   if(breaker_type.eq.1)then
		call exchange2d(endx,endy,overlap,pleft,pright,pbottom,
     -		ptop,LRstride,BTstride,comm2d,dims,t_breaking(:,:),
     -		dim)
	   elseif(breaker_type.eq.0)then
		call exchange2d(endx,endy,overlap,pleft,pright,pbottom,
     -		ptop,LRstride,BTstride,comm2d,dims,B_mat(:,:,4),
     -		dim)  
       endif 
	  endif


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  End Predictor Stage CCCCCCCCCCCCCCCCCCCCCCCC

C.............. update depth for landslide problem		
		if(sim_opt.eq.2) call dhdt_calc(4,nn+1) 

CCCCCCCCCCCCCCCCCCCCCCCCCCCC  Begin Iterative Corrector Stage  CCCCCCCCCCCC      

		maxerr=100      ! initialize maximum error in first
		maxerr2=100      ! initialize maximum error in first
		count_min=0      ! initialize number of iterations to zero
          do while((maxerr.gt.ep.and.maxerr2.gt.ep*100).or.
     -                        count_min.lt.min_itr)            ! start iterative loop

			count_min=count_min+1  ! iteration number counter

		    if(sh_mov.eq.0.and.count_min.le.4) call move_shoreline(n+1)  ! adjust values of bl_hor_wall (wet-dry boolean)

C...........  Store previous iteration values
C...........      in '_iter' for error calculation
			call store_previous_iter(4)

C.............. all the work
		
		  do cur_level=1,num_levels
C.............Enforce boundary conditions 
			call bc(n+1,nn+1)
			if(left_wvmk.eq.1) call bc_wvmk(n+1,nn+1)

		  enddo	
		  			                        
C.........subroutine 'FV_var_grp' groups the physical variables into the 
C........... finite volume based form of Kim & Lynett  		
		        call FV_var_grp(n+1,nn+1)		    
			
C.......exchange boundary variables - domain decomposition
cZZZZZZZZZZ
		  do cur_level=1,num_levels
cZZZZZZZZZZ
			call exchange2d(endx,endy,overlap,pleft,pright,pbottom,
     -		 ptop,LRstride,BTstride,comm2d,dims,zeta(:,:,4,cur_level),
     -		 dim)

			call exchange2d(endx,endy,
     -		 overlap,pleft,pright,pbottom,
     -		 ptop,LRstride,BTstride,comm2d,dims,u(:,:,4,cur_level),
     -		 dim)
			if(dim.eq.2)then
			 call exchange2d(endx,endy,
     -		  overlap,pleft,pright,pbottom,
     -		  ptop,LRstride,BTstride,comm2d,dims,v(:,:,4,cur_level),
     -		  dim)
			endif
cZZZZZZZZZZ
		  enddo
cZZZZZZZZZZ

		  if(wave_breaking.eq.1)then
		   if(breaker_type.eq.1)then
		    call exchange2d(endx,endy,overlap,pleft,pright,pbottom,
     -		ptop,LRstride,BTstride,comm2d,dims,t_breaking(:,:),
     -		dim)
	       elseif(breaker_type.eq.0)then
		    call exchange2d(endx,endy,overlap,pleft,pright,pbottom,
     -		ptop,LRstride,BTstride,comm2d,dims,B_mat(:,:,4),
     -		dim)  
           endif 
	      endif

C...............compute maximum error 
			call calc_corrector_error(4)
						  
		enddo
CCCCCCCCCCCCCCCCCCCCCCCCCCCC  End Iterative Corrector Loop CCCCCCCCCCCCCCCCCCCC


C...............Internal source, if added post-calculation (Type 2)
		if(int_src.eq.2.or.plus_tide.eq.1) call internal_source_type2(4,nn+1)

C................Determine if simulation has crashed 

		call calc_overflow(4)
		if(overflow.eq.1)then
			print*,nn,'SIMULATION ERROR - OVERFLOW'
			goto 50
		endif


C................Calc max Courant #	
		call calc_Courant_number(4)

C..............Calculate Maximum Free Surface Height, Mean Water Level,
C..............Mean Current, and Average Wave Hieght
		call calc_maxs_means(4,nn+1)
                
C..............Write data to file
		if(mod(nn,a).eq.0)then
              call write_surfaces(4,nn+1)
	    endif

		if(mod(nn,3).eq.0)then  ! write every third value to time series
			call write_timeseries(4,nn+1)
		endif

C..............shift matrices back one time level for next calculation
		call shift_matrices_back(4)  

      enddo
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      end time step loop CCCCCCCCCCCCCCCCCCCCCCCCCCC

C..............Write maxs and means data to file
	call write_maxs_means

	if(myrank.eq.0)then
		write(91,*) 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
		write(91,*) 'Simulation Completed!'
		write(91,*) 'Periods included in ave value computations=',
     -				(t(nn)-t_q_ss)/nt_per
	endif

50    continue

C.............Close all time series files (there can be so many that sometimes need to close out in code)
      do cur_ts=1,num_ts
		if(ind_ts(cur_ts,1).eq.1)then
			close(10119+cur_ts)
		endif
	enddo
	call mpi_barrier(mpi_comm_world,ierr)

C.............Perform spectral analysis on recorded time series
	if(spec_pp.eq.1.and.myrank.eq.0)then
		call spectral_analysis(num_ts,endt,spec_ts,end_t*0.99,dt)
	endif
	call mpi_barrier(mpi_comm_world,ierr) !wait until ts analysis is finished to close program

	call mpi_type_free(LRstride,ierr)
	call mpi_type_free(BTstride,ierr)
	call mpi_finalize(ierr)

      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  end program CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


    