      subroutine load_inputs

      use mainvar_module
	include 'mpif.h'
      real res_choice
	
            open(8,file='sim_set.dat',status='old')
            read(8,*) sim_opt                  ! Type of simulation,      1=surface wave evolution
            read(8,*) dim                        ! # of dimensions 1=1D, 2=2D
            read(8,*) deriv_order_ind   ! Order to which higher-order dispersive terms are finite-differenced: 1=all terms to at least O(dx**4)2=all terms to at least O(dx**2) recommended
            read(8,*) nonlin_ind            ! Include full nonlinear effects? 1=no, 2=yes, 0=linear simulation
            read(8,*) disp_prop                  ! Use Nwogu's dispersion 'depth'=1, or depth-average=2 (only for weakly nonlinear case)
            read(8,*) conv                        ! Use 'proper' convection=1, or Wei&Kirby's convection=2
            read(8,*) wave_type                  ! if wave_type=1, use input parameters below. Value
            read(8,*) depth                  ! water depth in meter
            read(8,*) wave_hgt                  ! solitary wave hieght in meter
            read(8,*) depth                        ! water depth in meter
            read(8,*) x0                        ! Initial location of crest of solitary wave
            read(8,*) mk2                        ! for cnoidal wave - square of the modulus of the 
            read(8,*) ramp                        ! # of waves to ramp
            read(8,*) dom_wave                  ! # of full (non-ramped) waves in initial domain
            read(8,*) end_wave                  ! # of waves to be created
            read(8,*) inc_ang                  ! Oblique angle of incident wave
            read(8,*) slide_type            ! Type of landslide, only for sim_opt=2.  Slide_type=1
            read(8,*) ts                        ! Start time of landslide - must be >0
            read(8,*) per            ! Time scale of lanslide
            read(8,*) te                  ! End time of landslide 
            read(8,*) filt                        ! Filter entire domain=1, at increment filt_int
            read(8,*) filt_int                  ! # of time steps between filtering
            read(8,*) filtsld_int            ! # of time steps between filtering over lanslide area (only for sim_opt=2)
            read(8,*) screen_output            ! print data to screen if =1
            read(8,*) display                  ! Increment to display information on screen
            read(8,*) itr                  ! Maximum # of iterations for Corrector Loop (will display
            read(8,*) min_itr            ! Minimum # of iterations for Corrector Loop      
            read(8,*) itr_or            ! Iteration after which over-relation should be used
            read(8,*) o                  ! Over-relaxtion coefficient
            read(8,*) bc_1                  ! Left side domain boundary condition
            read(8,*) bc_2                  ! Right side domain boundary condition
            read(8,*) bc_3                  ! Bottom side domain boundary condition
            read(8,*) bc_4                  ! Top side domain boundary condition
            read(8,*) spng_1            ! Use sponge layer on left side boundary when=1
            read(8,*) spng_2            ! Use sponge layer on right side boundary when=1
            read(8,*) spng_3            ! Use sponge layer on bottom side boundary when=1
            read(8,*) spng_4            ! Use sponge layer on top side boundary when=1      
            read(8,*) load_topo                  ! if load_topo=2, then use bathymetry data files created
            read(8,*) num_nodes            ! number of nodes required to create bottom profile, including first and last
            read(8,*) x_1
            read(8,*) h_1
            read(8,*) x_2
            read(8,*) h_2
            read(8,*) x_3
            read(8,*) h_3
            read(8,*) x_4
            read(8,*) h_4
            read(8,*) x_5
            read(8,*) h_5
            read(8,*) x_6
            read(8,*) h_6
            read(8,*) x_7
            read(8,*) h_7
            read(8,*) x_8
            read(8,*) h_8
            read(8,*) x_9
            read(8,*) h_9
            read(8,*) chan_width          ! width in meter of channel (ignore for 1D simulation)
            read(8,*) end_t                              ! time in seconds to run simulation
            read(8,*) writ_inc                        ! time interval to save to file, i.e. if a=10,
            read(8,*) pts_wvl            
            read(8,*) courant
            read(8,*) slide_node      
            read(8,*) L
            read(8,*) per_sld
            read(8,*) slope_reg
            read(8,*) sponge_width
            read(8,*) cdamp1
            read(8,*) cdamp2
            read(8,*) rotate
            read(8,*) smooth
            read(8,*) y1
            read(8,*) y2
            read(8,*) depth_top
            read(8,*) depth_bottom
            read(8,*) smooth_top
            read(8,*) radius_coef
            read(8,*) amp_top
            read(8,*) int_src
            read(8,*) wave_hgt_2
            read(8,*) bottom_fric
            read(8,*) f_BF       
            read(8,*) wave_breaking
            read(8,*) ep      
            read(8,*) ep2
            read(8,*) L_2
            read(8,*) swash_filt
            read(8,*) smooth_bathy
            read(8,*) use_av
            read(8,*) spec_type
            read(8,*) slope_ang
            read(8,*) gamma
            read(8,*) Cd
            read(8,*) Cm
            read(8,*) Frd
            read(8,*) num_ts
            read(8,*) offshore_cur
            read(8,*) FR_cur
            read(8,*) bf_ratio
            read(8,*) sh_mov 
            read(8,*) cutoff
		  read(8,*) num_levels 
		  read(8,*) aspect_ratio
					if(aspect_ratio.lt.1e-6) aspect_ratio=1.

		  read(8,*) is_oreint
		  read(8,*) upwind_baseline
		  read(8,*) visc_coef
		  read(8,*) rotationality
		  read(8,*) decomp_type
		  read(8,*) dims(1)
		  read(8,*) dims(2)
		  read(8,*) spec_pp
		  read(8,*) spec_ts
		  read(8,*) numerical_scheme
		  read(8,*) limiter_on
		  read(8,*) breaker_type
		  read(8,*) bf_type
		  read(8,*) Ch
		  read(8,*) backscatter 
		  read(8,*) CB 
		  read(8,*) ihvor
		  read(8,*) elder_length  
		  read(8,*) Ch_length
		  read(8,*) plus_tide
		  read(8,*) swash_d
		  read(8,*) left_wvmk
		  read(8,*) upwind_shore
		  read(8,*) step_shore
		  
          close(8)

		  if(batch.eq.2)then
           open(8,file='config_GUI.dat',status='old')
           read(8,*) res_choice
           read(8,*) f_BF  
           read(8,*) end_t 
           close(8)
          
           if(INT(res_choice).eq.0)then
               pts_wvl=100
           elseif(INT(res_choice).eq.1)then 
               pts_wvl=200
           elseif(INT(res_choice).eq.2)then 
               pts_wvl=400
			   if(dim.eq.2) pts_wvl=300
           else 
               pts_wvl=600
           endif
          endif
		              
		  if(num_ts.gt.0)then
			open(96,file='ts_locations.dat',status='old')
			do ii=1,num_ts
                  read(96,*) x_ts(ii),y_ts(ii)
			enddo                  
			close(96)

c			open(96,file='depth_input.dat',status='old')
c			read(96,*) depth
c			close(96)

		  endif

      return

      end

