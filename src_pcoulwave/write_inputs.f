      subroutine write_inputs

      use mainvar_module
	include 'mpif.h'

            open(8,file='sim_set.dat',status='unknown')
            write(8,*) sim_opt                  ! Type of simulation,      1=surface wave evolution
            write(8,*) dim                        ! # of dimensions 1=1D, 2=2D
            write(8,*) deriv_order_ind        ! Order to which higher-order dispersive terms are finite-differenced: 1=all terms to at least O(dx**4)2=all terms to at least O(dx**2) recommended
            write(8,*) nonlin_ind            ! Include full nonlinear effects? 1=no, 2=yes, 0=linear simulation
            write(8,*) disp_prop                  ! Use Nwogu's dispersion 'depth'=1, or depth-average=2 (only for weakly nonlinear case)
            write(8,*) conv                        ! Use 'proper' convection=1, or Wei&Kirby's convection=2
            write(8,*) wave_type                  ! if wave_type=1, use input parameters below. Value
            write(8,*) depth                  ! water depth in meter
            write(8,*) wave_hgt            ! amplitude of landslide displacement
            write(8,*) depth                        ! water depth in meter
            write(8,*) x0                        ! Initial location of crest of solitary wave
            write(8,*) mk2                        ! for cnoidal wave - square of the modulus of the 
            write(8,*) ramp                        ! # of waves to ramp
            write(8,*) dom_wave                  ! # of full (non-ramped) waves in initial domain
            write(8,*) end_wave                  ! # of waves to be created
            write(8,*) inc_ang                  ! Oblique angle of incident wave
            write(8,*) slide_type            ! Type of landslide, only for sim_opt=2.  Slide_type=1
            write(8,*) ts                        ! Start time of landslide - must be >0
            write(8,*) per            ! Time scale of lanslide
            write(8,*) te                  ! End time of landslide 
            write(8,*) filt                        ! Filter entire domain=1, at increment filt_int
            write(8,*) filt_int                  ! # of time steps between filtering
            write(8,*) filtsld_int            ! # of time steps between filtering over lanslide area (only for sim_opt=2)
            write(8,*) screen_output            ! print data to screen if =1
            write(8,*) display                  ! Increment to display information on screen
            write(8,*) itr                  ! Maximum # of iterations for Corrector Loop (will display
            write(8,*) min_itr            ! Minimum # of iterations for Corrector Loop      
            write(8,*) itr_or            ! Iteration after which over-relation should be used
            write(8,*) o                  ! Over-relaxtion coefficient
            write(8,*) bc_1                  ! Left side domain boundary condition
            write(8,*) bc_2                  ! Right side domain boundary condition
            write(8,*) bc_3                  ! Bottom side domain boundary condition
            write(8,*) bc_4                  ! Top side domain boundary condition
            write(8,*) spng_1            ! Use sponge layer on left side boundary when=1
            write(8,*) spng_2            ! Use sponge layer on right side boundary when=1
            write(8,*) spng_3            ! Use sponge layer on bottom side boundary when=1
            write(8,*) spng_4            ! Use sponge layer on top side boundary when=1      
            write(8,*) load_topo                  ! if load_topo=2, then use bathymetry data files created
            write(8,*) num_nodes            ! number of nodes required to create bottom profile, including first and last
            write(8,*) x_1
            write(8,*) h_1
            write(8,*) x_2
            write(8,*) h_2
            write(8,*) x_3
            write(8,*) h_3
            write(8,*) x_4
            write(8,*) h_4
            write(8,*) x_5
            write(8,*) h_5
            write(8,*) x_6
            write(8,*) h_6
            write(8,*) x_7
            write(8,*) h_7
            write(8,*) x_8
            write(8,*) h_8
            write(8,*) x_9
            write(8,*) h_9
            write(8,*) chan_width          ! width in meter of channel (ignore for 1D simulation)
            write(8,*) end_t                              ! time in seconds to run simulation
            write(8,*) writ_inc                        !
            write(8,*) pts_wvl            
            write(8,*) courant
            write(8,*) slide_node
            write(8,*) L      
            write(8,*) per_sld
            write(8,*) slope_reg
            write(8,*) sponge_width
            write(8,*) cdamp1
            write(8,*) cdamp2
            write(8,*) rotate
            write(8,*) smooth
            write(8,*) y1
            write(8,*) y2
            write(8,*) depth_top
            write(8,*) depth_bottom
            write(8,*) smooth_top
            write(8,*) radius_coef
            write(8,*) amp_top
            write(8,*) int_src
            write(8,*) wave_hgt_2
            write(8,*) bottom_fric
            write(8,*) f_BF       
            write(8,*) wave_breaking
            write(8,*) ep      
            write(8,*) ep2
            write(8,*) L_2
            write(8,*) swash_filt
            write(8,*) smooth_bathy
            write(8,*) use_av
            write(8,*) spec_type
            write(8,*) slope_ang
            write(8,*) gamma
            write(8,*) Cd
            write(8,*) Cm
            write(8,*) Frd
            write(8,*) num_ts
            write(8,*) offshore_cur
            write(8,*) FR_cur
            write(8,*) bf_ratio
            write(8,*) sh_mov 
            write(8,*) cutoff
		  write(8,*) num_levels 
		  write(8,*) aspect_ratio
		  write(8,*) is_oreint
		  write(8,*) upwind_baseline
		  write(8,*) visc_coef
		  write(8,*) rotationality
		  write(8,*) decomp_type
		  write(8,*) dims(1)
		  write(8,*) dims(2)
		  write(8,*) spec_pp
		  write(8,*) spec_ts
		  write(8,*) numerical_scheme
		  write(8,*) limiter_on
		  write(8,*) breaker_type
		  write(8,*) bf_type
		  write(8,*) Ch
		  write(8,*) backscatter 
		  write(8,*) CB 
		  write(8,*) ihvor  
		  write(8,*) elder_length  
		  write(8,*) Ch_length
		  write(8,*) plus_tide
		  write(8,*) swash_d
		  write(8,*) left_wvmk
		  write(8,*) upwind_shore
		  write(8,*) step_shore
			int_src=1

            do i=1,20
                  write(8,*) 0.0 ! Dummy space for future parameter additions
            enddo
            close(8) 

            if(num_ts.gt.0)then
             open(96,file='ts_locations.dat',status='unknown')
             do i=1,num_ts
                  write(96,*) x_ts(i),y_ts(i)
             enddo                  
             close(96)
            endif
      return

      end

