      subroutine user_interface

      use mainvar_module
	include 'mpif.h'


      open(99,file='batch.dec',status='unknown')
      read(99,*) batch

c		batch=2 ! for OSU wavemaker	  

c      batch=0
c      batch_time_limit=10000000

	call load_inputs

	if(batch.ne.0)then
		screen_output=0 
	endif
	
	  if(batch.eq.1)goto 111   ! skip UI

		PRINT*,'********************************************************'
		PRINT*,'********************************************************'
		Print*,'                        pCOULWAVE'
		Print*,' Copyright 2011 by Patrick Lynett, Texas A&M University'
		PRINT*,'Modeling Wave Generation, Evolution, and Interaction with'
		PRINT*,'    Depth-Integrated, Dispersive Long Wave Equations'
		PRINT*,'               Parallel MPI-based COULWAVE'
		PRINT*,' '
		PRINT*,'   direct all comments and feedback to plynett@tamu.edu'
		PRINT*,' '
		PRINT*,'                      Contributors:'
		PRINT*,'      Khairil Irfan Sitanggang - MPI Parallelization'
		PRINT*,'       Dae-Hong Kim - Finite Volume/Riemanns Solver'
		PRINT*,'********************************************************'
		PRINT*,'********************************************************'     
		PRINT*,' '
		PRINT*,' ' 
		
	  if(batch.eq.2)goto 111   ! skip UI
	    
      if(batch.eq.0)then

		PRINT*,' '
		PRINT*,' Loading simulation parameters from data file, sim_set'
		print*,'      you will be able to change the numerical'
		print*,'          parameters through menu choices.'

      endif


      print*,'******* Parallel/Domain Decomposition Setup **********'
      choice=1
      do while(choice.ne.0)
		print*,'Number of Processors Selected: ', nprocs
		
		if(decomp_type.eq.0)then
			dims(1)=0
			dims(2)=0
			call mpi_dims_create(nprocs,2,dims,ierr )
			print*,'ID 1: Automatic division of grids used'
		else
			print*,'ID 1: Manual division of grids used'
		endif

		print*,'Number of divisions in x-dir: ', dims(1)
		print*,'Number of divisions in y-dir: ', dims(2)
		if(dims(1)*dims(2).ne.nprocs)then
			print*,'!!!! Number of sub-domains, ',dims(1)*dims(2),
     - 'not equal to number of procs, ',nprocs
		endif
          print*,'Are the above choices OK? -
     -    Enter ID# to change or 0 for OK'
          read*, choice
          print*,' '

	    if(choice.eq.1)then
			print*,'ID 0: Automatic domain division'
			print*,'ID 1: Input domain divisions manually'
			print*,'Input ID >>'
			read*, decomp_type

			if(decomp_type.eq.0)then
				dims(1)=0
				dims(2)=0
				call mpi_dims_create(nprocs,2,dims,ierr )
			else
				print*,'Total number of domains
     - should be equal to nprocs= ',nprocs
				print*,'Enter number of equal divisions in x-dir'
				read*, dims(1)
				dims(2)=nprocs/dims(1)
			endif
		endif
	enddo

      print*,'********* Current Simulation setup *********************'
      choice=1
      do while(choice.ne.0)
      
         if(sim_opt.eq.1) print*,'ID 1: Surface wave evolution'
         if(sim_opt.eq.2) print*,'ID 1: 
     -Wave generation by submarine slide'
      
         if(dim.eq.1) print*,'ID 2: 1D simulation'
         if(dim.eq.2) print*,'ID 2: 2D simulation'

         if(nonlin_ind.eq.0) print*,'ID 3: Linear Simulation'
         if(nonlin_ind.eq.1) print*,'ID 3: Weakly Nonlinear Simulation'
         if(nonlin_ind.eq.2) print*,'ID 3: Fully Nonlinear Simulation'

         if(disp_prop.eq.1) PRINT*,'ID 4: Arbitrary level approximation'
         if(disp_prop.eq.2) PRINT*,'ID 4: Depth averaged approximation'
         if(disp_prop.eq.3) PRINT*,'ID 4: Shallow water wave equations'

	   num_levels=max(num_levels,1)
         if(num_levels.eq.1) PRINT*,'ID 40: One-layer model'
         if(num_levels.eq.2) PRINT*,'ID 40: Two-layer model'
         if(num_levels.eq.3) PRINT*,'ID 40: Three-layer model'
         if(num_levels.eq.4) PRINT*,'ID 40: Four-layer model'

         if(rotationality.eq.0) PRINT*,'ID 41: Irrotational Model'
         if(rotationality.eq.1) PRINT*,'ID 41: Weakly Rotational Model'

	   if(numerical_scheme.eq.0) PRINT*,'ID 42: Finite Difference Solver'
	   if(numerical_scheme.eq.1)then
			PRINT*,'ID 42: Finite Volume Solver'
			if(limiter_on.eq.1)then
				PRINT*,'   ID 43: Using flux limiter'
			else
				PRINT*,'   ID 43: Not using flux limiter'
			endif
		endif			
         if(sim_opt.eq.1)then
         if(wave_type.eq.1.or.wave_type.eq.0) 
     -                      PRINT*,'ID 5: Solitary Wave Evolution'
         if(wave_type.eq.2) PRINT*,'ID 5: Solitary Wave with a/h=0.58'
         if(wave_type.eq.3) PRINT*,'ID 5: Solitary Wave with a/h=0.42'
         if(wave_type.eq.4) PRINT*,'ID 5: Cnoidal Wave Evolution'
         if(wave_type.eq.5) PRINT*,'ID 5: Sine Wave Evolution'
         if(wave_type.eq.6)then
           PRINT*,'ID 5: Wave Spectrum\Input Time Series Evolution'
           if(spec_type.eq.0)then
              PRINT*,'ID 6:    Spectrum Data Files'
              PRINT*,'      created using spectrum.m in Matlab'
           endif
           if(spec_type.eq.1)then
              PRINT*,'ID 6:    Spectrum Data Files'
              PRINT*,'      created using MANOLO spec2D.exe'
           endif
           if(spec_type.eq.2)then
              PRINT*,'ID 6:    Waves Generated using free'
              PRINT*,'    surface time series data in eta_in.dat'
           endif
         endif
         if(wave_type.eq.7) PRINT*,'ID 5: Submarine slump hot start'
         if(wave_type.eq.8)then		 
              PRINT*,'ID 8:    Left boundary wavemaker'
              PRINT*,'   created using create_wm_signal.m in Matlab'
			  left_wvmk=1
		 endif
         endif
         print*,' '
      print*,'Are the above choices OK? -
     - Enter ID# to change or 0 for OK'
         read*, choice
         print*,' '

         if(choice.eq.1)then
           PRINT*,'ID=1:Surface wave evolution'
           PRINT*,'ID=2:Wave generation by submarine landslide'
           PRINT*,'Input ID >>'
           read*, sim_opt
         elseif(choice.eq.2)then
           PRINT*,'ID=1: 1D simulation'
           PRINT*,'ID=2: 2D simulation'
           PRINT*,'Input ID >>'
           read*, dim
         elseif(choice.eq.3)then
      PRINT*,'ID=0: Linear Simulation'
      PRINT*,'ID=1: Weakly Nonlinear Simulation'
      PRINT*,'ID=2: Fully Nonlinear Simulation'
      PRINT*,'Input ID >>'
      read*, nonlin_ind
         if(nonlin_ind.eq.2)then
           disp_prop=1
         endif

      elseif(choice.eq.4)then
         if(nonlin_ind.eq.2)then
         PRINT*,'ID=1: Use arbitary level approximation'
         PRINT*,'ID=3: Use shallow water wave equations'
         PRINT*,'Input ID >>'
         read*, disp_prop
         else
         PRINT*,'ID=1: Use arbitrary level approximation'
         PRINT*,'ID=2: Use depth averaged approximation'
         PRINT*,'ID=3: Use shallow water wave equations'
         PRINT*,'Input ID >>'
         read*, disp_prop
         endif
      elseif(choice.eq.40)then
		PRINT*,'ID=1: Use one-layer model (Extended Boussinesq)'
		PRINT*,'ID=2: Use two-layer model'
c		PRINT*,'ID=3: Use three-layer model'
c		PRINT*,'ID=4: Use four-layer model'
		PRINT*,'Input ID >>'
		read*, num_levels
      elseif(choice.eq.41)then
		PRINT*,'ID=0: Use irrotational (potential) model'
		PRINT*,'ID=1: Use rotational (vertical) model'
		PRINT*,'Input ID >>'
		read*, rotationality
      elseif(choice.eq.42)then
		PRINT*,'ID=0: Use finite difference solver'
		PRINT*,'ID=1: Use finite volume solver'
		PRINT*,'   FV solver only works for:'
		PRINT*,'      *Fully Nonlinear Simulation'
		PRINT*,'      *One-layer model'
		PRINT*,'      *Weakly Rotational Model'
		PRINT*,'   Compared to FD, FV solver will be:'
		PRINT*,'      *More dissipative (numerical dissipation)'
		PRINT*,'      *More stable'
		PRINT*,'      *Computationally slower'
		PRINT*,'Input ID >>'
		read*, numerical_scheme
		if(numerical_scheme.eq.1)then
			nonlin_ind=2
			num_levels=1
			rotationality=1
		endif
      elseif(choice.eq.43)then
		PRINT*,'ID=0: Do not use limiter in FV scheme'
		PRINT*,'ID=1: Use limiter in FV scheme'
		PRINT*,'   Using limiter will introduce some:'
		PRINT*,'      numerical dissipation, particularly'
		PRINT*,'      near steep fronts.  Limiter will'
		PRINT*,'      also increase model stability.'
		PRINT*,'Input ID >>'
		read*, limiter_on
      elseif(choice.eq.5)then
         PRINT*,'************** Type of Input Wave *****************'
         PRINT*,'ID=1: Solitary Wave'
c         PRINT*,'ID=2: Solitary Wave with a/h=0.58'
c         PRINT*,'ID=3: Solitary Wave with a/h=0.42'
c         PRINT*,'ID=4: Cnoidal Waves'
         PRINT*,'ID=5: Sine Waves'
         PRINT*,'ID=6: Wave Spectrum\Input Time series'
c         PRINT*,'ID=7: Submarine Slump Hot Start'
         PRINT*,'ID=8: Left Boundary Wavemaker'
         PRINT*,'Input ID >>'
         read*, wave_type
         PRINT*,'******************************************************'
      elseif(choice.eq.6)then
         PRINT*,'ID 0: Load Spectrum Data Files'
         PRINT*,'      created using spectrum.m in Matlab'
c         PRINT*,'ID 1: Load Spectrum Data Files'
c         PRINT*,'      created using MANOLO spec2D.exe'
         PRINT*,'ID 2: Waves Generated using free'
         PRINT*,'      surface time series data in eta_in.dat'
         PRINT*,'Input ID >>'
         read*, spec_type
      endif

      enddo

      print*,' '
      print*,'****************************************************'
      print*,' '


      choice=1
      do while(choice.ne.0)

      if(sim_opt.eq.1)then
         if(wave_type.eq.1) PRINT*,'Solitary Wave Evolution'
         if(wave_type.eq.2) PRINT*,
     - 'Solitary Wave with a/h=0.58 Evolution'
         if(wave_type.eq.3) PRINT*,
     - 'Solitary Wave with a/h=0.42 Evolution'
         if(wave_type.eq.4) PRINT*,'Cnoidal Wave Evolution'
         if(wave_type.eq.5) PRINT*,'Sine Wave Evolution'
         if(wave_type.eq.6)then
           PRINT*,'Wave Spectrum Evolution'
           if(spec_type.eq.0)then
              PRINT*,'   Spectrum Data Files'
              PRINT*,'   created using spectrum.m in Matlab'
           endif
           if(spec_type.eq.1)then
              PRINT*,'   Spectrum Data Files'
              PRINT*,'   created using MANOLO spec2D.exe'
           endif
           if(spec_type.eq.2)then
              PRINT*,'   Waves Generated using free'
              PRINT*,'   surface time series data in eta_in.dat'
           endif
         endif
         if(wave_type.eq.7) PRINT*,'Submarine slump hot start'
         if(wave_type.eq.8)then		 
              PRINT*,'ID 8:    Left boundary wavemaker'
              PRINT*,'   created using create_wm_signal.m in Matlab'
		 endif
         if(wave_type.le.5)then
           PRINT*,'ID 1. Wave height:(m)  ',wave_hgt
         elseif(wave_type.eq.6)then
		 if(spec_type.eq.2)then
              PRINT*,'ID 1. Wave height mod factor (near 1):  ',wave_hgt
		 else
              PRINT*,'ID 1. Significant wave height (approx): ',wave_hgt				
		 endif
		 int_src=1
         elseif(wave_type.eq.7)then
           PRINT*,'ID 1. Hot start depression:(m)  ',wave_hgt
           PRINT*,'ID 111. Elevation to depression ratio  ',bf_ratio
         endif

         PRINT*,'ID 2. Initial/characteristic depth:(m)  ',depth

         if(offshore_cur.eq.0)then
           PRINT*,'ID 200. No current sent through a boundary'
         else
           PRINT*,'ID 200. Sending in current through boundary along'
           if(offshore_cur.eq.1) print*,'    x=0 (left boundary)'   
           if(offshore_cur.eq.2) print*,'    x=ENDX (right boundary)'
           if(offshore_cur.eq.3) print*,'    y=0 (bottom boundary)'
           if(offshore_cur.eq.4) print*,'    x=ENDY (top boundary)'
           print*,'Current time series must reside in current_in.dat' 
         endif             
         if(wave_type.le.3)then
           print*,'ID 3. Initial location of crest of soliton (m)  ',x0
         elseif(wave_type.eq.7)then
           print*,'ID 3. Loc of centerpoint of the intial cond (m)  ',x0
           print*,'ID 8. Lengthscale of initial condition (m)  ', L
         else
           if(int_src.eq.0.and.left_wvmk.eq.0)then
           print*,'ID 4. # of wavelengths in initial domain  ',dom_wave
         print*,'ID 5. Total # of wavelengths to be created  ',end_wave
           print*,'ID 6. Total # of wavelengths to be ramped  ',ramp
           endif

           if(wave_type.eq.4)then
              print*,'ID 7. The modulus for cnoidal waves  ',mk2
              if(int_src.eq.0)then
                print*,'ID 71. Not using internal source.'
              else
                print*,'ID 71. Using internal source.'    
			  int_src=2          
                print*,'ID 72. Location of source wavemaker: ',x0 
			  if(is_oreint.eq.0)then
				print*,'ID 73. Source line normal to x-axis' 
			  else
				print*,'ID 73. Source line normal to y-axis' 
			  endif		                 
              endif

                  mk=sqrt(mk2)      ! mk is the modulus of the elliptic integral of the first kind
                  call cnoidal_rf(0.,1.-mk2,1.,CK)
	            call cnoidal_rd(0.,1.-mk2,1.,CE)
                  CE=CK-(1./3.)*mk2*CE
                  yt=(CE/CK-1+mk2)/mk2
                  yt2=(CE/CK*(mk2-2.)+2.*(1-mk2))/(4.*mk2**2)
                  rp1=(-CE/CK*(4.*CE/CK+5.*mk2-8.)+2.*(1.-mk2)*
     -                        (mk2-2.))/(4.*mk2**2)
                  r2=(2.*CE/CK+3.*mk2-4.)/(4.*mk2)
      
                  L=sqrt(16.*depth**3./(3.*wave_hgt))*mk*CK
                  c=SQRT(9.80665*depth)*(1.0+wave_hgt/depth*(2.0/mk2-
     -                  3.0*CE/mk2/CK))

              print*,'Cnoidal wave wavelength based on '
			print*,'        modulus & wave height (m)  ', L
              print*,'Cnoidal wave period based on '
			print*,'        modulus & wave height (s)  ', L/c

           elseif(wave_type.eq.5)then
c              if(int_src.eq.0)then
c                print*,'ID 71. Not using internal source.'
c              else
                print*,'Using internal source.'              
                print*,'ID 72. Location of source wavemaker: ',x0 
			  if(is_oreint.eq.0)then
				print*,'ID 73. Source line normal to x-axis' 
			  else
				print*,'ID 73. Source line normal to y-axis' 
			  endif		                 
c              endif

              print*,'ID 8. Sine wave wavelength (m)  ', L
              print*,'      Wave period based on wavelength (s)' 
                if(current.eq.1.and.int_src.ne.1)then
              print*, 1/(Fr_cur*sqrt(9.81*depth)/L+
     - sqrt(9.81/(2*3.1415*L)*tanh(2.*3.14*depth/L)))
                else
              print*, sqrt(2*3.14*L/9.81*1./tanh(2.*3.14*depth/L))
                endif

c              if(int_src.eq.1)then
                print*,'ID 9. 2nd harmonic wave height (m)  ', 
     - wave_hgt_2
                print*,'ID 91. 2nd harmonic wavelength (m)  ', L_2
                if(L_2.gt.eps) then
                print*,'      Wave period based on wavelength (s)', 
     - sqrt(2*3.14*L_2/9.81*1./tanh(2.*3.14*depth/L_2))
                endif
c              endif
           elseif(wave_type.eq.6)then
              print*,'ID 72. Location of source wavemaker: ',x0 
			  if(is_oreint.eq.0)then
				print*,'ID 73. Source line normal to x-axis' 
			  else
				print*,'ID 73. Source line normal to y-axis' 
			  endif				     
           endif
              
         endif
           
         if(dim.eq.2.and.wave_type.le.5)then
           print*,'ID 10. Incident angle of waves  ',inc_ang
         endif

      elseif(sim_opt.eq.2)then
         if(slide_type.eq.1) PRINT*,'Rotational/Translational Slide'
         if(slide_type.eq.2) PRINT*,'Sine Wave Transformation'
         if(slide_type.eq.3) PRINT*,'Hammocks bed thrust experiments'
         if(slide_type.eq.4) PRINT*,'Watts moving block experiments'
         if(slide_type.eq.6) PRINT*,'Elliptic slice slide'
         if(slide_type.eq.7) PRINT*,'Moving sin**2 disturbance'

         if(slide_type.le.2)then
           PRINT*,'ID 11. Depth at midpoint of the slide slope  ',depth
            PRINT*,'ID 12. Maximum change in depth due to the slide  ',
     - wave_hgt
            PRINT*,'ID 13. Slide duration: (in seconds)  ' , per_sld

           if(load_topo.eq.2)then
              print*,'ID 14. Rotate bathymetry 90 degrees, yes=1',rotate
              print*,'ID 15. Smooth between two y-transects, 
     - yes=1 (set=1 unless testing)', smooth
           endif

           if(smooth.eq.1)then
              print*,'ID 16. Y locations (m) to smooth between',y1,y2

              print*,'ID 17. Depth at top of slide slope (m)',depth_top

              print*,'ID 18. Depth at bottom of slide 
     - slope (m)',depth_bottom

              print*,'ID 19. Distance towards shoreline to 
     - smooth bathymetry (m)',smooth_top
           endif
      
           if(slide_type.eq.1)then
              print*, 'ID 20. Length of lip on top of slide (m):',
     - amp_top

              print*, 'ID 21. Radius coefficieny for 
     - rotational slide',radius_coef
           endif

         elseif(slide_type.eq.3)then
           print*,'NEED TO GET HAMMOCK PAPER'
         elseif(slide_type.eq.4)then
           print*,'NEED TO GET WATTS PAPER'
         elseif(slide_type.eq.6)then

c           load_topo=1
c           num_nodes=2
c           x_1=0
c           h_1=-0.25
c           x_2=10.
c           h_2=x_2*tan(slope_ang*3.1415/180)+h_1
c           spng_2=1

           PRINT*,'ID 22. Depth of water above the' 
           print*,'midpoint of the slide mass:', depth
           
            PRINT*,'ID 23. Maximum hieght of the slide mass'  
           print*,'normal to the slide slope:', wave_hgt

            PRINT*,'ID 24. Inclination angle of the slope
     - (degrees)',slope_ang

            PRINT*,'ID 25. Density ratio between the slide mass 
     - and water',gamma

            PRINT*,'ID 26. Drag coefficient', Cd

            PRINT*,'ID 27. Added-mass coefficient', Cm

		  if(aspect_ratio.lt.1e-6) aspect_ratio=1.

            if(dim.eq.2) PRINT*,'ID 32. Slide L/W', aspect_ratio

         elseif(slide_type.eq.7)then

           PRINT*,'ID 28. Constant water depth:', depth
           
            PRINT*,'ID 29. Maximum hieght of the mass:', wave_hgt

            PRINT*,'ID 30. Froude number of the disturbance:', Frd
         
           PRINT*,'ID 31. Length of the disturbance (in depths):', L

         endif

      endif

      print*,' '
      print*,'Are the above choices OK? -
     - Enter ID# to change or 0 for OK'
      read*, choice
      print*,' '

      if(choice.eq.1)then
         if(wave_type.le.5)then
           PRINT*,'Input wave height:(in meters)'      
           read*, wave_hgt
	   elseif(wave_type.eq.6)then
		 if(spec_type.eq.2)then
			PRINT*,'Input wave height mod factor (near 1):  '    
			read*, wave_hgt
		 else
			PRINT*,'Input significant wave height (approx):  '    
			read*, wave_hgt			
		 endif
         elseif(wave_type.eq.7)then
           PRINT*,'Input hot start depression:(m)  '
           read*, wave_hgt
           PRINT*,'Input elevation to depression ratio  '
           read*, bf_ratio
         endif
      elseif(choice.eq.2)then
         PRINT*,'Input initial/characteristic water depth:(in meters)'   
         read*, depth
      elseif(choice.eq.200)then
         PRINT*,'Input current through a lateral boundary?'
		 PRINT*,'ID=0. No input current'     
		 PRINT*,'ID=1. Input current along x=0 (left) boundary'  
		 PRINT*,'ID=2. Input current along x=ENDX (right) boundary'  
		 PRINT*,'ID=3. Input current along y=0 (bottom) boundary'  
		 PRINT*,'ID=4. Input current along y=ENDY (top) boundary'  
         PRINT*,'Input ID >>'       
         read*, offshore_cur
      elseif(choice.eq.3)then
         if(wave_type.le.3)then
           print*,'Input initial location of crest of soliton (meters)'
           read*, x0
         elseif(wave_type.eq.7)then
           print*,'Input loc of centerpoint of the intial cond (m)  '
           read*, x0
         endif
      elseif(choice.eq.4)then
         print*,'Enter number of wavelengths to be in initial domain'
         print*,'      - must be at least one'
         read*, dom_wave
      elseif(choice.eq.5)then
         print*,'Enter total number of wavelengths to be created'
         print*,'during entire simulation'
         read*, end_wave
      elseif(choice.eq.6)then
         print*,'Enter total number of wavelengths to be ramped to '
         print*,'eliminate free surface discontinuity - 1 recommended'
         read*, ramp
      elseif(choice.eq.7)then
         print*,'For cnoidal waves - enter the modulus of the'
         print*,'elliptic integral of the first kind'
         read*, mk2
      elseif(choice.eq.71)then
         print*,'Sine waves can be effectively modeled using'
         print*,'an internal source wavemaker.'
         PRINT*,'ID=0: Do not use internal source wavemaker'
         PRINT*,'ID=1: Use internal source wavemaker'
         PRINT*,'Input ID >>'
         read*, int_src
      elseif(choice.eq.72)then
         print*,'Input location of internal source (meters)-'
         if(inc_ang.ne.0)then
           print*,'This is the normal distance from the'
           print*,' origin to the source location.'
         endif
         read*, x0
	elseif(choice.eq.73)then
		print*,'Source orientation: ID=0  normal to x-axis'
		print*,'                    ID=1  normal to y-axis'
		print*,'Input ID >>'
		read*,is_oreint
      elseif(choice.eq.8)then
         if(wave_type.eq.5.or.wave_type.eq.4)then
           print*,'For sine/cnoidal waves - enter the wavelength (m)'
           read*, L
         elseif(wave_type.eq.7)then
           print*,'Input lengthscale of initial condition (m)  '
           read*, L
         endif
      elseif(choice.eq.9)then
         print*,'Enter second harmonic wave height:'
         read*, wave_hgt_2
      elseif(choice.eq.91)then
         print*,'Enter second harmonic wavelength:'
         read*, L_2
      elseif(choice.eq.10)then
         print*,'Enter incident angle of waves - 0 = waves traveling'
         print*,'in positive x-direction (degrees)'
         read*,inc_ang
      elseif(choice.eq.11)then      
         PRINT*,'Input the approximate depth at the' 
         print*,'midpoint of the slide slope:(in meters)'      
         read*, depth
      elseif(choice.eq.12)then
          PRINT*,'Input the maximum change in water depth'  
         print*,'due to the submarine slide:(in meters)'      
         read*, wave_hgt
      elseif(choice.eq.13)then
          PRINT*,'Input the slide duration: (in seconds)'      
         read*, per_sld
      elseif(choice.eq.14)then
         print*,'Rotate bathymetry 90 degrees, yes=1'
         read*, rotate
      elseif(choice.eq.15)then
         print*,'Smooth between two y-transects, 
     - yes=1 (set=1 unless testing)'
         read*, smooth
      elseif(choice.eq.16)then
         print*,'Y locations (m) to smooth between'
         read*, y1
         read*, y2
      elseif(choice.eq.17)then
         print*,'Depth at top of slide slope (m)'
         read*, depth_top
      elseif(choice.eq.18)then
         print*,'Depth at bottom of slide slope (m)'
         read*, depth_bottom
       elseif(choice.eq.19)then
         print*,'Distance towards shoreline to smooth bathymetry (m)'
         read*, smooth_top
      elseif(choice.eq.20)then
         print*, 'Length of lip on top of slide (m):'
         read*,amp_top
      elseif(choice.eq.21)then
         print*, 'Radius coefficient for rotational slide'
         read*,radius_coef
      elseif(choice.eq.22)then
         PRINT*,'Input the depth of water above the' 
         print*,'      midpoint of the slide mass:'      
         read*, depth
      elseif(choice.eq.23)then           
          PRINT*,'Input the maximum vertical hieght of the slide mass:'  
         read*, wave_hgt
      elseif(choice.eq.24)then
          PRINT*,'Input inclination angle of the slope (degrees)'  
         read*, slope_ang

         load_topo=1
         num_nodes=2
         x_1=0
         h_1=-1.
         x_2=5.
         h_2=x_2*tan(slope_ang*3.1415/180)+h_1

      elseif(choice.eq.25)then
          PRINT*,'Input the density ratio between
     - the slide mass and water:'
         read*, gamma
      elseif(choice.eq.26)then
          PRINT*,'Input the drag coefficient'
         read*, Cd
      elseif(choice.eq.27)then
          PRINT*,'Input the added-mass coefficient'
         read*, Cm
      elseif(choice.eq.32)then
          PRINT*,'Input the slide aspect ratio, L/W'
         read*, aspect_ratio
      elseif(choice.eq.28)then
         PRINT*,'Input the constant water depth:'      
         read*, depth
           x_1=0
           h_1=depth
           x_2=90
           h_2=depth
      elseif(choice.eq.29)then           
          PRINT*,'Input the maximum vertical hieght of the disturbance:'  
         read*, wave_hgt
      elseif(choice.eq.30)then
          PRINT*,'Input the Froude number of the disturbance'  
         read*, FRD
      elseif(choice.eq.31)then
         PRINT*,'Input the length of the disturbance (in depths):'      
         read*, L
      endif      
      
      enddo


      print*,' '
      print*,'****************************************************'
      print*,' '


       choice=1
      do while(choice.ne.0)

      if(load_topo.eq.1) PRINT*,'ID 99. Bottom Profile using 
     - Location/Depth Nodes'
      if(load_topo.eq.2) print*,'ID 99.
     - Load Topographical Data from Files'

      if(load_topo.eq.1)then
         print*,'ID 1. # of nodes used to create depth profile  ',
     - num_nodes
         
         if(num_nodes.ge.2)then
           print*,'      X(m)      H(m) '
           print*,'ID 2. ',x_1,h_1
           print*,'ID 3. ',x_2,h_2
         endif

         if(num_nodes.ge.3)print*,'ID 4. ',x_3,h_3
         if(num_nodes.ge.4)print*,'ID 5. ',x_4,h_4
         if(num_nodes.ge.5)print*,'ID 6. ',x_5,h_5
         if(num_nodes.ge.6)print*,'ID 7. ',x_6,h_6
         if(num_nodes.ge.7)print*,'ID 8. ',x_7,h_7
         if(num_nodes.ge.8)print*,'ID 9. ',x_8,h_8
         if(num_nodes.ge.9)print*,'ID 10. ',x_9,h_9

         if(dim.eq.2)then
           print*,'ID 11. Domain width (y direction) (m): ',chan_width
         endif

         if(sim_opt.eq.2) then 
           print*,'ID 12. Node index of top of slide slope',slide_node
         endif
      elseif(load_topo.eq.2.and.sim_opt.eq.2)then
         print*,'ID 13. Approx length of the slide slope',slope_reg
      endif


      print*,' '
      print*,'Are the above choices OK? -
     - Enter ID# to change or 0 for OK'
      read*, choice
      print*,' '

      if(choice.eq.99)then
         PRINT*,'ID=1: Specify Profile by Giving Location/Depth Nodes'
        PRINT*,'ID=2: Load Topographical Data from Matlab Created Files'
         PRINT*,'Input ID >>'
         read*, load_topo
      elseif(choice.eq.1)then
         print*,'Enter the number of nodes required to create' 
         print*,'bottom profile - including first and last
     - (min 2, max 9):'
         read*, num_nodes
         
         if(num_nodes.ge.2)then
           print*,'First x-node is x=0 meters'
           x_1=0
           print*,'Enter depth of first node (m):'
           read*, h_1

           print*,'Enter x-location of second node (m):'
           read*,x_2
           print*,'Enter depth of second node (m):'
           read*, h_2
         endif

         if(num_nodes.ge.3)then
           print*,'Enter x-location of third node (m):'
           read*,x_3
           print*,'Enter depth of third node (m):'
           read*, h_3
         endif

         if(num_nodes.ge.4)then
           print*,'Enter x-location of fourth node (m):'
           read*,x_4
           print*,'Enter depth of fourth node (m):'
           read*, h_4
         endif

         if(num_nodes.ge.5)then
           print*,'Enter x-location of fifth node (m):'
           read*,x_5
           print*,'Enter depth of fifth node (m):'
           read*, h_5
         endif

         if(num_nodes.ge.6)then
           print*,'Enter x-location of sixth node (m):'
           read*,x_6
           print*,'Enter depth of sixth node (m):'
           read*, h_6
         endif

         if(num_nodes.ge.7)then
           print*,'Enter x-location of seventh node (m):'
           read*,x_7
           print*,'Enter depth of seventh node (m):'
           read*, h_7
         endif

         if(num_nodes.ge.8)then
           print*,'Enter x-location of eighth node (m):'
           read*,x_8
           print*,'Enter depth of eighth node (m):'
           read*, h_8
         endif

         if(num_nodes.ge.9)then
           print*,'Enter x-location of ninth node (m):'
           read*,x_9
           print*,'Enter depth of ninth node (m):'
           read*, h_9
         endif

      elseif(choice.eq.2)then
           print*,'First x-node is x=0 meters'
           x_1=0
           print*,'Enter depth of first node (m):'
           read*, h_1

           if(slide_type.eq.6)h_2=x_2*tan(slope_ang*3.1415/180)+h_1
      elseif(choice.eq.3)then
           print*,'Enter x-location of second node (m):'
           read*,x_2

      
           if(slide_type.eq.6)then 
              h_2=x_2*tan(slope_ang*3.1415/180)+h_1
           else
              print*,'Enter depth of second node (m):'
              read*, h_2
           endif
      elseif(choice.eq.4)then
           print*,'Enter x-location of third node (m):'
           read*,x_3
           print*,'Enter depth of third node (m):'
           read*, h_3
      elseif(choice.eq.5)then
           print*,'Enter x-location of fourth node (m):'
           read*,x_4
           print*,'Enter depth of fourth node (m):'
           read*, h_4
      elseif(choice.eq.6)then
           print*,'Enter x-location of fifth node (m):'
           read*,x_5
           print*,'Enter depth of fifth node (m):'
           read*, h_5
      elseif(choice.eq.7)then
           print*,'Enter x-location of sixth node (m):'
           read*,x_6
           print*,'Enter depth of sixth node (m):'
           read*, h_6
      elseif(choice.eq.8)then
           print*,'Enter x-location of seventh node (m):'
           read*,x_7
           print*,'Enter depth of seventh node (m):'
           read*, h_7
      elseif(choice.eq.9)then
           print*,'Enter x-location of eighth node (m):'
           read*,x_8
           print*,'Enter depth of eighth node (m):'
           read*, h_8
      elseif(choice.eq.10)then
           print*,'Enter x-location of ninth node (m):'
           read*,x_9
           print*,'Enter depth of ninth node (m):'
           read*, h_9
      elseif(choice.eq.11)then      
         print*,'Enter domain width (y direction) (m):'
         read*,chan_width
      elseif(choice.eq.12)then      
         print*,'Enter the node index of the top of the slide slope'
         read*,slide_node
      elseif(choice.eq.13)then
         print*,'Enter the approx length of the slide slope'
         read*,slope_reg
      endif      
      
      enddo


      print*,' '
      print*,'****************************************************'
      print*,' '


       choice=1
      do while(choice.ne.0)

      print*,'ID 1. Simulation time in seconds  ',end_t
      print*,'ID 2. Time increment to write to file (s)  ',writ_inc

      print*,'ID 3. Number of grid points per wavelength  ',pts_wvl
      print*,'ID 4. Courant number = dx/dt/c_o  ',courant

c      PRINT*,'************** Boundary Conditions ******************'
c      PRINT*,'1 = Solid-Reflective Wall'
c      PRINT*,'2 = Input Wave - sending a wave through boundary'
c      PRINT*,'4 = Period - left/right or top/bottom must both use'
c      if(dim.eq.1)then
c         print*,'      Left Wall    Right Wall '
c         print*,'      ID 5      ID 6    '
c         print*, bc_1,bc_2
c      elseif(dim.eq.2)then
c      print*,'      Left Wall    Right Wall    Top Wall     Bottom Wall'
c      print*,'      ID 5      ID 6      ID 7      ID 8   '
c         print*, bc_1,bc_2,bc_4,bc_3
c      endif

      PRINT*,'************** Sponge Layer Absorbers ******************'
      PRINT*,'0 = Do Not Use Sponge Layer'
      PRINT*,'1 = Use Sponge Layer'
      if(dim.eq.1)then
         print*,'      Left Wall    Right Wall '
         print*,'      ID 9      ID 10    '
         print*, spng_1,spng_2
      elseif(dim.eq.2)then
      print*,'      Left Wall    Right Wall    Top Wall     Bottom Wall'
      print*,'      ID 9      ID 10      ID 11      ID 12   '
         print*, spng_1,spng_2,spng_4,spng_3
      endif


      if(screen_output.eq.1)then
         PRINT*,'ID 16. Display screen data at time step interval ',
     - display
      else
         PRINT*,'ID 16. Do not display screen data '
      endif

      print*,' '
      print*,'Are the above choices OK? -
     - Enter ID# to change or 0 for OK'
      read*, choice
      print*,' '

      if(choice.eq.1)then
         print*,'Enter the physical simulation time in seconds'
         read*,end_t
      elseif(choice.eq.2)then      
         print*,'Enter the time increment at which spatial snapshots'
         print*,'of free surface should be written to file (s)'
         read*,writ_inc
      elseif(choice.eq.3)then      
         print*,'Enter the number of grid points per wavelength'
         read*,pts_wvl
      elseif(choice.eq.4)then      
         print*,'Enter the courant number = dx/dt/c_o'
         read*,courant
      elseif(choice.eq.5)then      
         PRINT*,'Input ID for left wall >>'
         if(wave_type.ge.4) print*, 'Should be 2 as cnoidal is used'
         read*, bc_1
      elseif(choice.eq.6)then      
         PRINT*,'Input ID for right wall >>'
         read*, bc_2
      elseif(choice.eq.7)then      
         PRINT*,'Input ID for top wall >>'
         read*, bc_4
      elseif(choice.eq.8)then      
         PRINT*,'Input ID for bottom wall >>'
         read*, bc_3
      elseif(choice.eq.9)then      
         PRINT*,'Input ID for left wall >>'
         read*, spng_1
      elseif(choice.eq.10)then      
         PRINT*,'Input ID for right wall >>'
         read*, spng_2
      elseif(choice.eq.11)then      
         PRINT*,'Input ID for top wall >>'
         read*, spng_4
      elseif(choice.eq.12)then      
         PRINT*,'Input ID for bottom wall >>'
         read*, spng_3
      elseif(choice.eq.16)then
         print*,'Time step interval to print
     - select information to screen'
         PRINT*,'ID=1: Display output data to screen'
         PRINT*,'ID=2: Do not display - use for batch jobs'
         PRINT*,'Input ID >>'
         read*, screen_output

         if(screen_output.eq.1)then
           print*,'Time step interval to print select info to screen'
           read*, display
         else
           display=1
         endif
      endif

      enddo


      print*,' '
      print*,'*************** TIME SERIES OUTPUT ****************'
      print*,' '

      choice=1
      do while(choice.ne.0)

      print*,'ID 1. Number of time series to write to file',num_ts
	if(num_ts.gt.50)then
	 print*,'Number of time series large, edit ts_locations.dat to change'
	else
       do i=1,num_ts
         if(i.eq.1) print*,'          X     ','    Y     '
         print*,'ID ',i+1,' ',x_ts(i),y_ts(i)
       enddo
	endif
	if(spec_pp.eq.1)then
	print*,'ID 1000. Performing spectral analysis on recorded time series'
	print*,'ID 1001. Time at which to start processing time series:',
     - spec_ts
	else
		spec_pp=0
	print*,'ID 1000. Not performing spectral analysis on time series'
	endif

      print*,' '
      print*,'Are the above choices OK? -
     - Enter ID# to change or 0 for OK'
      read*, choice
      print*,' '

      if(choice.eq.1)then
         print*,'Write time series of z,u,v at
     - specific locations to file?'
         print*,'Enter 0 for no, or enter the number of locations at'
         print*,'which you would like time series for: (max 1000) '
         read*,num_ts
      elseif(choice.ge.2.and.choice.lt.1000)then
         print*,'Spatial location for ID',choice
         print*,'Enter the x-location of time series:'
         read*,x_ts(choice-1)
         if(dim.eq.1)then
           y_ts(choice-1)=0
         else
           print*,'Enter the y-location of time series:'
           read*,y_ts(choice-1)
         endif
	elseif(choice.eq.1000)then
         print*,'Perform spectral analysis on recorded time series?'
         print*,'Enter 0 for no, or 1 for yes'
         read*,spec_pp
	elseif(choice.eq.1001)then
         print*,'Enter time at which to start processing time series:'
         read*,spec_ts		   		
      endif                  
      enddo


      print*,' '
      print*,'*************** PARAMETERIZATIONS ****************'
      print*,' '

      choice=1
      do while(choice.ne.0)

      if(wave_breaking.eq.0)then
         print*,'ID 1. Wave Breaking Model not implemented.'
      elseif(wave_breaking.eq.1)then
         print*,'ID 1. Wave Breaking Model implemented.'
         if(breaker_type.eq.0)then
			print*,'ID 11. Using transport-based breaking model.'
         elseif(breaker_type.eq.1)then
			print*,'ID 11. Using Kennedy et al type breaking model.'
		 endif
      endif
	
	if(numerical_scheme.eq.0)then 
      print*,'ID 2. Fraction of upwinded differences composing'
	print*,'      convective terms: ', upwind_baseline
	endif

      if(bottom_fric.eq.0)then
         print*,'ID 3. Bottom friction ignored.'
      elseif(bottom_fric.eq.1)then
         print*,'ID 3. Bottom friction included.'
	   if(bf_type.eq.0)then 
         print*,'ID 31. Using roughness height model, with friction'
	     print*,'    factor determined via an approx Moody diagram.'
         print*,'ID 4. Roughness height, ks (m): ', f_BF
       elseif(bf_type.eq.1)then
         print*,'ID 31. Using Mannings roughness parameter.'
         print*,'ID 4. Mannings n (s/m^1/3): ', f_BF
	   elseif(bf_type.eq.2)then 
         print*,'ID 31. Using spatially constant friction factor.'
         print*,'ID 4. Constant friction factor, f: ', f_BF
	   endif
	   print*,'ID 41. Coefficient for vertical eddy'
	   print*,'       viscosity (Elders model, 0.0667): ',Ch
	   print*,'ID 5. Coefficient for subgrid horizontal eddy'
	   print*,'       viscosity (Smag form, 0.08-0.2): ',visc_coef
	   if(backscatter.eq.0)then
		print*,'ID 6. Not using backscatter model.'
	   elseif(backscatter.eq.1)then
		print*,'ID 6. Using backscatter model.'
	    print*,'ID 7. Backscatter coefficient (default:70 in 2HD): ', CB
	   endif
	   if(ihvor.eq.0)then
		print*,'ID 8. Ignore high-order viscous flux terms.'
	   elseif(ihvor.eq.1)then
		print*,'ID 8. Including high-order viscous flux terms.'
		print*,'ID 81. Minimum water depth required for local'
		print*,'       use of high-order viscous terms (m): ', Ch_length
	   endif	   
	endif
      print*,' '
      print*,'Are the above choices OK? -
     - Enter ID# to change or 0 for OK'
      read*, choice
      print*,' '

      if(choice.eq.1)then      
         PRINT*,'An eddy viscosity type wave breaking parameterization'
         PRINT*,'can be utilized.'
         PRINT*,'ID=0: Do Not Use Wave Breaking Model'
         PRINT*,'ID=1: Use Wave Breaking Model'
         PRINT*,'Input ID >>'
         read*, wave_breaking
      elseif(choice.eq.11)then      
         print*,' Choose between a new (unpublished) transport-based'
	   print*,' breaking model (ID=0) or the traditional Kennedy et al'
	   print*,' expliciy eddy viscosity model (ID=1).'
         read*, breaker_type
      elseif(choice.eq.2)then      
         print*,' Input fraction of upwinded differences composing'
	   print*,' convective terms, must be between 0 (no upwinding)'
	   print*,' and 1 (full upwinding).'
         read*, upwind_baseline
      elseif(choice.eq.3)then      
         PRINT*,'Bottom friction effects can be included, modeled'
         PRINT*,'using the quadratice law.'
         PRINT*,'ID=0: Do Not Include Bottom Friction'
         PRINT*,'ID=1: Include Bottom Friction'
         PRINT*,'Input ID >>'
         read*, bottom_fric
      elseif(choice.eq.31)then      
         PRINT*,'Bottom friction effects can be modeled, using:'
         print*,'ID 0. Roughness height model, with friction factor'
	     print*,'       determined via a digitized approx Moody diagram.'
         print*,'ID 1. Mannings roughness parameter.'
         print*,'ID 2. Spatially constant friction factor.'
		 print*,'NOTE: If using a moving shoreline with energetic or'
		 print*,'breaking waves, use model ID=0 for stable results.' 
         PRINT*,'Input ID >>'
         read*, bf_type
      elseif(choice.eq.4)then
	   if(bf_type.eq.1)then      
         PRINT*,'A Mannings n must be inputted. This value'
         PRINT*,'is typically in the 10-3 to
     - 10-1 range. Enter coefficent:'
	   elseif(bf_type.eq.0)then
         PRINT*,'A characteristic roughness height must be'
         PRINT*,'inputted. Enter height (m):'
	   elseif(bf_type.eq.2)then
         PRINT*,'A constant friction factor must be'
         PRINT*,'inputted. Enter value:'
	   endif
         read*, f_BF
      elseif(choice.eq.41)then      
         print*,' Input coefficient for vertical eddy'
	    print*,' viscosity (Elders model, default: 0.0667):'
         read*, Ch
      elseif(choice.eq.5)then      
         print*,' Input coefficient for subgrid eddy viscosity'
	     print*,' (should be small < 1, default: 0.2)'
         read*, visc_coef
      elseif(choice.eq.6)then      
         PRINT*,'An empirical backscatter model (BSM) can be inlcuded.'
         PRINT*,'The BSM will provide reasonable velocity flucuation'
         PRINT*,'predictions, but at the cost of increased simulation'
         PRINT*,'time and decreased simulation stability.'
         PRINT*,'ID=0: Do Not Use BSM'
         PRINT*,'ID=1: Use BSM'
         PRINT*,'Input ID >>'
         read*, backscatter   
      elseif(choice.eq.7)then      
         print*,' Input coefficient for backscatter model, '
	     print*,' (default:70 in 2HD). '
         read*, CB  
      elseif(choice.eq.8)then      
         PRINT*,'High order terms due to viscous effects can be'
         PRINT*,'included, following Kim et al. (2009).'
         PRINT*,'ID=0: Do Not Include Viscous Vorticty terms'
         PRINT*,'ID=1: Include Viscous Vorticty terms'
         PRINT*,'Input ID >>'
         read*, ihvor   
      elseif(choice.eq.81)then      
         print*,'The high-order viscous terms are often unstable'
         PRINT*,'in very shallow water depths, due to their '
         PRINT*,'proportionality to 1/H^2.  For stable simulations'
         PRINT*,'with these terms included, it is ususally required'
         PRINT*,'to specify a "turn-off" total water depth.  If the'
	     print*,'total water depth is less than this value at a '
	     print*,'particular grid point at a particular time, these '
         PRINT*,'high-order viscous terms will not be calculated at '
         PRINT*,'that grid location at that time. A value of zero '
         PRINT*,'implies that the high-order viscous terms are '
         PRINT*,'included at all locations at all times.'
         PRINT*,' Input minimum water depth required for local'
		 print*,'       use of high-order viscous terms (m): '
         read*, Ch_length
      endif
      enddo

      print*,' '
      print*,'****************************************************'
      print*,' '


      choice=0 !1 skip filtering options
      do while(choice.ne.0)

c      if(filt.eq.1)then
c         print*,'ID 1. Filter entire domain ',filt_int,'
c     - times per period'
c      else
c         print*,'ID 1. Do not filter entire domain'
c      endif
c      
c      if(sim_opt.eq.2)then
c      print*,'ID 2. Filter over slide ',filtsld_int,' times per period'
c      endif
c
c      print*,' '
c      print*,'Are the above choices OK? - Enter ID# to change or
c     - 0 for OK'
c      read*, choice
c      print*,' '

      if(choice.eq.1)then      
         Print*,'1 = Filter entire domain using 9-point filter'
         Print*,'2 = Do not filter entire domain'
         PRINT*,'Input  >>'
         read*, filt

         if(filt.eq.1)then
           print*,'Number of times to filter entire domain
     - per wave period: '
           read*,filt_int
         endif
      elseif(choice.eq.2)then      
         print*,'Number of times to filter over slide region
     - per slide period'
         read*,filtsld_int
      endif
      enddo


      print*,' '
      print*,'****************** DEFAULT VALUES *******************'
      print*,' '


      choice=1
      do while(choice.ne.0)
      
c      print*,'ID 1. Finite difference order of dispersive terms',
c     -      deriv_order_ind
      Print*,'ID 2. Width of sponge layer, in wavelengths',sponge_width
      print*,'ID 3. Corrector stage convergence error',ep
c      print*,'ID 4. Second corrector convergence error',ep2
      Print*,'ID 5. Max # of allowable iterations in corrector loop',itr
      print*,'ID 6. Min # of iterations in corrector loop',min_itr
c      print*,'ID 7. Iteration at which to apply over-relaxation ',itr_or
c      print*,'ID 8. Relaxtion coeficient',o
      if(smooth_bathy.eq.1)then
      print*,'ID 9. Smoothing depth profile using 4-point filter.'
      else
      print*,'ID 9. Not smoothing depth profile using 4-point filter.'
      endif
c      if(use_av.eq.1)then
c      print*,'ID 10. Use Visual Fortran Array Viewer to examine free '
c      print*,'      surface as the program runs (Windows w/Array '
c      print*,'      Viewer only)'
c      else
c      print*,'ID 10. Do not use Visual Fortran Array Viewer'
c      endif
      print*,'ID 11. First sponge layer coef ',cdamp1
      print*,'ID 12. Second sponge layer coef',cdamp2
      print*,'ID 13. Shoreline can move? 0=Yes, 1=No', sh_mov
	  if(swash_d.lt.1e-10) swash_d=0.01
      print*,'ID 14. Min water depth for moving boundary (m)', swash_d
	  if(upwind_shore.eq.0)then  
        print*,'ID 15. Use upwinding for flux terms at first wet cell' 
	  else
        print*,'ID 15. Use linear averaging for flux terms at first wet' 
	  endif

c	print*,'       Approx min allowable depth (m)= ', depth/cutoff
      print*,' '
      print*,'Are the above choices OK?-Enter ID# to change or 0 for OK'
      read*, choice
      print*,' '

      if(choice.eq.1)then      
         print*,'Enter finite difference order of dispersive terms:'
         read*,deriv_order_ind
      elseif(choice.eq.2)then      
         print*,'Enter width of sponge layer, in wavelengths:'
         read*,sponge_width
      elseif(choice.eq.3)then      
         print*,'Enter corrector stage convergence error:'
         read*,ep
         ep2=ep
      elseif(choice.eq.4)then      
         print*,'Enter second corrector convergence error:'
         read*,ep2
      elseif(choice.eq.5)then      
         print*,'Enter max # of allowable iterations in corrector loop'
         read*,itr
      elseif(choice.eq.6)then      
         print*,'Enter min # of iterations in corrector loop'
         read*,min_itr
      elseif(choice.eq.7)then      
         print*,'Enter iteration at which to apply over-relaxation'
         read*,itr_or
      elseif(choice.eq.8)then      
         print*,'Enter relaxtion coeficient'
         read*,o
      elseif(choice.eq.9)then
         print*,'0 = Do not smooth depth profile using 4-point filter'
         print*,'1 = Smooth depth profile using 4-point filter'
         PRINT*,'Input  >>'
         read*, smooth_bathy
      elseif(choice.eq.10)then
         print*,'0 = Do not use Visual Fortran Array Viewer'
         print*,'1 = Use Visual Fortran Array Viewer to examine free '
         print*,'      surface as the program runs (Windows w/Array'
         print*,'      Viewer only)'
         PRINT*,'Input  >>'
         read*, use_av
      elseif(choice.eq.11)then      
         print*,'Enter first sponge layer coef'
         read*,cdamp1
      elseif(choice.eq.12)then      
         print*,'Enter second sponge layer coef'
         read*,cdamp2
      elseif(choice.eq.13)then      
         print*,'Shoreline can move? 0=Yes, 1=No'
         read*,sh_mov
      elseif(choice.eq.14)then      
         print*,'Enter min allowable water depth.  This is'
         print*,'      the vertical accuracy of the shoreline'
         print*,'      (runup) motion.'
         read*,swash_d
      elseif(choice.eq.15)then      
         print*,'0 = Use upwinding for flux terms at first wet cell'
         print*,'1 = Use linear averaging for flux at first wet cll'
         read*,upwind_shore
      endif
      enddo

111	continue

C Write all simulation parameters back to inputs files
            if(dim.eq.1)then
                  inc_ang=0

                  spng_3=0
                  spng_4=0
                  bc_3=1
                  bc_4=1
            endif

            if(wave_type.lt.4) int_src=0
            if(wave_type.ge.5.and.int_src.ge.1) int_src=1
            if(wave_type.eq.6) int_src=1
            if(wave_type.eq.6.and.spec_type.eq.2)  int_src=2

			
			plus_tide=0

	call write_inputs

      return

      end

