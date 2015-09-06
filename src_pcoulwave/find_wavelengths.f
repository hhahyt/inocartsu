      subroutine find_wavelengths

      use mainvar_module


	if(wave_type.eq.5.or.wave_type.eq.4.or.wave_type.eq.1)then
		freq_file=2
		dir_file=2				
				
		ALLOCATE(amp_spec(freq_file,dir_file),
     -				per_spec(freq_file,dir_file),
     -				L_spec(freq_file,dir_file),
     -				beta_src(freq_file,dir_file),
     -				I_src(freq_file,dir_file),
     -				D_src(freq_file,dir_file),
     -				shift_spec(freq_file,dir_file),
     -				inc_ang_spec(freq_file,dir_file),
     -				cosA(freq_file,dir_file),
     -				sinA(freq_file,dir_file),
     -				D_src_nl(freq_file,dir_file))	

      elseif(wave_type.eq.6.or.wave_type.eq.8)then
		if(spec_type.eq.2.or.wave_type.eq.8)then
 
 			freq_file=2
			dir_file=2				
				
			ALLOCATE(amp_spec(freq_file,dir_file),
     -				per_spec(freq_file,dir_file),
     -				L_spec(freq_file,dir_file),
     -				beta_src(freq_file,dir_file),
     -				I_src(freq_file,dir_file),
     -				D_src(freq_file,dir_file),
     -				shift_spec(freq_file,dir_file),
     -				inc_ang_spec(freq_file,dir_file),
     -				cosA(freq_file,dir_file),
     -				sinA(freq_file,dir_file),
     -				D_src_nl(freq_file,dir_file))	

                  L=10.*depth
				per=sqrt(2*3.14*L/9.81*1./tanh(2.*3.14*depth/L)) 

                  k_1=2.*3.14/L            ! wavenumber
                  w_1=2.*3.14/per            ! frequency

                  L_min=L
                  per_min=per
                  L_max=L
                  per_max=per
                  
                  num_freq=1
                  num_theta=1 

            elseif(spec_type.eq.0)then
                  open(93,file='spectrum.dat',status='old')            
                  read(93,*) tmp1,tmp2,tmp
                  read(93,*) wave_hgt,per,inc_ang
                  num_freq=nint(tmp1)
                  num_theta=nint(tmp2)  
								
				freq_file=num_freq
				dir_file=num_theta				
				
				ALLOCATE(amp_spec(freq_file,dir_file),
     -				per_spec(freq_file,dir_file),
     -				L_spec(freq_file,dir_file),
     -				beta_src(freq_file,dir_file),
     -				I_src(freq_file,dir_file),
     -				D_src(freq_file,dir_file),
     -				shift_spec(freq_file,dir_file),
     -				inc_ang_spec(freq_file,dir_file),
     -				cosA(freq_file,dir_file),
     -				sinA(freq_file,dir_file),
     -				D_src_nl(freq_file,dir_file))								
				          
                  L=2.*depth
                  step=L/10.
                  shift=0
                  do while(step.gt.L/100000)
                        tmp=sqrt(2*3.14*L/9.81*1./tanh(2.*3.14*depth/L))
                        if(tmp.gt.per)then
                              L=L-step
                              if(shift.eq.1)then
                                 shift=0
                                 step=step/10
                              else
                                 shift=-1
                              endif
                        else
                              L=L+step
                              if(shift.eq.-1)then
                                 shift=0
                                 step=step/10
                              else
                                 shift=1
                              endif
                        endif
                  enddo

                  k_1=2.*3.14/L            ! wavenumber
                  w_1=2.*3.14/per            ! frequency
                  L_max=0.
                  per_max=0.
                  L_min=1.e10
                  per_min=1.e10
                  wave_hgt=0.
                  do i=1,num_freq
                    do j=1,num_theta
                        read(93,*) amp_spec(i,j),per_spec(i,j),
     -					inc_ang_spec(i,j)
                        wave_hgt=wave_hgt+amp_spec(i,j)
                        L_spec(i,j)=2.*depth
                        step=L_spec(i,j)/10.
                        shift=0

					  if(j.eq.1)then
						    inc_ang=inc_ang_spec(i,j)
					  endif

                        do while(step.gt.L_spec(i,j)/100000)
                              tmp=sqrt(2*3.14*L_spec(i,j)/9.81*
     -                        1./tanh(2.*3.14*depth/L_spec(i,j)))
                              if(tmp.gt.per_spec(i,j))then
                                    L_spec(i,j)=L_spec(i,j)-step
                                    if(shift.eq.1)then
                                       shift=0
                                       step=step/10
                                    else
                                       shift=-1
                                    endif
                              else
                                    L_spec(i,j)=L_spec(i,j)+step
                                    if(shift.eq.-1)then
                                       shift=0
                                       step=step/10
                                    else
                                       shift=1
                                    endif
                              endif
                        enddo
c                        if(L_spec(i,j).gt.L_max)L_max=L_spec(i,j)
c                        if(per_spec(i,j).gt.per_max)per_max=per_spec(i,j)
c                        if(L_spec(i,j).lt.L_min)L_min=L_spec(i,j)
c                        if(per_spec(i,j).lt.per_min)per_min=per_spec(i,j)
                    enddo
                  enddo
                  wave_hgt=wave_hgt/2.
                  L_min=L
                  per_min=per
                  L_max=L
                  per_max=per
            elseif(spec_type.eq.1)then
                  write(flou4, '(a4,a6)') clave , 'in.dat'
                  open(93,file=flou4,status='old')      
                  read(93,*) num_freq,num_theta

				freq_file=num_freq
				dir_file=num_theta				
				
				ALLOCATE(amp_spec(freq_file,dir_file),
     -				per_spec(freq_file,dir_file),
     -				L_spec(freq_file,dir_file),
     -				beta_src(freq_file,dir_file),
     -				I_src(freq_file,dir_file),
     -				D_src(freq_file,dir_file),
     -				shift_spec(freq_file,dir_file),
     -				inc_ang_spec(freq_file,dir_file),
     -				cosA(freq_file,dir_file),
     -				sinA(freq_file,dir_file),
     -				D_src_nl(freq_file,dir_file))	

                  L_max=0.
                  per_max=0.
                  L_min=1.e10
                  per_min=1.e10
                  do i=1,num_freq
                        read(93,*) tmp1
                        do j=1,num_theta
                              read(93,*) amp_spec(i,j),inc_ang_spec(i,j)
                              per_spec(i,j)=tmp1
                              L_spec(i,j)=2.*depth

							if(num_theta.eq.1)then
								inc_ang=inc_ang_spec(i,j)
							endif

                              step=L_spec(i,j)/10.
                              shift=0
                              do while(step.gt.L_spec(i,j)/100000)
                                    tmp=sqrt(2*3.14*L_spec(i,j)/9.81*
     -                              1./tanh(2.*3.14*depth/L_spec(i,j)))
                                    if(tmp.gt.per_spec(i,j))then
                                          L_spec(i,j)=L_spec(i,j)-step
                                          if(shift.eq.1)then
                                             shift=0
                                             step=step/10
                                          else
                                             shift=-1
                                          endif
                                    else
                                          L_spec(i,j)=L_spec(i,j)+step
                                          if(shift.eq.-1)then
                                             shift=0
                                             step=step/10
                                          else
                                             shift=1
                                          endif
                                    endif
                              enddo
                              if(L_spec(i,j).gt.L_max)L_max=L_spec(i,j)
                      if(per_spec(i,j).gt.per_max)per_max=per_spec(i,j)
                              if(L_spec(i,j).lt.L_min)L_min=L_spec(i,j)
                      if(per_spec(i,j).lt.per_min)per_min=per_spec(i,j)
                        enddo
                  enddo

                  read(93,*) per,wave_hgt,tmp1,tmp2
                  wave_hgt=wave_hgt/2.
                  per=1./per
                  L=2.*depth
                  step=L/10.
                  shift=0
                  do while(step.gt.L/100000)
                        tmp=sqrt(2*3.14*L/9.81*1./tanh(2.*3.14*depth/L))
                        if(tmp.gt.per)then
                              L=L-step
                              if(shift.eq.1)then
                                 shift=0
                                 step=step/10
                              else
                                 shift=-1
                              endif
                        else
                              L=L+step
                              if(shift.eq.-1)then
                                 shift=0
                                 step=step/10
                              else
                                 shift=1
                              endif
                        endif
                  enddo
                  k_1=2.*3.14/L            ! wavenumber
                  w_1=2.*3.14/per            ! frequency
            endif
      endif


      if(sim_opt.eq.1)then      
            if(wave_type.le.3)then
                  wave_hgt_corr=1.
            else
                  wave_hgt_corr=1.      
            endif
      elseif(sim_opt.eq.2)then
            wave_hgt_corr=5.
            if(slide_type.eq.4)wave_hgt=.1
      endif

      alpha=wave_hgt/depth      ! nonlinearity parameter

      ep2=1e-6
      ep3=1e-2

      if(sim_opt.eq.1)then
            if(wave_type.le.3)then
                  a_abs=abs(alpha)
                  k_1=sqrt(3.*a_abs/(4.*(1.+a_abs)))      ! non-dimensional wave-number for solitary wave
                  k_dim=k_1/depth
                  c=sqrt(1.+a_abs)                              ! non-dimensional wave speed for solitary wave
                  co=sqrt(9.81*depth)                              ! linear shallow water wave speed
                  per=10/(k_1*c)*depth/co                        ! period of solitary wave
                  L=per*c*co/2.
                  L_min=L

				num_freq=1
				num_theta=1
				L_2=1e6
			    per_2=1e6
				k_2=2.*3.14/L_2            ! wavenumber
				w_2=2.*3.14/per_2            ! frequency
            elseif(wave_type.eq.4)then
C********************************** CNOIDAL WAVE CALCS  ************************

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
                  c=SQRT(9.81*depth*(1.0+wave_hgt/depth*(-1.0+2.0/mk2-
     -                  3.0*CE/mk2/CK)))
                  per=L/c
                  co=sqrt(9.81*depth)                              ! linear shallow water wave speed
                  L_min=L
			
				k_1=2.*3.14/L            ! wavenumber
				w_1=2.*3.14/per            ! frequency

				num_freq=1
				num_theta=1
				L_2=1e6
			    per_2=1e6
				k_2=2.*3.14/L_2            ! wavenumber
				w_2=2.*3.14/per_2            ! frequency

            elseif(wave_type.eq.5)then
C********************************** SIN WAVE CALCS  ************************
				if(current.eq.1)then
					per=1/(Fr_cur*sqrt(9.81*depth)/L+
     -				sqrt(9.81/(2*3.1415*L)*tanh(2.*3.14*depth/L)))            
				else
					per=sqrt(2*3.14*L/9.81*1./tanh(2.*3.14*depth/L))  ! calculated period
				endif
            
				k_1=2.*3.14/L            ! wavenumber
				w_1=2.*3.14/per            ! frequency

				num_freq=2
				num_theta=1
				if(abs(wave_hgt_2).lt.depth/1.e5)then
					L_2=1e6
					num_freq=1
				endif

				L_min=L
				if(current.eq.1)then
					per_2=1/(Fr_cur*sqrt(9.81*depth)/L_2+
     -               sqrt(9.81/(2*3.1415*L_2)*tanh(2.*3.14*depth/L_2)))            
				else                              
					per_2=sqrt(2*3.14*L_2/9.81*1./tanh(2.*3.14*depth/L_2))  ! calculated period
				endif
				 k_2=2.*3.14/L_2            ! wavenumber
				w_2=2.*3.14/per_2            ! frequency
            
		        co=sqrt(9.81*depth)                              ! linear shallow water wave speed
                  c=co
            elseif(wave_type.eq.6)then
C********************************** SIN WAVE CALCS  ************************

				k_1=2.*3.14/L            ! wavenumber
				w_1=2.*3.14/per            ! frequency

		        co=sqrt(9.81*depth)                              ! linear shallow water wave speed
                  c=co
            elseif(wave_type.eq.7.or.wave_type.eq.8)then
C********************************** HOT START FOR LANDSLIDE  ************************
                  co=sqrt(9.81*depth)                              ! linear shallow water wave speed
                  per=L/co
                  L_min=L
            endif
      endif


      return

      end

