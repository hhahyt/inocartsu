C*******************************************************************************
C......  Internal source generator
      subroutine internal_source_type2(n_loc,nn_loc)
      use mainvar_module
	integer n_loc,nn_loc
	
	if(nn_loc.eq.4.and.wave_type.eq.6.and.spec_type.eq.2)then
		call load_eta_in
	elseif(nn_loc.eq.4.and.plus_tide.eq.1)then
		call load_eta_in
	endif

	if(plus_tide.ne.1)then  ! turn off this loop for added tide...

       do j=1,endy
        do i=1,endx
			  tmp2=min(1.,max(0.,h(i,j,n_loc)/depth))**2.

			  if(is_oreint.eq.1)then
				x_c=x(i)-end_x_t/2.
				y_c=y(j)-x0
			  else
c				x_c=x(i)- x0 !1900  !special for IHNC
c				y_c=y(j)- end_y_t/2. !4000

				x_c=x(i)-x0
				y_c=y(j)-end_y_t/2.
			  endif

				tmp=0

				do freq=1,num_freq
					do theta=1,num_theta

						if(nn_loc.le.4)then

							cosA(freq,theta)=
     -							cos(inc_ang_spec(freq,theta))
							sinA(freq,theta)=
     -							sin(inc_ang_spec(freq,theta))

							wA(freq,theta)=2.*3.1415/per_spec(freq,theta)

						endif

						loc=abs(((x_c- !nint(real(theta/6))*L_spec(freq,theta)   !to spread multi dir spec
     -						+L_spec(freq,theta)*
     -						sinA(freq,theta))*cosA(freq,theta)+
     -						(y_c+L_spec(freq,theta)*
     -						sinA(freq,theta))*sinA(freq,theta)))


						w_tmp=wA(freq,theta)
						k_tmp=2.*3.1415/L_spec(freq,theta)
						
						if(wave_type.eq.1)then
							zc=0.
							call solit(uc,vc,zc,c*co,co,t(nn_loc),
     -                           0.,0.,0.,
     -                           depth,alp,inc_ang,
     -                           wave_type,wave_hgt,
     -                           L,bf_ratio,
     -						   bet,cur_level,num_levels)


							tmp=tmp+
     @							D_src(freq,theta)*
     @							exp(-beta_src(freq,theta)*loc**2.)*
     @							zc/wave_hgt

						elseif(wave_type.lt.5)then
							call cnoidal(wave_hgt/depth,mk2,
     -                           0.,0.,t(nn_loc),zc,uc,vc,yt,
     -                           yt2,rp1,r2,CK,depth,
     -                           L,per,co,inc_ang,
     -                           wave_type,bet(cur_level)*depth)

							tmp=tmp+
     @							D_src(freq,theta)*
     @							exp(-beta_src(freq,theta)*loc**2.)*
     @							2.*zc/wave_hgt

						elseif(wave_type.eq.6.and.spec_type.eq.2)then

							tmp=tmp+
     @							D_src(freq,theta)*
     @							exp(-0.1*beta_src(freq,theta)*loc**2.)*  ! coef is 0.1 for USGS/large scale tsunami sims, 10 otherwise
     @							2.*eta_in(nn_loc)*(1+eta_in(nn_loc)/depth)

						elseif(spec_type.eq.2)then
											
											
!							if(loc/(1.0*L_spec(freq,theta)).le.0.75) ! y>1600 for 17th st runs 1,2 1240 for run 3, 1200 for Orleans, 1400 for london, 3300 for IHNC, x<700 MRGOb2d
     							tmp=tmp+
     -							   D_src(freq,theta)*
     -							   exp(-beta_src(freq,theta)*loc**2.)*
     -							   cos(-w_tmp*(t(nn_loc)+
     -							   shift_spec(freq,theta)))*tmp2*
     -                               (min(1.,t(nn_loc)/(2.*per)))
						endif
					enddo
				enddo

				zeta(i,j,n_loc,1)=zeta(i,j,n_loc,1)+tmp
     
     
     
		   enddo
		  enddo
		 endif

      return
      end      ! end subroutine
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

 
