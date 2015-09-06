      subroutine write_surfaces(n_loc,nn_loc)
      use mainvar_module
	integer n_loc,nn_loc
	real u_DA, v_DA, w_DA

	if(nn_loc.eq.3)then
		  a=max(nint(writ_inc/dt),1)
		  if(myrank.eq.0)then 
				write(90,*) a, num_levels, nprocs
				close(90)
		  endif
	endif


	if(wave_type.eq.8)then
	 if(dim.eq.1)then
		if(myrank.eq.0)then 
			write(9,*) t(nn_loc), wm_d(nn_loc,1)
		endif            
	 else
		if(myrank.eq.0)then 
			write(9,FMT='(E14.7,$)') t(nn_loc), (wm_d(nn_loc,j), j=1,29)
			write(9,FMT='(E14.7)') wm_d(nn_loc,30)
		endif 
	 endif
	else
	 write(9,*) t(nn_loc)
	endif
		 
	write(UNIT=7000+myrank) bl_hor_wall(
     -				startx_sponge:endx_sponge,
     -			starty_sponge:endy_sponge),
     -						  nut(
     -				startx_sponge:endx_sponge,
     -				starty_sponge:endy_sponge,n_loc),
     -						  h(
     -				startx_sponge:endx_sponge,
     -				starty_sponge:endy_sponge,n_loc)

cZZZZZZZZZZ
			do cur_level=1,num_levels
cZZZZZZZZZZ

				 do j=2,endy-1
					do i=2,endx-1
					    call calc_vel_DA(n_loc,u_DA,v_DA,w_DA)

                          u_z(i,j)=u_DA   
                          v_z(i,j)=v_DA
					enddo
				 enddo


                 write(UNIT=3000+myrank) zeta(
     -				startx_sponge:endx_sponge,
     -				starty_sponge:endy_sponge,n_loc,cur_level)

                 write(UNIT=6000+myrank) u(
     -				startx_sponge:endx_sponge,
     -				starty_sponge:endy_sponge,n_loc,cur_level),
     -				v(startx_sponge:endx_sponge,
     -				starty_sponge:endy_sponge,n_loc,cur_level)

			   if(dim.eq.2)then
				 do j=2,endy-1
					do i=2,endx-1
	
						z=zeta(i,j,4,cur_level)
						call calc_vel_z(4,z,u_z1,v_z1,w_z1)

                          u_z(i,j)=u_z1   
                          v_z(i,j)=v_z1
					enddo
				 enddo

				 do j=3,endy-2
					do i=3,endx-2
						dvdx(i,j,cur_level)=(v_z(i+1,j)
     -							-v_z(i-1,j))/(2.*dx)

						dudy(i,j,cur_level)=(u_z(i,j+1)
     -							-u_z(i,j-1))/(2.*dy)
					enddo
				 enddo


				 do j=1,endy
					do i=1,endx
						vort(i,j)=dvdx(i,j,cur_level)
     -							 -dudy(i,j,cur_level)
					enddo
				 enddo
				endif

				write(UNIT=8000+myrank) vort(
     -				startx_sponge:endx_sponge,
     -			starty_sponge:endy_sponge),u_z(
     -				startx_sponge:endx_sponge,
     -			starty_sponge:endy_sponge),v_z(
     -				startx_sponge:endx_sponge,
     -			starty_sponge:endy_sponge)

cZZZZZZZZZZ
			enddo
cZZZZZZZZZZ

      return
      end      ! end subroutine
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

 
