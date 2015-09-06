      subroutine write_timeseries(n_loc,nn_loc)
      use mainvar_module
	integer n_loc,nn_loc
	real u_DA,v_DA,w_DA,flux_x

cZZZZZZZZZZ
         do cur_level=1,1
cZZZZZZZZZZ
            do cur_ts=1,num_ts
			if(ind_ts(cur_ts,1).eq.1)then
	
				i=ind_ts(cur_ts,2)
				j=ind_ts(cur_ts,3)

			  if(nn_loc.eq.3)then
				write(10119+cur_ts,*) x(i),y(j),h(i,j,n_loc), myrank
			  else

                  if(bl_hor_wall(
     -             ind_ts(cur_ts,2),ind_ts(cur_ts,3)).ne.99)then
					call calc_vel_DA(n_loc,u_DA,v_DA,w_DA)

					flux_x=u_DA*(zeta(i,j,n_loc,1)+h(i,j,n_loc))

c					write(10119+cur_ts,*) t(nn_loc+1),
c     -                  zeta(i,j,n_loc,1),u(i,j,n_loc,1),flux_x

					write(10119+cur_ts,*) t(nn_loc+1),
     -                  zeta(i,j,n_loc,1),u_DA,v_DA

c					write(10119+cur_ts,*) t(nn_loc+1),
c     -                  zeta(i,j,n_loc,1),u(i,j,n_loc,1),v(i,j,n_loc,1)

                  else
					write(10119+cur_ts,*) t(nn_loc+1),
     -                  -h(i,j,n_loc),0.,0.
                  endif

			  endif
			endif
            enddo
cZZZZZZZZZZ
         enddo
cZZZZZZZZZZ

      return
      end      ! end subroutine
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

 
