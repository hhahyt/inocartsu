C*******************************************************************************
C......  Calculates total mass (difference from still water condition) in numerical domain
      subroutine calc_Courant_number(n_loc)
      use mainvar_module, only: Cr_max,j,js1,je1,i,endx,Cr_c,nut,
     -		u,v,dx,dt,screen_output,endy,overlap,bl_hor_wall,myrank
	integer n_loc
	real Diff_c,Diff_max

		   Cr_max=0.
		   do j=overlap+1,endy-overlap
				do i=overlap+1,endx-overlap
					if(bl_hor_wall(i,j).le.50)then
						Cr_c=sqrt(u(i,j,n_loc,1)**2.+
     -						      v(i,j,n_loc,1)**2.)*dt/dx
						if(Cr_c.gt.Cr_max)Cr_max=Cr_c

						Diff_c=nut(i,j,n_loc)*dt**2./dx
						if(Diff_c.gt.Diff_max)Diff_max=Diff_c

					endif
				enddo
		   enddo

		  if(Cr_max.gt.0.5.and.screen_output.eq.1)then
       print*,'!!!!!! Courant Stability exceeded !!!!!!!!'
       print*,'!! Local Max Cr: ', Cr_max, ' on subgrid rank: ', myrank
       print*,'!! Local Max Diff: ', Diff_max, ' on subgrid rank: ', myrank
          endif
      return
      end      ! end subroutine
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

 
