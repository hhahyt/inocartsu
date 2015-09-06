      subroutine bl_define

      use mainvar_module

		do ii=1,endx
			do jj=1,endy
				bl_hor_wall(ii,jj)=0
				level(ii,jj)=0
				bl_x_wall(ii,jj)=0
				bl_y_wall(ii,jj)=0
			enddo
		enddo

		if(bc_1.eq.1.or.bc_1.eq.2.or.bc_1.eq.4)then
			ii=1+overlap
			do jj=1,endy
                  bl_x_wall(ii,jj)=1
                  bl_x_wall(ii+1,jj)=3
			enddo
		elseif(bc_1.eq.3)then
			ii=1+overlap
			do jj=1,endy
                  bl_x_wall(ii,jj)=1
                  bl_x_wall(ii+1,jj)=3
			enddo
		endif

		if(bc_2.eq.1.or.bc_2.eq.2.or.bc_2.eq.4)then
			ii=endx-overlap
			do jj=1,endy
                  bl_x_wall(ii,jj)=2
                  bl_x_wall(ii-1,jj)=4
			enddo
		elseif(bc_2.eq.3)then
			ii=endx-overlap
			do jj=1,endy
                  bl_x_wall(ii,jj)=2
                  bl_x_wall(ii-1,jj)=4
			enddo
		endif

		if(dim.eq.2)then
			if(bc_3.eq.1.or.bc_3.eq.2.or.bc_3.eq.4)then
                  jj=1+overlap
                  do ii=1,endx
                        bl_y_wall(ii,jj)=1
                        bl_y_wall(ii,jj+1)=3
                  enddo
              endif

			if(bc_4.eq.1.or.bc_4.eq.2.or.bc_4.eq.3.or.bc_4.eq.4)then
                  jj=endy-overlap
                  do ii=1,endx
                        bl_y_wall(ii,jj)=2
                        bl_y_wall(ii,jj-1)=4
                  enddo
			endif
		endif


      return

      end

