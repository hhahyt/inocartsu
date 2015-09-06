      subroutine write_maxs_means
      use mainvar_module


      do i=1,endx
		do j=1,endy
			mean_z(i,j,1)=mean_z(i,j,1)/
     -              max(1.,(mean_z(i,j,3)-mean_z(i,j,4)))
			mean_hgt(i,j,1)=mean_hgt(i,j,1)/max(1.,mean_hgt(i,j,4))
			cur_x(i,j,1)=cur_x(i,j,1)/
     -               max(1.,(cur_x(i,j,3)-cur_x(i,j,4)))
			cur_y(i,j,1)=cur_y(i,j,1)/
     -               max(1.,(cur_y(i,j,3)-cur_y(i,j,4)))
		enddo
	enddo

      write(UNIT=9000+myrank) max_fs(startx_sponge:endx_sponge,
     -				starty_sponge:endy_sponge),
     -						max_vel(startx_sponge:endx_sponge,
     -				starty_sponge:endy_sponge),
     -						mean_z(startx_sponge:endx_sponge,
     -				starty_sponge:endy_sponge,1),
     -						mean_hgt(startx_sponge:endx_sponge,
     -				starty_sponge:endy_sponge,1),
     -						cur_x(startx_sponge:endx_sponge,
     -				starty_sponge:endy_sponge,1),
     -						cur_y(startx_sponge:endx_sponge,
     -				starty_sponge:endy_sponge,1)

      return
      end      ! end subroutine
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

 
