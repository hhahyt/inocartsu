C*******************************************************************************
C......  Calculates total mass (difference from still water condition) in numerical domain
      subroutine find_max_zeta(mx_zeta_loc,n_loc)
      use mainvar_module, only:i,j,n,overlap,endy,endx,zeta,
     -			bl_hor_wall,h
	integer n_loc
      real mx_zeta_loc,H_total

	mx_zeta_loc=0.
	do j=overlap+1,endy-overlap
		do i=overlap+1,endx-overlap
			if(abs(zeta(i,j,n_loc,1)).gt.abs(mx_zeta_loc).and.
     -						bl_hor_wall(i,j).le.50)
     -            mx_zeta_loc=abs(zeta(i,j,n_loc,1))
			
c			H_total=h(i,j,n_loc)+zeta(i,j,n_loc,1)

c			if(abs(H_total).lt.abs(mx_zeta_loc))
c     -					mx_zeta_loc=H_total

          enddo
      enddo

      return
      end      ! end subroutine
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

 
