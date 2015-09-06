
      subroutine store_previous_iter(n_loc)
      use mainvar_module, only: i,j,overlap,zeta,u,v,zeta1_iter,
     -	u1_iter,u_iter,v_iter,endx,endy
	integer n_loc

	do j=1+overlap,endy-overlap                       
		do i=1+overlap,endx-overlap
			zeta1_iter(i,j)=abs(zeta(i,j,n_loc,1))

			u1_iter(i,j)=sqrt(u(i,j,n_loc,1)**2.+
     -                    v(i,j,n_loc,1)**2.)

			u_iter(i,j)=u(i,j,n_loc,1)

			v_iter(i,j)=v(i,j,n_loc,1)
		enddo
	enddo

      return
      end      ! end subroutine
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

 
