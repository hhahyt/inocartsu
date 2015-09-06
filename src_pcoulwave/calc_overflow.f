C*******************************************************************************
      subroutine calc_overflow(n_loc)
      use mainvar_module, only:i,j,js1,je1,endx,overflow,zeta
	integer n_loc

	overflow=0
	do j=js1,je1
		do i=1,endx
			if(isnan(zeta(i,j,n_loc,1)))then
				overflow=1  !isnanf for PGI compilers, isnan for CVF
				return
			endif
		enddo
	enddo

      return
      end      ! end subroutine
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

 
