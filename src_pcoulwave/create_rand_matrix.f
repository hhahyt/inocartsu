C*******************************************************************************
C........    Creates matrix of random numbers for backscatter model
      subroutine create_rand_matrix
	  use mainvar_module, only:overlap,endx,endy,myrank
	  use FV_var_module, only:rand
      real tmp_rn
      
      integer i,j
	 
	  do j=1,myrank+1
          call RANDOM_NUMBER(tmp_rn) ! shift random # sequences for different grids
      enddo    

	  do j=1+overlap,endy-overlap
	    do i=1+overlap,endx-overlap
          call RANDOM_NUMBER(tmp_rn)
	      rand(i,j,1)=1.-2.*tmp_rn
          call RANDOM_NUMBER(tmp_rn)
	      rand(i,j,2)=1.-2.*tmp_rn
	    enddo
	  enddo

      return 

      end      ! end subroutine
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


