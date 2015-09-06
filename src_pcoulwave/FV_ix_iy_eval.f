      subroutine FV_ix_iy_eval(n_loc,nn_cur)
	use mainvar_module, only:dx,dy,nx,ny,grav,dim,bl_hor_wall,
     -	level,cutoff,endy,endx,overlap,bc_1,bc_2,bc_3,bc_4,
     -	pleft,pright,pbottom,ptop,LRstride,BTstride,comm2d,dims
	use FV_var_module

	integer i,j,n_loc,nn_cur

      do j=2,ny-1
       do i=2,nx-1
	      id(i,j)=0
          ix(i,j)=0
          ix2(i,j)=0  
          iy(i,j)=0  
          iy2(i,j)=0  
       enddo
      enddo

	if(dim.eq.1)then
       do j=1+overlap,endy-overlap
        do i=1+overlap,endx-overlap
		if(bl_hor_wall(i,j).eq.99)then
			 id(i,j)=9
		elseif(bl_hor_wall(i+1,j).eq.99)then
			 id(i,j)=2
		elseif(bl_hor_wall(i-1,j).eq.99)then
			 id(i,j)=1
		endif	
        enddo
       enddo

	 call exchange2d_int(endx,endy,overlap,pleft,pright,pbottom,
     -		ptop,LRstride,BTstride,comm2d,dims,id,dim)  ! exhange shoreline locations

! find ix, ix2
      do j=1,ny
       do i=3,nx-1
		if(id(i,j).eq.2)then
			 ix(i,j)=4
			 ix(i-1,j)=3
c			 ix(i-2,j)=6
			 ix2(i,j)=4
c			 ix2(i-1,j)=3
		elseif(id(i,j).eq.1)then
			 ix(i-1,j)=1
			 ix(i,j)=2
c			 ix(i+1,j)=5
			 ix2(i-1,j)=1
c			 ix2(i,j)=2
		endif	
        enddo
       enddo

      do j=1,ny
       do i=3,nx-1
		if(id(i,j).eq.9)then
			 ix(i,j)=9
			 ix2(i,j)=9
		endif	
        enddo
       enddo

	elseif(dim.eq.2)then
! find id's first
      do j=2,ny-1
       do i=2,nx-1
		if(bl_hor_wall(i,j).eq.99)then
			 id(i,j)=9
		elseif(bl_hor_wall(i+1,j).eq.99)then
			if(bl_hor_wall(i,j+1).eq.99)then
				id(i,j)=24
			elseif(bl_hor_wall(i,j-1).eq.99)then
				id(i,j)=23
			else
				id(i,j)=2
		    endif
		elseif(bl_hor_wall(i-1,j).eq.99)then
			if(bl_hor_wall(i,j+1).eq.99)then
				id(i,j)=14
			elseif(bl_hor_wall(i,j-1).eq.99)then
				id(i,j)=13
			else
				id(i,j)=1
		    endif
		elseif(bl_hor_wall(i,j+1).eq.99)then
			id(i,j)=4
		elseif(bl_hor_wall(i,j-1).eq.99)then
			id(i,j)=3
		endif	
        enddo
       enddo

	call exchange2d_int(endx,endy,overlap,pleft,pright,pbottom,
     -		ptop,LRstride,BTstride,comm2d,dims,id,dim)  ! exhange shoreline locations

! find ix, ix2
      do j=1,ny
       do i=3,nx-1
		if(id(i,j).eq.2.or.id(i,j).eq.23.or.id(i,j).eq.24)then
			 ix(i,j)=4
			 ix(i-1,j)=3
c			 ix(i-2,j)=6
			 ix2(i,j)=4
c			 ix2(i-1,j)=3
		elseif((id(i,j).eq.1.or.id(i,j).eq.13.or.id(i,j).eq.14))then
			 ix(i-1,j)=1
			 ix(i,j)=2
c			 ix(i+1,j)=5
			 ix2(i-1,j)=1
c			 ix2(i,j)=2
		endif	
        enddo
       enddo

! find iy, iy2
      do j=3,ny-1
       do i=1,nx
		if(id(i,j).eq.4.or.id(i,j).eq.14.or.id(i,j).eq.24)then
			 iy(i,j)=4
			 iy(i,j-1)=3
c			 iy(i,j-2)=6
			 iy2(i,j)=4
c			 iy2(i,j-1)=3
		elseif((id(i,j).eq.3.or.id(i,j).eq.13.or.id(i,j).eq.23))then
			 iy(i,j-1)=1
			 iy(i,j)=2
c			 iy(i,j+1)=5
			 iy2(i,j-1)=1
c			 iy2(i,j)=2
		endif	
        enddo
       enddo

      do j=1,ny
       do i=1,nx
		if(id(i,j).eq.9)then
			 ix(i,j)=9
			 ix2(i,j)=9
			 iy(i,j)=9
			 iy2(i,j)=9
		endif	
        enddo
       enddo

	endif                                                            
      return
      end      ! end subroutine
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC












    