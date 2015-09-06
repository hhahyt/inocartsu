	subroutine decomp2d(dims,coords,nx,ny,sx,ex,sy,ey,
     -		bc_1,bc_2,bc_3,bc_4,overlap,sx_X,ex_X,sy_X,ey_X,
     -        ho_glob)

	implicit none
	integer nx,ny,dims(2),coords(2),sx,ex,sy,ey,
     -		bc_1,bc_2,bc_3,bc_4,overlap,sx_X,ex_X,sy_X,ey_X,
     -		i,j,cur_coord
	real ho_glob(nx,ny),sum_x(nx),sum_y(ny),sum,sum_cur,sum_div

	if(dims(1).gt.1)then
		do i=1,nx
			sum_x(i)=0
			do j=1,ny
				if(ho_glob(i,j).gt.0)then
					sum_x(i)=sum_x(i)+1.
				else
					sum_x(i)=sum_x(i)+1.
				endif
			enddo
			sum_x(i)=sum_x(i)/nx
		enddo

		sum=0.
		do i=1,nx
			sum=sum+sum_x(i)
		enddo

		sum_div=sum/dims(1)

		cur_coord=0
		if(coords(1).eq.0)then
			sx=1
		endif

		sum_cur=0.
		do i=1,nx
			sum_cur=sum_cur+sum_x(i)

			if(sum_cur.ge.sum_div)then
				if(coords(1).eq.cur_coord)then
					ex=i
				endif
				if(coords(1).eq.cur_coord+1)then
					sx=i+1
				endif
				sum_cur=0.
				cur_coord=cur_coord+1
			endif
		enddo

		if(coords(1).eq.dims(1)-1)then
			ex=nx
		endif

	else
		sx=1
		ex=nx
	endif



	if(dims(2).gt.1)then
		do j=1,ny
			sum_y(j)=0.
			do i=1,nx
				if(ho_glob(i,j).gt.0)then
					sum_y(j)=sum_y(j)+1.
				else
					sum_y(j)=sum_y(j)+1.
				endif
			enddo
			sum_y(j)=sum_y(j)/ny
		enddo

		sum=0.
		do j=1,ny
			sum=sum+sum_y(j)
		enddo

		sum_div=sum/dims(2)

		cur_coord=0
		if(coords(2).eq.0)then
			sy=1
		endif

		sum_cur=0.
		do j=1,ny
			sum_cur=sum_cur+sum_y(j)

			if(sum_cur.ge.sum_div)then
				if(coords(2).eq.cur_coord)then
					ey=j
				endif
				if(coords(2).eq.cur_coord+1)then
					sy=j+1
				endif
				sum_cur=0
				cur_coord=cur_coord+1
			endif
		enddo

		if(coords(2).eq.dims(2)-1)then
			ey=ny
		endif

	else
		sy=1
		ey=ny
	endif



	if(sx.eq.1.and.ex.eq.nx)then
		sx_X=sx
		ex_X=ex
		bc_1=bc_1
		bc_2=bc_2
	elseif(sx.eq.1)then
		sx_X=sx
		ex_X=ex+overlap
		bc_1=bc_1
		bc_2=99
	elseif(ex.eq.nx)then
		sx_X=sx-overlap
		ex_X=ex
		bc_1=99
		bc_2=bc_2
	else
		sx_X=sx-overlap
		ex_X=ex+overlap
		bc_1=99
		bc_2=99
	endif
						
	if(sy.eq.1.and.ey.eq.ny)then
		sy_X=sy
		ey_X=ey
		bc_3=bc_3
		bc_4=bc_4
	elseif(sy.eq.1)then
		sy_X=sy
		ey_X=ey+overlap
		bc_3=bc_3
		bc_4=99
	elseif(ey.eq.ny)then
		sy_X=sy-overlap
		ey_X=ey
		bc_3=99
		bc_4=bc_4
	else
		sy_X=sy-overlap
		ey_X=ey+overlap
		bc_3=99
		bc_4=99
	endif


	return
	end
