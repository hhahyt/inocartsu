	subroutine exchange2d_int(endx,endy,overlap,pleft,pright,pbottom,
     -		ptop,LRstride,BTstride,comm2d,dims,var,dim)

	implicit none
	include 'mpif.h'
	integer ierr,endx,endy,overlap,sx,ex,sy,ey,dim,
     -		pleft,pright,pbottom,ptop,
     -		status(mpi_status_size,4),request(4),LRstride,BTstride,
     -		comm2d,dims(2),var(endx,endy)

c	! left and right exchange

	sx=1+overlap
	ex=endx-overlap

	sy=1+overlap
	ey=endy-overlap

	if(dims(2).gt.1.and.dim.eq.2)then
		call mpi_irecv(var(1,ey+1),1,BTstride,ptop,
     -			10,comm2d,request(1),ierr)
		call mpi_irecv(var(1,sy-overlap),1,BTstride,pbottom,
     -			20,comm2d,request(2),ierr)
		call mpi_isend(var(1,sy),1,BTstride,pbottom,
     -			10,comm2d,request(3),ierr)
		call mpi_isend(var(1,ey-overlap+1),1,BTstride,ptop,
     -			20,comm2d,request(4),ierr)
		call mpi_waitall(4,request,status,ierr)
	endif

c	! call mpi_barrier(mpi_comm_world,ierr)

	if(dims(1).gt.1)then
		call mpi_irecv(var(ex+1,1),1,LRstride,pright,
     -			10,comm2d,request(1),ierr)
		call mpi_irecv(var(sx-overlap,1),1,LRstride,pleft,
     -			20,comm2d,request(2),ierr)
		call mpi_isend(var(sx,1),1,LRstride,pleft,
     -			10,comm2d,request(3),ierr)
		call mpi_isend(var(ex-overlap+1,1),1,LRstride,pright,
     -			20,comm2d,request(4),ierr)
		call mpi_waitall(4,request,status,ierr)
	endif

c	!call mpi_barrier(mpi_comm_world,ierr)



	return
	end

