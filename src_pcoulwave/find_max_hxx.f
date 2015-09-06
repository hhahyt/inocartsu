      subroutine find_max_hxx

      use mainvar_module
	include 'mpif.h'
	
	hxx_max=0.
	hx_max=0.
	hx=hx*0.
	hy=hy*0.
	hxx=hxx*0.
	hyy=hyy*0.
	hxy=hxy*0.
	
		do i=overlap+1,endx-overlap
			do j=overlap+1,endy-overlap

				hx(i,j)=(-ho(i+2,j)+
     -					8.*ho(i+1,j)-
     -					8.*ho(i-1,j)+
     -					ho(i-2,j))/(12.*dx)

				hxx(i,j)=(-ho(i+2,j)+
     -					16.*ho(i+1,j)-
     -					30.*ho(i,j)+
     -					16.*ho(i-1,j)-
     -					ho(i-2,j))/(12.*dx**2)
                
				if(dim.eq.2)then
C............ d1v and d2v are used to convert VV to v
					hy(i,j)=(-ho(i,j+2)+
     -					8.*ho(i,j+1)-
     -					8.*ho(i,j-1)+
     -					ho(i,j-2))/(12.*dy)

					hyy(i,j)=(-ho(i,j+2)+
     -					16.*ho(i,j+1)-
     -					30.*ho(i,j)+
     -					16.*ho(i,j-1)-
     -					ho(i,j-2))/(12.*dy**2)

                      hxy(i,j)=(ho(i+1,j+1)+
     -                      ho(i-1,j-1)-
     -                      ho(i+1,j-1)-
     -                      ho(i-1,j+1))/(4.*dy*dx)

				endif

							if(bl_y_wall(i,j).eq.1.or.
     -							bl_y_wall(i,j).eq.2)then
     
									hxy(i,j)=0
							endif

							if(bl_x_wall(i,j).eq.1.or.
     -							bl_x_wall(i,j).eq.2)then
     
									hxy(i,j)=0
							endif
                  hxx_max=max(hxx_max,abs(hxx(i,j)),
     -                             abs(hyy(i,j)))
                  hx_max=max(hx_max,abs(hx(i,j)),
     -                             abs(hy(i,j)))
                  
            enddo
      enddo      

	call mpi_barrier(mpi_comm_world,ierr)

	call MPI_ALLREDUCE(hx_max, tmp, 1,MPI_REAL, 
     -					MPI_MAX, mpi_comm_world,ierr) 

	hx_max=tmp

	call MPI_ALLREDUCE(hxx_max, tmp, 1,MPI_REAL, 
     -					MPI_MAX, mpi_comm_world,ierr) 

	hxx_max=tmp

	if(screen_output.eq.1.and.myrank.eq.0)then
		print*,'Max slope of seafloor ',hx_max
		print*,'Max curvature of seafloor ',hxx_max
	endif

      return

      end

