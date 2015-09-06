      subroutine calc_corrector_error(n_loc)
      use mainvar_module
	include 'mpif.h'

	integer n_loc

      maxerr=0.
	maxerr2=0.
	sum_u=0.
	sum_u_dif=0.
	sum_z=wave_hgt/real(endx*endy)
	sum_z_dif=0.
	sum_v=0.
	sum_v_dif=0.

      cur_level=1
	do j=overlap+1,endy-overlap
		do i=wvmk_loc_i(j),endx-overlap
			if(bl_hor_wall(i,j).lt.99)then
				erru2=0.
				errv2=0.
				errz2=0.

				z=abs(zeta(i,j,n_loc,1))
				co_c=sqrt(9.81*abs(h(i,j,n_loc)+zeta(i,j,n_loc,1)))
				u_vec=sqrt(u(i,j,n_loc,1)**2.+v(i,j,n_loc,1)**2.)
				
				if(u_vec/co_c.lt.1.e-5.or.u_vec.lt.1.e-5)then
						erru2=0.
				else
						erru2=abs(u_vec-u1_iter(i,j))/u_vec !local error
				endif

c				sum_u_dif=sum_u_dif+abs(u_vec-u1_iter(i,j))  !global error
c				sum_u=sum_u+u_vec  !global sum, for relative error calc

				if(abs(z/h(i,j,n_loc)).lt.1.e-4.or.z.lt.1e-4)then
					errz2=0.
				sum_z_dif=sum_z_dif !global error
				sum_z=sum_z !global sum, for relative error calc
				else
					errz2=abs((z-zeta1_iter(i,j))/z) !local error
				sum_z_dif=sum_z_dif+abs(z-zeta1_iter(i,j)) !global error
				sum_z=sum_z+z !global sum, for relative error calc
				endif

				err2=max(abs(erru2),abs(errz2))
				
				if(err2.gt.maxerr2) then
					maxerr2=err2
				endif

			endif
		enddo
	enddo   

	call mpi_barrier(mpi_comm_world,ierr)

	call MPI_ALLREDUCE(sum_u, tmp, 1,MPI_REAL, 
     -					MPI_SUM, mpi_comm_world,ierr) 
						    
	sum_u=tmp

	call MPI_ALLREDUCE(sum_u_dif, tmp, 1,MPI_REAL, 
     -					MPI_SUM, mpi_comm_world,ierr) 
						    
	sum_u_dif=tmp

	call MPI_ALLREDUCE(sum_z, tmp, 1,MPI_REAL, 
     -					MPI_SUM, mpi_comm_world,ierr) 
						    
	sum_z=tmp

	call MPI_ALLREDUCE(sum_z_dif, tmp, 1,MPI_REAL, 
     -					MPI_SUM, mpi_comm_world,ierr) 
						    
	sum_z_dif=tmp
                              
	if(sum_u.lt.ep)then
		erru=sum_u_dif
	else
		erru=sum_u_dif/sum_u
	endif

	if(sum_z.lt.ep)then
		errz=sum_z_dif
	else
		errz=sum_z_dif/sum_z
	endif

	maxerr=max(abs(erru),abs(errz),1e-10)

	call MPI_ALLREDUCE(maxerr2, tmp2, 1,MPI_REAL, 
     -					MPI_MAX, mpi_comm_world,ierr) 

	maxerr2=tmp2


CCCCCCCCCCCCCCCCCCCCCCCCCCCCC end error calculation CCCCCCCCCCCCCCCCCCCCCCCCCCC
C ..............if maximum number of iterations is exceeded, artifically end loop
	if(count_min.gt.itr)then
		if(screen_output.eq.1.and.myrank.eq.0)then
c			print*, 'Maximum number of iterations reached'
c			print*,'Fractional difference between
c     & last two iterations: ', maxerr
c			print*,'Maximum local fractional difference between
c     & last two iterations: ', maxerr2
		endif
		maxerr=0.
		maxerr2=0.
	endif 

      return
      end      ! end subroutine
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

 
