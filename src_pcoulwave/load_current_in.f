      subroutine load_current_in

      use mainvar_module
	include 'mpif.h'

	 do cur_rank=0,nprocs-1
	  if(myrank.eq.cur_rank)then
		open(999,file='current_in.dat',status='old')
		do i=1,1000000
			read(999,*) tmp,tmp2
			if(tmp2.gt.1e6) goto 13
			time_cur(i)=tmp
			curf(i)=tmp2
		enddo
13		continue
		close(999)
	  endif
	  call mpi_barrier(mpi_comm_world,ierr)
	 enddo

		do ii=1,tt
			do jj=1,i-1
				if(t(ii).ge.time_cur(jj).and.
     -			   t(ii).lt.time_cur(jj+1))then
					current_in(ii)=(curf(jj+1)-curf(jj))/
     -					(time_cur(jj+1)-time_cur(jj))*
     -					(t(ii)-time_cur(jj))+curf(jj)
					goto 14

				endif
			enddo
14		continue
		enddo

		Deallocate(time_cur,curf)

      return

      end

