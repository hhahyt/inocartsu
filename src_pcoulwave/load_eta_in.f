      subroutine load_eta_in

      use mainvar_module
	include 'mpif.h'

	 do cur_rank=0,nprocs-1
	  if(myrank.eq.cur_rank)then
		open(999,file='eta_in.dat',status='old')
		do i=1,1000000
			read(999,*) tmp,tmp2
			if(tmp2.gt.1e6) goto 13
			time_bm(i)=tmp
			fse(i)=tmp2
		enddo
13		continue
		close(999)
	  endif
	  call mpi_barrier(mpi_comm_world,ierr)
	 enddo

		do ii=1,tt
			do jj=1,i-1
				if(t(ii).ge.time_bm(jj).and.
     -			   t(ii).lt.time_bm(jj+1))then
					eta_in(ii)=(fse(jj+1)-fse(jj))/
     -					(time_bm(jj+1)-time_bm(jj))*
     -					(t(ii)-time_bm(jj))+fse(jj)
					goto 14

				endif
			enddo
14		continue
		enddo

		Deallocate(time_bm,fse)

      return

      end

