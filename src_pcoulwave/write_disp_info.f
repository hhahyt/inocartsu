      subroutine write_disp_info(nn_loc)

      use mainvar_module
cc	include 'mpif.h'
	integer nn_loc

C    START OF TIMED SEQUENCE
	if(nn_loc.eq.display.or.nn_loc.eq.3)then
		cpu_0 = SECNDS(0.0)
		cpu_c = SECNDS(cpu_0)      ! DELTA gives the elapsed time
cc		cpu_dt_inc = real(MPI_Wtime())
		call cpu_time(cpu_dt_inc)
	endif
 
	cpu_last = cpu_c     ! DELTA gives the elapsed time
	cpu_c = SECNDS(cpu_0)      ! DELTA gives the elapsed time
cc	cpu_dt = real(MPI_Wtime())  ! gives cpu time
	call cpu_time(cpu_dt)

	count_cpu=count_cpu+1
	cpu_per_step=cpu_c-cpu_last
c	sum_per_step=sum_per_step+cpu_per_step
	sum_per_step=cpu_dt-cpu_dt_inc
	ave_per_step=sum_per_step/real(count_cpu)/real(display)
	hrs=int(ave_per_step*real(endt-1-nn_loc)/60./60.)
	minutes=int((ave_per_step*
     -		real(endt-1-nn_loc)/60./60.-real(hrs))*60.)
	seconds=-hrs*60*60+int((ave_per_step*
     -		real(endt-1-nn_loc)/60.-real(minutes))*60.)

	if(screen_output.eq.1)then
C...........Write to screen 
c            write(*,*) 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
            write(*,*) 'Current time (sec)',t(nn_loc),' out of', end_t
            write(*,*) 'Current step',nn_loc,' out of', endt-1
c            write(*,*) 'Elapsed wall clock time (sec)=',cpu_c
           write(*,*) 'Elapsed wall clock time (sec)', cpu_dt-cpu_dt_inc
            write(*,*) 'Number of iterations in last step=',count_min

	      write(*,*) 'Max free surface=',mx_zeta_g
c	      write(*,*) 'Relative mass in domain=',mass_g/init_mass_g
c            write(*,*) 'Average wall clock time per time step (sec)', 
c     -                    ave_per_step
            write(*,*) 'Estimated time until completion (hrs:min:sec)'
            write(*,*) hrs,minutes,seconds	
		  write(*,*) '-----------------------------------------------'
	endif

C...........Write to file
	write(91,*) 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
	write(91,*) 'Step',nn_loc,' out of', endt-1
	write(91,*) 'Elapsed wall clock time (sec)',cpu_dt-cpu_dt_inc
	write(91,*) 'Number of iterations in last step=',count_min
	write(91,*) 'Average Wall Clock time per time step (sec)', 
     -                    ave_per_step
	write(91,*) 'Estimated time until completion (hrs:min:sec)'
	write(91,*) hrs,minutes,seconds
	write(91,*) 'Max free surface=',mx_zeta_g
	write(91,*) 'Relative mass in domain=',mass_g/init_mass_g
	write(91,*) '----------------------------------------------'


      return
      end      ! end subroutine
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

 
