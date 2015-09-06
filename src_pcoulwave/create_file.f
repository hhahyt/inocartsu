      subroutine create_file(myrank)

C********************  USED FILE MARKERS *****************************

C      8 -- sim_set.dat, simulation setup parameters
C      9 -- time.dat, times that snapshots are written to file
C      56 -- x_topo.dat, input data file containing loaded bathymetry x-axis values
C      57 -- y_topo.dat, input data file containing loaded bathymetry y-axis values
C      58 -- f_topo.dat, input data file containing loaded bathymetry depth values
C      59 -- size_topo.dat, input data file containing dimensions of the loaded depth matrix
C      90 -- s.dat, time step increment that snapshots are written to file
C      91 -- prog.dat, information about each time step
C      92 -- overflows.dat, information about any overflows which may have occured
C      93 -- spectrum.dat, input data file when running a simulation with an amplitude specturm
C      94 -- grid_info.dat, nested grid information
C      95 -- slide_param.dat, landslide information
C      96 -- ts_locations.dat, time series locations
C      99 -- batch.dec, decision variables required to run program in batch on UNIX workstations
C      120-179      -- ts_grp.dat, time series of zeta, u, and v, for up tp 60 locations
C      myrank+1000      -- xlocXXX.dat, x axis values 
C      myrank+2000      -- ylocXXX.dat, y axis values 
C      myrank+3000      -- zetaXXX.dat, snapshot zeta values 
C      myrank+4000      -- dpthXXX.dat, depth matrices     
C      myrank+5000      -- psixXXX.dat, psix, psiy snapshots 
C      myrank+6000      -- veloXXX.dat, u,v velocity snapshots 
C      myrank+7000      -- blvsXXX.dat, shoreline location and wave breaking snapshots 
C      myrank+8000      -- vortXXX.dat, vorticity snapshots
C      myrank+9000      -- mxmnXXX.dat, maximum free surface elevation in time 


      integer myrank
	character*11 filename
	character*3 ranknum

	if(myrank.eq.0)then
		open(9,file='time.dat',status='unknown')
		open(90,file='s.dat',status='unknown')
		open(91,file='prog.dat',status='unknown')
		open(92,file='overflows.dat',status='unknown')
	endif

	filename = 'xlocXXX.dat'  ! x file
	write(ranknum,'(i3.3)')myrank
	filename(5:7)=ranknum
	open(myrank+1000,file=filename,status='unknown')

	filename = 'ylocXXX.dat'  ! y file
	write(ranknum,'(i3.3)')myrank
	filename(5:7)=ranknum
	open(myrank+2000,file=filename,status='unknown')

	filename = 'zetaXXX.dat'  ! free surface file
	write(ranknum,'(i3.3)')myrank
	filename(5:7)=ranknum
	open(myrank+3000,file=filename,form='unformatted',status='unknown')

	filename = 'dpthXXX.dat'   ! initial depth file
	write(ranknum,'(i3.3)')myrank
	filename(5:7)=ranknum
	open(myrank+4000,file=filename,form='unformatted',status='unknown')
	
	filename = 'psixXXX.dat'  ! velocity file
	write(ranknum,'(i3.3)')myrank
	filename(5:7)=ranknum
	open(myrank+5000,file=filename,form='unformatted',status='unknown')

	filename = 'veloXXX.dat'  ! velocity file
	write(ranknum,'(i3.3)')myrank
	filename(5:7)=ranknum
	open(myrank+6000,file=filename,form='unformatted',status='unknown')

	filename = 'blvsXXX.dat'   ! bl_visc file
	write(ranknum,'(i3.3)')myrank
	filename(5:7)=ranknum
	open(myrank+7000,file=filename,form='unformatted',status='unknown')

	filename = 'vortXXX.dat'   ! vorticity file
	write(ranknum,'(i3.3)')myrank
	filename(5:7)=ranknum
	open(myrank+8000,file=filename,form='unformatted',status='unknown')

	filename = 'mxmnXXX.dat'   ! max,min,mean file
	write(ranknum,'(i3.3)')myrank
	filename(5:7)=ranknum
	open(myrank+9000,file=filename,form='unformatted',status='unknown')

c	filename = 'sediXXX.dat'   ! current file
c	write(ranknum,'(i3.3)')myrank
c	filename(5:7)=ranknum
c	open(myrank+10000,file=filename,form='unformatted',status='unknown')

c	filename = 'fluxXXX.dat'   ! current file
c	write(ranknum,'(i3.3)')myrank
c	filename(5:7)=ranknum
c	open(myrank+11000,file=filename,form='unformatted',status='unknown')

c	time series files created ~line 3108 in coulwave_MPI


      return
      end


