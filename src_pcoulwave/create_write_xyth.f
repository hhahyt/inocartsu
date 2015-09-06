      subroutine create_write_xyth

      use mainvar_module

	if(bc_1.eq.1)then
		startx_sponge=1+overlap
	else
		startx_sponge=1+overlap-1
	endif
		 
	if(bc_2.eq.1)then
		endx_sponge=endx-overlap
	else
		endx_sponge=endx-overlap+1
	endif


	if(bc_3.eq.1)then
		starty_sponge=1+overlap
	else
		starty_sponge=1+overlap-1
	endif
		 
	if(bc_4.eq.1)then
		endy_sponge=endy-overlap
	else
		endy_sponge=endy-overlap+1
	endif


      do i=startx_sponge,endx_sponge
            write(1000+myrank,*) x(i)      
      enddo

      close(1000+myrank)

      do j=starty_sponge,endy_sponge
            write(2000+myrank,*) y(j)      
      enddo
      close(2000+myrank)


	t(1)=0.
      do n=2,endt
            t(n)=t(n-1)+dt
      enddo

      do n=1,4
		do i=1,endx
			do j=1,endy
				h(i,j,n)=ho(i,j)
			enddo
		enddo
      enddo

	write(UNIT=4000+myrank) ho(
     -		startx_sponge:endx_sponge,
     -		starty_sponge:endy_sponge)


      close(4000+myrank)



      return

      end

