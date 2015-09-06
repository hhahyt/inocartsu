C*******************************************************************************
C........    Groups the physical variables into the 
C........... conservative forms given in Wei & Kirby
      subroutine global_dims

	use mainvar_module

	real xt_1,xt_m,yt_1,yt_n

C************ Sponge Layer Characteristics
      if(wave_type.eq.6)then
            sponge_width=sponge_width*L_max
            sigma=2*3.14159/(per_max)            !freqency to be generated
      else
            sponge_width=sponge_width*L
            sigma=2*3.14159/(per)            !freqency to be generated
      endif


      if(spng_1.eq.1.and.dim.eq.1.and.load_topo.eq.1)then
            do i=1,num_ts
                  x_ts(i)=x_ts(i)+sponge_width
            enddo
      endif

      if(spng_3.eq.1.and.dim.eq.1.and.load_topo.eq.1)then
            do i=1,num_ts
                  y_ts(i)=y_ts(i)+sponge_width
            enddo
      endif


      if(load_topo.eq.1)then
 
			if(spng_1.eq.1)then
                  if(num_nodes.eq.8)then
                        x_9=x_8+sponge_width
                        h_9=h_8
                  endif
                  if(num_nodes.ge.7)then
                        x_8=x_7+sponge_width
                        h_8=h_7
                  endif
                   if(num_nodes.ge.6)then
                        x_7=x_6+sponge_width
                        h_7=h_6
                  endif
                   if(num_nodes.ge.5)then
                        x_6=x_5+sponge_width
                        h_6=h_5
                  endif
                   if(num_nodes.ge.4)then
                        x_5=x_4+sponge_width
                        h_5=h_4
                  endif
                   if(num_nodes.ge.3)then
                        x_4=x_3+sponge_width
                        h_4=h_3
                  endif
                  x_3=x_2+sponge_width
                  h_3=h_2
                  
                  x_2=x_1+sponge_width
                  h_2=h_1
                  x_1=x_1
                  h_1=h_1
                  num_nodes=num_nodes+1
                  slide_node=slide_node+1

                  if(slide_node.eq.1) slope_reg=x_2-x_1
                  if(slide_node.eq.2) slope_reg=x_3-x_2
                  if(slide_node.eq.3) slope_reg=x_4-x_3
                  if(slide_node.eq.4) slope_reg=x_5-x_4
                  if(slide_node.eq.5) slope_reg=x_6-x_5
                  if(slide_node.eq.6) slope_reg=x_7-x_6
                  if(slide_node.eq.7) slope_reg=x_8-x_7
                  if(slide_node.eq.8) slope_reg=x_9-x_8

                  if(slide_node.eq.1) init_reg=0
                  if(slide_node.eq.2) init_reg=x_2-x_1
                  if(slide_node.eq.3) init_reg=x_3-x_1
                  if(slide_node.eq.4) init_reg=x_4-x_1
                  if(slide_node.eq.5) init_reg=x_5-x_1
                  if(slide_node.eq.6) init_reg=x_6-x_1
                  if(slide_node.eq.7) init_reg=x_7-x_1
                  if(slide_node.eq.8) init_reg=x_8-x_1

                  x0=x0+sponge_width
			endif

			if(spng_2.eq.1)then
                  if(num_nodes.eq.8)then
                        x_9=x_8+sponge_width
                        h_9=h_8
                  endif
                  if(num_nodes.eq.7)then
                        x_8=x_7+sponge_width
                        h_8=h_7
                  endif
                   if(num_nodes.eq.6)then
                        x_7=x_6+sponge_width
                        h_7=h_6
                  endif
                   if(num_nodes.eq.5)then
                        x_6=x_5+sponge_width
                        h_6=h_5
                  endif
                   if(num_nodes.eq.4)then
                        x_5=x_4+sponge_width
                        h_5=h_4
                  endif
                   if(num_nodes.eq.3)then
                        x_4=x_3+sponge_width
                        h_4=h_3
                  endif
                  if(num_nodes.eq.2)then
                        x_3=x_2+sponge_width
                        h_3=h_2
					  if(slide_type.eq.6.and.sim_opt.eq.2)then
						h_3=h_2+sponge_width*(h_2-h_1)/(x_2-x_1)
					endif

                  endif
                  num_nodes=num_nodes+1
			endif


			if(num_nodes.eq.2) end_x=x_2
			if(num_nodes.eq.3) end_x=x_3
			if(num_nodes.eq.4) end_x=x_4
			if(num_nodes.eq.5) end_x=x_5
			if(num_nodes.eq.6) end_x=x_6
			if(num_nodes.eq.7) end_x=x_7
			if(num_nodes.eq.8) end_x=x_8
			if(num_nodes.eq.9) end_x=x_9

			if(dim.eq.1)then
                  end_y=dy                                      ! total length of y-domain
			elseif(dim.eq.2)then
                  if(chan_width.le.dy*2*overlap) 
     -					chan_width=dy*2*overlap+2*dy

                  if(spng_3.eq.1)chan_width=chan_width+sponge_width
                  if(spng_4.eq.1)chan_width=chan_width+sponge_width

                  end_y=chan_width

			endif


			endx=nint(end_x/dx)+2*overlap                  ! total # of grid points in x-direction
			
			endy=nint(end_y/dy)+2*overlap                  ! total # of grid points in y-direction

			m=1
			n=1

      elseif(load_topo.eq.2)then

		open(56,file='x_topo.dat',status='old') 
		open(57,file='y_topo.dat',status='old')
		open(58,file='f_topo.dat',status='old')
		open(59,file='size_topo.dat',status='old')

		read(59,*) tmp
		m=nint(tmp)
		read(59,*) tmp
		n=nint(tmp)

		do i=1,m
			read(56,*) tmp
			if(i.eq.1)xt_1=tmp
			if(i.eq.m)xt_m=tmp
		enddo

		do i=1,n
			read(57,*) tmp
			if(i.eq.1)yt_1=tmp
			if(i.eq.n)yt_n=tmp
		enddo
 
		close(56)
		close(57)
		close(58)
		close(59)
C....... INTERPOLATE FOR NEW GRID
  
		mm=nint((xt_m-xt_1)/dx)+1+2*overlap

		if(n.eq.1)then
			nn=1+2*overlap
		else
			nn=nint((yt_n-yt_1)/dy)+1+2*overlap
		endif

		endx=mm
		endy=nn  
		
c		if(spng_1.eq.1)then
c              x0=x0+sponge_width
c		endif		
		x0=max(x0,dx*real(1+overlap))
		    
      endif

	endx_glob=endx
	endy_glob=endy

	end_x_t=end_x
	end_y_t=end_y

      return
      end      ! end subroutine
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


