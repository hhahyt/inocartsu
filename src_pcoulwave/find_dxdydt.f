      subroutine find_dxdydt

      use mainvar_module

C.....Grid spaces
      if(sim_opt.eq.1)then
            dx=L_min/pts_wvl
      elseif(sim_opt.eq.2.and.slide_type.lt.3.or.slide_type.ge.6)then
            if(load_topo.eq.1.and.slide_type.lt.3)then
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
            elseif(load_topo.eq.1.and.slide_type.eq.6)then
                  slope_reg=1.   !FOR COSTAS SLIDE 1./wave_hgt !
            elseif(load_topo.eq.1.and.slide_type.eq.7)then
                  slope_reg=L
            endif

            co=sqrt(9.81*depth)                              ! linear shallow water wave speed
            L=slope_reg
            per=L/co
            dx=slope_reg/pts_wvl                        ! delta x

      elseif(sim_opt.eq.2.and.slide_type.eq.3.or.slide_type.eq.5)then
            co=sqrt(9.81*depth)                              ! linear shallow water wave speed
            dx=.01
c            depth/4
            L=3.
            per=L/co
      elseif(sim_opt.eq.2.and.slide_type.eq.4)then
            depth=0.1
            co=sqrt(9.81*depth)                              ! linear shallow water wave speed
            dx=1./pts_wvl
      endif

      dy=dx                                           ! delta y
      dt=min(dx,dy)/co*courant                              ! delta t

c	open(96,file='dxdydt.dat',status='old')
c	read(96,*) dx, dy, dt
c	close(96)

	endt=nint(end_t/dt)+2                  ! total # of time steps

      return

      end

