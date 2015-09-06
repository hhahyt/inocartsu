      subroutine add_sponge(n_loc,nn_loc)

      use mainvar_module
      integer n_loc,nn_loc

	cur_level=1
      do j=1+overlap,endy-overlap
       do i=1+overlap,endx-overlap

        if(bl_hor_wall(i,j).le.50.or.level(i,j).le.1)then

C%%%%%%%%%%%%%%%%%%%  SPONGE LAYER  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        w1x=0.
        w2x=0.
        w1y=0.
        w2y=0.
	  w1x_zeta=0.
	  w1y_zeta=0.
        dampfx=0.
        dampfy=0.
        x_rel=0
        tmp=1.
        if(spng_1.eq.1.or.spng_2.eq.1.or.spng_3.eq.1.or.spng_4.eq.1)then
          if(dim.eq.1)then
            if((spng_1.eq.1.and.x(i).le.sponge_width).or.
     -          (spng_2.eq.1.and.x(i).
     -          ge.end_x_t-sponge_width))then

              dudyy=0

              dvdxx=0


              dzdyy=0

              dzdxx=(zeta(i+1,j,n_loc,cur_level)-
     -              2.*zeta(i,j,n_loc,cur_level)+
     -              zeta(i-1,j,n_loc,cur_level))/(dx**2)

              if(x(i).le.sponge_width)then
                x_rel=x(i)-sponge_width      
                x_l=-sponge_width
                dampfx=(exp((x_rel/x_l)**2.)-1.)/(exp(1.)-1.)     
                           
              elseif(x(i).ge.end_x_t-sponge_width)then
                x_rel=x(i)-(end_x_t-sponge_width)
                x_r=end_x_t-(end_x_t-sponge_width)
                dampfx=(exp((x_rel/x_r)**2.)-1.)/(exp(1.)-1.) 
              endif

              w1x=cdamp1*sigma*dampfx
              w2x=cdamp2*sigma*dampfx
              
			  if(x(i).le.sponge_width.and.offshore_cur.eq.1)then
			     w1x=0.		  
			  elseif(x(i).ge.end_x_t-sponge_width.and.offshore_cur.eq.2)then
			     w1x=0.
			  endif
			  
			  w1x_zeta=w1x
              w1y_zeta=0.

            endif
          elseif(dim.eq.2)then      
            if((spng_1.eq.1.and.x(i).le.sponge_width).or.
     -           (spng_2.eq.1.and.x(i).
     -            ge.end_x_t-sponge_width).or.
     -           (spng_3.eq.1.and.y(j).le.sponge_width).or.
     -           (spng_4.eq.1.and.y(j).ge.
     -            end_y_t-sponge_width))then

              dudyy=(u(i,j+1,n_loc,cur_level)-
     -              2.*u(i,j,n_loc,cur_level)+
     -              u(i,j-1,n_loc,cur_level))/(dy**2)

              dvdxx=(v(i+1,j,n_loc,cur_level)-
     -              2.*v(i,j,n_loc,cur_level)+
     -              v(i-1,j,n_loc,cur_level))/(dx**2)

              dzdyy=(zeta(i,j+1,n_loc,cur_level)-
     -              2.*zeta(i,j,n_loc,cur_level)+
     -              zeta(i,j-1,n_loc,cur_level))/(dy**2)

              dzdxx=(zeta(i+1,j,n_loc,cur_level)-
     -              2.*zeta(i,j,n_loc,cur_level)+
     -              zeta(i-1,j,n_loc,cur_level))/(dx**2)

              if(x(i).le.sponge_width.and.spng_1.eq.1)then
                x_rel=x(i)-sponge_width      
                x_l=-sponge_width
                dampfx=(exp((x_rel/x_l)**2.)-1.)/(exp(1.)-1.)   
                tmp=(1-dampfx)  	          
              elseif(x(i).ge.end_x_t-sponge_width.
     -            and.spng_2.eq.1)then
                x_rel=x(i)-(end_x_t-sponge_width)
                x_r=end_x_t-(end_x_t-sponge_width)
                dampfx=(exp((x_rel/x_r)**2.)-1.)/(exp(1.)-1.)
                tmp=(1-dampfx)       
              endif

			  if(y(j).le.sponge_width.
     -            and.spng_3.eq.1)then
                x_rel=y(j)-sponge_width      
                x_l=-sponge_width
                dampfy=(exp((x_rel/x_l)**2.)-1.)/(exp(1.)-1.) 
                tmp=(1-dampfy)    
              elseif(y(j).ge.end_y_t-sponge_width.
     -            and.spng_4.eq.1)then
                x_rel=y(j)-(end_y_t-sponge_width)
                x_r=end_y_t-(end_y_t-sponge_width)
                dampfy=(exp((x_rel/x_r)**2.)-1.)/(exp(1.)-1.)
                tmp=(1-dampfy)        
              endif

              w1x=cdamp1*sigma*dampfx/dim
              w2x=cdamp2*sigma*dampfx/dim
              w1y=cdamp1*sigma*dampfy/dim
              w2y=cdamp2*sigma*dampfy/dim
			  w1x_zeta=w1x
              w1y_zeta=0.
              
              
			  if(x(i).le.sponge_width.and.offshore_cur.eq.1)then
			     w1x=0.	
			     w1y=0.		  	  
			  elseif(x(i).ge.end_x_t-sponge_width.and.offshore_cur.eq.2)then
			     w1x=0.	 
			     w1y=0.		  
			  endif

			  if(y(j).le.sponge_width.and.offshore_cur.eq.3)then  
			     w1x=0.	
			     w1y=0.		  
			  elseif(y(j).ge.end_y_t-sponge_width.and.offshore_cur.eq.4)then  
			     w1x=0.	
			     w1y=0.
			  endif
			  

            endif
          endif
        endif
          
        param_G=(-(w1y)*v(i,j,n_loc,cur_level)+
     -            (w2y)*(dvdxx+dvdyy(i,j,cur_level)+
     -              dvdxy(i,j,cur_level)))
        param_F=(-(w1x)*u(i,j,n_loc,cur_level)+
     -            (w2x)*(dudxx(i,j,cur_level)+dudyy+
     -              dudxy(i,j,cur_level)))
        param_E=-(w1x_zeta+w1y_zeta)*zeta(i,j,n_loc,cur_level)

	  tmp=1.
	  if(numerical_scheme.eq.1) 
     -		tmp=h(i,j,n_loc)+zeta(i,j,n_loc,1)

        E(i,j,n_loc,cur_level)=E(i,j,n_loc,cur_level)+param_E
        F(i,j,n_loc,cur_level)=F(i,j,n_loc,cur_level)+param_F*tmp
        G(i,j,n_loc,cur_level)=G(i,j,n_loc,cur_level)+param_G*tmp

	  endif
	 enddo
	enddo

      return
      end      ! end subroutine
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


