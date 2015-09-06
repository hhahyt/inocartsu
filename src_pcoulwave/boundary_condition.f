      subroutine bc(n_loc,nn_cur)
      use mainvar_module
      integer n_loc,nn_cur,nn1,nn_max,per_shift,ii_bc,jj_bc
	real u_interp1,u_interp2,u_interp3,u_interp4,
     -            u_interp5,u_interp6,u_interp7,u_interp8,
     -            v_interp1,v_interp2,v_interp3,v_interp4,
     -            v_interp5,v_interp6,v_interp7,v_interp8,
     -	u_1,v_1,zeta_1,u_2,v_2,zeta_2,d_interp,H_interp

      nn_max=2*overlap

      do jj=1,endy
       do ii=1,endx
          v_iter(ii,jj)=v(ii,jj,n_loc,cur_level)
          u_iter(ii,jj)=u(ii,jj,n_loc,cur_level)
          zeta_iter(ii,jj)=zeta(ii,jj,n_loc,cur_level)
          if(level(ii,jj).lt.90) level(ii,jj)=0
c          if(bl_hor_wall(ii,jj).eq.99)then
c            v(ii,jj,n_loc,cur_level)=0
c            u(ii,jj,n_loc,cur_level)=0
c            zeta(ii,jj,n_loc,cur_level)=0
c          endif            
       enddo
      enddo


      if(dim.eq.1)then
       do jj=1+overlap,endy-overlap
          do ii=3,endx-2
          if(bl_hor_wall(ii,jj).eq.99)then

            int_count=0
            z_interp1=0
            z_interp2=0

            u_interp1=0
            u_interp2=0
      
            if(bl_hor_wall(ii+1,jj).eq.0)then

c                u_interp1=u(ii+1,jj,n_loc,cur_level)


                int_count=int_count+1
            endif


            if(bl_hor_wall(ii-1,jj).eq.0)then

c                u_interp2=u(ii-1,jj,n_loc,cur_level)

                int_count=int_count+1
            endif
      
            if(int_count.gt.0)then

c                u(ii,jj,n_loc,cur_level)=
c     -				(u_interp1+u_interp2)/int_count/3.
       
                level(ii,jj)=1
            endif
          endif
       enddo
      enddo


      do nn1=1,nn_max
       total_count=0
       do jj=1+overlap,endy-overlap
          do ii=3,endx-2
            if(bl_hor_wall(ii,jj).eq.99.and.
     -                  level(ii,jj).eq.0)then
			  int_count=0
			  z_interp1=0
			  z_interp2=0

			  u_interp1=0
			  u_interp2=0

                if(level(ii+1,jj).eq.nn1)then

c				  u_interp1=u(ii+1,jj,n_loc,cur_level)

				  int_count=int_count+1
                endif


                if(level(ii-1,jj).eq.nn1)then

c				  u_interp2=u(ii-1,jj,n_loc,cur_level)

				  int_count=int_count+1
                endif
      
                if(int_count.gt.0)then
                
c                  u(ii,jj,n_loc,cur_level)=
c     -                  (u_interp1+u_interp2)/int_count
       
                  level(ii,jj)=nn1+1
                endif

                total_count=total_count+int_count

            endif
          enddo
       enddo
       if(total_count.eq.0)goto 50
      enddo

	elseif(dim.eq.2)then
       do jj=3,endy-2
        do ii=3,endx-2
          if(bl_hor_wall(ii,jj).eq.99)then

            int_count=0

            u_interp1=0.
            u_interp2=0.
            u_interp3=0.
            u_interp4=0.

            v_interp1=0.
            v_interp2=0.
            v_interp3=0.
            v_interp4=0.
      
            if(bl_hor_wall(ii+1,jj).le.50)then

c                u_interp1=u(ii+1,jj,n_loc,cur_level)

c                v_interp1=v(ii+1,jj,n_loc,cur_level)

                int_count=int_count+1
            endif


            if(bl_hor_wall(ii-1,jj).le.50)then

c                u_interp2=u(ii-1,jj,n_loc,cur_level)

c                v_interp2=v(ii-1,jj,n_loc,cur_level)

                int_count=int_count+1
            endif
      
            if(bl_hor_wall(ii,jj+1).le.50)then

c                u_interp3=u(ii,jj+1,n_loc,cur_level)

c                v_interp3=v(ii,jj+1,n_loc,cur_level)

                int_count=int_count+1
            endif   

            if(bl_hor_wall(ii,jj-1).le.50)then
c                u_interp4=u(ii,jj-1,n_loc,cur_level)

c                v_interp4=v(ii,jj-1,n_loc,cur_level)

                int_count=int_count+1
            endif

            if(int_count.gt.0)then

c                u(ii,jj,n_loc,cur_level)=(u_interp1+u_interp2+
c     -                u_interp3+u_interp4)/int_count
c
c                v(ii,jj,n_loc,cur_level)=(v_interp1+v_interp2+
c     -                v_interp3+v_interp4)/int_count

                level(ii,jj)=1
            endif
          endif
       enddo
      enddo


      do nn1=1,nn_max
       total_count=0
       do jj=3,endy-2
        do ii=3,endx-2
            if(bl_hor_wall(ii,jj).eq.99.and.
     -                  level(ii,jj).eq.0)then
                int_count=0.

c			  u_interp1=0.
c			  u_interp2=0.
c			  u_interp3=0.
c			  u_interp4=0.
c
c			  v_interp1=0.
c			  v_interp2=0.
c			  v_interp3=0.
c			  v_interp4=0.

                if(level(ii+1,jj).eq.nn1)then

c                  u_interp1=u(ii+1,jj,n_loc,cur_level)

c                  v_interp1=v(ii+1,jj,n_loc,cur_level)

                  int_count=int_count+1
                endif


                if(level(ii-1,jj).eq.nn1)then

c                  u_interp2=u(ii-1,jj,n_loc,cur_level)

c                  v_interp2=v(ii-1,jj,n_loc,cur_level)

                  int_count=int_count+1
                endif
      
                if(level(ii,jj+1).eq.nn1)then

c                  u_interp3=u(ii,jj+1,n_loc,cur_level)

c                  v_interp3=v(ii,jj+1,n_loc,cur_level)

                  int_count=int_count+1
                endif

                if(level(ii,jj-1).eq.nn1)then

c                  u_interp4=u(ii,jj-1,n_loc,cur_level)

c                  v_interp4=v(ii,jj-1,n_loc,cur_level)

                  int_count=int_count+1
                endif  
      
                if(int_count.gt.0)then

c                  u(ii,jj,n_loc,cur_level)=(u_interp1+u_interp2+
c     -                  u_interp3+u_interp4)/int_count
c
c                  v(ii,jj,n_loc,cur_level)=(v_interp1+v_interp2+
c     -                  v_interp3+v_interp4)/int_count
       
                  level(ii,jj)=nn1+1
                endif
                total_count=total_count+int_count
            endif
          enddo
       enddo
       if(total_count.eq.0)goto 50
      enddo

	endif


50      continue


      do jj=2,endy-1
       do ii=2,endx-1
         if(bl_hor_wall(ii,jj).eq.99.and.
     -	level(ii,jj).eq.0)then
            if(level(ii,jj).lt.50) level(ii,jj)=49
         endif
       enddo
      enddo


      do j=1,endy            
       do i=1,endx
          zeta_iter(i,j)=zeta(i,j,n_loc,cur_level)
          u_iter(i,j)=u(i,j,n_loc,cur_level)
          v_iter(i,j)=v(i,j,n_loc,cur_level)
       enddo
      enddo

	per_shift=1  !need to change in E_bc as well

	 do jj=1,endy
	   wvmk_loc_i(jj)=overlap+1
	enddo

      do jj=1,endy
       do ii=1,endx
          if(bl_x_wall(ii,jj).eq.1)then
		  if(bc_1.eq.1)then

			u(ii,jj,n_loc,cur_level)=0
			do ii_bc=1,overlap
				zeta(ii-ii_bc,jj,n_loc,cur_level)=
     -                  zeta(ii+ii_bc,jj,n_loc,cur_level)
				u(ii-ii_bc,jj,n_loc,cur_level)=
     -                  -u(ii+ii_bc,jj,n_loc,cur_level)
				h(ii-ii_bc,jj,n_loc)=h(ii+ii_bc,jj,n_loc)

				if(dim.eq.2)then
					v(ii-ii_bc,jj,n_loc,cur_level)=
     -                  v(ii+ii_bc,jj,n_loc,cur_level)
				endif
			enddo

		  elseif(bc_1.eq.2)then

				call solit(u_2,v_2,zeta_2,c,co,t(nn_cur),x0,x(ii),
     -				y(jj),ho,alp,inc_ang,wave_type,
     -                  wave_hgt,L,bf_ratio,
     -				  bet,cur_level,num_levels)
				
				u(ii,jj,n_loc,cur_level)=u_2

				do ii_bc=1,overlap
					call solit(u_1,v_1,zeta_1,c,co,t(nn_cur),x0,
     -					x(ii-ii_bc),y(jj),ho,alp,inc_ang,wave_type,
     -					wave_hgt,L,bf_ratio,
     -					bet,cur_level,num_levels)

					zeta(ii-ii_bc,jj,n_loc,cur_level)=zeta_1
					u(ii-ii_bc,jj,n_loc,cur_level)=u_1
					if(dim.eq.2)then
						v(ii-ii_bc,jj,n_loc,cur_level)=v_1
					endif
				enddo
		  elseif(bc_1.eq.4)then

			u(ii,jj,n_loc,cur_level)=
     -                 u(endx-overlap-per_shift,jj,n_loc,cur_level)
			do ii_bc=1,overlap
		     zeta(ii-ii_bc,jj,n_loc,cur_level)=
     -                 zeta(endx-overlap-per_shift-ii_bc,jj,
     -                 n_loc,cur_level)
				u(ii-ii_bc,jj,n_loc,cur_level)=
     -                    u(endx-overlap-per_shift-ii_bc,jj,
     -                 n_loc,cur_level)
				h(ii-ii_bc,jj,n_loc)=
     -					h(endx-overlap-per_shift-ii_bc,jj,n_loc)

				if(dim.eq.2)then
					v(ii-ii_bc,jj,n_loc,cur_level)=
     -                    v(endx-overlap-per_shift-ii_bc,jj,
     -                 n_loc,cur_level)
				endif
			enddo


		  endif
          elseif(bl_x_wall(ii,jj).eq.2)then
		  if(bc_2.eq.1)then

			u(ii,jj,n_loc,cur_level)=0

			do ii_bc=1,overlap
				zeta(ii+ii_bc,jj,n_loc,cur_level)=
     -                  zeta(ii-ii_bc,jj,n_loc,cur_level)
				u(ii+ii_bc,jj,n_loc,cur_level)=
     -                  -u(ii-ii_bc,jj,n_loc,cur_level)
				h(ii+ii_bc,jj,n_loc)=h(ii-ii_bc,jj,n_loc)

				if(dim.eq.2)then
					v(ii+ii_bc,jj,n_loc,cur_level)=
     -                  v(ii-ii_bc,jj,n_loc,cur_level)
				endif
			enddo
		  elseif(bc_2.eq.4)then

			u(ii,jj,n_loc,cur_level)=
     -                 u(overlap+1+per_shift,jj,n_loc,cur_level)

			do ii_bc=1,overlap
		     zeta(ii+ii_bc,jj,n_loc,cur_level)=
     -                  zeta(overlap+1+per_shift+ii_bc,jj,
     -                 n_loc,cur_level)
				u(ii+ii_bc,jj,n_loc,cur_level)=
     -                     u(overlap+1+per_shift+ii_bc,jj,
     -                 n_loc,cur_level)
				h(ii+ii_bc,jj,n_loc)=
     -					 h(overlap+1+per_shift+ii_bc,jj,n_loc)

				if(dim.eq.2)then
					v(ii+ii_bc,jj,n_loc,cur_level)=
     -                  v(overlap+1+per_shift+ii_bc,jj,n_loc,cur_level)
				endif
			enddo

		  endif
          endif

          if(dim.eq.2)then
			if(bl_y_wall(ii,jj).eq.1)then
			  if(bc_3.eq.4)then

				v(ii,jj,n_loc,cur_level)=
     -                 v(ii,endy-overlap-per_shift,n_loc,cur_level)
				
				do jj_bc=1,overlap
				  zeta(ii,jj-jj_bc,n_loc,cur_level)=
     -                 zeta(ii,endy-overlap-per_shift-jj_bc,
     -                 n_loc,cur_level)
				     u(ii,jj-jj_bc,n_loc,cur_level)=
     -                    u(ii,endy-overlap-per_shift-jj_bc,
     -                 n_loc,cur_level)
				     h(ii,jj-jj_bc,n_loc)=
     -					h(ii,endy-overlap-per_shift-jj_bc,n_loc)
				     v(ii,jj-jj_bc,n_loc,cur_level)=
     -                    v(ii,endy-overlap-per_shift-jj_bc,
     -                 n_loc,cur_level)
				enddo

			  else

				v(ii,jj,n_loc,cur_level)=0

				do jj_bc=1,overlap  
					zeta(ii,jj-jj_bc,n_loc,cur_level)=
     -                  zeta(ii,jj+jj_bc,n_loc,cur_level)
					v(ii,jj-jj_bc,n_loc,cur_level)=
     -                  -v(ii,jj+jj_bc,n_loc,cur_level)
					u(ii,jj-jj_bc,n_loc,cur_level)=
     -					u(ii,jj+jj_bc,n_loc,cur_level)
					h(ii,jj-jj_bc,n_loc)=h(ii,jj+jj_bc,n_loc)
				enddo
			  endif

			elseif(bl_y_wall(ii,jj).eq.2)then
			  if(bc_4.eq.4)then

				v(ii,jj,n_loc,cur_level)=
     -                 v(ii,overlap+1+per_shift,n_loc,cur_level)
				
				do jj_bc=1,overlap
				  zeta(ii,jj+jj_bc,n_loc,cur_level)=
     -                 zeta(ii,overlap+1+per_shift+jj_bc,
     -                 n_loc,cur_level)
				     u(ii,jj+jj_bc,n_loc,cur_level)=
     -                    u(ii,overlap+1+per_shift+jj_bc,
     -                 n_loc,cur_level)
				     h(ii,jj+jj_bc,n_loc)=
     -					h(ii,overlap+1+per_shift+jj_bc,n_loc)
				     v(ii,jj+jj_bc,n_loc,cur_level)=
     -                    v(ii,overlap+1+per_shift+jj_bc,
     -                 n_loc,cur_level)
				enddo

			  else

				v(ii,jj,n_loc,cur_level)=0

				do jj_bc=1,overlap  
					zeta(ii,jj+jj_bc,n_loc,cur_level)=
     -					zeta(ii,jj-jj_bc,n_loc,cur_level)
					v(ii,jj+jj_bc,n_loc,cur_level)=
     -                  -v(ii,jj-jj_bc,n_loc,cur_level)
					u(ii,jj+jj_bc,n_loc,cur_level)=
     -                  u(ii,jj-jj_bc,n_loc,cur_level)
					h(ii,jj+jj_bc,n_loc)=h(ii,jj-jj_bc,n_loc)
				enddo
			  endif


			endif
		endif

       enddo
      enddo

C FIX CORNERS
      if(dim.eq.2)then
		jj=1+overlap
		if(bc_1.eq.1.and.bc_3.eq.1)then 
		  do ii=1,wvmk_loc_i(jj)-1
c			v(ii,jj,n_loc,cur_level)=0
			do jj_bc=1,overlap  
				zeta(ii,jj-jj_bc,n_loc,cur_level)=
     -                  zeta(ii,jj+jj_bc,n_loc,cur_level)
				v(ii,jj-jj_bc,n_loc,cur_level)=
     -                  -v(ii,jj+jj_bc,n_loc,cur_level)
				u(ii,jj-jj_bc,n_loc,cur_level)=
     -					u(ii,jj+jj_bc,n_loc,cur_level)
				h(ii,jj-jj_bc,n_loc)=h(ii,jj+jj_bc,n_loc)
			enddo				 	
		  enddo
		endif

		if(bc_2.eq.1.and.bc_3.eq.1)then 
		  do ii=endx-overlap,endx
c			v(ii,jj,n_loc,cur_level)=0
			do jj_bc=1,overlap  
				zeta(ii,jj-jj_bc,n_loc,cur_level)=
     -                  zeta(ii,jj+jj_bc,n_loc,cur_level)
				v(ii,jj-jj_bc,n_loc,cur_level)=
     -                  -v(ii,jj+jj_bc,n_loc,cur_level)
				u(ii,jj-jj_bc,n_loc,cur_level)=
     -					u(ii,jj+jj_bc,n_loc,cur_level)
				h(ii,jj-jj_bc,n_loc)=h(ii,jj+jj_bc,n_loc)
			enddo				 	
		  enddo
		endif

		jj=endy-overlap
		if(bc_1.eq.1.and.bc_4.eq.1)then 
		  do ii=1,wvmk_loc_i(jj)-1
c			v(ii,jj,n_loc,cur_level)=0
			do jj_bc=1,overlap  
				zeta(ii,jj+jj_bc,n_loc,cur_level)=
     -					zeta(ii,jj-jj_bc,n_loc,cur_level)
				v(ii,jj+jj_bc,n_loc,cur_level)=
     -                  -v(ii,jj-jj_bc,n_loc,cur_level)
				u(ii,jj+jj_bc,n_loc,cur_level)=
     -                  u(ii,jj-jj_bc,n_loc,cur_level)
				h(ii,jj+jj_bc,n_loc)=h(ii,jj-jj_bc,n_loc)
			enddo				 	
		  enddo
		endif

		if(bc_2.eq.1.and.bc_4.eq.1)then 
		  do ii=endx-overlap,endx
c			v(ii,jj,n_loc,cur_level)=0
			do jj_bc=1,overlap  
				zeta(ii,jj+jj_bc,n_loc,cur_level)=
     -					zeta(ii,jj-jj_bc,n_loc,cur_level)
				v(ii,jj+jj_bc,n_loc,cur_level)=
     -                  -v(ii,jj-jj_bc,n_loc,cur_level)
				u(ii,jj+jj_bc,n_loc,cur_level)=
     -                  u(ii,jj-jj_bc,n_loc,cur_level)
				h(ii,jj+jj_bc,n_loc)=h(ii,jj-jj_bc,n_loc)
			enddo				 	
		  enddo
		endif

      endif

      do j=1,endy            
       do i=1,endx
          zeta_iter(i,j)=zeta(i,j,n_loc,cur_level)
          u_iter(i,j)=u(i,j,n_loc,cur_level)
          v_iter(i,j)=v(i,j,n_loc,cur_level)
       enddo
      enddo

      if(dim.eq.3)then
        jj=1+overlap
        do ii=3,endx-2
            if(bl_hor_wall(ii,jj).eq.99.and.
     -            level(ii,jj).lt.40)then
                u(ii,jj,n_loc,cur_level)=1./6.*
     -                (u_iter(ii-2,jj)+2.*u_iter(ii-1,jj)+
     -              2.*u_iter(ii+1,jj)+u_iter(ii+2,jj)+0.*u_iter(ii,jj))

                zeta(ii,jj,n_loc,cur_level)=
     -                  1./6.*(zeta_iter(ii-2,jj)+
     -                  2.*zeta_iter(ii-1,jj)+
     -                  2.*zeta_iter(ii+1,jj)+
     -                  zeta_iter(ii+2,jj)+0.*zeta_iter(ii,jj))
            endif
       enddo

      elseif(dim.eq.4)then
        do jj=3,endy-2
          do ii=3,endx-2

            if(bl_hor_wall(ii,jj).eq.99.and.
     -            level(ii,jj).le.40)then

                i=ii
                j=jj

                u(ii,jj,n_loc,cur_level)=1./24.*(
     -                u_iter(i-2,j)+2.*u_iter(i-1,j)+
     -                u_iter(i+2,j)+2.*u_iter(i+1,j)+
     -                u_iter(i,j-2)+2.*u_iter(i,j-1)+
     -                u_iter(i,j+2)+2.*u_iter(i,j+1)+
     -                u_iter(i-2,j-2)+2.*u_iter(i-1,j-1)+
     -                u_iter(i+2,j+2)+2.*u_iter(i+1,j+1)+
     -                u_iter(i+2,j-2)+2.*u_iter(i+1,j-1)+
     -                u_iter(i-2,j+2)+2.*u_iter(i-1,j+1))


                v(ii,jj,n_loc,cur_level)=1./24.*(
     -                v_iter(i-2,j)+2.*v_iter(i-1,j)+
     -                v_iter(i+2,j)+2.*v_iter(i+1,j)+
     -                v_iter(i,j-2)+2.*v_iter(i,j-1)+
     -                v_iter(i,j+2)+2.*v_iter(i,j+1)+
     -                v_iter(i-2,j-2)+2.*v_iter(i-1,j-1)+
     -                v_iter(i+2,j+2)+2.*v_iter(i+1,j+1)+
     -                v_iter(i+2,j-2)+2.*v_iter(i+1,j-1)+
     -                v_iter(i-2,j+2)+2.*v_iter(i-1,j+1))


                zeta(ii,jj,n_loc,cur_level)=1./24.*(
     -                zeta_iter(i-2,j)+2.*zeta_iter(i-1,j)+
     -                zeta_iter(i+2,j)+2.*zeta_iter(i+1,j)+
     -                zeta_iter(i,j-2)+2.*zeta_iter(i,j-1)+
     -                zeta_iter(i,j+2)+2.*zeta_iter(i,j+1)+
     -                zeta_iter(i-2,j-2)+2.*zeta_iter(i-1,j-1)+
     -                zeta_iter(i+2,j+2)+2.*zeta_iter(i+1,j+1)+
     -                zeta_iter(i+2,j-2)+2.*zeta_iter(i+1,j-1)+
     -                zeta_iter(i-2,j+2)+2.*zeta_iter(i-1,j+1))

            endif
          enddo

       enddo
      endif      

c      do ii=2,endx-1
c        do jj=2,endy-1
c          if(bl_hor_wall(ii,jj).eq.99.and.level(ii,jj).ge.2
c     -          .and.zeta(ii,jj,n_loc,cur_level)+
c     -            h(ii,jj,n_loc).ge.cutoff_mat(i,j))then
c                zeta(ii,jj,n_loc,cur_level)=
c     -                -h(ii,jj,n_loc)          
c          endif
c       enddo
c      enddo

                                                            
      return
      end      ! end subroutine
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine bc_wvmk(n_loc,nn_cur)
      use mainvar_module
      integer n_loc,nn_cur,nn1,nn_max,per_shift,ii_bc,jj_bc,n_wv,wi,n_wi
	real u_interp1,u_interp2,u_interp3,u_interp4,
     -            u_interp5,u_interp6,u_interp7,u_interp8,
     -            v_interp1,v_interp2,v_interp3,v_interp4,
     -            v_interp5,v_interp6,v_interp7,v_interp8,
     -	u_1,v_1,zeta_1,u_2,v_2,zeta_2,d_interp,H_interp,wm_y(30),
     -   x_i,ui,hi,flux,wm_d_m,wm_d_p,dx1(endy),dx2(endy),u_int
      
      nn_max=2*overlap

! mod for wavemaker

      if(nn_cur.eq.1)then
        ALLOCATE(wm_d_c(endy),wm_u_c(endy))
		if(dim.eq.1)then
            open(888,file='wv_signal.dat',status='old')
            read(888,*) tmp1, tmp2
            n_wv=real(tmp1)
            
            ALLOCATE(wv_signal(n_wv,2))
            ALLOCATE(wm_d(tt,1))
            do i=1,n_wv
                read(888,*) wv_signal(i,1), wv_signal(i,2)
            enddo
            close(888)
		
		    do ii=1,tt
		        wm_d(ii,1)=wv_signal(n_wv,2)
			    do jj=1,n_wv-1
				  if(t(ii).ge.wv_signal(jj,1).and.
     -			   t(ii).lt.wv_signal(jj+1,1))then
					wm_d(ii,1)=(wv_signal(jj+1,2)-wv_signal(jj,2))/
     -					(wv_signal(jj+1,1)-wv_signal(jj,1))*
     -					(t(ii)-wv_signal(jj,1))+wv_signal(jj,2)
     
                    goto 77
			 	  endif
			    enddo
77              continue
		    enddo
		elseif(dim.eq.2)then
            open(888,file='wv_signal.dat',status='old')
            read(888,*) tmp1, tmp2, tmp2, tmp2, tmp2, 
     -			  tmp2, tmp2, tmp2, tmp2, tmp2, 
     -			  tmp2, tmp2, tmp2, tmp2, tmp2, tmp2,
     -			  tmp2, tmp2, tmp2, tmp2, tmp2, 
     -			  tmp2, tmp2, tmp2, tmp2, tmp2, tmp2, 
     -			  tmp2, tmp2, tmp2, tmp2
            n_wv=real(tmp1)
            
            ALLOCATE(wv_signal(n_wv,31))
            ALLOCATE(wm_d(tt,30))
            do i=1,n_wv
                read(888,*) wv_signal(i,1), wv_signal(i,2), 
     -			  wv_signal(i,3), wv_signal(i,4), 
     -			  wv_signal(i,5), wv_signal(i,6), wv_signal(i,7), 
     -			  wv_signal(i,8), wv_signal(i,9), 
     -			  wv_signal(i,10), wv_signal(i,11), wv_signal(i,12), 
     -			  wv_signal(i,13), wv_signal(i,14), 
     -			  wv_signal(i,15), wv_signal(i,16), wv_signal(i,17), 
     -			  wv_signal(i,18), wv_signal(i,19), 
     -			  wv_signal(i,20), wv_signal(i,21), wv_signal(i,22), 
     -			  wv_signal(i,23), wv_signal(i,24), 
     -			  wv_signal(i,25), wv_signal(i,26), wv_signal(i,27), 
     -			  wv_signal(i,28), wv_signal(i,29), 
     -			  wv_signal(i,30), wv_signal(i,31)
            enddo
            close(888)
		
			do wi=1,30
		      do ii=1,tt
		        wm_d(ii,wi)=wv_signal(n_wv,wi+1)
			    do jj=1,n_wv-1
				  if(t(ii).ge.wv_signal(jj,1).and.
     -			   t(ii).lt.wv_signal(jj+1,1))then
					wm_d(ii,wi)=(wv_signal(jj+1,wi+1)-wv_signal(jj,wi+1))/
     -					(wv_signal(jj+1,1)-wv_signal(jj,1))*
     -					(t(ii)-wv_signal(jj,1))+wv_signal(jj,wi+1)
     
                    goto 78
			 	  endif
			    enddo
78              continue
		      enddo
		   enddo
		endif
		Deallocate(wv_signal)            
      endif      

	if(dim.eq.1)then
		n_wi=1
	else
		n_wi=30
	endif
	
	do wi=1,n_wi
		wm_y(wi)=real(wi-1)*26.5/29.
	enddo

	do jj=1,endy
	  do wi=1,n_wi
	   if(dim.eq.1)then
		wm_d_c(jj)=wm_d(nn_cur,wi)
		wm_d_m=wm_d(max(1,nn_cur-1),wi)
		wm_d_p=wm_d(min(tt,nn_cur+1),wi)

		wm_u_c(jj)=(wm_d_p-wm_d_m)/(2.*dt)
	   else
		call interp1D_single(wm_y,wm_d(nn_cur,:),
     -					y(jj),wm_d_c(jj),n_wi,1)	
		call interp1D_single(wm_y,wm_d(max(1,nn_cur-1),:),
     -					y(jj),wm_d_m,n_wi,1)	
		call interp1D_single(wm_y,wm_d(min(tt,nn_cur+1),:),
     -					y(jj),wm_d_p,n_wi,1)	

		wm_u_c(jj)=(wm_d_p-wm_d_m)/(2.*dt)
	   endif
	  enddo
	 enddo

	  if(dim.eq.2)then
	   do jj=1,overlap+1
		 wm_d_c(jj)=wm_d_c(overlap+2)
		 wm_u_c(jj)=wm_u_c(overlap+2)
	   enddo 

	   do jj=endy-overlap,endy
	 	 wm_d_c(jj)=wm_d_c(endy-overlap-1)
	 	 wm_u_c(jj)=wm_u_c(endy-overlap-1)
	   enddo 
      endif

      do jj=1,endy
	   wvmk_loc_i(jj)=overlap+1
       do ii=5,endx-4
			if(x(ii).lt.wm_d_c(jj).and.x(ii+1).ge.wm_d_c(jj))then
				bl_x_wall(ii-1,jj)=0 
				bl_x_wall(ii,jj)=1
				bl_x_wall(ii+1,jj)=0
				wvmk_loc_i(jj)=ii+1
				dx1(jj)=wm_d_c(jj)-x(ii)
				dx2(jj)=x(ii+1)-wm_d_c(jj)
			endif
        enddo
      enddo 
	  
      do jj=1,endy
       do ii=1,endx
          if(bl_x_wall(ii,jj).eq.1)then
		   if(bc_1.eq.1)then
			imag_leng=ii-1
			do ii_bc=0,imag_leng
                	
                if(ii_bc.eq.0)then
!                    x_i=wm_d_c+dx1+dx*ii_bc
!                    call interp1D_single(x,
!     -                   zeta(:,jj,n_loc,cur_level),x_i,u_1,endx,1)
!                    zeta(ii-ii_bc,jj,n_loc,cur_level)=u_1   
                    zeta(ii-ii_bc,jj,n_loc,cur_level)=
     -                   zeta(ii-ii_bc+1,jj,n_loc,cur_level)
 !                   u(ii-ii_bc,jj,n_loc,cur_level)=wm_u_c
                    if(dx1(jj).le.dx/2.)then
                        u_1=wm_u_c(jj)
                        u_2=u(ii+1,jj,n_loc,cur_level)
                    
                        u_int=u_1+dx1(jj)*(u_2-u_1)/(dx2(jj))
                        u(ii-ii_bc,jj,n_loc,cur_level)=
     -			  2.*wm_u_c(jj)-u_int
                    else
                        u_1=wm_u_c(jj)
                        u_2=u(ii+2,jj,n_loc,cur_level)
                    
                        u_int=u_1+dx1(jj)*(u_2-u_1)/(dx+dx2(jj))
     
                        u(ii-ii_bc,jj,n_loc,cur_level)=
     -			  2.*wm_u_c(jj)-u_int
                    endif   
                else
                    x_i=wm_d_c(jj)+dx1(jj)+dx*ii_bc
                    call interp1D_single(x,
     -                   u(:,jj,n_loc,cur_level),x_i,u_1,endx,1)
                    u(ii-ii_bc,jj,n_loc,cur_level)=2.*wm_u_c(jj)-u_1 
                     
                    call interp1D_single(x,
     -                   zeta(:,jj,n_loc,cur_level),x_i,u_1,endx,1)
                    zeta(ii-ii_bc,jj,n_loc,cur_level)=u_1
 
                endif
                    
				h(ii-ii_bc,jj,n_loc)=h(ii+ii_bc+1,jj,n_loc)
                UU(ii-ii_bc,jj,n_loc,cur_level)=
     -                   u(ii-ii_bc,jj,n_loc,cur_level)*
     -                   (zeta(ii-ii_bc,jj,n_loc,cur_level)+
     -                   h(ii-ii_bc,jj,n_loc))
     
				if(dim.eq.2)then
					v(ii-ii_bc,jj,n_loc,cur_level)=
     -                  v(ii+ii_bc,jj,n_loc,cur_level)
				endif
			enddo
		   endif
		  endif
		 enddo
		enddo
		  
		jj=endy-overlap
		if(bc_1.eq.1.and.bc_4.eq.1)then 
		  do ii=1,wvmk_loc_i(jj)-1
c			v(ii,jj,n_loc,cur_level)=0
			do jj_bc=1,overlap  
				zeta(ii,jj+jj_bc,n_loc,cur_level)=
     -					zeta(ii,jj-jj_bc,n_loc,cur_level)
				v(ii,jj+jj_bc,n_loc,cur_level)=
     -                  -v(ii,jj-jj_bc,n_loc,cur_level)
				u(ii,jj+jj_bc,n_loc,cur_level)=
     -                  u(ii,jj-jj_bc,n_loc,cur_level)
				h(ii,jj+jj_bc,n_loc)=h(ii,jj-jj_bc,n_loc)
			enddo				 	
		  enddo
		endif


      do j=1,endy            
       do i=1,endx
          zeta_iter(i,j)=zeta(i,j,n_loc,cur_level)
          u_iter(i,j)=u(i,j,n_loc,cur_level)
          v_iter(i,j)=v(i,j,n_loc,cur_level)
       enddo
      enddo
                                                 
      return
      end      ! end subroutine
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      subroutine bc_E(n_loc)
      use mainvar_module
      integer n_loc,per_shift,ii_bc,jj_bc


	per_shift=1
      do jj=1+overlap,endy-overlap
       do ii=1+overlap,endx-overlap
          if(bl_x_wall(ii,jj).eq.1)then
		  
			if(bc_1.eq.4)then
				do ii_bc=1,overlap
				E(ii-ii_bc,jj,n_loc,cur_level)=
     -                 E(endx-overlap-per_shift-ii_bc,jj,
     -                 n_loc,cur_level)
				enddo
			else
				do ii_bc=1,overlap
					E(ii-ii_bc,jj,n_loc,cur_level)=
     -                  E(ii+ii_bc,jj,n_loc,cur_level)
				enddo
			endif

          elseif(bl_x_wall(ii,jj).eq.2)then

			if(bc_2.eq.4)then
				do ii_bc=1,overlap
				E(ii+ii_bc,jj,n_loc,cur_level)=
     -                  E(overlap+1+per_shift+ii_bc,jj,
     -                 n_loc,cur_level)
				enddo
			else
				do ii_bc=1,overlap
					E(ii+ii_bc,jj,n_loc,cur_level)=
     -                  E(ii-ii_bc,jj,n_loc,cur_level)
				enddo
			endif
          endif

          if(dim.eq.2)then
			if(bl_y_wall(ii,jj).eq.1)then
				do jj_bc=1,overlap  
					E(ii,jj-jj_bc,n_loc,cur_level)=
     -                  E(ii,jj+jj_bc,n_loc,cur_level)
				enddo

			elseif(bl_y_wall(ii,jj).eq.2)then
				do jj_bc=1,overlap  
					E(ii,jj+jj_bc,n_loc,cur_level)=
     -					E(ii,jj-jj_bc,n_loc,cur_level)
				enddo
			endif
		endif

       enddo
      enddo

C FIX CORNERS
      if(dim.eq.2)then
		jj=1+overlap
		if(bc_1.eq.1.and.bc_3.eq.1)then 
		  do ii=1,overlap
			do jj_bc=1,overlap  
				E(ii,jj-jj_bc,n_loc,cur_level)=
     -                  E(ii,jj+jj_bc,n_loc,cur_level)
			enddo				 	
		  enddo
		endif

		if(bc_2.eq.1.and.bc_3.eq.1)then 
		  do ii=endx-overlap,endx
			do jj_bc=1,overlap  
				E(ii,jj-jj_bc,n_loc,cur_level)=
     -                  E(ii,jj+jj_bc,n_loc,cur_level)
			enddo				 	
		  enddo
		endif

		jj=endy-overlap
		if(bc_1.eq.1.and.bc_4.eq.1)then 
		  do ii=1,overlap
			do jj_bc=1,overlap  
				E(ii,jj+jj_bc,n_loc,cur_level)=
     -					E(ii,jj-jj_bc,n_loc,cur_level)
		    enddo
	      enddo
		endif

		if(bc_2.eq.1.and.bc_4.eq.1)then 
		  do ii=endx-overlap,endx
			do jj_bc=1,overlap  
				E(ii,jj+jj_bc,n_loc,cur_level)=
     -					E(ii,jj-jj_bc,n_loc,cur_level)
			enddo
		  enddo
		endif

      endif

                                                            
      return
      end      ! end subroutine
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC












    