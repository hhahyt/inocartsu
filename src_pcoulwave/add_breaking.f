C*******************************************************************************
C........    Groups the physical variables into the 
C........... conservative forms given in Wei & Kirby
      subroutine add_breaking(n_loc,nn_loc)

      use mainvar_module
      integer n_loc,nn_loc
	real mov_fac,x_mov,y_mov,dzdt_max2,dzdt_max3,
     -	d_b_I,d_b_F,z1,dzdx_c,dzdy_c

      cur_level=1

C%%%%%%%%%%%%%%%%%%%% KENNEDY et al TYPE WAVE BREAKING %%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if(wave_breaking.eq.1)then
	if(n_loc.eq.3)then  ! allow breaking front to move max of one grid point per time step
      do j=js1,je1
        do i=1,endx
          if(nn_loc.le.4)then
            roller(i,j)=0
          endif
          if(bl_hor_wall(i,j).eq.99)then
            t_breaking(i,j)=0
          endif
	  tmp_mat(i,j)=t_breaking(i,j)
	  enddo
      enddo
	
	mov_fac=0.0

      do j=js2,je2
        do i=2,endx-1
           if(t_breaking(i,j).gt.dt.and.
     -              t(nn_loc)-t_breaking(i,j).gt.dt/10.)then

			x_mov=u(i,j,3,1)*dt/dx

			if(dim.eq.1)then
				if(x_mov.gt.mov_fac)then
					if(t_breaking(i+1,j).lt.dt)then
						tmp_mat(i+1,j)=
     -						t_breaking(i,j)
					endif
				elseif(x_mov.lt.-mov_fac)then
					if(t_breaking(i-1,j).lt.dt)then
						tmp_mat(i-1,j)=
     -						t_breaking(i,j)
					endif
				endif
			elseif(dim.eq.2)then
				y_mov=v(i,j,3,1)*dt/dy

				if(abs(x_mov).ge.abs(y_mov))then
					if(x_mov.gt.mov_fac)then
						if(t_breaking(i+1,j).lt.dt)then
							tmp_mat(i+1,j)=
     -						t_breaking(i,j)
						endif
						if(y_mov.gt.mov_fac)then
							if(t_breaking(i+1,j+1).lt.dt)then
								tmp_mat(i+1,j+1)=
     -								t_breaking(i,j)
							endif
						elseif(y_mov.lt.-mov_fac)then
							if(t_breaking(i+1,j-1).lt.dt)then
								tmp_mat(i+1,j-1)=
     -								t_breaking(i,j)
							endif
						endif
					elseif(x_mov.lt.-mov_fac)then
						if(t_breaking(i-1,j).lt.dt)then
							tmp_mat(i-1,j)=
     -							t_breaking(i,j)
						endif
						if(y_mov.gt.mov_fac)then
							if(t_breaking(i-1,j+1).lt.dt)then
								tmp_mat(i-1,j+1)=
     -								t_breaking(i,j)
							endif
						elseif(y_mov.lt.-mov_fac)then
							if(t_breaking(i-1,j-1).lt.dt)then
								tmp_mat(i-1,j-1)=
     -								t_breaking(i,j)
							endif
						endif
					endif 
				else
					if(y_mov.gt.mov_fac)then
						if(t_breaking(i,j+1).lt.dt)then
							tmp_mat(i,j+1)=
     -							t_breaking(i,j)
						endif
						if(x_mov.gt.mov_fac)then
							if(t_breaking(i+1,j+1).lt.dt)then
								tmp_mat(i+1,j+1)=
     -								t_breaking(i,j)
							endif
						elseif(x_mov.lt.mov_fac)then
							if(t_breaking(i-1,j+1).lt.dt)then
								tmp_mat(i-1,j+1)=
     -								t_breaking(i,j)
							endif
						endif
					elseif(y_mov.lt.mov_fac)then
						if(t_breaking(i,j-1).lt.dt)then
							tmp_mat(i,j-1)=
     -							t_breaking(i,j)
						endif
						if(x_mov.gt.mov_fac)then
							if(t_breaking(i+1,j-1).lt.dt)then
								tmp_mat(i+1,j-1)=
     -								t_breaking(i,j)
							endif
						elseif(x_mov.lt.mov_fac)then
							if(t_breaking(i-1,j-1).lt.dt)then
								tmp_mat(i-1,j-1)=
     -								t_breaking(i,j)
							endif
						endif
					endif
				endif
		    endif
		 endif
	   enddo
	enddo

      do j=js1,je1
        do i=1,endx
		 t_breaking(i,j)=tmp_mat(i,j)
	  enddo
      enddo	

      do j=js1,je1
        do i=1,endx
		 tmp_mat(i,j)=0.
		 eddy_visc(i,j)=0.
	  enddo
      enddo

	dzdt_max2=0.
	dzdt_max3=0.
      do j=js2,je2
        do i=2,endx-1
          if(bl_hor_wall(i,j).le.50)then
		  tmp=0

		  z=zeta(i,j,n_loc,1)

		  if(h(i,j,n_loc)+z.lt.cutoff_mat(i,j))then
			H_total=cutoff_mat(i,j)
		  else
			H_total=h(i,j,n_loc)+z
		  endif

		  co_c=sqrt(9.81*H_total)

		  dzdt=E(i,j,n_loc,1)

		  if(num_levels.eq.1)then
			  delta_breaking=6.5  ! FOR FD: 6.5 - 1L   !10 - 2L
			  T_star=7.*sqrt(abs(H_total/grav))  ! FOR FD: 7 - 1L  !10 - 2L
			  dzdt_I=0.55*co_c   ! FOR FD: 65 - 1L   !50 - 2L
			  dzdt_F=0.05*co_c   ! FOR FD: 08 - 1L   !05 - 2L

		  else
			  delta_breaking=10.  !6.5 - 1L   !10 - 2L
			  T_star=10.*sqrt(abs(H_total/grav))  !7 - 1L  !10 - 2L
			  dzdt_I=0.50*co_c   !65 - 1L   !50 - 2L
			  dzdt_F=0.05*co_c   !08 - 1L   !05 - 2L
		  endif
		  
		  dzdt_star=0.
		  d_b_I=delta_breaking
		  d_b_F=delta_breaking
            
		  if(t_breaking(i,j).lt.dt)then
              dzdt_star=dzdt_I
			delta_breaking=d_b_I
            elseif(t(nn_loc)-t_breaking(i,j).lt.T_star)then
              dzdt_star=dzdt_I+
     -              (t(nn_loc)-t_breaking(i,j))/
     -              T_star*(dzdt_F-dzdt_I)
              delta_breaking=d_b_I+
     -              (t(nn_loc)-t_breaking(i,j))/
     -              T_star*(d_b_F-d_b_I)
            elseif(t(nn_loc)-t_breaking(i,j).ge.T_star)then
              dzdt_star=dzdt_F
			delta_breaking=d_b_F
            endif

		  tmp=1.
		  if(int_src.ge.1)then
			x_c=x(i)-x0*cosA(1,1)
			y_c=y(j)-x0*sinA(1,1)

			loc=abs(x_c*cosA(1,1)+y_c*sinA(1,1))

			tmp=min(1.,loc/L)		  
		  endif

		  z=zeta(i,j,n_loc,1)
		  z1=z_alp(i,j,n_loc,1)
		  hc=h(i,j,n_loc)

		  if(dzdt.le.dzdt_star)then
			tmp_mat(i,j)=0.
			t_breaking(i,j)=0.
		  else
              if(t_breaking(i,j).lt.dt)then
                t_breaking(i,j)=t(nn_loc)
c				print*,'Breaking at x,y= ',x(i),y(j)
              endif

              if(dzdt.ge.2.*dzdt_star)then
                tmp_mat(i,j)=1.*delta_breaking
              else
                tmp_mat(i,j)=(dzdt/dzdt_star-1.)*delta_breaking
              endif

            endif

          endif
        enddo
      enddo


      do j=js2,je2
        do i=3,endx-2
          if(bl_hor_wall(i,j).le.50)then
		  if(h(i,j,n_loc)+zeta(i,j,n_loc,1).lt.cutoff_mat(i,j))then
			H_total=cutoff_mat(i,j)
		  else
			H_total=h(i,j,n_loc)+zeta(i,j,n_loc,1)
		  endif

		  B=tmp_mat(i,j)
		  dzdt=max(0.,E(i,j,n_loc,cur_level))

          eddy_visc(i,j)=abs(B*H_total*dzdt)
          
		  z=zeta(i,j,n_loc,1)
		  hc=h(i,j,n_loc)

c		  if(B.gt.0.0001)then  ! Breaker Velocity modification of Lynett, 2007, CE 
c
c			co_c=sqrt(9.81*H_total)
c			z1=z_alp(i,j,n_loc,1)
c
c			u_s=u(i,j,n_loc,1)
c     -			-0.5*(z**2.-z1**2.)*dSdx(i,j,1)-
c     -            (z-z1)*dTdx(i,j,1)  
c
c
c			if(dim.eq.2)then
c				v_s=v(i,j,n_loc,1)
c     -				-0.5*(z**2.-z1**2.)*dSdy(i,j,1)-
c     -				(z-z1)*dTdy(i,j,1) 
c			else
c				v_s=0.
c			endif
c
c			k_br(i,j)=5./H_total
c
c			z_br(i,j)=-hc
c     -				+log(exp(H_total*k_br(i,j))/H_total/k_br(i,j)
c     -				-1./H_total/k_br(i,j))/k_br(i,j)
c
c			u0_b(i,j,1)=(-dzdx(i,j,cur_level)/
c     -			sqrt(dzdx(i,j,cur_level)**2.+dzdy(i,j,cur_level)**2.)*
c     -					co_c-u_s)/(1.-exp(k_br(i,j)*(z_br(i,j)-z))) 
c
c			u0_b(i,j,2)=(-dzdy(i,j,cur_level)/
c     -			sqrt(dzdx(i,j,cur_level)**2.+dzdy(i,j,cur_level)**2.)*
c     -					co_c-v_s)/(1.-exp(k_br(i,j)*(z_br(i,j)-z)))  
c
c		  else
c			u0_b(i,j,1)=0.
c			u0_b(i,j,2)=0.
c			k_br(i,j)=1.
c			z_br(i,j)=0.
c		  endif

          endif
        enddo
      enddo


      do j=js1,je1
        do i=1,endx
		  tmp_mat(i,j)=eddy_visc(i,j)  ! limit for diffusion stability
        enddo
      enddo

	if(dim.eq.1)then
       do j=js2,je2
        do i=3,endx-2
		nut(i,j,n_loc)=0.5*tmp_mat(i,j)
     -		+0.25*tmp_mat(i+1,j)+0.25*tmp_mat(i-1,j)
        enddo
       enddo
	else
       do j=js2,je2
        do i=3,endx-2
		nut(i,j,n_loc)=0.5*tmp_mat(i,j)
     -				+0.125*tmp_mat(i+1,j)+0.125*tmp_mat(i-1,j)
     -				+0.125*tmp_mat(i,j+1)+0.125*tmp_mat(i,j-1)
        enddo
       enddo
	endif

	endif  	! if(n_loc.eq.3)then

      do j=js6,je6
        do i=overlap,endx-overlap+1
	   stress_xx(1,i,j)=nut(i,j,n_loc)*dhzudx(i,j,cur_level)
	   stress_xx(2,i,j)=nut(i,j,n_loc)*(dhzudy(i,j,cur_level)+
     -      dhzvdx(i,j,cur_level))
	   stress_xx(3,i,j)=nut(i,j,n_loc)*dhzvdy(i,j,cur_level)
	  enddo
	enddo


      do j=1+overlap,endy-overlap
       do i=1+overlap,endx-overlap
C%%%%%%%%%%%%%%%%%%%%%%%  WAVE BREAKING AND EDDY VISC DISSIPATION %%%%%%%%%%%%%%%%%%%%%%
        WB_x(i,j)=0.
        WB_y(i,j)=0.
        deriv_1x=0.
        deriv_1y=0.
        deriv_2x=0.
        deriv_2y=0.

	  if(h(i,j,n_loc)+zeta(i,j,n_loc,1).lt.cutoff_mat(i,j))then
          H_total=cutoff_mat(i,j)
        else
          H_total=h(i,j,n_loc)+zeta(i,j,n_loc,1)
        endif

        if(bl_hor_wall(i,j).le.50)then

	    deriv_1x=(stress_xx(1,i+1,j)-
     -            stress_xx(1,i-1,j))/(2.*dx)
	    if(dim.eq.2)then
	      deriv_2y=(stress_xx(2,i+1,j)-
     -              stress_xx(2,i-1,j))/(2.*dx)

	      deriv_1y=(stress_xx(3,i,j+1)-
     -              stress_xx(3,i,j-1))/(2.*dy)

	      deriv_2x=(stress_xx(2,i,j+1)-
     -              stress_xx(2,i,j-1))/(2.*dy)

            WB_y(i,j)=1./H_total*(deriv_1y+
     -            0.5*deriv_2y)
          endif   

          WB_x(i,j)=1./H_total*(deriv_1x+
     -          0.5*deriv_2x)

          param_G=WB_y(i,j)
          param_F=WB_x(i,j)

		tmp=H_total

		F(i,j,n_loc,cur_level)=F(i,j,n_loc,cur_level)+param_F*tmp
		G(i,j,n_loc,cur_level)=G(i,j,n_loc,cur_level)+param_G*tmp

          endif
        enddo
      enddo


      endif           !       if(wave_breaking.eq.1)then


      return 

      end      ! end subroutine
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


