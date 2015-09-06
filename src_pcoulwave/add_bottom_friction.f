      subroutine add_bottom_friction(n_loc,nn_loc)

      use mainvar_module
      integer n_loc,nn_loc
	real ub,vb,vel_BF,f_BF_sw


      if(bottom_fric.eq.1)then

	cur_level=1
      do j=1+overlap,endy-overlap
       do i=1+overlap,endx-overlap

        if(bl_hor_wall(i,j).le.50.or.level(i,j).le.1)then

C%%%%%%%%%%%%%%%%%%%%%%%% BOTTOM FRICTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        BF_x=0
        BF_y=0
        
        if((h(i,j,3)+zeta(i,j,3,1)).lt.cutoff_mat(i,j))then
          H_total=cutoff_mat(i,j)
        else
          H_total=h(i,j,3)+zeta(i,j,3,1)
        endif

		z=-h(i,j,n_loc)
		call calc_vel_z(n_loc,z,ub,vb,tmp)

		vel_BF=sqrt(ub**2.+vb**2.)

c Convert the given Mannings n to a dimensionless friction factor
          f_BF_sw=f_BF**2.*8.*9.81/H_total**0.3333
	    
c		if(z.ge.0.)then
c			if(hx(i,j)*ub.gt.0)then  ! downrush swash special treatment
c				tmp=exp(-50.*min(1.,
c     -			abs(H_total/depth)))*
c     -			abs(hx(i,j))
c				f_BF_sw=max(f_BF_sw,tmp**2.*8.*9.81/H_total**0.3333)
c			endif
c		endif

          BF_x=f_BF_sw/H_total*ub*vel_BF

          if(dim.eq.2)then
c			if(z.ge.0.)then
c				if(hy(i,j)*vb.gt.0)then  ! downrush swash special treatment
c					tmp=exp(-50.*min(1.,
c     -				abs(H_total/depth)))*
c     -				abs(hy(i,j))
c					f_BF_sw=max(f_BF_sw,tmp**2.*8.*9.81/H_total**0.3333)
c				endif
c			endif          
            BF_y=f_BF_sw/H_total*vb*vel_BF
          endif

          
        param_G=-BF_y
        param_F=-BF_x

	  tmp=H_total

        F(i,j,n_loc,cur_level)=F(i,j,n_loc,cur_level)+param_F*tmp
        G(i,j,n_loc,cur_level)=G(i,j,n_loc,cur_level)+param_G*tmp

	  endif
	 enddo
	enddo

	endif

      return
      end      ! end subroutine
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


c	!---------------- add bottom friction ---------------
c
c
c	if(ibotfric.eq.1)then
c	   do j=sy,ey
c		do i=sx,ex
c			 ua = u(i,j) - ((z(i,j)**2-z(i,j)*d(i,j)+d(i,j)**2)/6.0
c     -				-0.5*za(i,j)**2)*(uxx(i,j)+vxy(i,j)) 
c     -				- (0.5*(z(i,j)-d(i,j))-za(i,j))*
c     -				(huxx(i,j)+hvxy(i,j))
c			 va = v(i,j) - ((z(i,j)**2-z(i,j)*d(i,j)+d(i,j)**2)/6.0
c     -				-0.5*za(i,j)**2)*(uxy(i,j)+vyy(i,j)) 
c     -				- (0.5*(z(i,j)-d(i,j))-za(i,j))*
c     -				(huxy(i,j)+hvyy(i,j))
c		   Rh = max(100.0,sqrt(ua**2+va**2)*
c     -				abs(z(i,j)+d(i,j)-0.000)*1000000.)
c			 Cf = -1.8*log10(6.9/Rh+(ks/H(i,j)/3.7)**1.11)
c			 Cf = 0.25/Cf**2
c			 if (Rh.le.411.58) Cf = 0.25*64./Rh
c		   Rhx = Cf*ua*sqrt(va**2+ua**2)
c		   Rhy = Cf*va*sqrt(va**2+ua**2)
c		   F(i,j) = F(i,j)- Rhx
c		   G(i,j) = G(i,j)- Rhy
c		enddo
c	   enddo
c	endif 
