C*******************************************************************************
C......  Calculates total mass (difference from still water condition) in numerical domain
      subroutine calc_vel_dz(n_loc,z,u_z1,v_z1,w_z1)
      use mainvar_module, only:i,j,zeta,u,v,z_alp,cur_level,num_levels,
     -	cur_level_z,dSdx,dTdx,dSdy,dTdy,T_grp,S_grp,h,disp_prop,dim,
     -	bl_hor_wall
	integer n_loc,sum_m,ii,jj
      real z,u_z1,v_z1,w_z1

	sum_m=0
	do jj=j,j
		do ii=i,i
			sum_m=sum_m+bl_hor_wall(ii,jj)
		enddo
	enddo 
	
	if(sum_m.ge.99.or.disp_prop.eq.3)then
		u_z1=u(i,j,n_loc,1)
		v_z1=v(i,j,n_loc,1)
		w_z1=0.
	else
		if(num_levels.eq.1)then
			cur_level_z=1
		elseif(num_levels.eq.2)then
			if(z.gt.zeta(i,j,n_loc,2))then
				cur_level_z=1
			else
				cur_level_z=2
			endif
		elseif(num_levels.eq.3)then
			if(z.gt.zeta(i,j,n_loc,2))then
				cur_level_z=1
			elseif(z.gt.zeta(i,j,n_loc,3))then
				cur_level_z=2
			else
				cur_level_z=3
			endif	
		elseif(num_levels.eq.4)then
			if(z.gt.zeta(i,j,n_loc,2))then
				cur_level_z=1
			elseif(z.gt.zeta(i,j,n_loc,3))then
				cur_level_z=2
			elseif(z.gt.zeta(i,j,n_loc,4))then
				cur_level_z=3
			else
				cur_level_z=4
			endif							                    
		endif

		u_z1=-0.5*(2.*z)*dSdx(i,j,cur_level_z)-dTdx(i,j,cur_level_z)    

		if(dim.eq.1)then
			v_z1=0.
		else
			v_z1=-0.5*(2.*z)*dSdy(i,j,cur_level_z)-dTdy(i,j,cur_level_z)
c     -  +u0_b(i,j,2)*(exp(k_br(i,j)*(z-zeta(cur_grid,i,j,4,1)))-
c     -			 exp(k_br(i,j)*(z_br(i,j)-zeta(cur_grid,i,j,4,1))))   
		endif
	
		w_z1=-S_grp(i,j,cur_level_z)
	endif

      return
      end      ! end subroutine
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

 
