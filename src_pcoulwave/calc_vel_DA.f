C*******************************************************************************
C......  Calculates total mass (difference from still water condition) in numerical domain
      subroutine calc_vel_DA(n_loc,u_DA,v_DA,w_DA)
      use mainvar_module, only:i,j,zeta,u,v,z_alp,cur_level,num_levels,
     -	cur_level_z,dSdx,dTdx,dSdy,dTdy,T_grp,S_grp,h,disp_prop,dim,
     -	bl_hor_wall,cutoff_mat,psix,psiy
	integer n_loc,ii,jj,sum_m
      real u_DA,v_DA,w_DA,z

	sum_m=0
	do jj=j-1,j+1
		do ii=i-1,i+1
			sum_m=sum_m+bl_hor_wall(ii,jj)
		enddo
	enddo 
	
	z=zeta(i,j,n_loc,1)
      
      if((z+h(i,j,n_loc)).lt.1.e-8) sum_m=99
        
	if(sum_m.ge.99.or.disp_prop.eq.3.or.h(i,j,n_loc).
     -			le.cutoff_mat(i,j))then
		u_DA=u(i,j,n_loc,1)
		v_DA=v(i,j,n_loc,1)
		w_DA=0.
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
		u_DA=u(i,j,n_loc,1)+1./(z+h(i,j,n_loc))*
     -			( -dSdx(i,j,1)*
     -			(1./6.*(z**3+h(i,j,n_loc)**3.)
     -			-z_alp(i,j,n_loc,1)**2.
     -			*(z+h(i,j,n_loc)))
     -			-dTdx(i,j,1)*
     -			(1./2.*(z**2-h(i,j,n_loc)**2.)
     -			-z_alp(i,j,n_loc,1)
     -			*(z+h(i,j,n_loc))) ) 
     -          + psix(i,j,n_loc)*(z_alp(i,j,n_loc,1)**2*0.5
     -			-z_alp(i,j,n_loc,1)*z+
     -			(2*z**2-2*z*h(i,j,n_loc)-h(i,j,n_loc)**2)/6.0) 

		if(dim.eq.1)then
			v_DA=0.
		else
			v_DA=v(i,j,n_loc,1)+1./(z+h(i,j,n_loc))*
     -			( -dSdy(i,j,1)*
     -			(1./6.*(z**3+h(i,j,n_loc)**3.)
     -			-z_alp(i,j,n_loc,1)**2.
     -			*(z+h(i,j,n_loc)))
     -			-dTdy(i,j,1)*
     -			(1./2.*(z**2-h(i,j,n_loc)**2.)
     -			-z_alp(i,j,n_loc,1)
     -			*(z+h(i,j,n_loc))) )
     -          + psiy(i,j,n_loc)*(z_alp(i,j,n_loc,1)**2*0.5
     -			-z_alp(i,j,n_loc,1)*z+
     -			(2*z**2-2*z*h(i,j,n_loc)-h(i,j,n_loc)**2)/6.0) 
		endif
	
		w_DA=-S_grp(i,j,cur_level_z)*(z**2-h(i,j,n_loc)**2.)/2.
     -		 -T_grp(i,j,cur_level_z)*(z+h(i,j,n_loc))
	endif

      return
      end      ! end subroutine
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

 
