! z --> zi

	subroutine shoreflux(i,j,id,zx_4th,flux,z,u,v,d,dx)
	
      use mainvar_module, only:nx,ny,dim,overlap,upwind_shore,step_shore
	implicit none

	integer :: i,j,id

	real :: zx_4th(nx,ny),flux(nx,ny,3),d(nx,ny)
	real :: z(nx,ny),u(nx,ny),v(nx,ny)
	real :: dx,Hs_i,us_i,vs_i

c	id=1  shoreline with id=1, x-dir fluxes
c	id=2  shoreline with id=2, x-dir fluxes
c	id=3  low-order upwinding

!	upwind_shore=0 !=0 use upwinding for flux terms at wet side of first wet cell, use =1 for centered
	step_shore=0 !=0 use stair-stepped formulation for shoreline flooding, use =1 for peice-wise linear formulation

	if(id.eq.2)then
		
		zx_4th(i-1,j)= (z(i,j)-z(i-2,j))/dx/2.
		if(z(i,j).gt.z(i+1,j))then
			zx_4th(i,j)= (z(i+1,j)-z(i-1,j))/dx/2.
		else
			zx_4th(i,j)= (z(i,j)-z(i-1,j))/dx
		endif

		if(step_shore.eq.0)then			
			if(z(i,j).gt.-d(i+1,j))then 
				Hs_i=z(i,j)+min(d(i,j),d(i+1,j))
				us_i=u(i,j)
				vs_i=v(i,j)
			else
				Hs_i=0.
				us_i=0.
				vs_i=0.
			endif
		else
			if(z(i,j).gt.-(d(i+1,j)+d(i,j))/2.)then
				Hs_i=z(i,j)+(d(i+1,j)+d(i,j))/2.
				us_i=u(i,j)
				vs_i=v(i,j)
			else
				Hs_i=0.
				us_i=0.
				vs_i=0.
			endif
		endif

		flux(i,j,1)=Hs_i*us_i
		flux(i,j,2)=Hs_i*us_i*us_i ! + 0.5*grav*Hs_i**2.
		flux(i,j,3)=Hs_i*us_i*vs_i
	elseif(id.eq.1)then

		zx_4th(i+1,j)= (z(i+2,j)-z(i,j))/dx/2.
		if(z(i,j).gt.z(i-1,j))then
			zx_4th(i,j)= (z(i+1,j)-z(i-1,j))/dx/2.
		else
			zx_4th(i,j)= (z(i+1,j)-z(i,j))/dx
		endif

		if(step_shore.eq.0)then				
			if(z(i,j).gt.-d(i-1,j))then
				Hs_i=z(i,j)+min(d(i,j),d(i-1,j))
				us_i=u(i,j)
				vs_i=v(i,j)
			else
				Hs_i=0.
				us_i=0.
				vs_i=0.
			endif
		else
			if(z(i,j).gt.-(d(i-1,j)+d(i,j))/2.)then
				Hs_i=z(i,j)+(d(i-1,j)+d(i,j))/2.
				us_i=u(i,j)
				vs_i=v(i,j)
			else
				Hs_i=0.
				us_i=0.
				vs_i=0.
			endif
		endif	
			
		flux(i-1,j,1)=Hs_i*us_i
		flux(i-1,j,2)=Hs_i*us_i*us_i ! + 0.5*grav*Hs_i**2.
		flux(i-1,j,3)=Hs_i*us_i*vs_i
	elseif(id.eq.3)then
		if(upwind_shore.eq.0)then
			us_i=(u(i,j)+u(i+1,j))/2.
			
			if(us_i.ge.0.)then
				us_i=u(i,j)
				vs_i=v(i,j)
				Hs_i=z(i,j)+d(i,j)
			else
				us_i=u(i+1,j)
				vs_i=v(i+1,j)
				Hs_i=z(i+1,j)+d(i+1,j)
			endif

		 else
			us_i=(u(i,j)+u(i+1,j))/2.
			vs_i=(v(i,j)+v(i+1,j))/2.
			Hs_i=(z(i,j)+d(i,j)+z(i+1,j)+d(i+1,j))/2.
		 endif
		 flux(i,j,1)=Hs_i*us_i
		 flux(i,j,2)=Hs_i*us_i*us_i ! + 0.5*grav*Hs_i**2.
		 flux(i,j,3)=Hs_i*us_i*vs_i
	endif

	return
	end