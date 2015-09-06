! z --> zi

	subroutine shorefluy(i,j,id,zy_4th,fluy,z,u,v,d,dy)
	
      use mainvar_module, only:nx,ny,dim,overlap,upwind_shore,step_shore
	implicit none

	integer :: i,j,id

	real :: zy_4th(nx,ny),fluy(nx,ny,3),d(nx,ny)
	real :: z(nx,ny),u(nx,ny),v(nx,ny)
	real :: dy,Hs_i,us_i,vs_i

c	id=1  shoreline with id=3, y-dir fluxes
c	id=2  shoreline with id=4, y-dir fluxes
c	id=3  low-order upwinding

!	upwind_shore=0 !=0 use upwinding for flux terms at wet side of first wet cell, use =1 for centered
	step_shore=0 !=0 use stair-stepped formulation for shoreline flooding, use =1 for peice-wise linear formulation

	if(id.eq.2)then
			zy_4th(i,j-1)= (z(i,j)-z(i,j-2))/dy/2.
			
			if(z(i,j).gt.z(i,j+1))then
				zy_4th(i,j)= (z(i,j+1)-z(i,j-1))/dy/2.
			else
				zy_4th(i,j)= (z(i,j)-z(i,j-1))/dy
			endif
			
			if(z(i,j).gt.-d(i,j+1))then
				Hs_i=z(i,j)+min(d(i,j),d(i,j+1))
				us_i=u(i,j)
				vs_i=v(i,j)
			else
				Hs_i=0.
				us_i=0.
				vs_i=0.
			endif
	
			fluy(i,j,1)=Hs_i*vs_i
			fluy(i,j,2)=Hs_i*us_i*vs_i
			fluy(i,j,3)=Hs_i*vs_i*vs_i ! + 0.5*grav*Hs_i**2.
	elseif(id.eq.1)then
			zy_4th(i,j+1)= (z(i,j+2)-z(i,j))/dy/2.

			if(z(i,j).gt.z(i,j-1))then
				zy_4th(i,j)= (z(i,j+1)-z(i,j-1))/dy/2.
			else
				zy_4th(i,j)= (z(i,j+1)-z(i,j))/dy
			endif
			
			if(z(i,j).gt.-d(i,j-1))then
				Hs_i=z(i,j)+min(d(i,j),d(i,j-1))
				us_i=u(i,j)
				vs_i=v(i,j)
			else
				Hs_i=0.
				us_i=0.
				vs_i=0.
			endif
			
			fluy(i,j-1,1)=Hs_i*vs_i
			fluy(i,j-1,2)=Hs_i*us_i*vs_i 
			fluy(i,j-1,3)=Hs_i*vs_i*vs_i ! + 0.5*grav*Hs_i**2.
	elseif(id.eq.3)then

		if(upwind_shore.eq.0)then
			vs_i=(v(i,j)+v(i,j+1))/2.

			if(vs_i.ge.0.)then
				us_i=u(i,j)
				vs_i=v(i,j)
				Hs_i=z(i,j)+d(i,j)
			else
				us_i=u(i,j+1)
				vs_i=v(i,j+1)
				Hs_i=z(i,j+1)+d(i,j+1)
			endif
		 else
			us_i=(u(i,j)+u(i,j+1))/2.
			vs_i=(v(i,j)+v(i,j+1))/2.
			Hs_i=(z(i,j)+d(i,j)+z(i,j+1)+d(i,j+1))/2.
		 endif	

		fluy(i,j,1)=Hs_i*vs_i
		fluy(i,j,2)=Hs_i*us_i*vs_i
		fluy(i,j,3)=Hs_i*vs_i*vs_i ! + 0.5*grav*Hs_i**2.
	endif

	return
	end