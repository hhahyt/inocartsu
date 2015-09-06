C*******************************************************************************
C......  Calculates total mass (difference from still water condition) in numerical domain
      subroutine calc_mass(mass_loc,n_loc)
      use mainvar_module, only:dim,i,j,overlap,endx,endy,
     -			bl_hor_wall,h,zeta,dx,dy,level
	integer n_loc
      real mass_loc,h1_mass,z_loc,fact,mass_sw

      mass_loc=0.
	mass_sw=0.
      if(dim.eq.2)then
         do j=overlap+1,endy-overlap
            do i=overlap+1,endx-overlap
                h1_mass=h(i,j,n_loc)
			  z_loc=zeta(i,j,n_loc,1)
			  fact=1.0
			  if(i.eq.overlap+1.or.i.eq.endx-overlap) fact=fact/2.
			  if(j.eq.overlap+1.or.j.eq.endy-overlap) fact=fact/2.
			  if(z_loc+h1_mass.lt.0) fact=0.0

                mass_loc=mass_loc+(z_loc+h1_mass)*dx*dy*fact
                mass_sw=mass_sw+max(0.,h1_mass)*dx*dy*fact
            enddo
         enddo
      else
	   j=overlap+1
         do i=overlap+1,endx-overlap
         
                h1_mass=h(i,j,n_loc)
			  z_loc=zeta(i,j,n_loc,1)
			  fact=1.0
			  if(i.eq.overlap+1.or.i.eq.endx-overlap) fact=0.5
			  if(z_loc+h1_mass.lt.0) fact=0.0

                mass_loc=mass_loc+(z_loc+h1_mass)*dx*fact
                mass_sw=mass_sw+max(0.,h1_mass)*dx
                
         enddo
      endif

	mass_loc=mass_loc-mass_sw

      return
      end      ! end subroutine
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

 
