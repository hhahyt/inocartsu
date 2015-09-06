! z --> zi

	subroutine interp1D_single(x,z,xi,zi,n,ni)	  ! linear 1D matrix interpolation
	implicit none

	integer :: i,j,n,ni

	real :: x(n),z(n)
	real :: xi,zi,dz,dx


	do i=1,ni
		do j=1,n-1
			if(x(j).le.xi.and.x(j+1).gt.xi)then
				dz=z(j+1)-z(j)
				dx=x(j+1)-x(j)
				zi=z(j)+  dz/dx*(xi-x(j))
			endif
		enddo
		if(xi.ge.x(n))then
			zi=z(n)  
		elseif(xi.lt.x(1))then
			zi=z(1)
		endif
	enddo			

	return
	end