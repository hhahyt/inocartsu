! z --> zi

	subroutine interp1D(x,z,xi,zi,n,ni)	  ! linear 1D matrix interpolation
	implicit none

	integer :: i,j,n,ni

	real :: x(n),z(n)
	real :: xi(ni),zi(ni),dz,dx


	do i=1,ni
		do j=1,n-1
			if(x(j).le.xi(i).and.x(j+1).gt.xi(i))then
				dz=z(j+1)-z(j)
				dx=x(j+1)-x(j)
				zi(i)=z(j)+  dz/dx*(xi(i)-x(j))
			endif
		enddo
		if(xi(i).ge.x(n))then
			zi(i)=z(n)  
		elseif(xi(i).lt.x(1))then
			zi(i)=z(1)
		endif
	enddo			

	return
	end