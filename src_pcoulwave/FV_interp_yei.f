! z --> zj

	subroutine interp_yei(f,fi,opt,id)

      use mainvar_module, only:nx,ny,dim,overlap
	implicit none

	integer :: x3,x4,y3,y4
	integer :: i,j,opt
	integer :: id(nx,ny)

	real :: f(nx,ny),fi(nx,ny)
	real :: lo(ny),di(ny),up(ny)

	if(dim.eq.1)then
		y3=overlap+1
		y4=overlap+1
	else
		y3=2
		y4=ny-2
	endif
	x3=2
	x4=nx-2

	if (opt.eq.2) then ! compact interp.

	   do i = x3,x4
		do j = y3,y4
		   if (j.eq.y3) then
			  fi(i,j) = 1.5*f(i,j+1)
				lo(j) = 0.0
				di(j) = 1.0
				up(j) = 0.5
		   elseif (j.eq.y4) then
			  fi(i,j) = 1.5*f(i,j)
				lo(j) = 0.5
				di(j) = 1.0
				up(j) = 0.0
		   else
			  fi(i,j) = 0.75*(f(i,j)+f(i,j+1))
				lo(j) = 0.25
				di(j) = 1.0
				up(j) = 0.25
		   endif
		enddo
c		call tridag(lo,di,up,fi(i,:),fi(i,:),y4-y3+1)
		print *,'need to fix FVM tridiag interp_yei'
	   enddo

	elseif (opt.eq.1) then ! explicit interpolation

	  do j = y3,y4
	   do i = x3,x4
c		   if (id(i,j).eq.1) then
c			  fi(i,j)=(7.*(f(i,j+1)+f(i,j+1))-(f(i,j+2)+f(i,j+2)))/12.0
c			 elseif (id(i,j).eq.2) then
c			  fi(i,j) = (7.*(f(i,j)+f(i,j+1))-(f(i,j)+f(i,j+2)))/12.0
c			 elseif (id(i,j).eq.3) then
c			  fi(i,j) =(7.*(f(i,j)+f(i,j+1))-(f(i,j-1)+f(i,j+1)))/12.0
c		   elseif (id(i,j).eq.4) then
c			  fi(i,j) = (7.*(f(i,j)+f(i,j))-(f(i,j-1)+f(i,j-1)))/12.0
c		   elseif (id(i,j).eq.9) then
c			  fi(i,j) = 0.0
c		   else
			  fi(i,j) = (7*(f(i,j)+f(i,j+1))-(f(i,j-1)+f(i,j+2)))/12.0
c			 endif
		enddo
	   enddo

	elseif (opt.eq.3) then ! explicit interpolation

	  do j = y3,y4
	   do i = x3,x4
c		   if (id(i,j).eq.3.or.id(i,j).eq.13.or.id(i,j).eq.23) then
c			  fi(i,j) = (f(i,j+1)+f(i,j+1))/2.0
c		   elseif (id(i,j).eq.4.or.id(i,j).eq.14.or.id(i,j).eq.24)then
c			  fi(i,j) = (f(i,j)+f(i,j))/2.0
c		   elseif (id(i,j).eq.9) then
c			  fi(i,j) = 0.0
c		   else
			  fi(i,j) = (f(i,j)+f(i,j+1))/2.0
c			 endif
		enddo
	   enddo

	endif

	return
	end