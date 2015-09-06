! z --> zi

	subroutine interp_xei(f,fi,opt,id)


      use mainvar_module, only:nx,ny,dim,overlap
	implicit none

	integer :: x3,x4,y3,y4
	integer :: i,j,opt
	integer :: id(nx,ny)

	real :: f(nx,ny),fi(nx,ny)
	real :: lo(nx),di(nx),up(nx)

	if(dim.eq.1)then
		y3=overlap+1
		y4=overlap+1
	else
		y3=2
		y4=ny-2
	endif
	x3=2
	x4=nx-2

	if (opt.eq.2) then ! 2 = compact, 1 = explicit

	   do j = y3,y4
		do i = x3,x4
		   if (i.eq.x3) then
			  fi(i,j) = 1.5*f(i+1,j)
			  lo(i) = 0.0
				di(i) = 1.0
				up(i) = 0.5
		   elseif (i.eq.x4) then
			  fi(i,j) = 1.5*f(i,j)
				lo(i) = 0.5
				di(i) = 1.0
				up(i) = 0.0
		   else
			  fi(i,j) = 0.75*(f(i,j)+f(i+1,j))
				lo(i) = 0.25
				di(i) = 1.0
				up(i) = 0.25
		   endif
		enddo
c	   call tridag(lo,di,up,fi(:,j),fi(:,j),x4-x3+1)
		print *,'need to fix FVM tridiag interp_xei'
	   enddo

	elseif (opt.eq.1) then ! explicit interpolation

	   do j = y3,y4
		do i = x3,x4
c		   if (id(i,j).eq.1) then
c			  fi(i,j)=(7.*(f(i+1,j)+f(i+1,j))-(f(i+2,j)+f(i+2,j)))/12.0
c			 elseif (id(i,j).eq.2) then
c			  fi(i,j) = (7.*(f(i,j)+f(i+1,j))-(f(i,j)+f(i+2,j)))/12.0
c			 elseif (id(i,j).eq.3) then
c			  fi(i,j)=(7.*(f(i,j)+f(i+1,j))-(f(i-1,j)+f(i+1,j)))/12.0
c		   elseif (id(i,j).eq.4) then
c			  fi(i,j) = (7.*(f(i,j)+f(i,j))-(f(i-1,j)+f(i-1,j)))/12.0
c		   elseif (id(i,j).eq.9) then
c			  fi(i,j) = 0.0
c		   else
			  fi(i,j)=(7.*(f(i,j)+f(i+1,j))-(f(i-1,j)+f(i+2,j)))/12.0
c		   endif
		enddo
	   enddo

	elseif (opt.eq.3) then ! explicit interpolation

	   do j = y3,y4
		do i = x3,x4
c		   if (id(i,j).eq.1.or.id(i,j).eq.13.or.id(i,j).eq.14) then
c			  fi(i,j) = (f(i+1,j)+f(i+1,j))/2.
c		   elseif (id(i,j).eq.2.or.id(i,j).eq.23.or.id(i,j).eq.24) then
c			  fi(i,j) = (f(i,j)+f(i,j))/2.
c		   elseif (id(i,j).eq.9) then
c			  fi(i,j) = 0.0
c		   else
			  fi(i,j) = (f(i,j)+f(i+1,j))/2.
c		   endif
		enddo
	   enddo

	endif

	return
	end