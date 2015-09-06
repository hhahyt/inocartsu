C*******************************************************************************

      subroutine tridiag(a,b,c,d,endx,endy,u,v,flag)
      implicit none
      integer endx,endy,flag,i,j,ind
      real a(endx,endy),b(endx,endy),c(endx,endy),d(endx,endy),
     -		u(endx,endy),v(endx,endy)


	if(flag.eq.1)then
		do j=1,endy

C          decomposition
			do ind=2,endx
				a(ind,j)=a(ind,j)/b(ind-1,j)
				b(ind,j)=b(ind,j)-a(ind,j)*c(ind-1,j)
			enddo

C          forward substitution
			do ind=2,endx
				d(ind,j)=d(ind,j)-a(ind,j)*d(ind-1,j)
			enddo

C          back subsitution	
			u(endx,j)=d(endx,j)/b(endx,j)
			do ind=endx-1,1,-1
				u(ind,j)=(d(ind,j)-
     -				c(ind,j)*u(ind+1,j))/b(ind,j)
			enddo
		enddo
	elseif(flag.eq.2)then
		do i=1,endx

C          decomposition
			do ind=2,endy
				a(i,ind)=a(i,ind)/b(i,ind-1)
				b(i,ind)=b(i,ind)-a(i,ind)*c(i,ind-1)
			enddo

C          forward substitution
			do ind=2,endy
				d(i,ind)=d(i,ind)-a(i,ind)*d(i,ind-1)
			enddo

C          back subsitution
			v(i,endy)=d(i,endy)/b(i,endy)
			do ind=endy-1,1,-1
				v(i,ind)=(d(i,ind)-
     -				c(i,ind)*v(i,ind+1))/b(i,ind)

			enddo
		enddo
	endif
      

      return
      end      ! end subroutine
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
