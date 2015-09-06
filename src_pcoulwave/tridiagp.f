	subroutine tridagp(a,b,c,r,M,N,xrank,xcomm,xchunks,x)
	implicit none

	!=============================================
	! a       = lower diagonal
	! b       = diagonal
	! c       = upper diagonal
	! r       = right hand side   
	! M       = # of equations
	! N       = # of system of equations
	! xchunks = # of chunks in x-direction  
	! xrank   = rank of the new x-group 
	! xcomm   = new x-communicator
	! x       = solution
	! ipivot  : 0: no pivot  1: pivot
	!=============================================

	include 'mpif.h'

	integer ipivot
	parameter(ipivot=1) 

	integer M,N,ierr,xrank,xcomm,xchunks,MM,MMM,i,j,
     -		indx(2*xchunks-2),ii,k,l

	real tiny
	parameter(tiny=1.e-20) 

	real a(M,N),b(M,N),c(M,N),r(M,N),xlh(M,N),xuh(M,N),xr(M,N),
     -		OutData(4*N,2),reduc_abcd(4*N,2*xchunks),
     -		bet,denom,gam(M),coeffs(2*xchunks-2,N),
     -		uhcoeff(N),lhcoeff(N),x(M,N),
     -		reduc_a(2*xchunks-2,N),reduc_b(2*xchunks-2,N),
     -		reduc_c(2*xchunks-2,N),
     -		reduc_r(2*xchunks-2,N),
     -		d,al(2*xchunks-2,3),dum

c	! forward substitution

	do j=1,N
	   bet = b(1,j)
	   gam(1) = c(1,j)/bet
	   xlh(1,j) = r(1,j)/bet

	   do i=2,M
		denom = b(i,j)-a(i,j)*gam(i-1)
		gam(i) = c(i,j)/denom
		xlh(i,j) = (r(i,j)-a(i,j)*xlh(i-1,j))/denom
	   enddo

c	   ! back substitution

	   xr(M,j) = xlh(M,j)
	   xlh(M,j) = -gam(M)
	   xuh(M,j) = a(M,j)/b(M,j)

	   do i=M-1,1,-1
		xr(i,j) = xlh(i,j)-gam(i)*xr(i+1,j)
		xlh(i,j) = -gam(i)*xlh(i+1,j)
		denom = b(i,j)-c(i,j)*xuh(i+1,j)
		xuh(i,j) = a(i,j)/denom
	   enddo

c	   ! forward substitution

	   xuh(1,j) = -xuh(1,j)
	   do i=2,M
		xuh(i,j) = -xuh(i,j)*xuh(i-1,j)
	   enddo

	enddo

c	! write constributions of current processor into OutData:

	do i=1,N
	   OutData(4*i-3,1) = -1.0
	   OutData(4*i-2,1) = xuh(1,i)
	   OutData(4*i-1,1) = xlh(1,i)
	   OutData(4*i,1) = -xr(1,i)
	   OutData(4*i-3,2) = xuh(M,i)
	   OutData(4*i-2,2) = xlh(M,i)
	   OutData(4*i-1,2) = -1.0
	   OutData(4*i,2) = -xr(M,i)
	enddo

	call mpi_allgather(OutData(1,1),8*N,mpi_real,
     -		reduc_abcd(1,1),8*N,mpi_real,xcomm,ierr)

c	! solve reduced system

	MM = 2*xchunks-1
	MMM = 2*(xchunks-1)

	do i=1,N
	   reduc_a(:,i) = reduc_abcd(4*i-3,2:MM)
	   reduc_b(:,i) = reduc_abcd(4*i-2,2:MM)
	   reduc_c(:,i) = reduc_abcd(4*i-1,2:MM)
	   reduc_r(:,i) = reduc_abcd(4*i,2:MM)
	enddo

	do i=1,n
	   if(ipivot.eq.0)then
		bet=reduc_b(1,i)
		if(bet.eq.0.)pause 'tridag: rewrite equations'
		coeffs(1,i)=reduc_r(1,i)/bet
		do j=2,MMM
		   gam(j)=reduc_c(j-1,i)/bet
		   bet=reduc_b(j,i)-reduc_a(j,i)*gam(j)
		   if(bet.eq.0.)pause 'tridag failed'
		   coeffs(j,i)=(reduc_r(j,i)-reduc_a(j,i)*coeffs(j-1,i))/bet
		enddo
   
		do j=MMM-1,1,-1
		   coeffs(j,i)=coeffs(j,i)-gam(j+1)*coeffs(j+1,i)
		enddo
	   else
		reduc_a(1,i)=reduc_b(1,i)
		reduc_b(1,i)=reduc_c(1,i)
		reduc_c(1,i) = 0. 

		d=1.
		l = 1
		do k=1,mmm
		   dum= reduc_a(k,i)
		   ii=k
		   if(l.lt.mmm)l=l+1
		   do j=k+1,l
			  if(abs(reduc_a(j,i)).gt.abs(dum))then
				 dum=reduc_a(j,i)
				 ii=j
			  endif
		   enddo
   
		   indx(k)=ii
		   if(dum.eq.0.) reduc_a(k,i)=tiny
		   if(ii.ne.k)then
			  d=-d
			  dum=reduc_a(k,i)
			  reduc_a(k,i)=reduc_a(ii,i)
			  reduc_a(ii,i)=dum

			  dum=reduc_b(k,i)
			  reduc_b(k,i)=reduc_b(ii,i)
			  reduc_b(ii,i)=dum

			  dum=reduc_c(k,i)
			  reduc_c(k,i)=reduc_c(ii,i)
			  reduc_c(ii,i)=dum

		   endif
   
		   do ii=k+1,l
			  dum=reduc_a(ii,i)/reduc_a(k,i)
			  al(k,ii-k)=dum
			  reduc_a(ii,i)=reduc_b(ii,i)-dum*reduc_b(k,i)
			  reduc_b(ii,i)=reduc_c(ii,i)-dum*reduc_c(k,i)
			  reduc_c(ii,i) = 0. 

		   enddo
		enddo

		l = 1
		do k=1,mmm
		   ii=indx(k)
		   if(ii.ne.k)then
			  dum=reduc_r(k,i)
			  reduc_r(k,i)=reduc_r(ii,i)
			  reduc_r(ii,i)=dum
		   endif
   
		   if(l.lt.mmm)l=l+1

		   do ii=k+1,l
			  reduc_r(ii,i)=reduc_r(ii,i)-al(k,ii-k)*reduc_r(k,i)
		   enddo
		enddo

		l=1
		do ii=mmm,1,-1
		   dum=reduc_r(ii,i)

		   do k=2,l 
			  if(k.eq.2)then
				   dum=dum-reduc_b(ii,i)*coeffs(k+ii-1,i)
			  elseif(k.eq.3)then
				   dum=dum-reduc_c(ii,i)*coeffs(k+ii-1,i)
			  endif
		   enddo
		   coeffs(ii,i)=dum/reduc_a(ii,i)
		   if(l.lt.3) l=l+1
		enddo
	   endif
	enddo

c	! pick out the appropriate elements of coeffs

	j = 2*xrank

	if(xrank.ne.0)then
	   uhcoeff(:) = coeffs(j,:)
	else
	   uhcoeff(:) = 0.0
	endif
   
	if(xrank.ne.xchunks-1)then
	   lhcoeff(:) = coeffs(j+1,:)
	else
	   lhcoeff(:) = 0.0
	endif

c	! compute final solution

	do i=1,N
		x(:,i) = xr(:,i)+uhcoeff(i)*xuh(:,i)+lhcoeff(i)*xlh(:,i)
	enddo


	return
	end



















