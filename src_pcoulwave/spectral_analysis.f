C*******************************************************************************
C......  Calculates total mass in numerical domain
      subroutine spectral_analysis(num_ts,endt,ts,te,dt_sim)
      implicit none
      integer k,N,endt
	integer cur_ts,count,i,j,jj,ii,nfc,num_ts
      real ts,te,t_e,xloc,yloc,hloc,dum,tloc, zloc, uloc, fluxloc,
     -		tspec(endt),file_ts(endt,3),dt,dt_sim,
     -		mean(num_ts,3),Hs(num_ts,3),sum_eta,var
	character*12 filename
	character*4 ranknum

      real, ALLOCATABLE :: t(:),fr(:,:),fi(:,:),eta(:)
      complex, ALLOCATABLE :: fn(:)

c	print*,2, endt,te,ts,dt_sim, k,k

	call Log2( nint( (te-ts)/dt_sim), k)
	k=min(k,14)
	N=2**k

c	print*,'Log2 of processed time series length=',k

	ALLOCATE(t(N),fr(N,3),fi(N,3),eta(N) )
	ALLOCATE(fn(N))

      do cur_ts=1,num_ts
		print*,'Processing time series',cur_ts,'of',num_ts
		filename = 'tmsrXXXX.dat'
		write(ranknum,'(i4.4)') cur_ts
		filename(5:8)=ranknum
		open(10119+cur_ts,file=filename,status='unknown')

		count=0
		read(10119+cur_ts,*) xloc, yloc, hloc, dum
		do i=1,endt-4
			read(10119+cur_ts,*) tloc, zloc, uloc, fluxloc
			if(tloc.gt.ts)then
				count=count+1
				tspec(count)=tloc
				file_ts(count,1)=zloc						
				file_ts(count,2)=uloc
				file_ts(count,3)=fluxloc
			endif
 
			if(tloc.gt.te) goto 10

		enddo

10		continue

		close(10119+cur_ts)

		tspec=tspec-tspec(1)
		t_e=tspec(count)
		dt=t_e/(N-1)
		do i=1,N
			t(i)=(i-1)*dt
		enddo

		do j=1,N
			do jj=1,count-1
				if(t(j).ge.tspec(jj).and.t(j).lt.tspec(jj+1))then
					do i=1,3
						fr(j,i)=(file_ts(jj+1,i)-file_ts(jj,i))/
     -					(tspec(jj+1)-tspec(jj))*
     -					(t(j)-tspec(jj))+file_ts(jj,i)

						fi(j,i)=0.

					enddo
					goto 14						
				endif
			enddo
14			continue
		enddo

		do i=1,3
			mean(cur_ts,i)=sum(fr(:,i))/N
			fr(:,i)=fr(:,i)-mean(cur_ts,i)

c compute Fourier coefficients for eta
			call FFT(fr(:,i),fi(:,i),k,0,N)
			nfc=N/2-1
			do j=1,nfc
				fn(j)=cmplx(fr(j+1,i),fi(j+1,i))
			enddo
c estimate spectrum for u

			do ii=1,N
			   if(ii.ge.1.and.ii.le.nfc)then
					fr(ii+1,i)=real(fn(ii))
					fi(ii+1,i)=-aimag(fn(ii))
			   else
					if(ii.ge.nfc+1.and.ii.le.(N-nfc-1))then 
						fr(ii+1,i)=0.
						fi(ii+1,i)=0.
					else
						if(ii.ge.(N-nfc+1).and.ii.le.N)then
							fr(ii,i)=real(conjg(fn(n-ii+1)))
							fi(ii,i)=-aimag(conjg(fn(n-ii+1)))
						else
						
						endif
					endif
			   endif
			enddo
			call fft(fr(:,i),fi(:,i),k,1,N)
			do ii=1,N
				eta(ii)=fr(ii,i)
			enddo

			sum_eta=0.
			do ii=1,N
				sum_eta=sum_eta+(eta(ii))**2
			enddo
			var=sum_eta/N

c calculate significant properties
			Hs(cur_ts,i) = 4.004*sqrt(var)

		enddo
		print*,'    Significant wave height (m)=',Hs(cur_ts,1)

		if(cur_ts.eq.1)then
			open(20000,file='spectral_out.dat',
     -				form='unformatted',status='unknown') ! spectral output file
			write(UNIT=20000) num_ts-1
		endif

		write(UNIT=20000) xloc, yloc, hloc, 
     -			mean(cur_ts,:), Hs(cur_ts,:)

	enddo

	close(20000)

      return
      end      ! end subroutine
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE FFT(FR,FI,K,ICO,N)                                               
	INTEGER N, ICO, K, I, MR, NN, M, L, ISTEP, EL, J, AN
      REAL FR(N),FI(N),TR, TI, A, WR, WI                                                                                                           
	K=K
      IF(ICO.EQ.0)GO TO 10                                                      
      DO 8 I=1,N                                                                
    8 FI(I)=-FI(I)                                                              
   10 CONTINUE                                                                  
      MR=0                                                                      
      NN=N-1                                                                    
      DO 2 M=1,NN                                                               
      L=N                                                                       
    1 L=L/2                                                                     
      IF(MR+L.GT.NN)GO TO 1                                                     
      MR=MOD(MR,L)+L                                                            
      IF(MR.LE.M)GO TO 2                                                        
      TR=FR(M+1)                                                                
      FR(M+1)=FR(MR+1)                                                          
      FR(MR+1)=TR                                                               
      TI=FI(M+1)                                                                
      FI(M+1)=FI(MR+1)                                                          
      FI(MR+1)=TI                                                               
    2 CONTINUE                                                                  
      L=1                                                                       
    3 IF(L.GE.N) GO TO 7                                                        
      ISTEP=2*L                                                                 
      EL=L                                                                      
      DO 4 M=1,L                                                                
      A=3.1415926535*FLOAT(1-M)/EL                                              
      WR=COS(A)                                                                 
      WI=SIN(A)                                                                 
      DO 4 I=M,N,ISTEP                                                          
      J=I+L                                                                     
      IF(ICO.EQ.1)GO TO 11                                                      
      TR=WR*FR(J)-WI*FI(J)                                                      
      TI=WR*FI(J)+WI*FR(J)                                                      
      GO TO 12                                                                  
   11 TR=WR*FR(J)+WI*FI(J)                                                      
      TI=WR*FI(J)-WI*FR(J)   
   12 FR(J)=FR(I)-TR                                                            
      FI(J)=FI(I)-TI                                                            
      FR(I)=FR(I)+TR                                                            
    4 FI(I)=FI(I)+TI                                                            
      L=ISTEP                                                                   
      GO TO 3                                                                   
    7 CONTINUE                                                                  
      AN=N                                                                      
      IF(ICO.EQ.1)GO TO 6                                                       
      DO 5 I=1,N                                                                
      FR(I)=FR(I)/AN                                                            
    5 FI(I)=-FI(I)/AN                                                           
    6 RETURN                                                                    
      END   


      SUBROUTINE Log2(number,Logtwo)

*	       ====
*     determine 2_log of number
*     -------------------------
	integer number, Logtwo, n, k

      if(number.le.0)
     /   stop 'log2: invalid argument'
      n=number
      do 10 k=1,64
         n=n/2
         Logtwo=k-1
         if(n.eq.0) return
   10 continue
      stop 'log2: still looping'
      end
