C*******************************************************************************
      subroutine create_depth

      use mainvar_module

	integer ni
      real tmp_rn

	ALLOCATE(x_glob(endx_glob),y_glob(endy_glob),xs(endx_glob),
     -			xe(endx_glob),ys(endy_glob),ye(endy_glob),
     -			topo_tmp(endx_glob,endy_glob),
     -			ho_glob(endx_glob,endy_glob),
     -			x_topo(m+2*overlap),y_topo(n+2*overlap),
     -			topo(n+2*overlap,m+2*overlap))


	endx=endx_glob
	endy=endy_glob
	end_x=end_x_t
	end_y=end_y_t

	cutoff=depth/cutoff

      if(load_topo.eq.1)then
 
CCCCCCCCCCCCCCCC  Assemble and write to file time, x, and y matrices CCCCCCCCCC

			x_glob(1)=-dx*overlap 
			y_glob(1)=-dy*overlap 

			do j=2,endy
                  y_glob(j)=y_glob(j-1)+dy
			enddo

            
			do i=2,endx
                  x_glob(i)=x_glob(i-1)+dx
			enddo
   
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCC Assemble and write to file water depth matrix CCCCCCCCC

			node_1=nint((x_2-x_1)/dx+overlap)
			node_2=node_1+nint((x_3-x_2)/dx)
			node_3=node_2+nint((x_4-x_3)/dx)
			node_4=node_3+nint((x_5-x_4)/dx)
			node_5=node_4+nint((x_6-x_5)/dx)
			node_6=node_5+nint((x_7-x_6)/dx)
			node_7=node_6+nint((x_8-x_7)/dx)
			node_8=node_7+nint((x_9-x_8)/dx)


			if(num_nodes.eq.2) node_2=endx
			if(num_nodes.eq.3) node_3=endx
			if(num_nodes.eq.4) node_4=endx
			if(num_nodes.eq.5) node_5=endx
			if(num_nodes.eq.6) node_6=endx
			if(num_nodes.eq.7) node_7=endx
			if(num_nodes.eq.8) node_8=endx


			slope_1=0.
			if(abs(x_1-x_2).gt.1e-8)slope_1=(h_2-h_1)/(x_2-x_1)
			do i=1,node_1
                  do j=1,endy
                        ho_glob(i,j)=h_1+slope_1*dx*(i-2)      
                  enddo
			enddo

            
			if(num_nodes.ge.2)then
                  slope_2=(h_3-h_2)/(x_3-x_2)
                  do i=node_1+1,node_2
                        do j=1,endy
                              ho_glob(i,j)=h_2+slope_2*dx*(i-node_1)
                        enddo
                  enddo
			endif

			if(num_nodes.ge.3)then
                  slope_3=(h_4-h_3)/(x_4-x_3)
                  do i=node_2+1,node_3
                        do j=1,endy
                              ho_glob(i,j)=h_3+slope_3*dx*(i-node_2)
                        enddo
                  enddo
			endif

			if(num_nodes.ge.4)then
                  slope_4=(h_5-h_4)/(x_5-x_4)
                    do i=node_3+1,node_4
                        do j=1,endy
                              ho_glob(i,j)=h_4+slope_4*dx*(i-node_3)
                        enddo
                  enddo
			endif

			if(num_nodes.ge.5)then
                  slope_5=(h_6-h_5)/(x_6-x_5)
                  do i=node_4+1,node_5
                        do j=1,endy
                              ho_glob(i,j)=h_5+slope_5*dx*(i-node_4)
                        enddo
                  enddo
			endif

			if(num_nodes.ge.6)then 
                  slope_6=(h_7-h_6)/(x_7-x_6)
                  do i=node_5+1,node_6
                        do j=1,endy
                              ho_glob(i,j)=h_6+slope_6*dx*(i-node_5)
                        enddo
                  enddo
			endif

			if(num_nodes.ge.7)then 
                  slope_7=(h_8-h_7)/(x_8-x_7)
                  do i=node_6+1,node_7
                        do j=1,endy
                              ho_glob(i,j)=h_7+slope_7*dx*(i-node_6)
                        enddo
                  enddo
			endif

			if(num_nodes.ge.8)then
                  slope_8=(h_9-h_8)/(x_9-x_8)
                  do i=node_7+1,node_8
                        do j=1,endy
                              ho_glob(i,j)=h_8+slope_8*dx*(i-node_7)
                        enddo
                  enddo
			endif

	 elseif(load_topo.eq.2)then

		open(56,file='x_topo.dat',status='old') 
		open(57,file='y_topo.dat',status='old')
		open(58,file='f_topo.dat',status='old')
		open(59,file='size_topo.dat',status='old')

		read(59,*) tmp
		m=nint(tmp)
		read(59,*) tmp
		n=nint(tmp)

		if(dim.eq.1)then
          do i=1,overlap
            x_topo(i)=-(overlap-i+1)*dx
          enddo
        
		  do i=overlap+1,m+overlap
			read(56,*) x_topo(i)
		  enddo

          do i=m+overlap+1,m+2*overlap
            x_topo(i)=x_topo(i-1)+dx
          enddo

		  do i=1,n
			read(57,*) y_topo(i)
		  enddo

		  do i=1,n
			do j=overlap+1,m+overlap
				read(58,*) topo(i,j)
			enddo
		  enddo 

		  do i=1,n
			do j=1,overlap
				topo(i,j)=topo(i,overlap+1)
			enddo
		  enddo
		
		  do i=1,n
			do j=m+overlap+1,m+2*overlap
				topo(i,j)=topo(i,j-1)
			enddo
		  enddo

		else
		  do i=1,overlap
            x_topo(i)=-(overlap-i+1)*dx
          enddo
        
		  do i=overlap+1,m+overlap
			read(56,*) x_topo(i)
		  enddo

          do i=m+overlap+1,m+2*overlap
            x_topo(i)=x_topo(i-1)+dx
          enddo

          do i=1,overlap
            y_topo(i)=-(overlap-i+1)*dx
          enddo
        
		  do i=overlap+1,n+overlap
			read(57,*) y_topo(i)
		  enddo

          do i=n+overlap+1,n+2*overlap
            y_topo(i)=y_topo(i-1)+dx
          enddo

		  do i=overlap+1,n+overlap
			do j=overlap+1,m+overlap
				read(58,*) topo(i,j)
			enddo
		  enddo 

		  do i=overlap+1,n+overlap
			do j=1,overlap
				topo(i,j)=topo(i,overlap+1)
			enddo
		  enddo
		
		  do i=overlap+1,n+overlap
			do j=m+overlap+1,m+2*overlap
				topo(i,j)=topo(i,j-1)
			enddo
		  enddo

		  do j=1,m+2*overlap
			do i=1,overlap
				topo(i,j)=topo(overlap+1,j)
			enddo
		  enddo

		  do j=1,m+2*overlap
			do i=n+overlap+1,n+2*overlap
				topo(i,j)=topo(i-1,j)
			enddo
		  enddo

		endif

		close(56)
		close(57)
		close(58)
		close(59)

C....... INTERPOLATE FOR NEW GRID
  
		mm=nint((x_topo(m+2*overlap)-x_topo(1))/dx)+1

		if(n.eq.1)then
			nn=1+2*overlap
		else
			nn=nint((y_topo(n+2*overlap)-y_topo(1))/dy)+1
		endif

		endx=mm
		endy=nn

		do i=1,mm
			x_glob(i)=x_topo(1)+dx*(i-1)
		enddo

		do i=1,nn
			y_glob(i)=y_topo(1)+dy*(i-1)
		enddo
      
		end_x=x_glob(mm)
		end_y=y_glob(nn)

		do i=1,mm
			xs(i)=1
			xe(i)=m+2*overlap
			do ii=1,m+2*overlap
				if(x_topo(ii).le.x_glob(i).and.
     -			x_topo(ii).ge.x_topo(xs(i)))then
					xs(i)=ii
				endif
				if(x_topo(ii).ge.x_glob(i).and.
     -			x_topo(ii).le.x_topo(xe(i)))then
					xe(i)=ii
				endif
			enddo
		enddo


		if(n.eq.1)then
		 do i=1,mm	  
			j=1
			ys(j)=1
			ye(j)=1

			if(xe(i).eq.xs(i))then
			  topo_tmp(i,j)=topo(ys(j),xs(i))
			else
			  s=(topo(ye(j),xe(i))-topo(ye(j),xs(i)))/
     -            (x_topo(xe(i))-x_topo(xs(i)))
			  topo_tmp(i,j)=topo(ye(j),xs(i))+s*(x_glob(i)-x_topo(xs(i)))
			endif
			ho_glob(i,j)=-topo_tmp(i,j)
		 enddo

		else

		 do i=1,nn
			ys(i)=1
			ye(i)=n
			do ii=1,n
				if(y_topo(ii).le.y_glob(i).and.
     -			y_topo(ii).ge.y_topo(ys(i)))then
					ys(i)=ii
				endif
				if(y_topo(ii).ge.y_glob(i).and.
     -			y_topo(ii).le.y_topo(ye(i)))then
					ye(i)=ii
				endif
			enddo
		 enddo

		 do i=1,mm
		  do j=1,nn
			if(xe(i).eq.xs(i).and.ye(j).eq.ys(j))then
			  topo_tmp(i,j)=topo(ys(j),xs(i))
			elseif(ye(j).eq.ys(j))then
			  s=(topo(ye(j),xe(i))-topo(ye(j),xs(i)))/
     -            (x_topo(xe(i))-x_topo(xs(i)))
			topo_tmp(i,j)=topo(ye(j),xs(i))+s*(x_glob(i)-x_topo(xs(i)))
			elseif(xe(i).eq.xs(i))then
			  s=(topo(ye(j),xe(i))-topo(ys(j),xs(i)))/
     -           (y_topo(ye(j))-y_topo(ys(j)))
			topo_tmp(i,j)=topo(ys(j),xs(i))+s*(y_glob(j)-y_topo(ys(j)))
			else
			s=(topo(ye(j),xe(i))-topo(ye(j),xs(i)))/
     -       (x_topo(xe(i))-x_topo(xs(i)))
			t1=topo(ye(j),xs(i))+s*(x_glob(i)-x_topo(xs(i)))

			s=(topo(ys(j),xe(i))-topo(ys(j),xs(i)))/
     -           (x_topo(xe(i))-x_topo(xs(i)))
			  t2=topo(ys(j),xs(i))+s*(x_glob(i)-x_topo(xs(i)))

			s=(t1-t2)/(y_topo(ye(j))-y_topo(ys(j)))
			   t5=t2+s*(y_glob(j)-y_topo(ys(j))) 
  

			s=(topo(ye(j),xe(i))-topo(ys(j),xe(i)))/
     -           (y_topo(ye(j))-y_topo(ys(j)))
			  t3=topo(ys(j),xe(i))+s*(y_glob(j)-y_topo(ys(j)))

			s=(topo(ye(j),xs(i))-topo(ys(j),xs(i)))/
     -           (y_topo(ye(j))-y_topo(ys(j)))
			  t4=topo(ys(j),xs(i))+s*(y_glob(j)-y_topo(ys(j)))

			s=(t3-t4)/(x_topo(xe(i))-x_topo(xs(i)))
			  t6=t4+s*(x_glob(i)-x_topo(xs(i)))

			topo_tmp(i,j)=(t5+t6)/2.
			endif
			ho_glob(i,j)=-topo_tmp(i,j)
		  enddo
		 enddo
		endif

	 extrap=0

       if(extrap.eq.1)then
      do jj=1,nn
          do ii=1,mm
			if(ho_glob(ii,jj).lt.0)then
				level(ii,jj)=-999
			endif
		enddo
	enddo

      do jj=3,nn-2
          do ii=3,mm-2
          if(ho_glob(ii,jj).gt.-1.e-6)then

            int_count=0
            z_interp1=0
            z_interp2=0
            z_interp3=0
            z_interp4=0
            z_interp5=0
            z_interp6=0
            z_interp7=0
            z_interp8=0

      
            if(level(ii+1,jj).eq.-999)then
                z_interp1=2.*ho_glob(ii+1,jj)-
     -                ho_glob(ii+2,jj)

                int_count=int_count+1
            endif


            if(level(ii-1,jj).eq.-999)then
                z_interp2=2.*ho_glob(ii-1,jj)-
     -                ho_glob(ii-2,jj)

                int_count=int_count+1
            endif
      
            if(level(ii,jj+1).eq.-999)then
                z_interp3=ho_glob(ii,jj+1)-
     -                (ho_glob(ii,jj+2)-
     -                ho_glob(ii,jj+1))

                int_count=int_count+1
            endif


            if(level(ii,jj-1).eq.-999)then
                z_interp4=ho_glob(ii,jj-1)-
     -                (ho_glob(ii,jj-2)-
     -                ho_glob(ii,jj-1))

                int_count=int_count+1
            endif   

            if(level(ii+1,jj+1).eq.-999)then
                z_interp5=ho_glob(ii+1,jj+1)-
     -            (ho_glob(ii+2,jj+2)-
     -                ho_glob(ii+1,jj+1))

                int_count=int_count+1
            endif


            if(level(ii-1,jj-1).eq.-999)then
                z_interp6=ho_glob(ii-1,jj-1)-
     -            (ho_glob(ii-2,jj-2)-
     -                ho_glob(ii-1,jj-1))

                int_count=int_count+1
            endif
      
            if(level(ii-1,jj+1).eq.-999)then
                z_interp7=ho_glob(ii-1,jj+1)-
     -            (ho_glob(ii-2,jj+2)-
     -                ho_glob(ii-1,jj+1))

                int_count=int_count+1
            endif


            if(level(ii+1,jj-1).eq.-999)then
                z_interp8=ho_glob(ii+1,jj-1)-
     -                (ho_glob(ii+2,jj-2)-
     -                ho_glob(ii+1,jj-1))

                int_count=int_count+1
            endif   

            if(int_count.gt.0)then
                ho_glob(ii,jj)=min(1.,(z_interp1+z_interp2+
     -                z_interp3+z_interp4+z_interp5+z_interp6+
     -                z_interp7+z_interp8)/int_count)

                level(ii,jj)=1
            endif
          endif
       enddo
      enddo


	slope=1./50.
      do ni=1,10000
       total_count=0
      do jj=2,nn-1
          do ii=2,mm-1
            if(ho_glob(ii,jj).gt.-1.e-6.and.
     -                  level(ii,jj).eq.0)then
                int_count=0
                z_interp1=0
                z_interp2=0
                z_interp3=0
                z_interp4=0
                z_interp5=0
                z_interp6=0
                z_interp7=0
                z_interp8=0


                if(level(ii+1,jj).eq.ni)then
                  z_interp1=ho_glob(ii+1,jj)+dx*slope

                  int_count=int_count+1
                endif


                if(level(ii-1,jj).eq.ni)then
                  z_interp2=ho_glob(ii-1,jj)+dx*slope

                  int_count=int_count+1
                endif
      
                if(level(ii,jj+1).eq.ni)then
                  z_interp3=ho_glob(ii,jj+1)+dx*slope

                  int_count=int_count+1
                endif

                if(level(ii,jj-1).eq.ni)then
                  z_interp4=ho_glob(ii,jj-1)+dx*slope

                  int_count=int_count+1
                endif   

                if(level(ii+1,jj+1).eq.ni)then
                  z_interp5=ho_glob(ii+1,jj+1)+dx*slope

                  int_count=int_count+1
                endif


                if(level(ii-1,jj-1).eq.ni)then
                  z_interp6=ho_glob(ii-1,jj-1)+dx*slope

                  int_count=int_count+1
                endif
      
                if(level(ii-1,jj+1).eq.ni)then
                  z_interp7=ho_glob(ii-1,jj+1)+dx*slope

                  int_count=int_count+1
                endif


                if(level(ii+1,jj-1).eq.ni)then
                  z_interp8=ho_glob(ii+1,jj-1)+dx*slope

                  int_count=int_count+1
                endif   
      
                if(int_count.gt.0)then
                  ho_glob(ii,jj)=min(30.,(z_interp1+z_interp2+
     -                  z_interp3+z_interp4+z_interp5+z_interp6+
     -                  z_interp7+z_interp8)/int_count)
       
                  level(ii,jj)=ni+1
                endif
                total_count=total_count+int_count
            endif
          enddo
       enddo
       if(total_count.eq.0)goto 50
      enddo
      
50      continue

	endif


		if(dim.eq.1)then
			do i=1,mm
				do j=1,nn
					ho_glob(i,j)=ho_glob(i,1)
				enddo
			enddo
		endif
      
      endif


      if(smooth_bathy.eq.1)then
		count_min=0
          if(dim.eq.1)then
			do i=3,endx-2
				do j=1,endy
					topo_tmp(i,j)=1./6.*
     -                        (ho_glob(i-2,j)+2.*ho_glob(i-1,j)+
     -                         ho_glob(i+2,j)+2.*ho_glob(i+1,j))
				enddo
			enddo
                        
              max_hxx=10.
			do while(max_hxx.ge.0.1)
				max_hxx=0.
				do j=1,endy
					do i=3,endx-2
						ho_glob(i,j)=topo_tmp(i,j)
					enddo

                      do i=3,endx-2
						tmp=(-ho_glob(i+2,j)+
     -						16.*ho_glob(i+1,j)-
     -						30.*ho_glob(i,j)+
     -						16.*ho_glob(i-1,j)-
     -						ho_glob(i-2,j))/(12.*dx**2)

						if(abs(tmp).gt.0.1.and.topo_tmp(i,j).ge.0)then
							if(tmp.gt.0)then
									topo_tmp(i,j)=1./30.*
     -                                         (-0.1*12*dx**2.
     -										-ho_glob(i+2,j)+
     -										16.*ho_glob(i+1,j)+
     -										16.*ho_glob(i-1,j)-
     -										ho_glob(i-2,j))
							else
									topo_tmp(i,j)=1./30.*
     -                                         (0.1*12*dx**2.
     -										-ho_glob(i+2,j)+
     -										16.*ho_glob(i+1,j)+
     -										16.*ho_glob(i-1,j)-
     -										ho_glob(i-2,j))
							endif
							max_hxx=abs(tmp)
						endif
					enddo

					do i=3,endx-2
						ho_glob(i,j)=topo_tmp(i,j)
					enddo
				enddo

				count_min=count_min+1
				if(count_min.gt.2)then
					max_hxx=0
				endif
			enddo
		elseif(dim.eq.2)then
			do i=2,endx-1
				do j=2,endy-1
					if(i.ge.3.and.i.le.endx-2.and.
     -				j.ge.3.and.j.le.endy-2)then
                                    
						topo_tmp(i,j)=1./24.*(
     -                         ho_glob(i-2,j)+2.*ho_glob(i-1,j)+
     -                         ho_glob(i+2,j)+2.*ho_glob(i+1,j)+
     -                         ho_glob(i,j-2)+2.*ho_glob(i,j-1)+
     -                         ho_glob(i,j+2)+2.*ho_glob(i,j+1)+
     -                     ho_glob(i-2,j-2)+2.*ho_glob(i-1,j-1)+
     -                     ho_glob(i+2,j+2)+2.*ho_glob(i+1,j+1)+
     -                     ho_glob(i+2,j-2)+2.*ho_glob(i+1,j-1)+
     -                     ho_glob(i-2,j+2)+2.*ho_glob(i-1,j+1))
                       elseif(ho_glob(i,j).ge.0)then
						topo_tmp(i,j)=1./8.*(
     -							ho_glob(i-1,j)+ho_glob(i+1,j)+
     -                            ho_glob(i,j-1)+ho_glob(i,j+1)+
     -                            ho_glob(i-1,j-1)+ho_glob(i+1,j+1)+
     -                            ho_glob(i+1,j-1)+ho_glob(i-1,j+1))
                       else
						topo_tmp(i,j)=ho_glob(i,j)
                       endif
				enddo
			enddo

              max_hxx=10.
              do while(max_hxx.ge.0.1)
				max_hxx=0.
				do j=3,endy-2
					do i=3,endx-2
						ho_glob(i,j)=topo_tmp(i,j)
					enddo
				enddo

				do j=3,endy-2
					do i=3,endx-2
						tmp=abs((-ho_glob(i+2,j)+
     -                              16.*ho_glob(i+1,j)-
     -                              30.*ho_glob(i,j)+
     -                              16.*ho_glob(i-1,j)-
     -                              ho_glob(i-2,j))/(12.*dx**2))

						tmp2=abs((-ho_glob(i,j+2)+
     -                              16.*ho_glob(i,j+1)-
     -                              30.*ho_glob(i,j)+
     -                              16.*ho_glob(i,j-1)-
     -                              ho_glob(i,j-2))/(12.*dy**2))

						if(tmp.ge.0.1.or.tmp2.ge.0.1)then
							topo_tmp(i,j)=(1./10.*(2.*ho_glob(i,j-1)+
     -                        2.*ho_glob(i,j+1)+4.*ho_glob(i,j)+
     -                                          ho_glob(i,j-2)+
     -                                          ho_glob(i,j+2))+
     -                          1./10.*(2.*ho_glob(i-1,j)+
     -                        2.*ho_glob(i+1,j)+4.*ho_glob(i,j)+
     -                                          ho_glob(i-2,j)+
     -                                          ho_glob(i+2,j))+
     -                          1./10.*(2.*ho_glob(i-1,j-1)+
     -                      2.*ho_glob(i+1,j+1)+4.*ho_glob(i,j)+
     -                                          ho_glob(i-2,j-2)+
     -                                          ho_glob(i+2,j+2))+
     -                          1./10.*(2.*ho_glob(i-1,j+1)+
     -                      2.*ho_glob(i+1,j-1)+4.*ho_glob(i,j)+
     -                                          ho_glob(i-2,j+2)+
     -                                          ho_glob(i+2,j-2)))/4.

							max_hxx=max(tmp,tmp2)
						endif

					enddo
				enddo
				count_min=count_min+1
				if(count_min.gt.2)then
					max_hxx=0
				endif
				do j=3,endy-2
					do i=3,endx-2
						ho_glob(i,j)=topo_tmp(i,j)
					enddo
				enddo
			enddo
		endif
      endif

	max_h=0.
	do j=3,endy-2
		do i=3,endx-2
			if(ho_glob(i,j).gt.max_h) max_h=ho_glob(i,j)
		enddo
	enddo

	DEALLOCATE(topo_tmp,topo)

	end_x_t=end_x
	end_y_t=end_y

      return
      end      ! end subroutine
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC





C*******************************************************************************
      subroutine create_depth_wvmk

      use mainvar_module

	integer ni
      real ho_wv

	ALLOCATE(x_glob(endx_glob),y_glob(endy_glob),xs(endx_glob),
     -			xe(endx_glob),ys(endy_glob),ye(endy_glob),
     -			topo_tmp(endx_glob,endy_glob),
     -			ho_glob(endx_glob,endy_glob),
     -			x_topo(m+2*overlap),y_topo(n+2*overlap),
     -			topo(n+2*overlap,m+2*overlap))


	endx=endx_glob
	endy=endy_glob
	end_x=end_x_t
	end_y=end_y_t

	cutoff=depth/cutoff

      open(888,file='depth_at_wm.dat',status='old')
      read(888,*) ho_wv
      close(888)
         
      if(load_topo.eq.1)then
 
CCCCCCCCCCCCCCCC  Assemble and write to file time, x, and y matrices CCCCCCCCCC

			x_glob(1)=-dx*overlap 
			y_glob(1)=-dy*overlap 

			do j=2,endy
                  y_glob(j)=y_glob(j-1)+dy
			enddo

            
			do i=2,endx
                  x_glob(i)=x_glob(i-1)+dx
			enddo
   
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCC Assemble and write to file water depth matrix CCCCCCCCC

			node_1=nint((x_2-x_1)/dx+overlap)
			node_2=node_1+nint((x_3-x_2)/dx)
			node_3=node_2+nint((x_4-x_3)/dx)
			node_4=node_3+nint((x_5-x_4)/dx)
			node_5=node_4+nint((x_6-x_5)/dx)
			node_6=node_5+nint((x_7-x_6)/dx)
			node_7=node_6+nint((x_8-x_7)/dx)
			node_8=node_7+nint((x_9-x_8)/dx)


			if(num_nodes.eq.2) node_2=endx
			if(num_nodes.eq.3) node_3=endx
			if(num_nodes.eq.4) node_4=endx
			if(num_nodes.eq.5) node_5=endx
			if(num_nodes.eq.6) node_6=endx
			if(num_nodes.eq.7) node_7=endx
			if(num_nodes.eq.8) node_8=endx


			slope_1=0.
			if(abs(x_1-x_2).gt.1e-8)slope_1=(h_2-h_1)/(x_2-x_1)
			do i=1,node_1
                  do j=1,endy
                        ho_glob(i,j)=h_1+slope_1*dx*(i-2)      
                  enddo
			enddo

            
			if(num_nodes.ge.2)then
                  slope_2=(h_3-h_2)/(x_3-x_2)
                  do i=node_1+1,node_2
                        do j=1,endy
                              ho_glob(i,j)=h_2+slope_2*dx*(i-node_1)
                        enddo
                  enddo
			endif

			if(num_nodes.ge.3)then
                  slope_3=(h_4-h_3)/(x_4-x_3)
                  do i=node_2+1,node_3
                        do j=1,endy
                              ho_glob(i,j)=h_3+slope_3*dx*(i-node_2)
                        enddo
                  enddo
			endif

			if(num_nodes.ge.4)then
                  slope_4=(h_5-h_4)/(x_5-x_4)
                    do i=node_3+1,node_4
                        do j=1,endy
                              ho_glob(i,j)=h_4+slope_4*dx*(i-node_3)
                        enddo
                  enddo
			endif

			if(num_nodes.ge.5)then
                  slope_5=(h_6-h_5)/(x_6-x_5)
                  do i=node_4+1,node_5
                        do j=1,endy
                              ho_glob(i,j)=h_5+slope_5*dx*(i-node_4)
                        enddo
                  enddo
			endif

			if(num_nodes.ge.6)then 
                  slope_6=(h_7-h_6)/(x_7-x_6)
                  do i=node_5+1,node_6
                        do j=1,endy
                              ho_glob(i,j)=h_6+slope_6*dx*(i-node_5)
                        enddo
                  enddo
			endif

			if(num_nodes.ge.7)then 
                  slope_7=(h_8-h_7)/(x_8-x_7)
                  do i=node_6+1,node_7
                        do j=1,endy
                              ho_glob(i,j)=h_7+slope_7*dx*(i-node_6)
                        enddo
                  enddo
			endif

			if(num_nodes.ge.8)then
                  slope_8=(h_9-h_8)/(x_9-x_8)
                  do i=node_7+1,node_8
                        do j=1,endy
                              ho_glob(i,j)=h_8+slope_8*dx*(i-node_7)
                        enddo
                  enddo
			endif

      elseif(load_topo.eq.2)then

		open(56,file='x_topo.dat',status='old') 
		open(57,file='y_topo.dat',status='old')
		open(58,file='f_topo.dat',status='old')
		open(59,file='size_topo.dat',status='old')

		read(59,*) tmp
		m=nint(tmp)
		read(59,*) tmp
		n=nint(tmp)

		if(dim.eq.1)then
          do i=1,overlap
            x_topo(i)=-(overlap-i+1)*dx
          enddo
        
		  do i=overlap+1,m+overlap
			read(56,*) x_topo(i)
		  enddo

          do i=m+overlap+1,m+2*overlap
            x_topo(i)=x_topo(i-1)+dx
          enddo

		  do i=1,n
			read(57,*) y_topo(i)
		  enddo

		  do i=1,n
			do j=1,overlap
				topo(i,j)=-ho_wv
			enddo
		  enddo

		  do i=1,n
			do j=overlap+1,m+overlap
				read(58,*) topo(i,j)
				topo(i,j)=topo(i,j)-ho_wv
			enddo
		  enddo 
		
		  do i=1,n
			do j=m+overlap+1,m+2*overlap
				topo(i,j)=topo(i,j-1)
			enddo
		  enddo

		else
		  do i=1,overlap
            x_topo(i)=-(overlap-i+1)*dx
          enddo
        
		  do i=overlap+1,m+overlap
			read(56,*) x_topo(i)
		  enddo

          do i=m+overlap+1,m+2*overlap
            x_topo(i)=x_topo(i-1)+dx
          enddo

          do i=1,overlap
            y_topo(i)=-(overlap-i+1)*dx
          enddo
        
		  do i=overlap+1,n+overlap
			read(57,*) y_topo(i)
		  enddo

          do i=n+overlap+1,n+2*overlap
            y_topo(i)=y_topo(i-1)+dx
          enddo

		  do i=overlap+1,n+overlap
			do j=overlap+1,m+overlap
				read(58,*) topo(i,j)
				topo(i,j)=topo(i,j)-ho_wv
			enddo
		  enddo 

		  do i=1,n
			do j=1,overlap
				topo(i,j)=topo(i,overlap+1)
			enddo
		  enddo
		
		  do i=1,n
			do j=m+overlap+1,m+2*overlap
				topo(i,j)=topo(i,j-1)
			enddo
		  enddo

		  do j=1,m
			do i=1,overlap
				topo(i,j)=topo(overlap+1,j)
			enddo
		  enddo

		  do j=1,m
			do i=n+overlap+1,n+2*overlap
				topo(i,j)=topo(i-1,j)
			enddo
		  enddo

		endif

		close(56)
		close(57)
		close(58)
		close(59)

C....... INTERPOLATE FOR NEW GRID
  
		mm=nint((x_topo(m+2*overlap)-x_topo(1))/dx)+1

		if(n.eq.1)then
			nn=1+2*overlap
		else
			nn=nint((y_topo(n+2*overlap)-y_topo(1))/dy)+1
		endif

		endx=mm
		endy=nn

		do i=1,mm
			x_glob(i)=x_topo(1)+dx*(i-1)
		enddo

		do i=1,nn
			y_glob(i)=y_topo(1)+dy*(i-1)
		enddo
      
		end_x=x_glob(mm)
		end_y=y_glob(nn)

		do i=1,mm
			xs(i)=1
			xe(i)=m+2*overlap
			do ii=1,m+2*overlap
				if(x_topo(ii).le.x_glob(i).and.
     -			x_topo(ii).ge.x_topo(xs(i)))then
					xs(i)=ii
				endif
				if(x_topo(ii).ge.x_glob(i).and.
     -			x_topo(ii).le.x_topo(xe(i)))then
					xe(i)=ii
				endif
			enddo
		enddo


		if(n.eq.1)then
		 do i=1,mm	  
			j=1
			ys(j)=1
			ye(j)=1

			if(xe(i).eq.xs(i))then
			  topo_tmp(i,j)=topo(ys(j),xs(i))
			else
			  s=(topo(ye(j),xe(i))-topo(ye(j),xs(i)))/
     -            (x_topo(xe(i))-x_topo(xs(i)))
			  topo_tmp(i,j)=topo(ye(j),xs(i))+s*(x_glob(i)-x_topo(xs(i)))
			endif
			ho_glob(i,j)=-topo_tmp(i,j)
		 enddo

		else

		 do i=1,nn
			ys(i)=1
			ye(i)=n
			do ii=1,n
				if(y_topo(ii).le.y_glob(i).and.
     -			y_topo(ii).ge.y_topo(ys(i)))then
					ys(i)=ii
				endif
				if(y_topo(ii).ge.y_glob(i).and.
     -			y_topo(ii).le.y_topo(ye(i)))then
					ye(i)=ii
				endif
			enddo
		 enddo

		 do i=1,mm
		  do j=1,nn
			if(xe(i).eq.xs(i).and.ye(j).eq.ys(j))then
			  topo_tmp(i,j)=topo(ys(j),xs(i))
			elseif(ye(j).eq.ys(j))then
			  s=(topo(ye(j),xe(i))-topo(ye(j),xs(i)))/
     -            (x_topo(xe(i))-x_topo(xs(i)))
			topo_tmp(i,j)=topo(ye(j),xs(i))+s*(x_glob(i)-x_topo(xs(i)))
			elseif(xe(i).eq.xs(i))then
			  s=(topo(ye(j),xe(i))-topo(ys(j),xs(i)))/
     -           (y_topo(ye(j))-y_topo(ys(j)))
			topo_tmp(i,j)=topo(ys(j),xs(i))+s*(y_glob(j)-y_topo(ys(j)))
			else
			s=(topo(ye(j),xe(i))-topo(ye(j),xs(i)))/
     -       (x_topo(xe(i))-x_topo(xs(i)))
			t1=topo(ye(j),xs(i))+s*(x_glob(i)-x_topo(xs(i)))

			s=(topo(ys(j),xe(i))-topo(ys(j),xs(i)))/
     -           (x_topo(xe(i))-x_topo(xs(i)))
			  t2=topo(ys(j),xs(i))+s*(x_glob(i)-x_topo(xs(i)))

			s=(t1-t2)/(y_topo(ye(j))-y_topo(ys(j)))
			   t5=t2+s*(y_glob(j)-y_topo(ys(j))) 
  

			s=(topo(ye(j),xe(i))-topo(ys(j),xe(i)))/
     -           (y_topo(ye(j))-y_topo(ys(j)))
			  t3=topo(ys(j),xe(i))+s*(y_glob(j)-y_topo(ys(j)))

			s=(topo(ye(j),xs(i))-topo(ys(j),xs(i)))/
     -           (y_topo(ye(j))-y_topo(ys(j)))
			  t4=topo(ys(j),xs(i))+s*(y_glob(j)-y_topo(ys(j)))

			s=(t3-t4)/(x_topo(xe(i))-x_topo(xs(i)))
			  t6=t4+s*(x_glob(i)-x_topo(xs(i)))

			topo_tmp(i,j)=(t5+t6)/2.
			endif
			ho_glob(i,j)=-topo_tmp(i,j)
		  enddo
		 enddo
		endif

	 extrap=0

       if(extrap.eq.1)then
      do jj=1,nn
          do ii=1,mm
			if(ho_glob(ii,jj).lt.0)then
				level(ii,jj)=-999
			endif
		enddo
	enddo

      do jj=3,nn-2
          do ii=3,mm-2
          if(ho_glob(ii,jj).gt.-1.e-6)then

            int_count=0
            z_interp1=0
            z_interp2=0
            z_interp3=0
            z_interp4=0
            z_interp5=0
            z_interp6=0
            z_interp7=0
            z_interp8=0

      
            if(level(ii+1,jj).eq.-999)then
                z_interp1=2.*ho_glob(ii+1,jj)-
     -                ho_glob(ii+2,jj)

                int_count=int_count+1
            endif


            if(level(ii-1,jj).eq.-999)then
                z_interp2=2.*ho_glob(ii-1,jj)-
     -                ho_glob(ii-2,jj)

                int_count=int_count+1
            endif
      
            if(level(ii,jj+1).eq.-999)then
                z_interp3=ho_glob(ii,jj+1)-
     -                (ho_glob(ii,jj+2)-
     -                ho_glob(ii,jj+1))

                int_count=int_count+1
            endif


            if(level(ii,jj-1).eq.-999)then
                z_interp4=ho_glob(ii,jj-1)-
     -                (ho_glob(ii,jj-2)-
     -                ho_glob(ii,jj-1))

                int_count=int_count+1
            endif   

            if(level(ii+1,jj+1).eq.-999)then
                z_interp5=ho_glob(ii+1,jj+1)-
     -            (ho_glob(ii+2,jj+2)-
     -                ho_glob(ii+1,jj+1))

                int_count=int_count+1
            endif


            if(level(ii-1,jj-1).eq.-999)then
                z_interp6=ho_glob(ii-1,jj-1)-
     -            (ho_glob(ii-2,jj-2)-
     -                ho_glob(ii-1,jj-1))

                int_count=int_count+1
            endif
      
            if(level(ii-1,jj+1).eq.-999)then
                z_interp7=ho_glob(ii-1,jj+1)-
     -            (ho_glob(ii-2,jj+2)-
     -                ho_glob(ii-1,jj+1))

                int_count=int_count+1
            endif


            if(level(ii+1,jj-1).eq.-999)then
                z_interp8=ho_glob(ii+1,jj-1)-
     -                (ho_glob(ii+2,jj-2)-
     -                ho_glob(ii+1,jj-1))

                int_count=int_count+1
            endif   

            if(int_count.gt.0)then
                ho_glob(ii,jj)=min(1.,(z_interp1+z_interp2+
     -                z_interp3+z_interp4+z_interp5+z_interp6+
     -                z_interp7+z_interp8)/int_count)

                level(ii,jj)=1
            endif
          endif
       enddo
      enddo


	slope=1./50.
      do ni=1,10000
       total_count=0
      do jj=2,nn-1
          do ii=2,mm-1
            if(ho_glob(ii,jj).gt.-1.e-6.and.
     -                  level(ii,jj).eq.0)then
                int_count=0
                z_interp1=0
                z_interp2=0
                z_interp3=0
                z_interp4=0
                z_interp5=0
                z_interp6=0
                z_interp7=0
                z_interp8=0


                if(level(ii+1,jj).eq.ni)then
                  z_interp1=ho_glob(ii+1,jj)+dx*slope

                  int_count=int_count+1
                endif


                if(level(ii-1,jj).eq.ni)then
                  z_interp2=ho_glob(ii-1,jj)+dx*slope

                  int_count=int_count+1
                endif
      
                if(level(ii,jj+1).eq.ni)then
                  z_interp3=ho_glob(ii,jj+1)+dx*slope

                  int_count=int_count+1
                endif

                if(level(ii,jj-1).eq.ni)then
                  z_interp4=ho_glob(ii,jj-1)+dx*slope

                  int_count=int_count+1
                endif   

                if(level(ii+1,jj+1).eq.ni)then
                  z_interp5=ho_glob(ii+1,jj+1)+dx*slope

                  int_count=int_count+1
                endif


                if(level(ii-1,jj-1).eq.ni)then
                  z_interp6=ho_glob(ii-1,jj-1)+dx*slope

                  int_count=int_count+1
                endif
      
                if(level(ii-1,jj+1).eq.ni)then
                  z_interp7=ho_glob(ii-1,jj+1)+dx*slope

                  int_count=int_count+1
                endif


                if(level(ii+1,jj-1).eq.ni)then
                  z_interp8=ho_glob(ii+1,jj-1)+dx*slope

                  int_count=int_count+1
                endif   
      
                if(int_count.gt.0)then
                  ho_glob(ii,jj)=min(30.,(z_interp1+z_interp2+
     -                  z_interp3+z_interp4+z_interp5+z_interp6+
     -                  z_interp7+z_interp8)/int_count)
       
                  level(ii,jj)=ni+1
                endif
                total_count=total_count+int_count
            endif
          enddo
       enddo
       if(total_count.eq.0)goto 50
      enddo
      
50      continue

	endif


		if(dim.eq.1)then
			do i=1,mm
				do j=1,nn
					ho_glob(i,j)=ho_glob(i,1)
				enddo
			enddo
		endif
      
      endif


      if(smooth_bathy.eq.1)then
		count_min=0
          if(dim.eq.1)then
			do i=3,endx-2
				do j=1,endy
					topo_tmp(i,j)=1./6.*
     -                        (ho_glob(i-2,j)+2.*ho_glob(i-1,j)+
     -                         ho_glob(i+2,j)+2.*ho_glob(i+1,j))
				enddo
			enddo
                        
              max_hxx=10.
			do while(max_hxx.ge.0.1)
				max_hxx=0.
				do j=1,endy
					do i=3,endx-2
						ho_glob(i,j)=topo_tmp(i,j)
					enddo

                      do i=3,endx-2
						tmp=(-ho_glob(i+2,j)+
     -						16.*ho_glob(i+1,j)-
     -						30.*ho_glob(i,j)+
     -						16.*ho_glob(i-1,j)-
     -						ho_glob(i-2,j))/(12.*dx**2)

						if(abs(tmp).gt.0.1.and.topo_tmp(i,j).ge.0)then
							if(tmp.gt.0)then
									topo_tmp(i,j)=1./30.*
     -                                         (-0.1*12*dx**2.
     -										-ho_glob(i+2,j)+
     -										16.*ho_glob(i+1,j)+
     -										16.*ho_glob(i-1,j)-
     -										ho_glob(i-2,j))
							else
									topo_tmp(i,j)=1./30.*
     -                                         (0.1*12*dx**2.
     -										-ho_glob(i+2,j)+
     -										16.*ho_glob(i+1,j)+
     -										16.*ho_glob(i-1,j)-
     -										ho_glob(i-2,j))
							endif
							max_hxx=abs(tmp)
						endif
					enddo

					do i=3,endx-2
						ho_glob(i,j)=topo_tmp(i,j)
					enddo
				enddo

				count_min=count_min+1
				if(count_min.gt.2)then
					max_hxx=0
				endif
			enddo
		elseif(dim.eq.2)then
			do i=2,endx-1
				do j=2,endy-1
					if(i.ge.3.and.i.le.endx-2.and.
     -				j.ge.3.and.j.le.endy-2)then
                                    
						topo_tmp(i,j)=1./24.*(
     -                         ho_glob(i-2,j)+2.*ho_glob(i-1,j)+
     -                         ho_glob(i+2,j)+2.*ho_glob(i+1,j)+
     -                         ho_glob(i,j-2)+2.*ho_glob(i,j-1)+
     -                         ho_glob(i,j+2)+2.*ho_glob(i,j+1)+
     -                     ho_glob(i-2,j-2)+2.*ho_glob(i-1,j-1)+
     -                     ho_glob(i+2,j+2)+2.*ho_glob(i+1,j+1)+
     -                     ho_glob(i+2,j-2)+2.*ho_glob(i+1,j-1)+
     -                     ho_glob(i-2,j+2)+2.*ho_glob(i-1,j+1))
                       elseif(ho_glob(i,j).ge.0)then
						topo_tmp(i,j)=1./8.*(
     -							ho_glob(i-1,j)+ho_glob(i+1,j)+
     -                            ho_glob(i,j-1)+ho_glob(i,j+1)+
     -                            ho_glob(i-1,j-1)+ho_glob(i+1,j+1)+
     -                            ho_glob(i+1,j-1)+ho_glob(i-1,j+1))
                       else
						topo_tmp(i,j)=ho_glob(i,j)
                       endif
				enddo
			enddo

              max_hxx=10.
              do while(max_hxx.ge.0.1)
				max_hxx=0.
				do j=3,endy-2
					do i=3,endx-2
						ho_glob(i,j)=topo_tmp(i,j)
					enddo
				enddo

				do j=3,endy-2
					do i=3,endx-2
						tmp=abs((-ho_glob(i+2,j)+
     -                              16.*ho_glob(i+1,j)-
     -                              30.*ho_glob(i,j)+
     -                              16.*ho_glob(i-1,j)-
     -                              ho_glob(i-2,j))/(12.*dx**2))

						tmp2=abs((-ho_glob(i,j+2)+
     -                              16.*ho_glob(i,j+1)-
     -                              30.*ho_glob(i,j)+
     -                              16.*ho_glob(i,j-1)-
     -                              ho_glob(i,j-2))/(12.*dy**2))

						if(tmp.ge.0.1.or.tmp2.ge.0.1)then
							topo_tmp(i,j)=(1./10.*(2.*ho_glob(i,j-1)+
     -                        2.*ho_glob(i,j+1)+4.*ho_glob(i,j)+
     -                                          ho_glob(i,j-2)+
     -                                          ho_glob(i,j+2))+
     -                          1./10.*(2.*ho_glob(i-1,j)+
     -                        2.*ho_glob(i+1,j)+4.*ho_glob(i,j)+
     -                                          ho_glob(i-2,j)+
     -                                          ho_glob(i+2,j))+
     -                          1./10.*(2.*ho_glob(i-1,j-1)+
     -                      2.*ho_glob(i+1,j+1)+4.*ho_glob(i,j)+
     -                                          ho_glob(i-2,j-2)+
     -                                          ho_glob(i+2,j+2))+
     -                          1./10.*(2.*ho_glob(i-1,j+1)+
     -                      2.*ho_glob(i+1,j-1)+4.*ho_glob(i,j)+
     -                                          ho_glob(i-2,j+2)+
     -                                          ho_glob(i+2,j-2)))/4.

							max_hxx=max(tmp,tmp2)
						endif

					enddo
				enddo
				count_min=count_min+1
				if(count_min.gt.2)then
					max_hxx=0
				endif
				do j=3,endy-2
					do i=3,endx-2
						ho_glob(i,j)=topo_tmp(i,j)
					enddo
				enddo
			enddo
		endif
      endif

	max_h=0.
	do j=3,endy-2
		do i=3,endx-2
			if(ho_glob(i,j).gt.max_h) max_h=ho_glob(i,j)
		enddo
	enddo

	DEALLOCATE(topo_tmp,topo)

	end_x_t=end_x
	end_y_t=end_y

      return
      end      ! end subroutine
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


