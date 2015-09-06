      subroutine find_ts_indices

      use mainvar_module
	
		character*19 filenamedir

CCCCCCCCCCCCCCCCCCCCCCCCCC  FIND INDEX'S FOR TIME SERIES  CCCCCCCCCCCCCCCCCCCCCC


      do cur_ts=1,num_ts

                  if(dim.eq.1)then
                        ind_ts(cur_ts,1)=0   
                        ind_ts(cur_ts,2)=endx-1    
                        do i=overlap+1,endx-overlap-1
                              if(x(i).le.x_ts(cur_ts).and.
     -                              x(i+1).gt.x_ts(cur_ts))then

                                  ind_ts(cur_ts,1)=1

                                  if(abs(x(i)-x_ts(cur_ts)).lt.
     -                               abs(x(i+1)-x_ts(cur_ts)))then
                                              ind_ts(cur_ts,2)=i
                                              ind_ts(cur_ts,3)=overlap+1
                                  else
                                              ind_ts(cur_ts,2)=i+1
                                              ind_ts(cur_ts,3)=overlap+1
                                  endif
                              endif
                        enddo
                  else
                        ind_ts(cur_ts,1)=0   
                        ind_ts(cur_ts,2)=endx-1  
                        ind_ts(cur_ts,3)=endy-1      
					  					      
                        do i=overlap+1,endx-overlap-1          
                              if(x(i).le.x_ts(cur_ts).and.
     -                           x(i+1).gt.x_ts(cur_ts))then

                                  ind_ts(cur_ts,1)=2

                                  if(abs(x(i)-x_ts(cur_ts)).lt.
     -                               abs(x(i+1)-x_ts(cur_ts)))then
                                       ind_ts(cur_ts,2)=i
                                   else
                                       ind_ts(cur_ts,2)=i+1
                                   endif
                              endif
                        enddo            
                                                                  
                        do j=overlap+1,endy-overlap-1
                          if(y(j).le.y_ts(cur_ts).and.
     -                       y(j+1).gt.y_ts(cur_ts).and.
     -					   ind_ts(cur_ts,1).eq.2)then

                                 ind_ts(cur_ts,1)=1

                                 if(abs(y(j)-y_ts(cur_ts)).lt.
     -                              abs(y(j+1)-y_ts(cur_ts)))then
                                          ind_ts(cur_ts,3)=j
                                 else
                                          ind_ts(cur_ts,3)=j+1
                                 endif
                          endif
                        enddo
                  endif

		if(ind_ts(cur_ts,1).eq.1)then
			filename = 'tmsrXXXX.dat'
			write(ranknum,'(i4.4)') cur_ts
			filename(5:8)=ranknum
			
			filenamedir='./tmsr/' // filename
			
			open(10119+cur_ts,file=filenamedir,status='unknown')
		endif
	enddo
            
     

      return

      end
