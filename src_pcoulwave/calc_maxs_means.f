C*******************************************************************************
C......  Internal source generator
      subroutine calc_maxs_means(n_loc,nn_loc)
      use mainvar_module
	integer n_loc,nn_loc


				if(nn_loc.le.4)then   
				   do i=1,endx
						do j=1,endy
							  max_fs(i,j)=0
							  max_vel(i,j)=0
						enddo
				   enddo
				endif

                  cur_level=1
                   do j=js1,je1
					do i=1,endx
                              if(zeta(i,j,n_loc,1).
     -                    gt.max_fs(i,j).and.
     -                    bl_hor_wall(i,j).lt.50)
     -             max_fs(i,j)=zeta(i,j,n_loc,1)

						vel=u(i,j,n_loc,1)
c						sqrt(u(i,j,n+1,cur_level)**2.+
c     -                    v(i,j,n+1,cur_level)**2.)

                              if(vel.
     -                    gt.max_vel(i,j).and.
     -                    bl_hor_wall(i,j).lt.50)
     -             max_vel(i,j)=vel

						
                        enddo
                  enddo

                  t_q_ss=per*40.  !time the quasi steady-state has been reached
                  nt_per=nint(per/dt)
                  if(t(nn_loc).gt.t_q_ss)then
                      do j=js1,je1
                        do i=1,endx
C.........................Sum values over one wave period
                            if(ho(i,j).gt.0)then
                                if(mod(nn_loc,nt_per).eq.0)then
                                    if(mean_z(i,j,n_loc).gt.1.)then

                                          mean_z(i,j,1)=
     -                                          mean_z(i,j,1)+
     -                                          mean_z(i,j,2)
                                          
                                          mean_z(i,j,3)=real(nn_loc)

                                          mean_z(i,j,2)=
     -                                  zeta(i,j,n_loc,1)

                                          mean_hgt(i,j,1)=
     -                                         mean_hgt(i,j,1)+
     -                                         mean_hgt(i,j,2)-
     -                                          mean_hgt(i,j,3)
                                          mean_hgt(i,j,4)=
     -                                       mean_hgt(i,j,4)+1.      

                                          cur_x(i,j,1)=
     -                                          cur_x(i,j,1)+
     -                                          cur_x(i,j,2)
                                          
                                          cur_x(i,j,3)=real(nn_loc)

                                          cur_x(i,j,2)=
     -                                    u(i,j,n_loc,1)


                                          cur_y(i,j,1)=
     -                                          cur_y(i,j,1)+
     -                                          cur_y(i,j,2)
                                          
                                          cur_y(i,j,3)=real(nn_loc)

                                          cur_y(i,j,2)=
     -                                    v(i,j,n_loc,1)



                                elseif(mean_z(i,j,4).lt.1.)then
     -                                          
                                          mean_z(i,j,2)=
     -                                  zeta(i,j,n_loc,1)

                                          mean_z(i,j,4)=real(nn_loc)

                                          cur_x(i,j,2)=
     -                                     u(i,j,n_loc,1)

                                          cur_x(i,j,4)=real(nn_loc)

                                          cur_y(i,j,2)=
     -                                    v(i,j,n_loc,1)

                                          cur_y(i,j,4)=real(nn_loc)
                                    endif
                                
                                    mean_hgt(i,j,2)=0
                                    mean_hgt(i,j,3)=0

                                elseif(mean_z(i,j,4).gt.1.)then
                                
                                    mean_z(i,j,2)=
     -                                          mean_z(i,j,2)+
     -                                 zeta(i,j,n_loc,1)      
     
                                    cur_x(i,j,2)=
     -                                          cur_x(i,j,2)+
     -                                    u(i,j,n_loc,1)                                         

                                    cur_y(i,j,2)=
     -                                          cur_y(i,j,2)+
     -                                    v(i,j,n_loc,1)      
                                    
                                if(zeta(i,j,n_loc,1).gt.
     -                                   mean_hgt(i,j,2))then
                                          mean_hgt(i,j,2)=
     -                                  zeta(i,j,n_loc,1)
                             elseif(zeta(i,j,n_loc,1).lt.
     -                                   mean_hgt(i,j,3))then
                                          mean_hgt(i,j,3)=
     -                                  zeta(i,j,n_loc,1)
                                    endif
                                endif
                            endif
                        enddo
                  enddo
              endif  

      return
      end      ! end subroutine
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

 
