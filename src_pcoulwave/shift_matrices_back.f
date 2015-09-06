      subroutine shift_matrices_back(n_loc)
      use mainvar_module
	integer n_loc

cZZZZZZZZZZ
            do cur_level=1,num_levels
cZZZZZZZZZZ

                   do j=js1,je1
                        do i=1,endx
                              if(bl_hor_wall(i,j).eq.99)then
								tmp=zeta(i,j,n_loc,cur_level)+h(i,j,n_loc)

                                    UU(i,j,n_loc,cur_level)=
     -                                    u(i,j,n_loc,cur_level)*tmp
                                    VV(i,j,n_loc,cur_level)=
     -                                    v(i,j,n_loc,cur_level)*tmp
                              endif
                        enddo
                  enddo

            do j=js1,je1
              do i=1,endx
                     do ii=1,3
                        zeta(i,j,ii,cur_level)=
     -                        zeta(i,j,ii+1,cur_level)
                        u(i,j,ii,cur_level)=
     -                        u(i,j,ii+1,cur_level)
                        v(i,j,ii,cur_level)=
     -                        v(i,j,ii+1,cur_level)
                        E(i,j,ii,cur_level)=
     -                        E(i,j,ii+1,cur_level)
                        F1(i,j,ii,cur_level)=
     -                        F1(i,j,ii+1,cur_level)
                        G1(i,j,ii,cur_level)=
     -                        G1(i,j,ii+1,cur_level)
                        UU(i,j,ii,cur_level)=
     -                        UU(i,j,ii+1,cur_level)
                        VV(i,j,ii,cur_level)=
     -                        VV(i,j,ii+1,cur_level)
                        F(i,j,ii,cur_level)=
     -                        F(i,j,ii+1,cur_level)
                        G(i,j,ii,cur_level)=
     -                        G(i,j,ii+1,cur_level)
                        z_alp(i,j,ii,cur_level)=
     -                        z_alp(i,j,ii+1,cur_level)	

                        if(cur_level.eq.1)then
                           h(i,j,ii)=h(i,j,ii+1)
						   nut(i,j,ii)=nut(i,j,ii+1)	
						   B_mat(i,j,ii)=B_mat(i,j,ii+1)	
						   B_RHS(i,j,ii)=B_RHS(i,j,ii+1)
						   k_RHS(i,j,ii)=k_RHS(i,j,ii+1)
						   if(ihvor.eq.1)then
							psix(i,j,ii)=psix(i,j,ii+1)
							psiy(i,j,ii)=psiy(i,j,ii+1)
						   endif
                        endif

                        dz_alpdt(i,j,ii,cur_level)=
     -                           dz_alpdt(i,j,ii+1,cur_level)

                        if(cur_level.eq.1)then
                               dhdxt(i,j,ii)=dhdxt(i,j,ii+1)
                               dhdyt(i,j,ii)=dhdyt(i,j,ii+1)
                               dhdtt(i,j,ii)=dhdtt(i,j,ii+1)
                               dhdxtt(i,j,ii)=dhdxtt(i,j,ii+1)
                               dhdytt(i,j,ii)=dhdytt(i,j,ii+1)
                               hp(i,j,ii)=hp(i,j,ii+1)
                               dhdt(i,j,ii)=dhdt(i,j,ii+1)
                        endif
                     enddo
              enddo
           enddo
cZZZZZZ
            enddo
cZZZZZZ
      return
      end      ! end subroutine
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

 
