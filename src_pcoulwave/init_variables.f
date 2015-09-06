      subroutine init_variables

      use mainvar_module

            do n_tmp=1,4
                  do i=1,endx
                        do j=1,endy
cZZZZZZZZZZ
                              do cur_level=1,num_levels
cZZZZZZZZZZ
                                    zeta(i,j,n_tmp,1)=0
                                    u(i,j,n_tmp,cur_level)=0
                                    v(i,j,n_tmp,cur_level)=0
                                    E(i,j,n_tmp,cur_level)=0
                                    F1(i,j,n_tmp,cur_level)=0
                                    G1(i,j,n_tmp,cur_level)=0
                                    UU(i,j,n_tmp,cur_level)=0
                                    VV(i,j,n_tmp,cur_level)=0
                                    F(i,j,n_tmp,cur_level)=0
                                    G(i,j,n_tmp,cur_level)=0
                                    mean_hgt(i,j,n_tmp)=0
                                    mean_z(i,j,n_tmp)=0
                                    cur_x(i,j,n_tmp)=0
                                    cur_y(i,j,n_tmp)=0

cZZZZZZZZ
                              enddo
cZZZZZZZZ
                        enddo
                  enddo
            enddo

      return

      end
