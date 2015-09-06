      subroutine move_shoreline(n_loc)

      use mainvar_module
	integer n_loc
	real hx_loc,hy_loc,co_glob,ru,rv,fr_sh,lr_cut,c_sl,updown_vel,
     -						 u_d,v_d,hd_loc

            do jj=js2,je2            
                  do ii=2,endx-1
					hx_loc=min(abs(h(ii,jj,n)-h(ii-1,jj,n)),
     -					   abs(h(ii+1,jj,n)-h(ii,jj,n)) )
					hy_loc=min(abs(h(ii,jj,n)-h(ii,jj-1,n)),
     -					   abs(h(ii,jj+1,n)-h(ii,jj,n)) )

!					( hx(ii,jj)**2.+hy(ii,jj)**2. )

					if(upwind_shore.eq.0)then
						u_d=u(ii,jj,n_loc,1)*hx(ii,jj)
						v_d=v(ii,jj,n_loc,1)*hy(ii,jj)
						
						c_sh=sqrt(9.81*max(zeta(ii,jj,n_loc,1)+
     -						  h(ii,jj,n_loc),swash_d))

						fr_sh=(u_d+v_d)/c_sh

						if(fr_sh.le.0.)then
							lr_cut=1.
						else
							lr_cut=1.-min(1.,2.*fr_sh)
						endif

						if(n_loc.ge.3)then
							cutoff_mat(ii,jj)=swash_d*lr_cut+
     -						 (1.-lr_cut)*max(swash_d,hx_loc,hy_loc)
						else
							cutoff_mat(ii,jj)=swash_d
						endif


c						ru=abs(u(ii,jj,n_loc,1))/c_sh
c						rv=abs(v(ii,jj,n_loc,1))/c_sh
c						
c						fr_sh=1.-max(0.,min(1.,ru,rv))
c
c						lr_cut=max(swash_d,hx_loc,hy_loc)
c
c
cc						cutoff_mat(ii,jj)=lr_cut-fr_sh*(lr_cut-swash_d) ! if Fr number is large, use coarse vertical shoreline resolution
cc						cutoff_mat_d(ii,jj)=swash_d
c						if(n_loc.ge.3)then
c						  if(dim.eq.1)then
c							updown_vel=u(ii,jj,n_loc,1)*hx(ii,jj)
c								hd_loc=hx_loc
c						  else
c							u_d=u(ii,jj,n_loc,1)*hx(ii,jj)
c							v_d=v(ii,jj,n_loc,1)*hy(ii,jj)
c							updown_vel=sqrt(u_d**2.+v_d**2.)
c							hd_loc=sqrt(hx_loc**2.+hy_loc**2.)
c							if(u_d.gt.1e-6.or.v_d.gt.1e-6)then
c								updown_vel=updown_vel
c							else
c								updown_vel=-updown_vel
c							endif
c						  endif
c						 
c						  if(updown_vel.gt.0.)then  ! downrush
c							cutoff_mat(ii,jj)=(cutoff_mat(ii,jj)+
c     -						  max(swash_d,hd_loc))/2.
c						  else
c							cutoff_mat(ii,jj)=(cutoff_mat(ii,jj)+swash_d)/2.
c						  endif
c						else
c						  cutoff_mat(ii,jj)=swash_d
c						endif
c						cutoff_mat(ii,jj)=max(swash_d,hx_loc,hy_loc)
						cutoff_mat_d(ii,jj)=max(swash_d,hx_loc,hy_loc)*1.1
c						cutoff_mat(ii,jj)=swash_d
c						cutoff_mat_d(ii,jj)=swash_d
					else
						cutoff_mat(ii,jj)=swash_d
						cutoff_mat_d(ii,jj)=swash_d
					endif

				enddo
		  enddo

            do jj=js1,je1            
                  do ii=1,endx
                        mod_co=1.
                        if(bl_hor_wall(ii,jj).eq.60)then
                              bl_hor_wall(ii,jj)=60
                        elseif(zeta(ii,jj,n_loc,1)+
     -						  h(ii,jj,n_loc).le.cutoff_mat(ii,jj))then
                                    bl_hor_wall(ii,jj)=99
                        else
c							if(level(ii,jj).eq.1)then
c								if(level(ii-1,jj).eq.0.and.
c     -								u(ii-1,jj,n_loc,1).gt.0)then
c										bl_hor_wall(ii,jj)=0
c								elseif(level(ii+1,jj).eq.0.and.
c     -								u(ii+1,jj,n_loc,1).lt.0)then
c										bl_hor_wall(ii,jj)=0
c								elseif(level(ii,jj-1).eq.0.and.
c     -								v(ii,jj-1,n_loc,1).gt.0)then
c										bl_hor_wall(ii,jj)=0
c								elseif(level(ii,jj+1).eq.0.and.
c     -								v(ii,jj+1,n_loc,1).lt.0)then
c										bl_hor_wall(ii,jj)=0
c								else
c									bl_hor_wall(ii,jj)=99
c								endif
c                              else
								bl_hor_wall(ii,jj)=0
c							endif
                        endif
                  enddo
            enddo


         

      return

      end
