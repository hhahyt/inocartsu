      subroutine solit_cnoidal_ic

      use mainvar_module


	stem_wv=2
cZZZZZZZZZZ
            do cur_level=1,num_levels
cZZZZZZZZZZ            
            if(sim_opt.eq.1.and.int_src.eq.0)then
                  do i=1,endx
                        do j=1,endy
                              if(sim_opt.eq.1.and.
     -                           (wave_type.le.3.or.wave_type.eq.7))then
C..................subroutine 'solit' calculates free surface and velocity
C....................  based on position  (see procedure for detail)
                                    tmp=depth
                                    if(wave_type.eq.7)then
                                          depth=ho(i,j)
                                    endif
                                    
								  alp=0.5*bet(1)**2.+
     -									bet(1)

								 zc=0.
								 zc1=0.
								 uc=0.
								 uc1=0.
								 vc=0.
								 vc1=0.

c								 if(y(j).le.end_y/2.)
                                    call solit(uc,vc,zc,c*co,co,t(n),
     -                                    x0,x(i),y(j)-end_y/2.,
     -                                    depth,alp,inc_ang,
     -                                    wave_type,wave_hgt,
     -                                    L,bf_ratio,
     -									bet,cur_level,num_levels)


                                    depth=tmp

c                                    if(stem_wv.eq.2)then
c                                    if(y(j).ge.end_y/2.)
c     -                              call solit(uc1,vc1,zc1,c*co,co,t(n),
c     -                                    x0,x(i),y(j)-end_y/2.,
c     -                                    depth,alp,-inc_ang,
c     -                                    wave_type,wave_hgt,
c     -                                    L,bf_ratio,
c     -									bet,cur_level,num_levels)
c                                    endif
       if(cur_level.eq.1) zeta(i,j,n,cur_level)=zc+zc1
       if(cur_level.eq.1.and.num_levels.gt.1) 
     -                              zeta(i,j,n,cur_level)=
     -                                    (zc+zc1)*1.025

                                    v(i,j,n,cur_level)=vc+vc1
                                    u(i,j,n,cur_level)=uc+uc1

                              elseif(sim_opt.eq.1.and.
     -                                    (wave_type.le.5))then
C*****************If using cnoidal wave, calculate initial wave profile

                                    call cnoidal(wave_hgt/depth,mk2,
     -                                    x(i),y(j),
     -                                    t(n),zc,uc,vc,yt,
     -                                    yt2,rp1,r2,CK,depth,
     -                                    L,per,co,inc_ang,
     -                                  wave_type,bet(cur_level)*depth)
                                    dom_loc=dom_wave*L+c*co*dt*(n-1)
                                    ramp_loc=(dom_wave+ramp)*L+
     -                                    c*co*dt*(n-1)
                                    cur_loc=x(i)+c*co*dt*(n-1)
                                    
                                    shift_cur=0
CCCCCCCCCCCCCCCCCCCCCCCC special for Sangsoo
                                    tmp=ramp*L+c*co*dt*(n-1)
                                    shift_cur=nint(5*L/dx)
                                    if(cur_loc.lt.tmp)then
                                          zc=zc*exp(-2*3.1415*
     -                                    (tmp-cur_loc)/(ramp*L))
                                          uc=uc*exp(-2*3.1415*
     -                                    (tmp-cur_loc)/(ramp*L))
                                          vc=vc*exp(-2*3.1415*
     -                                    (tmp-cur_loc)/(ramp*L))
CCCCCCCCCCCCCCCCCCCCCCCC end special                                                
                                    elseif(cur_loc.lt.dom_loc)then
                                          zc=zc
                                          uc=uc
                                          vc=vc
                                    else
                                          zc=zc*exp(-2*3.1415*
     -                                    (cur_loc-dom_loc)/(ramp*L))
                                          uc=uc*exp(-2*3.1415*
     -                                    (cur_loc-dom_loc)/(ramp*L))
                                          vc=vc*exp(-2*3.1415*
     -                                    (cur_loc-dom_loc)/(ramp*L))      
                                          
                                    endif
      
c       if(cur_level.eq.1) zeta(i+shift_cur,j,n,cur_level)=zc
c                                u(i+shift_cur,j,n,cur_level)=uc
c                                v(i+shift_cur,j,n,cur_level)=vc

                                    call cnoidal(wave_hgt_2/depth,mk2,
     -                                    x(i),y(j),
     -                                    t(n),zc,uc,vc,yt,
     -                                    yt2,rp1,r2,CK,depth,
     -                                    L_2,per_2,co,inc_ang,
     -                                   wave_type,bet(cur_level)*depth)
                                    dom_loc=dom_wave*L+c*co*dt*(n-1)
                                    ramp_loc=(dom_wave+ramp)*L+
     -                                    c*co*dt*(n-1)
                                    cur_loc=x(i)+c*co*dt*(n-1)
                                    if(x(i).lt.dom_loc)then
                                          zc=zc
                                          uc=uc
                                          vc=vc
                                   elseif(x(i).lt.ramp_loc)then
                                          zc=zc*(ramp_loc-cur_loc)/
     -                                    (ramp*L)
                                          uc=uc*(ramp_loc-cur_loc)/
     -                                    (ramp*L)
                                          vc=vc*(ramp_loc-cur_loc)/
     -                                    (ramp*L)                              
                                    else
                                          zc=0
                                          uc=0
                                          vc=0
                                    endif

c                                    if(cur_level.eq.1.and.
c     -                                    wave_hgt_2/wave_hgt.gt.0.001) 
c     -                                 zeta(i,j,n,cur_level)=zc
c
c                                    v(i,j,n,cur_level)=
c     -                                    v(i,j,n,cur_level)+vc
c                                    u(i,j,n,cur_level)=
c     -                                    u(i,j,n,cur_level)+uc

                              endif            
                        enddo
                  enddo
            endif

cZZZZZZZZZZ
            enddo
cZZZZZZZZZZ


      return

      end
