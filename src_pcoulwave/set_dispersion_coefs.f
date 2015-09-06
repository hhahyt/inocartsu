      subroutine set_dispersion_coefs

      use mainvar_module

C...Use Nwogu's constants when disp_prop=1, or depth-average when disp_prop=2..

cZZZZZZZZZZ
      do cur_level=1,num_levels
cZZZZZZZZZZ
            if(disp_prop.eq.1)then
                  if(num_levels.eq.1)then

                        bet2(cur_level)=0.
                        bet(cur_level)=-0.531
                        xi(cur_level)=0. !1+bet(cur_level)
                        do i=1,nx
                              do j=1,ny
                                    a1(i,j,cur_level)=
     -                             bet(cur_level)**2./2.-1./6.
                                    a2(i,j,cur_level)=
     -                             bet(cur_level)+1./2.
                                    b1(i,j,cur_level)=
     -                             bet(cur_level)**2./2.
                                    b2(i,j,cur_level)=bet(cur_level)      
                              enddo
                        enddo


                  elseif(num_levels.eq.2)then
CCC Coef
CCC%%%% 4,4 Pade
CCCCCCC -.2483      -.4586      -.7413

CCC%%%% Opt to 2pi
CCCCCCC -.202            -.38      -.684

CCC%%%% Opt to 3pi
CCCCCCC -.173            -.327      -.646

CCC%%%% Opt to 4pi
CCCCCCC -.15            -.284      -.616

CCC%%%% Best Linear Opt
CCCCCCC -.168            -.317      -.63


           

                         do i=1,nx
                                do j=1,ny
c                                    if(ho(1,i,j).ge.0)then
c                                          bet(1)=-.2483
c                                          bet(2)=-.7413
c                                          bet2(1)=0
c                                          bet2(2)=-.4586

c                                          bet(1)=-.20
c                                          bet(2)=-.68
c                                          bet2(1)=0
c                                          bet2(2)=-.376

                                          bet(1)=-.127
                                          bet(2)=-.618
                                          bet2(1)=0
                                          bet2(2)=-.266

c                                          bet(1)=-.092
c                                          bet(2)=-.625
c                                          bet2(1)=0
c                                          bet2(2)=-.256


CCCCCCCCCCCCC Nonlinear Opt            

                                          xi(1)=0.3
                                          xi(2)=0.
                                          xi2(1)=0.
                                          xi2(2)=-0.2

c                                          xi(1)=-bet(1)
c                                          xi(2)=-bet(2)
c                                          xi2(1)=0
c                                          xi2(2)=-bet2(2)

c                                          xi(1)=0.031
c                                          xi(2)=-0.063
c                                          xi2(1)=0
c                                          xi2(2)=-.001

c                                          xi(1)=0.135
c                                          xi(2)=-.045
c                                          xi2(1)=0.
c                                          xi2(2)=0.

c                                          xi(1)=0.1
c                                          xi(2)=0
c                                          xi2(1)=0.05
c                                          xi2(2)=0

c                                    else
c                                          bet(1)=-1.
c                                          bet(2)=-1.
c                                          bet2(1)=0.
c                                          bet2(2)=-1.
c
c                                          xi(1)=0
c                                          xi(2)=0
c                                          xi2(1)=0
c                                          xi2(2)=0
c                                    endif

                              enddo                        
                        enddo


                  elseif(num_levels.eq.3)then

                         do i=1,nx
                                do j=1,ny
                                    if(ho(i,j).ge.0)then

c                                          bet(1)=-0.1008708
c                                          bet(2)=-0.3608673
c                                          bet(3)=-0.7888216
c
c                                          bet2(1)=0.
c                                          bet2(2)=-0.1888076
c                                          bet2(3)=-0.5511979

                                          bet(1)=-.105
                                          bet(2)=-.3675
                                          bet(3)=-.7875

                                          bet2(1)=0
                                          bet2(2)=-.1965
                                          bet2(3)=-.555


c                                          bet(1)=-.13145
c                                          bet(2)=-.42545
c                                          bet(3)=-.83693
c
c                                          bet2(1)=0
c                                          bet2(2)=-.24292
c                                          bet2(3)=-.62015

                                    else
                                          bet(1)=-1.
                                          bet(2)=-1.
                                          bet(3)=-1.
                                          bet2(1)=0.
                                          bet2(2)=-1.
                                          bet2(3)=-1.

                                    endif
                                    if(cur_level.eq.1)then
                                          a1(i,j,cur_level)=
     -                                           -bet2(cur_level+1)*
     -                                           bet(cur_level)**2./2.+
     -                                           bet2(cur_level+1)**2.*
     -                                           bet2(cur_level+1)/6.
                                          a2(i,j,cur_level)=
     -                                           -bet2(cur_level+1)*
     -                                           bet(cur_level)+
     -                                          bet2(cur_level+1)**2./2.
                                          b1(i,j,cur_level)=
     -                                           bet(cur_level)**2./2.-
     -                                           bet2(cur_level+1)*
     -                                           bet(cur_level)
                                          b2(i,j,cur_level)=0

                                          
                                          NL_a1(i,j,cur_level)=0
                                          NL_a2(i,j,cur_level)=0
                                          NL_b1(i,j,cur_level)=0
                                          NL_b2(i,j,cur_level)=0 


                                    elseif(cur_level.eq.2)then
                                          a1(i,j,cur_level)=
     -                                           (bet2(cur_level)+1.)*
     -                                           bet(cur_level)**2./2.-
     -                                           (bet2(cur_level)**2.*
     -                                           bet2(cur_level)+1)/6.
                                          a2(i,j,cur_level)=
     -                                           (bet2(cur_level)+1.)*
     -                                           bet(cur_level)-
     -                                       (bet2(cur_level)**2.-1.)/2.
                                          b1(i,j,cur_level)=
     -                                           bet(cur_level)**2./2.+
     -                                         bet2(cur_level)**2./2.-
     -                                         bet2(cur_level)*
     -                                          bet(cur_level-1)+
     -                                         bet2(cur_level+1)*
     -                                         bet(cur_level-1)-
     -                                         bet2(cur_level+1)*
     -                                         bet(cur_level)

                                          b2(i,j,cur_level)=0

                                          NL_a1(i,j,cur_level)=0
                                          NL_a2(i,j,cur_level)=0
                                          NL_b1(i,j,cur_level)=0
                                          NL_b2(i,j,cur_level)=0

                                    elseif(cur_level.eq.3)then
                                          a1(i,j,cur_level)=
     -                                           (bet2(cur_level)+1.)*
     -                                           bet(cur_level)**2./2.-
     -                                           (bet2(cur_level)**2.*
     -                                           bet2(cur_level)+1)/6.
                                          a2(i,j,cur_level)=
     -                                           (bet2(cur_level)+1.)*
     -                                           bet(cur_level)-
     -                                       (bet2(cur_level)**2.-1.)/2.
                                          b1(i,j,cur_level)=
     -                                           bet(cur_level)**2./2.+
     -                                           bet2(cur_level)**2./2.-
     -                                          bet2(cur_level)*
     -                                          bet(cur_level-1)
                                          b2(i,j,cur_level)=
     -                                           bet(cur_level)-
     -                                           bet(cur_level-1)

                                          NL_a1(i,j,cur_level)=0
                                          NL_a2(i,j,cur_level)=0
                                          NL_b1(i,j,cur_level)=0
                                          NL_b2(i,j,cur_level)=0
                                    endif
                              enddo                        
                        enddo

c                                          bet(1)=-0.1008708
c                                          bet(2)=-0.3608673
c                                          bet(3)=-0.7888216
c
c                                          bet2(1)=0.
c                                          bet2(2)=-0.1888076
c                                          bet2(3)=-0.5511979

                                          bet(1)=-.105
                                          bet(2)=-.3675
                                          bet(3)=-.7875

                                          bet2(1)=0
                                          bet2(2)=-.1965
                                          bet2(3)=-.555

c                        bet(1)=-.13145
c                        bet(2)=-.42545
c                        bet(3)=-.83693
c
c                        bet2(1)=0
c                        bet2(2)=-.24292
c                        bet2(3)=-.62015

                  elseif(num_levels.eq.4)then

                         do i=1,nx
                                do j=1,ny
                                    if(ho(i,j).ge.0)then

c                                          bet(1)=-.083
c                                          bet(2)=-.315
c                                          bet(3)=-.748
c
c                                          bet2(1)=0
c                                          bet2(2)=-.156
c                                          bet2(3)=-.494


                                          bet(1)=-.06
                                          bet(2)=-.22
                                          bet(3)=-.5165
                                          bet(4)=-.88

                                          bet2(1)=0
                                          bet2(2)=-.1125
                                          bet2(3)=-.34
                                          bet2(4)=-.70

                                          xi(1)=0
                                          xi(2)=0
                                          xi(3)=0
                                          xi(4)=0

                                          xi2(1)=0
                                          xi2(2)=0
                                          xi2(3)=0
                                          xi2(4)=0

                                    else
                                          bet(1)=-1.
                                          bet(2)=-1.
                                          bet(3)=-1.
                                          bet(4)=-1.

                                          bet2(1)=0.
                                          bet2(2)=-1.
                                          bet2(3)=-1.
                                          bet2(4)=-1.

                                    endif
                                    if(cur_level.eq.1)then
                                          a1(i,j,cur_level)=
     -                                           -bet2(cur_level+1)*
     -                                           bet(cur_level)**2./2.+
     -                                           bet2(cur_level+1)**2.*
     -                                           bet2(cur_level+1)/6.
                                          a2(i,j,cur_level)=
     -                                           -bet2(cur_level+1)*
     -                                            bet(cur_level)+
     -                                         bet2(cur_level+1)**2./2.
                                          b1(i,j,cur_level)=
     -                                           bet(cur_level)**2./2.-
     -                                           bet2(cur_level+1)*
     -                                           bet(cur_level)
                                          b2(i,j,cur_level)=0

                                          
                                          NL_a1(i,j,cur_level)=0
                                          NL_a2(i,j,cur_level)=0
                                          NL_b1(i,j,cur_level)=0
                                          NL_b2(i,j,cur_level)=0 

                                    elseif(cur_level.eq.num_levels)then
                                          a1(i,j,cur_level)=
     -                                           (bet2(cur_level)+1.)*
     -                                            bet(cur_level)**2./2.-
     -                                           (bet2(cur_level)**2.*
     -                                            bet2(cur_level)+1)/6.
                                          a2(i,j,cur_level)=
     -                                           (bet2(cur_level)+1.)*
     -                                           bet(cur_level)-
     -                                      (bet2(cur_level)**2.-1.)/2.
                                          b1(i,j,cur_level)=
     -                                           bet(cur_level)**2./2.+
     -                                         bet2(cur_level)**2./2.-
     -                                         bet2(cur_level)*
     -                                          bet(cur_level-1)
                                          b2(i,j,cur_level)=
     -                                           bet(cur_level)-
     -                                           bet(cur_level-1)

                                          NL_a1(i,j,cur_level)=0
                                          NL_a2(i,j,cur_level)=0
                                          NL_b1(i,j,cur_level)=0
                                          NL_b2(i,j,cur_level)=0

                                    else
                                          a1(i,j,cur_level)=
     -                                           (bet2(cur_level)+1.)*
     -                                           bet(cur_level)**2./2.-
     -                                           (bet2(cur_level)**2.*
     -                                           bet2(cur_level)+1)/6.
                                          a2(i,j,cur_level)=
     -                                           (bet2(cur_level)+1.)*
     -                                           bet(cur_level)-
     -                                       (bet2(cur_level)**2.-1.)/2.
                                          b1(i,j,cur_level)=
     -                                           bet(cur_level)**2./2.+
     -                                          bet2(cur_level)**2./2.-
     -                                          bet2(cur_level)*
     -                                           bet(cur_level-1)+
     -                                          bet2(cur_level+1)*
     -                                          bet(cur_level-1)-
     -                                          bet2(cur_level+1)*
     -                                          bet(cur_level)

                                          b2(i,j,cur_level)=0

                                          NL_a1(i,j,cur_level)=0
                                          NL_a2(i,j,cur_level)=0
                                          NL_b1(i,j,cur_level)=0
                                          NL_b2(i,j,cur_level)=0

                                    endif
                              enddo                        
                        enddo

c                        bet(1)=-.083
c                        bet(2)=-.315
c                        bet(3)=-.748
c
c                        bet2(1)=0
c                        bet2(2)=-.156
c                        bet2(3)=-.494


                                          bet(1)=-.06
                                          bet(2)=-.22
                                          bet(3)=-.5165
                                          bet(4)=-.88

                                          bet2(1)=0
                                          bet2(2)=-.1125
                                          bet2(3)=-.34
                                          bet2(4)=-.70

                  endif
              alp=0.5*bet(1)**2.+bet(1)
            elseif(disp_prop.eq.2)then
                  alp=-.33333
                  bet(cur_level)=-0.5
                  bet2(cur_level)=0.
                  do i=1,nx
                        do j=1,ny
                              a1(i,j,cur_level)=0
                              a2(i,j,cur_level)=0
                              b1(i,j,cur_level)=1./6.
                              b2(i,j,cur_level)=-1./2.
                        enddo
                  enddo
            endif


                  do n=1,4
                            do i=1,endx
                              do j=1,endy
                                    h(i,j,n)=ho(i,j)
                                    hp(i,j,n)=ho(i,j)
                                      z_alp(i,j,n,cur_level)=
     -                                  bet(cur_level)*h(i,j,n)
                                      zeta(i,j,n,cur_level)=
     -                                 bet2(cur_level)*h(i,j,n)
                                      zeta(i,j,n,num_levels+1)=
     -                                             -h(i,j,n)
                               enddo
                        enddo
                  enddo

cZZZZZZZ
      enddo
cZZZZZZZ



cZZZZZZZZZZ
      do cur_level=1,num_levels
cZZZZZZZZZZ
      do n=1,4
               do i=1,endx
                  do j=1,endy
                        ho(i,j)=h(i,j,n)
                          z_alp(i,j,n,cur_level)=
     -                        bet(cur_level)*h(i,j,n)
                          zeta(i,j,n,cur_level)=
     -                        bet2(cur_level)*h(i,j,n)
                          zeta(i,j,n,num_levels+1)=
     -                        -h(i,j,n)
                   enddo
            enddo
      enddo
cZZZZZZZZZZ
      enddo
cZZZZZZZZZZ




      do i=2,endx-1
            do j=2,endy-1
                  if(bl_hor_wall(i,j).le.50)then

cZZZZZZZZZZ
                  do cur_level=1,num_levels
cZZZZZZZZZZ
                       d1u(i,j,cur_level)=ho(i,j)**2.*
     -                              (b1(i,j,cur_level)+
     -                             b2(i,j,cur_level))/(dx**2.)
                       d2u(i,j,cur_level)=(1+b2(i,j,cur_level)*
     -                              ho(i,j)*hxx(i,j))-
     -                              2.*ho(i,j)**2.
     -                             *(b1(i,j,cur_level)+
     -                              b2(i,j,cur_level))/(dx**2.)
                       d3u(i,j,cur_level)=2.*b2(i,j,cur_level)*
     -                              ho(i,j)*
     -                             hx(i,j)/(2.*dx)

C............ d1v and d2v are used to convert VV to v
                       d1v(i,j,cur_level)=ho(i,j)**2.*
     -                              (b1(i,j,cur_level)+
     -                             b2(i,j,cur_level))/(dy**2.)
                       d2v(i,j,cur_level)=(1+b2(i,j,cur_level)*
     -                              ho(i,j)*hyy(i,j))-
     -                              2.*ho(i,j)**2.*
     -                             (b1(i,j,cur_level)+
     -                              b2(i,j,cur_level))/(dy**2.)
                       d3v(i,j,cur_level)=2.*b2(i,j,cur_level)*
     -                              ho(i,j)*
     -                             hy(i,j)/(2.*dy)
c                  endif
cZZZZZZZZZ
                  enddo
cZZZZZZZZZ
              endif
            enddo
      enddo      









      return

      end

