      subroutine set_internal_source_coefs

      use mainvar_module

	  
	if(int_src.eq.1.and.num_levels.ge.3) int_src=2

	if(current.eq.1)then
	
            mod_ang=abs(inc_ang-90*nint(inc_ang/90))*3.1415/180.
            inc_ang_cur=inc_ang*3.1415/180.
		  x0_cur=x0

		  L_spec_cur=20*depth

            wid=L_spec_cur
            beta_src_cur=80./(wid)**2.

            I_src_cur=sqrt(3.14/beta_src_cur)

            amp_cur=depth*Fr_cur

            D_src_cur=2.*amp_cur*cos(mod_ang)*
     -                              sqrt(9.81*depth)/(I_src_cur)

	endif


c	open(60,file='tide.dat',status='old') ! TEMP FOR IPET-MRGO SURGE
c	read(60,*) D_src_nl(1,1)
c	close(60)

	D_src_nl(1,1)=0.0

c	open(60,file='f_bf.dat',status='old') ! TEMP FOR IPET-MRGO SURGE
c	read(60,*) f_BF
c	close(60)


	depth_tmp=depth
	depth=D_src_nl(1,1)+depth

	D_src_nl(1,1)=D_src_nl(1,1)*
     -     sqrt(1+0.5*D_src_nl(1,1)/(D_src_nl(1,1)+depth_tmp))
c     /
c     -	 (depth_tmp/(D_src_nl(1,1)+depth_tmp))**0.12

	if(D_src_nl(1,1).gt.1e-6)then
            I_src_nl=sqrt(3.14*depth_tmp**2.)

            D_src_nl(1,1)=2.*D_src_nl(1,1)*cos(mod_ang)*
     -                         sqrt(9.81*depth_tmp)/(I_src_nl)
	endif


      if(int_src.eq.1)then
        if(wave_type.eq.5)then
            mod_ang=abs(inc_ang-90*nint(inc_ang/90))*3.1415/180.
            inc_ang=inc_ang*3.1415/180.

            L=2.*depth
            step=L/10.
            shift=0
            do while(step.gt.L/100000)
                        tmp=sqrt(2*3.14*L/9.81*1./tanh(2.*3.14*depth/L))
                        if(tmp.gt.per)then
                              L=L-step
                              if(shift.eq.1)then
                                 shift=0
                                 step=step/10
                              else
                                 shift=-1
                              endif
                        else
                              L=L+step
                              if(shift.eq.-1)then
                                 shift=0
                                 step=step/10
                              else
                                 shift=1
                              endif
                        endif
            enddo

            L_2=2.*depth
            step=L_2/10.
            shift=0
            do while(step.gt.L_2/100000)
                    tmp=sqrt(2*3.14*L_2/9.81*1./tanh(2.*3.14*depth/L_2))
                        if(tmp.gt.per_2)then
                              L_2=L_2-step
                              if(shift.eq.1)then
                                 shift=0
                                 step=step/10
                              else
                                 shift=-1
                              endif
                        else
                              L_2=L_2+step
                              if(shift.eq.-1)then
                                 shift=0
                                 step=step/10
                              else
                                 shift=1
                              endif
                        endif
            enddo

            wid=max(40*dx,1.2*L)
            beta_src(1,1)=80./(wid)**2.

            wid=max(40*dx,1.2*L_2)
            beta_src(2,1)=80./(wid)**2.

            k_1=2*3.15/L
            L_spec(1,1)=L
            per_spec(1,1)=per
            amp_spec(1,1)=wave_hgt
            inc_ang_spec(1,1)=inc_ang
		  shift_spec(1,1)=0.

            k_2=2*3.15/L_2
            L_spec(2,1)=L_2
            per_spec(2,1)=per_2
            amp_spec(2,1)=wave_hgt_2
            inc_ang_spec(2,1)=inc_ang
		  shift_spec(2,1)=0.

            I_src(1,1)=sqrt(3.14/beta_src(1,1))*
     -                 exp(-k_1**2./4./beta_src(1,1))
            I_src(2,1)=sqrt(3.14/beta_src(2,1))*
     -                 exp(-k_2**2./4./beta_src(2,1))

            if(num_levels.eq.1)then
                  D_src(1,1)=2.*wave_hgt*cos(mod_ang)*
     -            (w_1**2.-(1./3.+alp)*
     -                 9.81*k_1**2.*k_1**2.*depth**2.*depth)/
     -            (w_1*k_1*I_src(1,1)*(1-alp*(k_1*depth)**2.))


                  D_src(2,1)=2.*wave_hgt_2*cos(mod_ang)*
     -            (w_2**2.-(1./3.+alp)*
     -                 9.81*k_2**2.*k_2**2.*depth**2.*depth)/
     -            (w_2*k_2*I_src(2,1)*(1-alp*(k_2*depth)**2.))


            elseif(num_levels.eq.2)then

              d1=-bet2(2)
              d2=1+bet2(2)
              d3=1./6.*(-2.*bet2(2)**2*bet2(2)+
     -                  6.*bet(1)*bet2(2)**2-3.*bet(1)**2*bet2(2))
              d4=1./6.*(2.*bet2(2)**2*bet2(2)-
     -                  6.*bet(1)*bet2(2)**2-6.*bet(1)*bet2(2)+
     -                  3.*bet(2)**2*bet2(2)+6.*bet(2)*bet2(2)+
     -                  3.*bet(2)**2+6.*bet(2)+2)
              d5=bet(1)**2/2.-bet(1)*bet2(2)
              d6=bet(1)*bet2(2)+bet(1)
              d7=-bet(1)**2./2.+bet(1)*bet2(2) -bet2(2)**2./2.
              d8=-bet(1)*bet2(2)+bet(2)**2./2.-
     -                  bet(1)+bet(2)+bet2(2)**2./2.

              Num2=d2*d7-d1*d8-d3-d4
              Num3=d3*d8-d4*d7
              Den1=d8+d5+d6
              Den2=-d5*d8+d6*d7

			w_1=sqrt(k_1**2.*9.81*depth*(1+k_1**2.*depth**2.*Num2+
     -				Num3*k_1**2.*k_1**2.*depth**2.*depth**2.)/
     -				(1-k_1**2.*depth**2.*Den1-
     -				Den2*k_1**2.*k_1**2.*depth**2.*depth**2.))


c           D_src(1,1)=.5*wave_hgt*w_1/k_1/I_src(1,1)

              D_src(1,1)=2.*wave_hgt*cos(mod_ang)*
     -            (w_1**2.+Num2*9.81*k_1**2.*k_1**2.*depth**2.*depth+
     -            Num3*9.81*k_1**2.*k_1**2.*
     -            k_1**2.*depth**2.*depth**2.*depth)/
     -            (w_1*k_1*I_src(1,1)*(1-Den1*k_1**2.*depth**2.-
     -            2./sqrt(-bet(2))*Den2*k_1**2.*
     -                 k_1**2.*depth**2.*depth**2.))


              D_src(2,1)=2.*wave_hgt_2*cos(mod_ang)*
     -            (w_2**2.+Num2*9.81*k_2**2.*k_2**2.*depth**2.*depth+
     -            Num3*9.81*k_2**2.*k_2**2.*
     -            k_2**2.*depth**2.*depth**2.*depth)/
     -            (w_2*k_2*I_src(2,1)*(1-Den1*k_2**2.*depth**2.-
     -            2./sqrt(-bet(2))*Den2*k_2**2.*
     -                 k_2**2.*depth**2.*depth**2.))

				
            elseif(num_levels.eq.3)then

                  Num2=(((3 * bet(3)**2 + 6 * bet(3) + 3 * bet2(3)**2
     -                        - 6 * bet(2) * bet2(3) + 3 *
     -                        bet(2)**2 + 3 * bet2(2)**2
     -                        - 6 * bet(1) * bet2(2) +
     -                        3 * bet(1)**2 + 2))/6.)
                  Num3= -(((bet2(3)**3 * bet(3)**2 + 3 * bet2(3)**2 *
     -                    bet(3)**2 - 6 * bet(2) * bet2(3) * bet(3)**2 +
     -                    3 * bet(2)**2 * bet(3)**2 + 3 * bet2(2)**2 *
     -                    bet(3)**2 - 6 * bet(1) * bet2(2) * bet(3)**2 +
     -                  3 * bet(1)**2 * bet(3)**2 + 2 * bet2(3)**3 *
     -                        bet(3) + 6 * bet2(3)**2 * bet(3) -
     -                 12 * bet(2) * bet2(3) * bet(3) + 6 * bet(2)**2
     -                        * bet(3) + 6 * bet2(2)**2 * bet(3) -
     -                  12 * bet(1) * bet2(2) * bet(3) + 6 * bet(1)**2 *
     -                         bet(3) - bet(2)**2 * bet2(3)**3 -
     -                        2 * bet(2) * bet2(3)**3 + 3 * bet2(2)**2 *
     -                 bet2(3)**2 - 6 * bet(1) * bet2(2) * bet2(3)**2 +
     -                  3 * bet(1)**2 * bet2(3)**2 + 2 * bet2(3)**2 - 6
     -                        * bet2(2)**2 * bet(2) * bet2(3) +
     -                  12 * bet(1) * bet2(2) * bet(2) * bet2(3) - 6 *
     -                        bet(1)**2 * bet(2) * bet2(3) -
     -                  4 * bet(2) * bet2(3) + bet2(2)**3 * bet(2)**2 +
     -                         3 * bet2(2)**2 * bet(2)**2 -
     -                  6 * bet(1) * bet2(2) * bet(2)**2 + 3 * bet(1)**2
     -                         * bet(2)**2 + 2 * bet(2)**2 +
     -                  2 * bet2(2)**3 * bet(2) - bet(1)**2 * bet2(2)**3
     -                         - 2 * bet(1) * bet2(2)**3 +
     -                  2 * bet2(2)**2 - 4 * bet(1) * bet2(2) + 2 *
     -                        bet(1)**2))/12.)
                  Num4=(((bet2(2)**2 * bet2(3)**3
     -                        * bet(3)**2 - 2 * bet(1)
     -                        * bet2(2) * bet2(3)**3 * bet(3)**2 +
     -                        bet(1)**2 * bet2(3)**3 * bet(3)**2 + 3 *
     -                         bet2(2)**2 * bet2(3)**2 * bet(3)**2 -
     -                        6 * bet(1) * bet2(2) *
     -                        bet2(3)**2 * bet(3)**2
     -                        + 3 * bet(1)**2 * bet2(3)**2 * bet(3)**2 -
     -                        2 * bet2(2)**3 * bet(2)
     -                        * bet2(3) * bet(3)**2 -
     -                         6 * bet2(2)**2 * bet(2)
     -                        * bet2(3) * bet(3)**2 +
     -                        12 * bet(1) * bet2(2) * bet(2) * bet2(3) *
     -                        bet(3)**2 - 6 * bet(1)**2 *
     -                         bet(2) * bet2(3) * bet(3)**2 +
     -                     2 * bet(1) * bet2(2)**3 * bet2(3) * bet(3)**2
     -                         + bet2(2)**3 * bet(2)**2 * bet(3)**2 +
     -                        3 * bet2(2)**2 * bet(2)**2 * bet(3)**2 - 6
     -                     * bet(1) * bet2(2) * bet(2)**2 * bet(3)**2 +
     -                3 * bet(1)**2 * bet(2)**2 * bet(3)**2 - bet(1)**2
     -                        * bet2(2)**3 * bet(3)**2 +
     -                2 * bet2(2)**2 * bet2(3)**3 * bet(3) - 4 * bet(1)
     -                        * bet2(2) * bet2(3)**3 * bet(3) +
     -                        2 * bet(1)**2 * bet2(3)**3 * bet(3) + 6 *
     -                         bet2(2)**2 * bet2(3)**2 * bet(3) -
     -                       12 * bet(1) * bet2(2) * bet2(3)**2 * bet(3)
     -                        + 6 * bet(1)**2 * bet2(3)**2 * bet(3) -
     -                     4 * bet2(2)**3 * bet(2) * bet2(3) * bet(3) -
     -                    12 * bet2(2)**2 * bet(2) * bet2(3) * bet(3) +
     -                 24 * bet(1) * bet2(2) * bet(2) * bet2(3) * bet(3)
     -                   - 12 * bet(1)**2 * bet(2) * bet2(3) * bet(3) +
     -                     4 * bet(1) * bet2(2)**3 * bet2(3) * bet(3) +
     -                         2 * bet2(2)**3 * bet(2)**2 * bet(3) +
     -                        6 * bet2(2)**2 * bet(2)**2 * bet(3) - 12 *
     -                         bet(1) * bet2(2) * bet(2)**2 * bet(3) +
     -                        6 * bet(1)**2 * bet(2)**2 * bet(3) - 2 *
     -                        bet(1)**2 * bet2(2)**3 * bet(3) -
     -                        bet2(2)**2 * bet(2)**2 * bet2(3)**3 + 2 *
     -                       bet(1) * bet2(2) * bet(2)**2 * bet2(3)**3 -
     -                        bet(1)**2 * bet(2)**2 * bet2(3)**3 - 2 *
     -                        bet2(2)**2 * bet(2) * bet2(3)**3 +
     -                        4 * bet(1) * bet2(2) * bet(2) * bet2(3)**3
     -                        - 2 * bet(1)**2 * bet(2) * bet2(3)**3 +
     -                        bet2(2)**3 * bet(2)**2 * bet2(3)**2 + 2 *
     -                        bet2(2)**3 * bet(2) * bet2(3)**2 -
     -                        bet(1)**2 * bet2(2)**3 * bet2(3)**2 - 2 *
     -                         bet(1) * bet2(2)**3 * bet2(3)**2 +
     -                        2 * bet2(2)**2 * bet2(3)**2 - 4 * bet(1)
     -                        * bet2(2) * bet2(3)**2 +
     -                        2 * bet(1)**2 * bet2(3)**2 - 2 * bet(1)
     -                        * bet2(2)**3 * bet(2)**2 * bet2(3) +
     -                        2 * bet(1)**2 * bet2(2)**3 * bet(2) *
     -                     bet2(3) - 4 * bet2(2)**2 * bet(2) * bet2(3) +
     -                        8 * bet(1) * bet2(2) * bet(2) * bet2(3) -
     -                        4 * bet(1)**2 * bet(2) * bet2(3) -
     -                        2 * bet(1) * bet2(2)**3 * bet(2)**2 + 2
     -                        * bet2(2)**2 * bet(2)**2 -
     -                        4 * bet(1) * bet2(2) * bet(2)**2 + 2 *
     -                        bet(1)**2 * bet(2)**2 +
     -                        2 * bet(1)**2 * bet2(2)**3 * bet(2)))/24.)


                   Den1= - (((bet(3)**2 + 2 * bet(3) + bet2(3)**2
     -                        - 2 * bet(2) * bet2(3) +
     -                         bet(2)**2 + bet2(2)**2 - 2 * bet(1) *
     -                         bet2(2) + bet(1)**2))/2.)
                   Den2= - (((2 * bet(2) * bet2(3) * bet(3)**2 -
     -                         bet(2)**2 * bet(3)**2 -
     -                        bet2(2)**2 * bet(3)**2 + 2 * bet(1) *
     -                        bet2(2) * bet(3)**2 -
     -                         bet(1)**2 * bet(3)**2 + 4 * bet(2) *
     -                        bet2(3) * bet(3) -
     -                        2 * bet(2)**2 * bet(3) - 2 * bet2(2)**2
     -                        * bet(3) + 4 * bet(1) * bet2(2) * bet(3) -
     -                        2 * bet(1)**2 * bet(3) - bet(2)**2 *
     -                        bet2(3)**2 - 2 * bet(2) * bet2(3)**2 -
     -                        bet2(2)**2 * bet2(3)**2 + 2 * bet(1) *
     -                   bet2(2) * bet2(3)**2 - bet(1)**2 * bet2(3)**2 +
     -                        2 * bet2(2)**2 * bet(2) * bet2(3) - 4 *
     -                         bet(1) * bet2(2) * bet(2) * bet2(3) +
     -                        2 * bet(1)**2 * bet(2) * bet2(3) + 2 *
     -                        bet(1) * bet2(2) * bet(2)**2 -
     -                         bet(1)**2 * bet(2)**2 + 2 * bet2(2)**2
     -                         * bet(2) - bet(1)**2 * bet2(2)**2 -
     -                         2 * bet(1) * bet2(2)**2))/4.)
                   Den3= - ((bet(1) * (4 * bet2(2) * bet(2) *
     -                         bet2(3) * bet(3)**2 -
     -                        2 * bet(1) * bet(2) * bet2(3) * bet(3)**2
     -                         - 2 * bet2(2)**2 * bet2(3) * bet(3)**2 -
     -                        2 * bet2(2) * bet(2)**2 * bet(3)**2 +
     -                        bet(1) * bet(2)**2 * bet(3)**2 +
     -                        bet(1) * bet2(2)**2 * bet(3)**2 + 8 *
     -                        bet2(2) * bet(2) * bet2(3) * bet(3) -
     -                        4 * bet(1) * bet(2) * bet2(3) * bet(3)
     -                         - 4 * bet2(2)**2 * bet2(3) * bet(3) -
     -                        4 * bet2(2) * bet(2)**2 * bet(3) + 2 *
     -                         bet(1) * bet(2)**2 * bet(3) +
     -                        2 * bet(1) * bet2(2)**2 * bet(3) - 2 *
     -                         bet2(2) * bet(2)**2 * bet2(3)**2 +
     -                        bet(1) * bet(2)**2 * bet2(3)**2 - 4 *
     -                        bet2(2) * bet(2) * bet2(3)**2 +
     -                        2 * bet(1) * bet(2) * bet2(3)**2 +
     -                         bet(1) * bet2(2)**2 * bet2(3)**2 +
     -                        2 * bet2(2)**2 * bet2(3)**2 + 2 *
     -                        bet2(2)**2 * bet(2)**2 * bet2(3) -
     -                        2 * bet(1) * bet2(2)**2 * bet(2) *
     -                        bet2(3) + 2 * bet2(2)**2 * bet(2)**2 -
     -                        2 * bet(1) * bet2(2)**2 * bet(2)))/8.)


              D_src(1,1)=2.*wave_hgt*cos(mod_ang)*
     -            (w_1**2.+Num2*9.81*k_1**2.*k_1**2.*depth**2.*depth+
     -            Num3*9.81*k_1**2.*k_1**2.*
     -            k_1**2.*depth**2.*depth**2.*depth+
     -            Num4*9.81*k_1**2.*k_1**2.*
     -            k_1**2.*k_1**2.*depth**2.*depth**2.*depth**2.*depth)/
     -           (w_1*k_1*I_src(1,1)*(1+(-bet(2))*Den1*k_1**2.*depth**2.
     -            +Den2*k_1**2.*k_1**2.*depth**2.*depth**2.+
     -             Den3*k_1**2.*k_1**2.*k_1**2.*
     -              depth**2.*depth**2.*depth**2.))



              D_src(2,1)=2.*wave_hgt_2*cos(mod_ang)*
     -            (w_2**2.+Num2*9.81*k_2**2.*k_2**2.*depth**2.*depth+
     -            Num3*9.81*k_2**2.*k_2**2.*
     -            k_2**2.*depth**2.*depth**2.*depth+
     -            Num4*9.81*k_2**2.*k_2**2.*
     -            k_2**2.*k_2**2.*depth**2.*depth**2.*depth**2.*depth)/
     -           (w_2*k_2*I_src(2,1)*(1+(-bet(2))*Den1*k_2**2.*depth**2.
     -            +Den2*k_2**2.*k_2**2.*depth**2.*depth**2.+
     -             Den3*k_2**2.*k_2**2.*k_2**2.*
     -              depth**2.*depth**2.*depth**2.))


            elseif(num_levels.eq.4)then

              D_src(1,1)=-wave_hgt*sqrt(beta_src(1,1))*w_1/
     -                  (sqrt(3.1415)*k_1)

              D_src(2,1)=-wave_hgt_2*sqrt(beta_src(2,1))*w_2/
     -                  (sqrt(3.1415)*k_2)

            endif
            shift_spec(1,1)=0
            shift_spec(2,1)=0
        elseif(wave_type.eq.6)then
		if(num_levels.eq.1)then
            do i=1,num_freq
                  do j=1,num_theta

                        inc_ang_spec(i,j)=
     -                             inc_ang_spec(i,j)*3.1415/180.

					  if(is_oreint.eq.1)then
								inc_ang_spec(i,j)=
     -								inc_ang_spec(i,j)+3.1415/2.
					  endif

                        mod_ang=abs(inc_ang_spec(i,j)-
     -                        1.57*nint(inc_ang_spec(i,j)/1.57))


                        L_spec(i,j)=2.*depth
                        step=L_spec(i,j)/10.
                        shift=0

					  if(j.eq.1)then
						    inc_ang=inc_ang_spec(i,j)
					  endif

                        do while(step.gt.L_spec(i,j)/100000)
                              tmp=sqrt(2*3.14*L_spec(i,j)/9.81*
     -                        1./tanh(2.*3.14*depth/L_spec(i,j)))
                              if(tmp.gt.per_spec(i,j))then
                                    L_spec(i,j)=L_spec(i,j)-step
                                    if(shift.eq.1)then
                                       shift=0
                                       step=step/10
                                    else
                                       shift=-1
                                    endif
                              else
                                    L_spec(i,j)=L_spec(i,j)+step
                                    if(shift.eq.-1)then
                                       shift=0
                                       step=step/10
                                    else
                                       shift=1
                                    endif
                              endif
                        enddo

                        k_tmp=2*3.1415/L_spec(i,j)
                        w_tmp=2*3.1415/per_spec(i,j)
                        beta_src(i,j)=80./(1.2*L_spec(i,j))**2.
                        I_src(i,j)=sqrt(3.14/beta_src(i,j))*
     -                        exp(-k_tmp**2./4./beta_src(i,j))
                        D_src(i,j)=4.*amp_spec(i,j)*cos(mod_ang)*
     -                        (w_tmp**2.-(1./3.+alp)*9.81*k_tmp**2.*
     -                        k_tmp**2.*depth**2.*depth)/
     -                        (w_tmp*k_tmp*I_src(i,j)*
     -                             (1-alp*(k_tmp*depth)**2.))


                       call RANDOM_NUMBER(shift_rn)
c			shift_rn=abs(real(i*j)*.384165-nint(real(i*j)*.384165))
                       shift_spec(i,j)=per_spec(i,j)*shift_rn
                  enddo
            enddo
		elseif(num_levels.eq.2)then

              d1=-bet2(2)
              d2=1+bet2(2)
              d3=1./6.*(-2.*bet2(2)**2*bet2(2)+
     -                  6.*bet(1)*bet2(2)**2-3.*bet(1)**2*bet2(2))
              d4=1./6.*(2.*bet2(2)**2*bet2(2)-
     -                  6.*bet(1)*bet2(2)**2-6.*bet(1)*bet2(2)+
     -                  3.*bet(2)**2*bet2(2)+6.*bet(2)*bet2(2)+
     -                  3.*bet(2)**2+6.*bet(2)+2)
              d5=bet(1)**2/2.-bet(1)*bet2(2)
              d6=bet(1)*bet2(2)+bet(1)
              d7=-bet(1)**2./2.+bet(1)*bet2(2) -bet2(2)**2./2.
              d8=-bet(1)*bet2(2)+bet(2)**2./2.-
     -                  bet(1)+bet(2)+bet2(2)**2./2.

              Num2=d2*d7-d1*d8-d3-d4
              Num3=d3*d8-d4*d7
              Den1=d8+d5+d6
              Den2=-d5*d8+d6*d7

			do i=1,num_freq
                  do j=1,num_theta

                        inc_ang_spec(i,j)=
     -                             inc_ang_spec(i,j)*3.1415/180.

					  if(is_oreint.eq.1)then
								inc_ang_spec(i,j)=
     -								inc_ang_spec(i,j)+3.1415/2.
					  endif

                        L_spec(i,j)=2.*depth
                        step=L_spec(i,j)/10.
                        shift=0

					  if(j.eq.1)then
						    inc_ang=inc_ang_spec(i,j)
					  endif

                        do while(step.gt.L_spec(i,j)/100000)
                              tmp=sqrt(2*3.14*L_spec(i,j)/9.81*
     -                        1./tanh(2.*3.14*depth/L_spec(i,j)))
                              if(tmp.gt.per_spec(i,j))then
                                    L_spec(i,j)=L_spec(i,j)-step
                                    if(shift.eq.1)then
                                       shift=0
                                       step=step/10
                                    else
                                       shift=-1
                                    endif
                              else
                                    L_spec(i,j)=L_spec(i,j)+step
                                    if(shift.eq.-1)then
                                       shift=0
                                       step=step/10
                                    else
                                       shift=1
                                    endif
                              endif
                        enddo

                        k_tmp=2*3.1415/L_spec(i,j)
                        w_tmp=2*3.1415/per_spec(i,j)
                        beta_src(i,j)=80./(1.2*L_spec(i,j))**2.
                        I_src(i,j)=sqrt(3.14/beta_src(i,j))*
     -                        exp(-k_tmp**2./4./beta_src(i,j))

					 D_src(i,j)=4.*amp_spec(i,j)*cos(mod_ang)*
     -					(w_tmp**2.+Num2*9.81*k_tmp**2.*
     -					k_tmp**2.*depth**2.*depth+
     -					Num3*9.81*k_tmp**2.*k_tmp**2.*
     -					k_tmp**2.*depth**2.*depth**2.*depth)/
     -					(w_tmp*k_tmp*I_src(i,j)*(1-
     -					Den1*k_tmp**2.*depth**2.-
     -					2./sqrt(-bet(2))*Den2*k_tmp**2.*
     -					 k_tmp**2.*depth**2.*depth**2.))


                       call RANDOM_NUMBER(shift_rn)
c			shift_rn=abs(real(i*j)*.384165-nint(real(i*j)*.384165))
                       shift_spec(i,j)=per_spec(i,j)*shift_rn


c						if(i.eq.1)then
c							shift_spec(i,j)=-0.25*3.1415
c						elseif(i.eq.3)then
c							shift_spec(i,j)=0.25*3.1415
c						else
c							shift_spec(i,j)=0.
c						endif

                  enddo
			enddo
		endif
        endif
      elseif(int_src.eq.2)then
        if(wave_type.eq.1)then
            mod_ang=abs(inc_ang-90*nint(inc_ang/90))*3.1415/180.
            inc_ang=inc_ang*3.1415/180.

		  beta=100.

		  k_1=2*3.1415/L
		  beta_src(1,1)=beta/L**2.

		  shift_spec(1,1)=0
		  D_src(1,1)=7.5*dt/per*
     -				2.*wave_hgt*
     -				0.5*(1.+2.*k_1*depth/sinh(2.*k_1*depth))

	  elseif(wave_type.eq.5.or.wave_type.eq.4)then
            mod_ang=abs(inc_ang-90*nint(inc_ang/90))*3.1415/180.
            inc_ang=inc_ang*3.1415/180.

		  if(wave_type.eq.4)then
			beta=400.   !400 plunging  100 spilling

			L_spec(1,1)=L
			per_spec(1,1)=per
			amp_spec(1,1)=wave_hgt/2.
			inc_ang_spec(1,1)=inc_ang
			beta_src(1,1)=beta/L**2.

			shift_spec(1,1)=0
			D_src(1,1)=11.5*dt/per_spec(1,1)*   !11.5 plunging 5 spilling
     -				2.*amp_spec(1,1)*
     -				0.5*(1.+2.*k_1*depth/sinh(2.*k_1*depth))
		  elseif(wave_type.eq.5)then

			beta=20.    !20

			L_spec(1,1)=L
			per_spec(1,1)=per
			amp_spec(1,1)=wave_hgt/2.
			inc_ang_spec(1,1)=inc_ang
			beta_src(1,1)=beta/L**2.

			shift_spec(1,1)=0

			D_src(1,1)=4.1357*dt/per_spec(1,1)*   !4.1357
     -				2.*amp_spec(1,1)*
     -				0.5*(1.+2.*k_1*depth/sinh(2.*k_1*depth))

			L_spec(2,1)=L_2
			per_spec(2,1)=per_2
			amp_spec(2,1)=wave_hgt_2/2.
			inc_ang_spec(2,1)=inc_ang
			beta_src(2,1)=beta/L_2**2.

			shift_spec(2,1)=0
			D_src(2,1)=4.1357*dt/per_spec(2,1)*
     -				2.*amp_spec(2,1)*
     -				0.5*(1.+2.*k_2*depth/sinh(2.*k_2*depth))
		   endif

        elseif(wave_type.eq.6)then

		if(spec_type.eq.2)then
			beta=20.

			L_spec(1,1)=L
			per_spec(1,1)=per
			amp_spec(1,1)=wave_hgt/2.
			inc_ang_spec(1,1)=inc_ang
			beta_src(1,1)=beta/L**2.

			shift_spec(1,1)=0
			D_src(1,1)=10.*dt/per_spec(1,1)*
     -				2.*amp_spec(1,1)*
     -				0.5*(1.+2.*k_1*depth/sinh(2.*k_1*depth))			
		else
            do i=1,num_freq
                  do j=1,num_theta
                        if(num_theta.eq.1)then
                              mod_ang=abs(inc_ang-90*
     -                        nint(inc_ang/90))*3.1415/180.
                              inc_ang_spec(i,j)=inc_ang*3.1415/180.
                        else
                              mod_ang=abs(inc_ang_spec(i,j)-
     -                        90*nint(inc_ang_spec(i,j)/90))*3.1415/180.

                              inc_ang_spec(i,j)=
     -                             inc_ang_spec(i,j)*3.1415/180.

							if(is_oreint.eq.1)then
								inc_ang_spec(i,j)=
     -								inc_ang_spec(i,j)+3.1415/2.
							endif
                        endif

                        L_spec(i,j)=2.*depth
                        step=L_spec(i,j)/10.
                        shift=0

					  if(j.eq.1)then
						    inc_ang=inc_ang_spec(i,j)
					  endif

                        do while(step.gt.L_spec(i,j)/100000)
                              tmp=sqrt(2*3.14*L_spec(i,j)/9.81*
     -                        1./tanh(2.*3.14*depth/L_spec(i,j)))
                              if(tmp.gt.per_spec(i,j))then
                                    L_spec(i,j)=L_spec(i,j)-step
                                    if(shift.eq.1)then
                                       shift=0
                                       step=step/10
                                    else
                                       shift=-1
                                    endif
                              else
                                    L_spec(i,j)=L_spec(i,j)+step
                                    if(shift.eq.-1)then
                                       shift=0
                                       step=step/10
                                    else
                                       shift=1
                                    endif
                              endif
                        enddo

                        k_tmp=2*3.1415/L_spec(i,j)
                        w_tmp=2*3.1415/per_spec(i,j)
					  beta=20.
				      beta_src(i,j)=beta/L_spec(i,j)**2.
                        D_src(i,j)=4.1357*dt/per_spec(i,j)*
     -					2.*amp_spec(i,j)*
     -					0.5*(1.+2.*k_tmp*depth/sinh(2.*k_tmp*depth))

                       call RANDOM_NUMBER(shift_rn)
c			shift_rn=abs(real(i*j)*.384165-nint(real(i*j)*.384165))
                       shift_spec(i,j)=per_spec(i,j)*shift_rn
                  enddo
            enddo
		endif
        endif
      endif


      if(plus_tide.eq.1)then
c			beta_tide=20./(10.*depth)**2. % for source_type=2
c			D_tide=10.*dt/per	


            beta_tide=80./(20.*depth)**2.

            I_src_cur=sqrt(3.14/beta_tide)

            D_tide=2.*sqrt(9.81*depth)/(I_src_cur)

      endif

	depth=depth_tmp

      return

      end
