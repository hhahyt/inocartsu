      subroutine dhdt_calc(n_loc,nn_loc)
      use mainvar_module
      integer n_loc,nn_loc,xind_top,xind_bottom,p,interp,c_ind,nnn,
     -	shape,x1_ind,y_mid_ind,n_xs,count_sl,mid_j
	real end_slope_2,y_scale(ny),amp_c,damp_dt,damp_dtt,
     -	dy_scaledy(ny),alp_3,sl,xsl,ang,Thick,steep,Cd_new,
     -	Cm_new,u_t,a_o,hor_len,x_o,ds,total,u_c,x_o_C,max_ellip,
     -	norm,ellip1,xo_offset,ellip2,L1,L2,Length,TH_max,
     -	ellip,du_tdt,du_tdtt,dbdt,dbdtt
	
	real xslip(endx),zslip(endx),headx,headz,toex,toez,cx,cz,
     -	dzh,dzt,cx0,err_sl,rad,ztop(endx),x1_sl,z1,pt(endx,2),zos(endx),
     -	dh0(endx),mx_rot,rad_bas,z_bas,dh0_proj(endx),xs_proj(endx),
     -	dw,ds_sl,rot,dxc,dzc,xrot(endx),zrot(endx),zosa(endx),
     -	dha(endx),dh_rot(endx),dh_proj(endx),dh(endx),ur(endt),
     -  time_scale,accel_phase,ac_t,proj,ang_sl(endx),xtemp(endx),
     -  ztemp(endx),int_ur,dh_max,alp_sl,x_m,z_m,slope_sl,zloc,
     -  wid_factor,midy,width,portion_to_rotate,dang_t
	
	
	slide_type=2		 ! 1=solid body motion, following JGR paper
						 ! 2=rotational - translation motion
	if(slide_type.eq.1)then

       ts=5      !start moving time
       ang=10*3.1415/180	 ! slope angle
       b=5.		! slide half-width in sliding plane
       Thick=0.1	!   max vertical slide thickness
       x_o=30.     ! initial location of slide centerpoint
       gamma=1.853	   ! slide material specific weight
       Cd=1.0		   ! slide drag coef
       Cm=1.0			! slide added mass coef
       
       steep=2. !Cd/cos(ang)  ! slide side steepeness factor


       u_t=sqrt(9.81*1.*3.1415*(gamma-1)/(2.*Cd)*sin(ang))
       a_o=9.81*(gamma-1)/(gamma+Cm)*sin(ang)
       t_o=u_t/a_o

       hor_len=b/2.*cos(ang)
  
       ds=0
          
       total=0.
       fs(nn_loc)=0.
       fa(nn_loc)=0.
       do jj=1,endy
          do i=1,endx
            if(h(i,jj,n_loc).gt.0)then
                fs(nn_loc)=fs(nn_loc)-h(i,jj,n_loc)+ho(i,jj)
            else
                fa(nn_loc)=fa(nn_loc)-h(i,jj,n_loc)+ho(i,jj)
            endif
          enddo
       enddo

       total=fs(nn_loc)+fa(nn_loc)
       if(total.gt.1e-20)then
          fs(nn_loc)=fs(nn_loc)/total
          fa(nn_loc)=fa(nn_loc)/total
       endif

       do nnn=1,nn_loc
          if(t(nnn).lt.ts)then
            ds=0
            u_c=0
          else
c          ds=(t(nn_loc)-ts)*u_t*tanh((t(nn_loc)-ts)/t_o)
            u_c=fs(nnn)*u_t*tanh((t(nnn)-ts)/t_o)+
     -                fa(nnn)*9.81*(t(nnn)-ts)*sin(ang)
            ds=ds+u_c*dt
          endif
       enddo



c		determine length
	 x_o_C=x(nint(endx/2.))  
       x_l=x_o_C-hor_len
       x_r=x_o_C+hor_len
	 max_ellip=0.
	 norm=1./((.5+.5*tanh((x_o_C-x_l)/(steep)))*
     -    (.5*(1.-tanh((x_o_C-x_r)/(steep*Cm/Cd)))))

	do i=1,endx-1
		ellip1=norm*(.5+.5*tanh((x(i)-x_l)/(1.*steep)))*
     -            (.5*(1.-tanh((x(i)-x_r)/(steep*Cm/Cd))))
     		if(ellip1.gt.max_ellip)then							
			max_ellip=ellip1
			xo_offset=x(nint(endx/2.))-x(i)
		endif
	enddo

	do i=1,endx-1
		ellip1=norm*(.5+.5*tanh((x(i)-x_l)/(1.*steep)))*
     -            (.5*(1.-tanh((x(i)-x_r)/(steep*Cm/Cd))))	
		ellip2=norm*(.5+.5*tanh((x(i+1)-x_l)/(1.*steep)))*
     -            (.5*(1.-tanh((x(i+1)-x_r)/(steep*Cm/Cd))))	
		if(ellip1.le.0.01*max_ellip.and.ellip2.gt.0.01*max_ellip)then
			L1=x(i)
		endif
		if(ellip1.gt.0.01*max_ellip.and.ellip2.le.0.01*max_ellip)then	
			L2=x(i+1)	
		endif
	enddo
	Length=L2-L1	

       x_o=x_o+ds*cos(ang)+xo_offset
       x_l=x_o-hor_len
       x_r=x_o+hor_len

	 Thick=Thick/((.5+.5*tanh((x_o-x_l)/(steep)))*
     -    (.5*(1.-tanh((x_o-x_r)/(steep*Cm/Cd)))))/max_ellip

       TH_max=0
       do jj=1,endy
          do i=1,endx
            ellip=Thick*(.5+.5*tanh((x(i)-x_l)/(1.*steep)))*
     -            (.5*(1.-tanh((x(i)-x_r)/(steep*Cm/Cd))))

            if(ellip.ge.TH_max) TH_max=ellip

            if(dim.eq.1)then
                y_scale(jj)=1.
            elseif(dim.eq.2)then
                yo=y(endy)/2.
                ao=Length/aspect_ratio/4.3
                y_scale(jj)=exp(-(y(jj)-yo)**2./ao**2.)
            endif

            h(i,jj,n_loc)=ho(i,jj)-ellip*y_scale(jj)

          enddo
       enddo
       
	elseif(slide_type.eq.2)then
       headx=1.8e4		! x coordinate of head
       toex=3.05e4		! x coordinate of toe
	   midy=5.75e4		! y coordinate of slide centerline
	   width=5.e4		! width of slide, in the y-direction
	   Thick=2000.		! max vertical slide thickness
	   ts=10.			! start motion time
	   portion_to_rotate=1.0  !fraction of rotational mass to rotate
							  ! = 1., head moves to toe	
							  ! = 0.5, head moves to halfpoint between toe and head
							  ! = 0., slide does not move, head stays at head
       time_scale=1000.	 ! time scale of sliding motion	  (secs)
       accel_phase=0.1; ! =0.5 -> equal acceleration and deacceleration phases
					   ! =0.4 -> ~40% of time_scale in acceleration and 60% in deacceleration	   
					   ! =0.3 -> ~30% of time_scale in acceleration and 70% in deacceleration	   
					   ! =0.0 -> ~20% acceleration and 80% deacceleration	
       
       
       ! determine velocity time history
       ac_t=time_scale*accel_phase

	   do i=1,endt
		if(t(i).le.ts)then
			ur(i)=0.
		else
			ur(i)=log(cosh(2.*3.1415/time_scale*(t(i)-ts)))*
     -            (0.5-0.5*tanh(2.*3.1415/time_scale*(t(i)-ac_t-ts)))
		endif   
	   enddo

	   int_ur=0;
	   do i=1,endt
		int_ur=int_ur+ur(i)
	   enddo

	   ur=ur/int_ur
       
       mid_j=overlap+1
       do j=1,endy-1
			if(y(j).le.midy.and.y(j+1).gt.midy)then
				mid_j=j
			endif
		enddo
		
		j=mid_j
		
       call interp1D_single(x,-ho(:,j),headx,headz,endx,1)		 
       call interp1D_single(x,-ho(:,j),toex,toez,endx,1)   	 
       
		rad=sqrt((headx-toex)**2. + (headz-toez)**2.)
		alp_sl=90.*3.1415/180.-60.*3.1415/180.-asin((headz-toez)/rad)
		cx=toex-rad*sin(alp_sl)
		cz=toez+rad*cos(alp_sl)

		x_m=(headx+toex)/2.
		z_m=(toez+headz)/2.

		slope_sl=(cz-z_m)/(cx-x_m)
	    n_xs=floor( (toex-headx)/dx )+1
	   
	    do i=1,n_xs
			xslip(i)=headx +real(i-1)*dx
	    enddo

		count_sl=0;
		err=1.
		do while(err.gt.1e-3)
			count_sl=count_sl+1
			dh_max=0.
			do i=1,n_xs
				zslip(i)=cz-sqrt( rad**2-(xslip(i)-cx)**2)
				call interp1D_single(x,-ho(:,j),xslip(i),zloc,endx,1)
				dh_max=max(dh_max,zloc-zslip(i))
			enddo
    
			err= abs(dh_max-Thick)
    
			if(dh_max.gt.Thick)then
				cz=cz+abs(slope_sl*err)
				cx=cx+abs(err)
			else
				cz=cz-abs(slope_sl*err)
				cx=cx-abs(err)
			endif
			rad=sqrt((headx-cx)**2 + (headz-cz)**2)
    
			if(count_sl.gt.1000)then
				print*, 'ERROR: Could create slide with given Thick'
				print*, 'Min possible Thick:', dh_max
				err=0.
			endif
		enddo	
	   	 
       call interp1D(x,-ho(:,j),xslip,ztop,endx,n_xs)	
	   do i=1,n_xs
		x1_sl=xslip(i)
		z1=ztop(i)
    
		pt(i,1)=sqrt( (x1_sl-cx)**2+(z1-cz)**2)
		pt(i,2)=atan( (z1-cz)/(x1_sl-cx) )   
		if(pt(i,2).lt.0.)then
			pt(i,2)=3.1415+pt(i,2)
		endif
	   enddo	
	  
	   do i=1,endx
		if(x(i).lt.xslip(1).or.x(i).gt.xslip(n_xs))then
			zos(i)=-ho(i,j)
		else
			call interp1D_single(xslip,zslip,x(i),zos(i),n_xs,1)
		endif
	   enddo
	   
	   dh0=-ho(:,j)-zos
	   
	   mx_rot=pt(n_xs,2);
	   do i=1,n_xs
		rad_bas=mx_rot-pt(i,2)
		z_bas=cz-sqrt(pt(i,1)**2-(toex-cx)**2)
		dh0_proj(i)=toez-z_bas
		xs_proj(i)=-rad_bas*rad+toex
	   enddo
       
       dang_t=abs(pt(1,2)-pt(n_xs,2))*portion_to_rotate
       rot=0.
       proj=0.
       
       do i=1,nn_loc
		   dw=ur(i)*dang_t		
		   rot=rot+dw
		   proj=proj+dw*rad
	   enddo
	   
	   ang_sl=pt(:,2)+rot
	   
	   do i=1,n_xs
        dxc=cos(ang_sl(i))*pt(i,1)
        dzc=sin(ang_sl(i))*pt(i,1)
        
        xrot(i)=cx-dxc
        zrot(i)=cz-dzc
	   enddo
	   
	   do j=1,endy
		do i=1,n_xs+3
			if(i.eq.1)then
				xtemp(i)=x(1)
				ztemp(i)=-ho(1,j)
			elseif(i.eq.2)then
				xtemp(i)=headx-.001	
				ztemp(i)=headz
			elseif(i.eq.n_xs+3)then
				xtemp(i)=x(endx)
				ztemp(i)=-ho(endx,j)
			else
				xtemp(i)=xrot(i-2)
				ztemp(i)=zrot(i-2)
			endif
		enddo		
	   
		call interp1D(xtemp,ztemp,x,zosa,n_xs+3,endx)
		dh_rot=-ho(:,j)-zosa
		   
		call interp1D(xs_proj+proj,dh0_proj,x,dh_proj,n_xs,endx)
		do i=1,endx
			if(x(i).lt.xrot(1))then
				dh(i)=-ho(i,j)-zos(i)
			elseif(x(i).lt.toex)then
				dh(i)=dh_rot(i)
			elseif(x(i).lt.xs_proj(n_xs)+proj)then
				dh(i)=min(0.,dh_proj(i))
			else
				dh(i)=0.
			endif
		enddo
       
		do i=1,endx
			if(dim.eq.1)then
				wid_factor=1.0
			else	
				wid_factor=exp(-4.*3.14*((y(j)-midy)/width)**4.);
			endif
			h_tmp(i,j)=dh(i)*wid_factor
		enddo

		do i=3,endx-2
            dh(i)=1./10.*(2.*h_tmp(i-1,j)+
     -                2.*h_tmp(i+1,j)+4.*h_tmp(i,j)+
     -                h_tmp(i-2,j)+
     -                h_tmp(i+2,j))
     
			h(i,j,n_loc)=-(-ho(i,j)-dh(i))
		enddo
       enddo
       
	endif

      if(dim.eq.1.and.n.ge.3)then
       j=overlap+1
       do i=3,endx-2
          dhdt(i,j,n_loc)=(h(i,j,n_loc)-h(i,j,n_loc-1))/(dt)
          dhdtt(i,j,n_loc)=(h(i,j,n_loc)-2.*h(i,j,n_loc-1)+
     -                h(i,j,n_loc-2))/(dt**2.)
       enddo

       do i=2,endx-1
          do cur_level=1,num_levels
            if(cur_level.gt.1)then
                zeta(i,j,n_loc,cur_level)=
     -                bet2(cur_level)*h(i,j,n_loc)                
            endif
             if(cur_level.eq.num_levels)then
                zeta(i,j,n_loc,num_levels+1)=
     -                -h(i,j,n_loc)      
            endif      

            z_alp(i,j,n_loc,cur_level)=
     -                bet(cur_level)*h(i,j,n_loc)

            dz_alpdt(i,j,n_loc,cur_level)=
     -                bet(cur_level)*dhdt(i,j,n_loc)
          enddo
       enddo

C............Calculate Derivates.....................
       do i=5,endx-4
          dhdxt(i,j,n_loc)=(-dhdt(i+2,j,n_loc)+8.*dhdt(i+1,j,n_loc)-
     -            8.*dhdt(i-1,j,n_loc)+dhdt(i-2,j,n_loc))/(12.*dx)

          dhdxtt(i,j,n_loc)=(-dhdtt(i+2,j,n_loc)+8.*dhdtt(i+1,j,n_loc)-
     -            8.*dhdtt(i-1,j,n_loc)+dhdtt(i-2,j,n_loc))/(12.*dx)
       enddo
      elseif(dim.eq.2.and.n.ge.3)then
       do j=1,endy
          do i=1,endx
            dhdt(i,j,n_loc)=(h(i,j,n_loc)-h(i,j,n_loc-1))/(dt)
            dhdtt(i,j,n_loc)=(h(i,j,n_loc)-2.*h(i,j,n_loc-1)+
     -                h(i,j,n_loc-2))/(dt**2.)
          enddo
       enddo

       do j=1,endy
          do i=1,endx
            do cur_level=1,num_levels
                if(cur_level.gt.1)then
                  zeta(i,j,n_loc,cur_level)=
     -                  bet2(cur_level)*h(i,j,n_loc)                
                endif
                if(cur_level.eq.num_levels)then
                  zeta(i,j,n_loc,num_levels+1)=
     -                  -h(i,j,n_loc)      
                endif      

                z_alp(i,j,n_loc,cur_level)=
     -                bet(cur_level)*h(i,j,n_loc)

                dz_alpdt(i,j,n_loc,cur_level)=
     -                bet(cur_level)*dhdt(i,j,n_loc)
            enddo
          enddo
       enddo

C............Calculate Derivates.....................
       do j=3,endy-2
          do i=3,endx-2
            dhdxt(i,j,n_loc)=(-dhdt(i+2,j,n_loc)+8.*dhdt(i+1,j,n_loc)-
     -            8.*dhdt(i-1,j,n_loc)+dhdt(i-2,j,n_loc))/(12.*dx)
           dhdxtt(i,j,n_loc)=(-dhdtt(i+2,j,n_loc)+8.*dhdtt(i+1,j,n_loc)-
     -            8.*dhdtt(i-1,j,n_loc)+dhdtt(i-2,j,n_loc))/(12.*dx)

            dhdyt(i,j,n_loc)=(-dhdt(i,j+2,n_loc)+8.*dhdt(i,j+1,n_loc)-
     -            8.*dhdt(i,j-1,n_loc)+dhdt(i,j-2,n_loc))/(12.*dy)
           dhdytt(i,j,n_loc)=(-dhdtt(i,j+2,n_loc)+8.*dhdtt(i,j+1,n_loc)-
     -            8.*dhdtt(i,j-1,n_loc)+dhdtt(i,j-2,n_loc))/(12.*dy)
          enddo      
       enddo      
      endif


      return
      end



 
