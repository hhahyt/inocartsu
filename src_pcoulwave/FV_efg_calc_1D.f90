!*******************************************************************************
!........    Groups the physical variables into the 
!........... conservative forms given in Wei & Kirby
subroutine FV_efg_calc_1D(z,d,u,E,F,F1,n_loc,nn_loc)
use mainvar_module, only:dx,dy,nx,ny,grav,dvdyy,dudxx,dvdxy,dudxy,dhzudx,dhzudy,dhzvdx,dhzvdy,dim,bl_hor_wall,level, &
 Ch,bf_type,disp_prop,cutoff,dt,S_grp,T_grp,dhdt,dSdx,dTdx,dhdxt,endx,endy,cutoff_mat,nut,B_mat,bf_type, & 
 f_int_src, theta, freq, num_theta, num_freq, depth, x_c,y_c,end_x_t, end_y_t, x0, cosA,sinA,wA,kA,inc_ang_spec,L_spec,per_spec, & 
 shift_spec,is_coef,is_sin,loc,x,y,D_src,beta_src,t,count_min,int_src,is_oreint,tmp,CB,backscatter,ihvor,psix, &
 elder_length,Ch_length,D_tide,beta_tide,eta_in,plus_tide,x_wvmk,u_wvmk,imag_leng,wm_d_c,wm_u_c,hx,numerical_scheme,cutoff_mat_d
use FV_var_module, only:sx,ex,sy,ey,idissip,ibotfric,b0,b1,cm,ks,alpha,ilim,ix,iy,ix2,iy2,id,limiter2,limiter3, &
 iriemann,icor,interp2,zx,zy,rand,BS,ish_mat
      
integer sum,ii,i,j,sum_0,ilim_c,n_loc,nn_loc,wvmk_right,sum_id
real Hs_i,us_i,BS_breaker
     
real :: u(nx,ny),z(nx,ny), d(nx,ny),za(nx,ny)
real :: Hu(nx,ny),H(nx,ny)
real :: zi(nx,ny), hi(nx,ny),ui(nx,ny)
real :: hui(nx,ny),H_i(nx,ny)
real :: phil(nx,ny),phir(nx,ny)
real :: ux(nx,ny),hux(nx,ny),Hcx(nx,ny)
real :: uixx(nx,ny),uix(nx,ny)
real :: huix(nx,ny),huixx(nx,ny)
real :: zai(nx,ny),Ei(nx,ny)
real :: uxx(nx,ny),dummy(nx,ny)
real :: huxx(nx,ny)
     
real :: dissipcoef
real :: vt(nx,ny),vti(nx,ny)
real :: tmpi1(nx,ny)
real :: tmpi2(nx,ny)
real :: tmp1,tmp2,tmp3,tmp4
real :: E(nx,ny),F(nx,ny),F1(nx,ny)
real :: temp1,temp2,temp3,temp4,temp5,temp6,temp7
real :: E1(nx,ny)
real :: H_l(nx,ny),H_r(nx,ny)
real :: Hu_l(nx,ny),Hu_r(nx,ny)
real :: H_b(nx,ny),H_t(nx,ny)
real :: Hu_b(nx,ny)
real :: Hu_t(nx,ny),zx_4th(nx,ny)
real :: conv_u,loc_freq(nx,ny)
real :: flux(nx,ny,3),fluy(nx,ny,3)
real :: dh0,dh1,dh2,dh3,dh4,src_amp
real :: dbh0,dbh1,dbh2,dbh3,dbh4,d3h,dsh1,dsh2,dsh3
real :: dshb1,dshb2,dshb3,dshb4,mx_dzdx
real :: EzSTi(nx,ny),TzS2i(nx,ny)
real :: uSxvSyi(nx,ny),uTxvTyi(nx,ny)

real :: ua(nx,ny),Rhx(nx,ny),vtv(nx,ny),flux_loc,ks_i,cm_i
real :: Rh, Cf(nx,ny),Cf_c, mx_theta, L_tanh,k_tmp,w_tmp
real :: dx1,dx2,z_wvmk,flux1_left,flux2_left,f_bc_add,C_loc,Cd_loc,Rwx
real minH, Hmin, phir1, phir2,phir1i,phir2i,omeg,dupw,dloc,delta,up_ratio
parameter(dissipcoef=1e-6)

dummy=1.
minH = cutoff*2.
Hmin = cutoff*2.

! conserved variables
H = d+z
Hu = H*u

! correction boolean matrix for near shoreline, SW simulations
	do j=sy,ey
		do i=sx,ex
			if(id(i,j).ne.0.or.disp_prop.eq.3.or.H(i,j).le.cutoff_mat_d(i,j))then
				ish_mat(i,j)=1
			else
				ish_mat(i,j)=0
			endif
		enddo
	enddo

!.....calculate derivatives and terms required for calculating E,F, and G



!if(ilim.eq.0)then

do j = sy,ey
	do i = sx-1,ex

!.....k = 1

		 dh0 = z(i-1,j)-z(i-2,j)
		 dh1 = z(i  ,j)-z(i-1,j)
		 dh2 = z(i+1,j)-z(i  ,j)
		 dh3 = z(i+2,j)-z(i+1,j)
		 dh4 = z(i+3,j)-z(i+2,j)

         ilim_c=ilim
		 mx_dzdx=max( abs(dh0),abs(dh1), abs(dh2), abs(dh3), abs(dh4) )/dx
		 if(mx_dzdx.gt.0.60) ilim_c=1

		 dbh0 = sign(1.0,dh0)*max(0.0,min(abs(dh0),b1*sign(1.0,dh0)*dh1,b1*sign(1.0,dh0)*dh2))
		 dbh1 = sign(1.0,dh1)*max(0.0,min(abs(dh1),b1*sign(1.0,dh1)*dh2,b1*sign(1.0,dh1)*dh0))
		 dbh2 = sign(1.0,dh2)*max(0.0,min(abs(dh2),b1*sign(1.0,dh2)*dh0,b1*sign(1.0,dh2)*dh1))

		 d3h =  dbh2-2*dbh1+dbh0
		 dsh1 = dh1-d3h/6

		 dbh1 = sign(1.0,dh1)*max(0.0,min(abs(dh1),b1*sign(1.0,dh1)*dh2,b1*sign(1.0,dh1)*dh3))
		 dbh2 = sign(1.0,dh2)*max(0.0,min(abs(dh2),b1*sign(1.0,dh2)*dh3,b1*sign(1.0,dh2)*dh1))
		 dbh3 = sign(1.0,dh3)*max(0.0,min(abs(dh3),b1*sign(1.0,dh3)*dh1,b1*sign(1.0,dh3)*dh2))

		 d3h =  dbh3-2*dbh2+dbh1
		 dsh2 = dh2-d3h/6

		 dbh2 = sign(1.0,dh2)*max(0.0,min(abs(dh2),b1*sign(1.0,dh2)*dh3,b1*sign(1.0,dh2)*dh4))
		 dbh3 = sign(1.0,dh3)*max(0.0,min(abs(dh3),b1*sign(1.0,dh3)*dh4,b1*sign(1.0,dh3)*dh2))
		 dbh4 = sign(1.0,dh4)*max(0.0,min(abs(dh4),b1*sign(1.0,dh4)*dh2,b1*sign(1.0,dh4)*dh3))

		 d3h =  dbh4-2*dbh3+dbh2
		 dsh3 = dh3-d3h/6

		 if (ilim_c.eq.0) then

		    phir1 = 1.
		    phir2 = 1.
		    phir1i = 1.
		    phir2i = 1.


	        H_l(i,j) = z(i,j)   + (phir1*dsh1+phir1i*2*dsh2)/6.0 + 0.5*(d(i,j)+d(i+1,j))
            H_r(i,j) = z(i+1,j) - (2*phir2*dsh2+phir2i*dsh3)/6.0 + 0.5*(d(i,j)+d(i+1,j))


         elseif (ilim_c.eq.1) then

	        dshb1 = sign(1.0,dsh1)*max(0.0,min(abs(dsh1),b0*dsh2*sign(1.0,dsh1)))
	        dshb2 = sign(1.0,dsh2)*max(0.0,min(abs(dsh2),b0*dsh1*sign(1.0,dsh2)))
	        dshb3 = sign(1.0,dsh2)*max(0.0,min(abs(dsh2),b0*dsh3*sign(1.0,dsh2)))	     
		    dshb4 = sign(1.0,dsh3)*max(0.0,min(abs(dsh3),b0*dsh2*sign(1.0,dsh3)))


	        H_l(i,j) = z(i,j)   + (dshb1+2*dshb2)/6.0 + 0.5*(d(i,j)+d(i+1,j))
            H_r(i,j) = z(i+1,j) - (2*dshb3+dshb4)/6.0 + 0.5*(d(i,j)+d(i+1,j))

		 endif

!.....k = 2
      

		 dh0 = Hu(i-1,j)-Hu(i-2,j)
		 dh1 = Hu(i  ,j)-Hu(i-1,j)
		 dh2 = Hu(i+1,j)-Hu(i  ,j)
		 dh3 = Hu(i+2,j)-Hu(i+1,j)
		 dh4 = Hu(i+3,j)-Hu(i+2,j)

		 dbh0 = sign(1.0,dh0)*max(0.0,min(abs(dh0),b1*sign(1.0,dh0)*dh1,b1*sign(1.0,dh0)*dh2))
		 dbh1 = sign(1.0,dh1)*max(0.0,min(abs(dh1),b1*sign(1.0,dh1)*dh2,b1*sign(1.0,dh1)*dh0))
		 dbh2 = sign(1.0,dh2)*max(0.0,min(abs(dh2),b1*sign(1.0,dh2)*dh0,b1*sign(1.0,dh2)*dh1))

		 d3h =  dbh2-2*dbh1+dbh0
		 dsh1 = dh1-d3h/6

		 dbh1 = sign(1.0,dh1)*max(0.0,min(abs(dh1),b1*sign(1.0,dh1)*dh2,b1*sign(1.0,dh1)*dh3))
		 dbh2 = sign(1.0,dh2)*max(0.0,min(abs(dh2),b1*sign(1.0,dh2)*dh3,b1*sign(1.0,dh2)*dh1))
		 dbh3 = sign(1.0,dh3)*max(0.0,min(abs(dh3),b1*sign(1.0,dh3)*dh1,b1*sign(1.0,dh3)*dh2))

		 d3h =  dbh3-2*dbh2+dbh1
		 dsh2 = dh2-d3h/6

		 dbh2 = sign(1.0,dh2)*max(0.0,min(abs(dh2),b1*sign(1.0,dh2)*dh3,b1*sign(1.0,dh2)*dh4))
		 dbh3 = sign(1.0,dh3)*max(0.0,min(abs(dh3),b1*sign(1.0,dh3)*dh4,b1*sign(1.0,dh3)*dh2))
		 dbh4 = sign(1.0,dh4)*max(0.0,min(abs(dh4),b1*sign(1.0,dh4)*dh2,b1*sign(1.0,dh4)*dh3))

		 d3h =  dbh4-2*dbh3+dbh2
		 dsh3 = dh3-d3h/6

!.... Extended variables

		 if (ilim_c.eq.0) then

		    phir1 = 1.
		    phir2 = 1.
		    phir1i = 1.
		    phir2i = 1.

	        Hu_l(i,j) = Hu(i,j)   + (phir1*dsh1+phir1i*2*dsh2)/6.0
	        Hu_r(i,j) = Hu(i+1,j) - (2*phir2*dsh2+phir2i*dsh3)/6.0

		 elseif (ilim_c.eq.1) then

	        dshb1 = sign(1.0,dsh1)*max(0.0,min(abs(dsh1),b0*dsh2*sign(1.0,dsh1)))
	        dshb2 = sign(1.0,dsh2)*max(0.0,min(abs(dsh2),b0*dsh1*sign(1.0,dsh2)))
	        dshb3 = sign(1.0,dsh2)*max(0.0,min(abs(dsh2),b0*dsh3*sign(1.0,dsh2)))	     
		    dshb4 = sign(1.0,dsh3)*max(0.0,min(abs(dsh3),b0*dsh2*sign(1.0,dsh3)))


	        Hu_l(i,j) = Hu(i,j)   + (dshb1+2*dshb2)/6.0
            Hu_r(i,j) = Hu(i+1,j) - (2*dshb3+dshb4)/6.0

		 endif

    enddo
enddo

   do j = sy,ey
	do i = sx,ex
	     Hcx(i,j) = 0.5*(H_r(i-1,j)+H_l(i,j))
      enddo
   enddo

!..... call approximate Riemann solvers
call Riemanns4 (H_l,H_r,Hu_l,Hu_r,dummy,dummy,H_b,H_t,Hu_b,Hu_t,dummy,dummy,flux,fluy,1,sx,ex,sy,ey,ix,iy,iriemann)

if(numerical_scheme.eq.0)then
! HO FD for leading order terms
!
    call interp_xei (H,H_i,1,ix)  ! fourth order interpolation
    call interp_xei (u,ui,1,ix)  ! fourth order interpolation
	
	do j=sy-1,ey
	   do i=sx-1,ex
	     if(id(i,j).ne.0)then
			flux(i,j,1)=H_i(i,j)*ui(i,j)
			flux(i,j,2)=H_i(i,j)*ui(i,j)*ui(i,j) ! + 0.5*grav*Hs_i**2.
		 endif
		enddo
	enddo
endif

    call interp_xei (z,zi,1,ix)  ! fourth order interpolation
	do j=sy,ey
	   do i=sx,ex
			zx_4th(i,j)  = (zi(i,j)-zi(i-1,j))/dx	
		enddo
	enddo

! correct near shoreline flux
	
	do j=sy,ey
	   do i=sx,ex
		if(id(i,j).eq.9)then
			flux(i,j,1)=0.
			flux(i,j,2)=0.
			flux(i,j,3)=0.

			zx_4th(i,j)=0.
		endif
	   enddo
	enddo

	do j=sy,ey
	   do i=sx,ex

		if(id(i,j).eq.2.or.id(i,j).eq.23.or.id(i,j).eq.24)then
			
			call shoreflux(i,j,2,zx_4th,flux,z,u,v,d,dx)  ! calcs flux at i,j for id=2

			call shoreflux(i-1,j,3,zx_4th,flux,z,u,v,d,dx)  ! calcs upwinded flux at i,j

		elseif(id(i,j).eq.1.or.id(i,j).eq.13.or.id(i,j).eq.14)then
			
			call shoreflux(i,j,1,zx_4th,flux,z,u,v,d,dx)  ! calcs flux at i-1,j for id=1

			call shoreflux(i,j,3,zx_4th,flux,z,u,v,d,dx)  ! calcs upwinded flux at i,j

		elseif(d(i,j).le.cutoff_mat(i,j).and.id(i,j).eq.0)then
			zx_4th(i,j)= (z(i+1,j)-z(i-1,j))/dx/2.			
		endif				

	  enddo
	enddo

	do j=sy,ey
	   do i=sx,ex
	    if(level(i,j).eq.0)then
			if(level(i-1,j).eq.1.and.level(i+1,j).eq.1)then  ! one wet alone, freeze it

				zx_4th(i-1,j)= 0.	
				zx_4th(i,j)= 0.
				zx_4th(i+1,j)= 0.	

				flux(i-1,j,1)=0.
				flux(i-1,j,2)=0.
				flux(i-1,j,3)=0.

				flux(i,j,1)=0.
				flux(i,j,2)=0.
				flux(i,j,3)=0.


			elseif(level(i-1,j).eq.1.and.level(i+2,j).eq.1)then  ! two wets alone, correct derivs and middle flux
				
				zx_4th(i-1,j)= 0.	
				if(z(i,j).gt.z(i-1,j))then
					zx_4th(i,j)= (z(i+1,j)-z(i-1,j))/dx/2.
				else
					zx_4th(i,j)= (z(i+1,j)-z(i,j))/dx
				endif

				if(z(i+1,j).gt.z(i+2,j))then
					zx_4th(i+1,j)= (z(i+2,j)-z(i,j))/dx/2.
				else
					zx_4th(i+1,j)= (z(i+1,j)-z(i,j))/dx
				endif
				zx_4th(i+2,j)= 0.
				
				call shoreflux(i,j,3,zx_4th,flux,z,u,v,d,dx)  ! calcs upwinded flux at i,j
								
			endif		
		endif
	  enddo
	enddo

! special boundary corrections
	do j=sy,ey
	   do i=sx,ex
	    if(level(i,j).eq.1)then
			if(level(i-1,j).eq.0.and.level(i+1,j).eq.0)then  ! one dry alone	
				zx_4th(i,j)= (z(i+1,j)-z(i-1,j))/dx/2.
			endif
		endif
	  enddo
	enddo
call interp_xei (z,zi,interp2,ix)
call interp_xei (d,hi,interp2,ix)
call interp_xoi (u,ui,interp2,ix) 
!za  = alpha*d
!zai = alpha*hi

za = -d+(1+alpha)*H
zai = -hi+(1+alpha)*(zi+hi)

hu = d*u

hui = hi*ui

!....

if (icor .eq.1) then
do j=sy,ey
   do i=sx,ex
      ux(i,j)  = (ui(i,j)-ui(i-1,j))/dx
   enddo
enddo

call interpolations(1,limiter2,id,z,phil,phir)
zx = (phir-phil)/dx

call interpolations(3,limiter2,id,ux,phil,phir)
uxx = (phir-phil)/dx

call interpolations(2,limiter2,id,hu,phil,phir)
hux = (phir-phil)/dx
call interpolations(3,limiter2,id,hux,phil,phir)
huxx = (phir-phil)/dx

else

do j=sy,ey
   do i=sx,ex
      zx(i,j)  = (zi(i,j)-zi(i-1,j))/dx
      ux(i,j)  = (ui(i,j)-ui(i-1,j))/dx
   enddo
enddo

call interp_xei (ux,tmpi1,interp2,ix)

do j=sy,ey
   do i=sx,ex
      uxx(i,j)  = (tmpi1(i,j)-tmpi1(i-1,j))/dx
      hux(i,j)  = (hui(i,j)-hui(i-1,j))/dx
   enddo
enddo

call interp_xei (hux,tmpi1,interp2,ix)

do j=sy,ey
   do i=sx,ex
      huxx(i,j)  = (tmpi1(i,j)-tmpi1(i-1,j))/dx
   enddo
enddo


endif


! add for external source layer subroutine
	dudxx(:,:,1)=uxx


! add for velocity profile calculations
do j=1,endy
	do i=1,endx
		S_grp(i,j,1)=ux(i,j)
		T_grp(i,j,1)=hux(i,j) !+dhdt(i,j,n_loc)
		dSdx(i,j,1)=uxx(i,j)
        dTdx(i,j,1)=huxx(i,j) !+dhdxt(i,j,n_loc)
	enddo
enddo




!.....nonlinear terms at x-dir interface
do j=sy,ey
   do i=sx-1,ex
!	     uix(i,j) = (15*(u(i+1,j)-u(i,j))-(u(i+2,j)-u(i-1,j)))/12./dx
!	     huix(i,j) = (15*(hu(i+1,j)-hu(i,j))-(hu(i+2,j)-hu(i-1,j)))/12./dx
!		 uixx(i,j) = (u(i+2,j)+u(i-1,j)-u(i+1,j)-u(i,j))/dx/dx*0.5
!		 huixx(i,j) = (hu(i+2,j)+hu(i-1,j)-hu(i+1,j)-hu(i,j))/dx/dx*0.5	 
	     uix(i,j) = (u(i+1,j)-u(i,j))/dx
	     huix(i,j) = (hu(i+1,j)-hu(i,j))/dx
		 uixx(i,j) = (u(i+2,j)+u(i-1,j)-u(i+1,j)-u(i,j))/dx/dx*0.5
		 huixx(i,j) = (hu(i+2,j)+hu(i-1,j)-hu(i+1,j)-hu(i,j))/dx/dx*0.5
   enddo
enddo


!...............................................................................

if(ibotfric.eq.1 .or. idissip.eq.1 )then
   do j=sy,ey
      do i=sx,ex
		 Hs_i=max(cutoff_mat(i,j),H(i,j))
         
		 if(ish_mat(i,j).eq.0)then
			ua(i,j) = u(i,j)  - ((z(i,j)**2-z(i,j)*d(i,j)+d(i,j)**2)/6.0-0.5*za(i,j)**2)*(uxx(i,j)) &
                 - (0.5*(z(i,j)-d(i,j))-za(i,j))*(huxx(i,j))  !&
                 ! + psix(i,j,n_loc)*(za(i,j)**2*0.5-za(i,j)*z(i,j)+(2*z(i,j)**2-2*z(i,j)*d(i,j)-d(i,j)**2)/6.0)      ! depth-averaged velocity 
         else
			ua(i,j) = u(i,j) 
		 endif

		 Cf(i,j)=0.
		 if(bf_type.eq.0)then
			Rh = max(100.0,sqrt(ua(i,j)**2)*Hs_i/dissipcoef)
			Cf_c = -1.8*log10(6.9/Rh+(ks/Hs_i/3.7)**1.11)  ! Moody approx
			Cf(i,j) = 0.25/Cf_c**2
			if (Rh.le.411.58) Cf(i,j) = 0.25*64./Rh
		 elseif(bf_type.eq.1)then		 
		    Cf(i,j)=ks**2.*9.81/Hs_i**0.3333  ! Mannings
		 elseif(bf_type.eq.2)then		 
			Cf(i,j)=ks  ! constant friction
         endif

		 if(id(i,j).ne.0.or.Hs_i.lt.cutoff_mat_d(i,j))then
			if(ua(i,j)*hx(i,j).ge.0)then
				Cf(i,j) = max(Cf(i,j),0.04)  ! Puleo & Holland downrush friction factor
			else
				Cf(i,j) = max(Cf(i,j),0.01)  ! Puleo & Holland uprush friction factor
			endif
		 endif

         Rhx(i,j) = Cf(i,j)*ua(i,j)*sqrt(ua(i,j)**2)
         
      enddo
   enddo

endif   


!--------- add dissipation --------------
if (idissip.eq.1) then
   do j = sy,ey
	do i = sx,ex
		    Hs_i=max(cutoff_mat(i,j),H(i,j))
            vtv(i,j) = Ch*Hs_i*sqrt(Cf(i,j)*(ua(i,j)**2)) + dissipcoef
			   
			cm_i=cm
			if(id(i,j).ne.0.or.Hs_i.le.cutoff_mat_d(i,j)) cm_i=1.0

            vt(i,j) =  cm*dx*dy*sqrt(4.*ux(i,j)**2) + dissipcoef 
      enddo
   enddo
   call interp_xei (vt,vti,interp2,ix)
endif

if (ihvor.eq.1) then
   do j = sy,ey
      do i = sx,ex
		    if(H(i,j).lt.Ch_length)then
				psix(i,j,n_loc)=0.
			else
				psix(i,j,n_loc) = Rhx(i,j)/(vtv(i,j)*H(i,j))
			endif
      enddo
   enddo
else
	do j=2,endy-2
	   do i=2,endx-2
			psix(i,j,n_loc) = Rhx(i,j)/(vtv(i,j)*H(i,j))
      enddo
   enddo      
endif

! Internal source (type 1) wave generation
	
	do j = sy,ey
	do i = sx,ex
			  E(i,j)=0.
			  F(i,j)=0.
		enddo
	enddo


	if(int_src.eq.1.or.plus_tide.eq.1)then  !^^ Internal source type 1
	 if(n_loc.eq.4.and.count_min.eq.1)then !&&
	   do j = sy,ey
		do i = sx,ex

		  f_int_src(i,j)=0.
		  
          if(bl_hor_wall(i,j).le.50)then !%%
          
!%%%%%%%%%%%%%%%%%  INTERNAL SOURCE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          tmp=0.

		  tmp=min(1.,max(0.,d(i,j)/depth))**2.

		  if(is_oreint.eq.1)then
				x_c=x(i)-end_x_t/2.
				y_c=y(j)-x0
		  else
				x_c=x(i)-x0 
				y_c=y(j)-end_y_t/2.
		  endif

          do theta=1,num_theta
            loc_freq(i,j)=0.
			do freq=1,num_freq
			
				  loc=abs(((x_c-L_spec(freq,theta)*& 
                         sinA(freq,theta))*cosA(freq,theta)+& 
                       (y_c+L_spec(freq,theta)*& 
                         sinA(freq,theta))*sinA(freq,theta)))
                         
                 if(nn_loc.le.4)then !@@

                  cosA(freq,theta)=cos(inc_ang_spec(freq,theta))
                  sinA(freq,theta)=sin(inc_ang_spec(freq,theta))

                  wA(freq,theta)=2.*3.1415/per_spec(freq,theta)  
                  kA(freq,theta)=wA(freq,theta)/tanh(depth*2.*3.1415/L_spec(freq,theta))  ! this IS NOT WAVENUMBER, used by velocity source            
                  								
                 endif !@@
                 
			      if(loc/L_spec(freq,theta).le.0.75)then !$$
				   
				   is_coef=0.5*D_src(freq,theta)*exp(-beta_src(freq,theta)*loc**2)	

				   is_sin=sin(-wA(freq,theta)*(t(nn_loc)+shift_spec(freq,theta)))
    		
				   src_amp=is_coef*is_sin*tmp*kA(freq,theta)
				   
				   f_int_src(i,j)=f_int_src(i,j)+src_amp
                    
			     endif !$$
			     
			   enddo
			   
               
            enddo
             if(plus_tide.eq.1)then
					loc=x_c
				    is_coef=0.5*D_tide*exp(-beta_tide*loc**2)	

				    is_sin=eta_in(nn_loc)
    		
				    src_amp=is_coef*is_sin*tmp
				
					f_int_src(i,j)=f_int_src(i,j)-src_amp*sqrt(9.81/depth)  ! minus for this particular setup only

			endif           
            
		endif !%%
		
       enddo
     enddo
     
     endif !&&
     
	 do j = sy-1,ey
		do i = sx-1,ex
				  F(i,j) = F(i,j)+f_int_src(i,j)*H(i,j)
        enddo
     enddo
	endif !^^

!-------------start calculating E, F, and G

  !----------calculate E------------
   do j=sy,ey
      do i=sx-1,ex
         temp2 = 1./6.*(zi(i,j)**3+hi(i,j)**3)-0.5*zai(i,j)**2*(zi(i,j)+hi(i,j))
	     temp3 = 0.5*(zi(i,j)**2-hi(i,j)**2)-zai(i,j)*(zi(i,j)+hi(i,j))
	     E1(i,j) = temp2*(uixx(i,j))+temp3*(huixx(i,j))
      enddo
   enddo


   do j=sy,ey
      do i=sx,ex
	     E(i,j) = E(i,j)+(E1(i,j)-E1(i-1,j))/dx -(flux(i,j,1)-flux(i-1,j,1))/dx !- dhdt(i,j,n_loc)
      enddo
   enddo

   if (ihvor.eq.1) then
      call interp_xoi (psix(:,:,n_loc),tmpi1,interp2,ix)
	  tmpi2 = tmpi1*(zi+hi)*(zai*zai*0.5-zai*zi+(2*zi*zi-2*zi*hi-hi*hi)/6.0)
	  do j = sy,ey
		do i = sx,ex
		    E(i,j) = E(i,j) - (tmpi2(i,j)-tmpi2(i-1,j))/dx
		 enddo
      enddo
   endif

   !----------caclulate F,F1------------

   call interp_xei (E,Ei,interp2,ix)

   do j=sy,ey
      do i=sx-1,ex
         EzSTi(i,j) = Ei(i,j)*(zi(i,j)*(uix(i,j))+(huix(i,j)))
	     TzS2i(i,j) = (zi(i,j)*(uix(i,j)))**2 &
                    + 2.0*zi(i,j)*(uix(i,j))*(huix(i,j)) &
                    + (huix(i,j))**2
      enddo
   enddo
   do j=sy,ey
      do i=sx-1,ex
	     uSxvSyi(i,j) = 0.5*(zai(i,j)**2-zi(i,j)**2)*(ui(i,j)*(uixx(i,j))) 
	     uTxvTyi(i,j) = (zai(i,j)-zi(i,j))*(ui(i,j)*(huixx(i,j)))
      enddo
   enddo	   		   
   
   do j=sy,ey
      do i=sx,ex
		 temp1 = -(flux(i,j,2)-flux(i-1,j,2))/dx+u(i,j)*((E1(i,j)-E1(i-1,j))/dx)
	     temp3 = -(EzSTi(i,j)-EzSTi(i-1,j))/dx*H(i,j) &
                    + E(i,j)*(0.5*(za(i,j)**2-z(i,j)**2)*dSdx(i,j,1)+(za(i,j)-z(i,j))*dTdx(i,j,1) &
                              -zx(i,j)*(z(i,j)*S_grp(i,j,1)+T_grp(i,j,1)))
         temp4 = -(uSxvSyi(i,j)-uSxvSyi(i-1,j))/dx*H(i,j)
	     temp5 = -(uTxvTyi(i,j)-uTxvTyi(i-1,j))/dx*H(i,j)
	     temp6 = -0.5*(TzS2i(i,j)-TzS2i(i-1,j))/dx*H(i,j)
         temp7 = 0.
		 dhzudx(i,j,1)=(flux(i,j,1)-flux(i-1,j,1))/dx	  
	     F(i,j) = F(i,j) + (temp1+temp3+temp4+temp5+temp6+temp7*H(i,j)) -grav*H(i,j)*zx_4th(i,j) !+ grav*Hcx(i,j)*(d(i+1,j)-d(i-1,j))/dx*0.5
		 F1(i,j) = 0.
	  enddo
   enddo

	
  if(ihvor.eq.1)then

   call interp_xoi (psix(:,:,n_loc),tmpi1,interp2,ix) 

   do j = sy,ey
		do i = sx,ex
	     tmp1 = (ui(i,j)*tmpi1(i,j)*zi(i,j)-ui(i-1,j)*tmpi1(i-1,j)*zi(i-1,j))/dx 
	     tmp2 = (ui(i,j)*tmpi1(i,j)-ui(i-1,j)*tmpi1(i-1,j))/dx
	     tmp3 = (ui(i,j)*tmpi1(i,j)*(0.5*zai(i,j)*zai(i,j)-zi(i,j)*zai(i,j)) &
		             -ui(i-1,j)*tmpi1(i-1,j)*(0.5*zai(i-1,j)*zai(i-1,j)-zi(i-1,j)*zai(i-1,j)))/dx 
		 tmp4 = psix(i,j,n_loc)*((z(i,j)**2.+z(i,j)*d(i,j)-2*d(i,j)**2.)*(ux(i,j))/6.0+0.5*H(i,j)*(hux(i,j)))
		 
		 F(i,j) = F(i,j) - tmp1*(z(i,j)-d(i,j))*0.5*H(i,j) &
						 + tmp2*(z(i,j)**2.-z(i,j)*d(i,j)+d(i,j)**2.)/6.0*H(i,j)  &
						 - tmp3*H(i,j) &
						 + tmp4*H(i,j)
      enddo
   enddo
   
  endif

! correct for near shoreline, SW simulations
	do j=sy,ey
		do i=sx,ex
			if(ish_mat(i,j).eq.1)then
				E(i,j) = -dhzudx(i,j,1) !- dhdt(i,j,n_loc)
				F(i,j)=   - (flux(i,j,2)-flux(i-1,j,2))/dx - grav*H(i,j)*zx_4th(i,j)   !+ grav*Hcx(i,j)*(d(i+1,j)-d(i-1,j))/dx*0.5
				F1(i,j)=0.
			endif
		enddo
	enddo

	if(ibotfric.eq.1)then  ! add bottom stress
		do j=sy,ey
			do i=sx,ex
		      F(i,j)=F(i,j)-Rhx(i,j)
			enddo
		enddo
	endif

	windstress=0
	if(windstress.eq.1)then  ! add wind stress
		U_10=44.7 ! ~100 mph
		do j=sy,ey
			do i=sx,ex
			  Hs_i=max(cutoff_mat(i,j),d(i,j))
			  C_loc=sqrt(9.81*Hs_i)
			  Cd_loc=0.001*(0.2+18.*abs(zx_4th(i,j))+0.065*U_10)
			  Rwx=1.25/1025.*Cd_loc*abs(U_10-C_loc)*(U_10-C_loc)
		      F(i,j)=F(i,j)+Rwx
!		      if(z(i,j).gt.0.) F(i,j)=F(i,j)+Rwx  ! on crest only, only for wavy flows
			enddo
		enddo
	endif
	
	if (idissip.eq.1) then  ! add dissipation
	 do j=sy,ey      
      do i=sx,ex
         temp1 = (vti(i,j)*uix(i,j)-vti(i-1,j)*uix(i-1,j))/dx
         temp2 = 0.
         temp3 = vtv(i,j)*(uix(i,j) - uix(i-1,j))/dx
         temp4 = 0.
         F(i,j) = F(i,j) + (2.*temp1+temp2-1*temp3+temp4)*H(i,j)
      enddo
	 enddo
	endif  ! for dissipation terms

! Backscatter model
	if(backscatter.eq.1)then
		if(count_min.eq.1) call create_rand_matrix
		do j=sy,ey
			do i=sx,ex
			  BS = CB*sqrt(ua(i,j)**2)*sqrt(dissipcoef*sqrt(Cf(i,j))/dt)
!		      Hs_i=max(cutoff_mat(i,j),H(i,j))
!			  BS_breaker = CB*nut(i,j,n_loc)/Hs_i*sqrt(dissipcoef/dt)  ! need to calibrate leading coef, can't assume it to be the same as for bottom turbulence
		      F(i,j)=F(i,j)+BS*rand(i,j,1) !+ BS_breaker*rand(i,j,1)
       enddo
     enddo	
	endif	
	     
return
end
