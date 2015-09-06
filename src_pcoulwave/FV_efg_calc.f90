!*******************************************************************************
!........    Groups the physical variables into the 
!........... conservative forms given in Wei & Kirby
subroutine FV_efg_calc(z,d,u,v,E,F,F1,G,G1,n_loc,nn_loc)
use mainvar_module, only:dx,dy,nx,ny,grav,dvdyy,dudxx,dvdxy,&
dudxy,dhzudx,dhzudy,dhzvdx,dhzvdy,dim,bl_hor_wall,level, & 
 disp_prop,cutoff,bc_1,bc_2,bc_3,bc_4,overlap,endx,endy,&
dt,S_grp,T_grp,dhdt,dSdx,dTdx,dhdxt,dSdy,dTdy,dhdyt,cutoff_mat,nut,B_mat, & 
 f_int_src, f_int_src_v,theta, freq, num_theta, num_freq,&
 depth, x_c,y_c,end_x_t, end_y_t, x0, cosA,sinA,inc_ang_spec,L_spec,per_spec, & 
 shift_spec,is_coef,is_sin,loc,x,y,D_src,beta_src,t,&
count_min,int_src,is_oreint,tmp,CB,backscatter,ihvor,Ch,bf_type,psix,psiy,&
 elder_length,Ch_length,wA,kA,D_tide,beta_tide,&
eta_in,plus_tide,hx,hy,numerical_scheme,cutoff_mat_d
use FV_var_module, only:sx,ex,sy,ey,idissip,ibotfric,b0,b1,&
cm,ks,alpha,ilim,ix,iy,ix2,iy2,id,limiter2,limiter3,iriemann,icor, &
 interp2,zx,zy,rand,BS,ish_mat
      
integer sum,ii,jj,i,j,sum_0,n_loc,ilim_c,nn_loc,sum_id
     
real :: u(nx,ny),v(nx,ny),z(nx,ny), d(nx,ny),za(nx,ny)
real :: Hu(nx,ny),Hv(nx,ny),H(nx,ny)
real :: zi(nx,ny),zj(nx,ny)
real :: hi(nx,ny),hj(nx,ny)
real :: H_i(nx,ny),H_j(nx,ny)
real :: ui(nx,ny),uj(nx,ny)
real :: vi(nx,ny),vj(nx,ny)
real :: hui(nx,ny),huj(nx,ny)
real :: hvi(nx,ny),hvj(nx,ny)
real :: phil(nx,ny),phir(nx,ny)
real :: ux(nx,ny),vx(nx,ny)
real :: uy(nx,ny),vy(nx,ny)
real :: hux(nx,ny),hvy(nx,ny)
real :: Hcx(nx,ny),Hcy(nx,ny)
real :: uixx(nx,ny),uixy(nx,ny)
real :: vixy(nx,ny),viyy(nx,ny)
real :: uix(nx,ny),viy(nx,ny)
real :: huix(nx,ny),hviy(nx,ny)
real :: huixx(nx,ny),huixy(nx,ny)
real :: hvixy(nx,ny),hviyy(nx,ny)
real :: zai(nx,ny),Ei(nx,ny)
real :: vix(nx,ny),hvix(nx,ny)
real :: ujxx(nx,ny),ujxy(nx,ny)
real :: vjyy(nx,ny),vjxy(nx,ny)
real :: hujxx(nx,ny),hujxy(nx,ny)
real :: hvjyy(nx,ny),hvjxy(nx,ny)
real :: hujx(nx,ny),hvjy(nx,ny)
real :: vjy(nx,ny),ujx(nx,ny)
real :: ujy(nx,ny),hujy(nx,ny)
real :: zaj(nx,ny),Ej(nx,ny)
real :: uxy(nx,ny),vxy(nx,ny)
real :: huxy(nx,ny),hvxy(nx,ny)
real :: cor(nx,ny),dim_coef
real :: uxx(nx,ny),vyy(nx,ny),zx_4th(nx,ny),zy_4th(nx,ny)
real :: huxx(nx,ny),hvyy(nx,ny)
     
real :: uyy(nx,ny),vxx(nx,ny),dissipcoef,u_DA,v_DA,w_DA,k_tmp,w_tmp
real :: vt(nx,ny),vti(nx,ny),vtj(nx,ny),src_amp
real :: tmpi1(nx,ny),tmpj1(nx,ny)
real :: tmpi2(nx,ny),tmpj2(nx,ny)
real :: tmp1(nx,ny),tmp2(nx,ny)
real :: E(nx,ny),F(nx,ny),F1(nx,ny)
real :: G(nx,ny),G1(nx,ny)
real :: temp1,temp2,temp3,temp4,temp5,temp6,temp7
real :: E1(nx,ny),E2(nx,ny)
real :: H_l(nx,ny),H_r(nx,ny)
real :: Hu_l(nx,ny),Hu_r(nx,ny)
real :: Hv_l(nx,ny),Hv_r(nx,ny)
real :: H_b(nx,ny),H_t(nx,ny)
real :: Hu_b(nx,ny),Hs_i,us_i,vs_i
real :: Hu_t(nx,ny),Hv_b(nx,ny)
real :: Hv_t(nx,ny),conv_u,conv_v
real :: flux(nx,ny,3),fluy(nx,ny,3)
real :: dh0,dh1,dh2,dh3,dh4
real :: dbh0,dbh1,dbh2,dbh3,dbh4,d3h,dsh1,dsh2,dsh3
real :: dshb1,dshb2,dshb3,dshb4
real :: EzSTi(nx,ny),TzS2i(nx,ny)
real :: EzSTj(nx,ny),TzS2j(nx,ny)
real :: uSxvSyi(nx,ny),uTxvTyi(nx,ny)
real :: uSxvSyj(nx,ny),uTxvTyj(nx,ny)
real :: huxzy(nx,ny),hvyzx(nx,ny),mx_dzdx

real :: ua(nx,ny),va(nx,ny),Rhx(nx,ny),Rhy(nx,ny),vtv(nx,ny),xc,conv_reduc(nx)
real :: Rh, Cf_c,Cf(nx,ny), ds_i, zs_i, L_tanh,mx_theta,BS_breaker,ks_i,cm_i

real minH, Hmin, phir1, phir2,phir1i,phir2i,omeg,dupw,dloc,delta,up_ratio
integer k,wn,wnx(20),wny(20)
parameter(dissipcoef=1.e-6)

minH = dx/10.
Hmin = dx/10.

! conserved variables
H = d+z
Hu = H*u
Hv = H*v

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


do j = sy,ey
	do i = sx-1,ex

!.....k = 1
         ilim_c=ilim

		 dh0 = z(i-1,j)-z(i-2,j)
		 dh1 = z(i  ,j)-z(i-1,j)
		 dh2 = z(i+1,j)-z(i  ,j)
		 dh3 = z(i+2,j)-z(i+1,j)
		 dh4 = z(i+3,j)-z(i+2,j)

		 mx_dzdx=max( abs(dh0),abs(dh1), abs(dh2), abs(dh3), abs(dh4) )/dx
		 if(mx_dzdx.gt.0.60) ilim_c=1

		 dbh0 = sign(1.0,dh0)*max(0.0,min(abs(dh0),b1*&
sign(1.0,dh0)*dh1,b1*sign(1.0,dh0)*dh2))
		 dbh1 = sign(1.0,dh1)*max(0.0,min(abs(dh1),b1*&
sign(1.0,dh1)*dh2,b1*sign(1.0,dh1)*dh0))
		 dbh2 = sign(1.0,dh2)*max(0.0,min(abs(dh2),b1*&
sign(1.0,dh2)*dh0,b1*sign(1.0,dh2)*dh1))

		 d3h =  dbh2-2*dbh1+dbh0
		 dsh1 = dh1-d3h/6

		 dbh1 = sign(1.0,dh1)*max(0.0,min(abs(dh1),b1*&
sign(1.0,dh1)*dh2,b1*sign(1.0,dh1)*dh3))
		 dbh2 = sign(1.0,dh2)*max(0.0,min(abs(dh2),b1*&
sign(1.0,dh2)*dh3,b1*sign(1.0,dh2)*dh1))
		 dbh3 = sign(1.0,dh3)*max(0.0,min(abs(dh3),b1*&
sign(1.0,dh3)*dh1,b1*sign(1.0,dh3)*dh2))

		 d3h =  dbh3-2*dbh2+dbh1
		 dsh2 = dh2-d3h/6

		 dbh2 = sign(1.0,dh2)*max(0.0,min(abs(dh2),b1*&
sign(1.0,dh2)*dh3,b1*sign(1.0,dh2)*dh4))
		 dbh3 = sign(1.0,dh3)*max(0.0,min(abs(dh3),b1*&
sign(1.0,dh3)*dh4,b1*sign(1.0,dh3)*dh2))
		 dbh4 = sign(1.0,dh4)*max(0.0,min(abs(dh4),b1*&
sign(1.0,dh4)*dh2,b1*sign(1.0,dh4)*dh3))

		 d3h =  dbh4-2*dbh3+dbh2
		 dsh3 = dh3-d3h/6

		 if (ilim_c.eq.0) then

		    phir1 = 1.
		    phir2 = 1.
		    phir1i = 1.
		    phir2i = 1.

	        H_l(i,j) = z(i,j)   + (phir1*dsh1+phir1i*2*dsh2)&
/6.0 + 0.5*(d(i,j)+d(i+1,j))
            H_r(i,j) = z(i+1,j) - (2*phir2*dsh2+phir2i*dsh3)&
/6.0 + 0.5*(d(i,j)+d(i+1,j))

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

		 dbh0 = sign(1.0,dh0)*max(0.0,min(abs(dh0),b1*&
sign(1.0,dh0)*dh1,b1*sign(1.0,dh0)*dh2))
		 dbh1 = sign(1.0,dh1)*max(0.0,min(abs(dh1),b1*&
sign(1.0,dh1)*dh2,b1*sign(1.0,dh1)*dh0))
		 dbh2 = sign(1.0,dh2)*max(0.0,min(abs(dh2),b1*&
sign(1.0,dh2)*dh0,b1*sign(1.0,dh2)*dh1))

		 d3h =  dbh2-2*dbh1+dbh0
		 dsh1 = dh1-d3h/6

		 dbh1 = sign(1.0,dh1)*max(0.0,min(abs(dh1),b1*&
sign(1.0,dh1)*dh2,b1*sign(1.0,dh1)*dh3))
		 dbh2 = sign(1.0,dh2)*max(0.0,min(abs(dh2),b1*&
sign(1.0,dh2)*dh3,b1*sign(1.0,dh2)*dh1))
		 dbh3 = sign(1.0,dh3)*max(0.0,min(abs(dh3),b1*&
sign(1.0,dh3)*dh1,b1*sign(1.0,dh3)*dh2))

		 d3h =  dbh3-2*dbh2+dbh1
		 dsh2 = dh2-d3h/6

		 dbh2 = sign(1.0,dh2)*max(0.0,min(abs(dh2),b1*&
sign(1.0,dh2)*dh3,b1*sign(1.0,dh2)*dh4))
		 dbh3 = sign(1.0,dh3)*max(0.0,min(abs(dh3),b1*&
sign(1.0,dh3)*dh4,b1*sign(1.0,dh3)*dh2))
		 dbh4 = sign(1.0,dh4)*max(0.0,min(abs(dh4),b1*&
sign(1.0,dh4)*dh2,b1*sign(1.0,dh4)*dh3))

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

!.....k = 3
      
		 dh0 = Hv(i-1,j)-Hv(i-2,j)
		 dh1 = Hv(i  ,j)-Hv(i-1,j)
		 dh2 = Hv(i+1,j)-Hv(i  ,j)
		 dh3 = Hv(i+2,j)-Hv(i+1,j)
		 dh4 = Hv(i+3,j)-Hv(i+2,j)

		 dbh0 = sign(1.0,dh0)*max(0.0,min(abs(dh0),b1*&
sign(1.0,dh0)*dh1,b1*sign(1.0,dh0)*dh2))
		 dbh1 = sign(1.0,dh1)*max(0.0,min(abs(dh1),b1*&
sign(1.0,dh1)*dh2,b1*sign(1.0,dh1)*dh0))
		 dbh2 = sign(1.0,dh2)*max(0.0,min(abs(dh2),b1*&
sign(1.0,dh2)*dh0,b1*sign(1.0,dh2)*dh1))

		 d3h =  dbh2-2*dbh1+dbh0
		 dsh1 = dh1-d3h/6

		 dbh1 = sign(1.0,dh1)*max(0.0,min(abs(dh1),b1*&
sign(1.0,dh1)*dh2,b1*sign(1.0,dh1)*dh3))
		 dbh2 = sign(1.0,dh2)*max(0.0,min(abs(dh2),b1*&
sign(1.0,dh2)*dh3,b1*sign(1.0,dh2)*dh1))
		 dbh3 = sign(1.0,dh3)*max(0.0,min(abs(dh3),b1*&
sign(1.0,dh3)*dh1,b1*sign(1.0,dh3)*dh2))

		 d3h =  dbh3-2*dbh2+dbh1
		 dsh2 = dh2-d3h/6

		 dbh2 = sign(1.0,dh2)*max(0.0,min(abs(dh2),b1*&
sign(1.0,dh2)*dh3,b1*sign(1.0,dh2)*dh4))
		 dbh3 = sign(1.0,dh3)*max(0.0,min(abs(dh3),b1*&
sign(1.0,dh3)*dh4,b1*sign(1.0,dh3)*dh2))
		 dbh4 = sign(1.0,dh4)*max(0.0,min(abs(dh4),b1*&
sign(1.0,dh4)*dh2,b1*sign(1.0,dh4)*dh3))

		 d3h =  dbh4-2*dbh3+dbh2
		 dsh3 = dh3-d3h/6

 	     if (ilim_c.eq.0) then

		    phir1 = 1.
            phir2 = 1.
		    phir1i = 1.
		    phir2i = 1.

 	        Hv_l(i,j) = Hv(i,j)   + (phir1*dsh1+phir1i*2*dsh2)/6.0
	        Hv_r(i,j) = Hv(i+1,j) - (2*phir2*dsh2+phir2i*dsh3)/6.0

		 elseif (ilim_c.eq.1) then

	        dshb1 = sign(1.0,dsh1)*max(0.0,min(abs(dsh1),b0*&
dsh2*sign(1.0,dsh1)))
	        dshb2 = sign(1.0,dsh2)*max(0.0,min(abs(dsh2),b0*&
dsh1*sign(1.0,dsh2)))
	        dshb3 = sign(1.0,dsh2)*max(0.0,min(abs(dsh2),b0*&
dsh3*sign(1.0,dsh2)))	     
		dshb4 = sign(1.0,dsh3)*max(0.0,min(abs(dsh3),b0*&
dsh2*sign(1.0,dsh3)))

	        Hv_l(i,j) = Hv(i,j)   + (dshb1+2*dshb2)/6.0
            Hv_r(i,j) = Hv(i+1,j) - (2*dshb3+dshb4)/6.0

		 endif

    enddo
enddo

do j = sy-1,ey
	do i = sx,ex

!.....k = 4
         ilim_c=ilim
      
		 dh0 = z(i,j-1)-z(i,j-2)
		 dh1 = z(i  ,j)-z(i,j-1)
		 dh2 = z(i,j+1)-z(i  ,j)
		 dh3 = z(i,j+2)-z(i,j+1)
		 dh4 = z(i,j+3)-z(i,j+2)

		 mx_dzdx=max( abs(dh0),abs(dh1), abs(dh2), abs(dh3), abs(dh4) )/dy
		 if(mx_dzdx.gt.0.60) ilim_c=1

		 dbh0 = sign(1.0,dh0)*max(0.0,min(abs(dh0),b1*&
sign(1.0,dh0)*dh1,b1*sign(1.0,dh0)*dh2))
		 dbh1 = sign(1.0,dh1)*max(0.0,min(abs(dh1),b1*&
sign(1.0,dh1)*dh2,b1*sign(1.0,dh1)*dh0))
		 dbh2 = sign(1.0,dh2)*max(0.0,min(abs(dh2),b1*&
sign(1.0,dh2)*dh0,b1*sign(1.0,dh2)*dh1))

		 d3h =  dbh2-2*dbh1+dbh0
		 dsh1 = dh1-d3h/6

		 dbh1 = sign(1.0,dh1)*max(0.0,min(abs(dh1),b1*&
sign(1.0,dh1)*dh2,b1*sign(1.0,dh1)*dh3))
		 dbh2 = sign(1.0,dh2)*max(0.0,min(abs(dh2),b1*&
sign(1.0,dh2)*dh3,b1*sign(1.0,dh2)*dh1))
		 dbh3 = sign(1.0,dh3)*max(0.0,min(abs(dh3),b1*&
sign(1.0,dh3)*dh1,b1*sign(1.0,dh3)*dh2))

		 d3h =  dbh3-2*dbh2+dbh1
		 dsh2 = dh2-d3h/6

		 dbh2 = sign(1.0,dh2)*max(0.0,min(abs(dh2),b1*&
sign(1.0,dh2)*dh3,b1*sign(1.0,dh2)*dh4))
		 dbh3 = sign(1.0,dh3)*max(0.0,min(abs(dh3),b1*&
sign(1.0,dh3)*dh4,b1*sign(1.0,dh3)*dh2))
		 dbh4 = sign(1.0,dh4)*max(0.0,min(abs(dh4),b1*&
sign(1.0,dh4)*dh2,b1*sign(1.0,dh4)*dh3))

		 d3h =  dbh4-2*dbh3+dbh2
		 dsh3 = dh3-d3h/6

		 if (ilim_c.eq.0) then

		    phir1 = 1.
		    phir2 = 1.
		    phir1i = 1.
		    phir2i = 1.

!.... Extended variables

	        H_b(i,j) = z(i,j)   + (phir1*dsh1+phir1i*2*dsh2)&
/6.0 + 0.5*(d(i,j)+d(i,j+1))
            H_t(i,j) = z(i,j+1) - (2*phir2*dsh2+phir2i*dsh3)&
/6.0 + 0.5*(d(i,j)+d(i,j+1))


		 elseif (ilim_c.eq.1) then

	        dshb1 = sign(1.0,dsh1)*max(0.0,min(abs(dsh1),b0*&
dsh2*sign(1.0,dsh1)))
	        dshb2 = sign(1.0,dsh2)*max(0.0,min(abs(dsh2),b0*&
dsh1*sign(1.0,dsh2)))
	        dshb3 = sign(1.0,dsh2)*max(0.0,min(abs(dsh2),b0*&
dsh3*sign(1.0,dsh2)))	     
		    dshb4 = sign(1.0,dsh3)*max(0.0,min(abs(dsh3),b0*&
dsh2*sign(1.0,dsh3)))

	        H_b(i,j) = z(i,j)   + (dshb1+2*dshb2)&
/6.0 + 0.5*(d(i,j)+d(i,j+1))
            H_t(i,j) = z(i,j+1) - (2*dshb3+dshb4)&
/6.0 + 0.5*(d(i,j)+d(i,j+1))


		 endif

!.....k = 5
      
		 dh0 = Hu(i,j-1)-Hu(i,j-2)
		 dh1 = Hu(i  ,j)-Hu(i,j-1)
		 dh2 = Hu(i,j+1)-Hu(i  ,j)
		 dh3 = Hu(i,j+2)-Hu(i,j+1)
		 dh4 = Hu(i,j+3)-Hu(i,j+2)


		 dbh0 = sign(1.0,dh0)*max(0.0,min(abs(dh0),b1*&
sign(1.0,dh0)*dh1,b1*sign(1.0,dh0)*dh2))
		 dbh1 = sign(1.0,dh1)*max(0.0,min(abs(dh1),b1*&
sign(1.0,dh1)*dh2,b1*sign(1.0,dh1)*dh0))
		 dbh2 = sign(1.0,dh2)*max(0.0,min(abs(dh2),b1*&
sign(1.0,dh2)*dh0,b1*sign(1.0,dh2)*dh1))

		 d3h =  dbh2-2*dbh1+dbh0
		 dsh1 = dh1-d3h/6

		 dbh1 = sign(1.0,dh1)*max(0.0,min(abs(dh1),b1*&
sign(1.0,dh1)*dh2,b1*sign(1.0,dh1)*dh3))
		 dbh2 = sign(1.0,dh2)*max(0.0,min(abs(dh2),b1*&
sign(1.0,dh2)*dh3,b1*sign(1.0,dh2)*dh1))
		 dbh3 = sign(1.0,dh3)*max(0.0,min(abs(dh3),b1*&
sign(1.0,dh3)*dh1,b1*sign(1.0,dh3)*dh2))

		 d3h =  dbh3-2*dbh2+dbh1
		 dsh2 = dh2-d3h/6

		 dbh2 = sign(1.0,dh2)*max(0.0,min(abs(dh2),b1*&
sign(1.0,dh2)*dh3,b1*sign(1.0,dh2)*dh4))
		 dbh3 = sign(1.0,dh3)*max(0.0,min(abs(dh3),b1*&
sign(1.0,dh3)*dh4,b1*sign(1.0,dh3)*dh2))
		 dbh4 = sign(1.0,dh4)*max(0.0,min(abs(dh4),b1*&
sign(1.0,dh4)*dh2,b1*sign(1.0,dh4)*dh3))

		 d3h =  dbh4-2*dbh3+dbh2
		 dsh3 = dh3-d3h/6

		 if (ilim_c.eq.0) then

		    phir1 = 1.
		    phir2 = 1.
		    phir1i = 1.
		    phir2i = 1.
 
	        Hu_b(i,j) = Hu(i,j)   + (phir1*dsh1+phir1i*2*dsh2)/6.0
	        Hu_t(i,j) = Hu(i,j+1) - (2*phir2*dsh2+phir2i*dsh3)/6.0

		 elseif (ilim_c.eq.1) then

	        dshb1 = sign(1.0,dsh1)*max(0.0,min(abs(dsh1),b0*&
dsh2*sign(1.0,dsh1)))
	        dshb2 = sign(1.0,dsh2)*max(0.0,min(abs(dsh2),b0*&
dsh1*sign(1.0,dsh2)))
	        dshb3 = sign(1.0,dsh2)*max(0.0,min(abs(dsh2),b0*&
dsh3*sign(1.0,dsh2)))	     
		dshb4 = sign(1.0,dsh3)*max(0.0,min(abs(dsh3),b0*&
dsh2*sign(1.0,dsh3)))

	        Hu_b(i,j) = Hu(i,j)   + (dshb1+2*dshb2)/6.0
            Hu_t(i,j) = Hu(i,j+1) - (2*dshb3+dshb4)/6.0

		 endif

!.....k = 6
      
		 dh0 = Hv(i,j-1)-Hv(i,j-2)
		 dh1 = Hv(i  ,j)-Hv(i,j-1)
		 dh2 = Hv(i,j+1)-Hv(i  ,j)
		 dh3 = Hv(i,j+2)-Hv(i,j+1)
		 dh4 = Hv(i,j+3)-Hv(i,j+2)

		 dbh0 = sign(1.0,dh0)*max(0.0,min(abs(dh0),b1*&
sign(1.0,dh0)*dh1,b1*sign(1.0,dh0)*dh2))
		 dbh1 = sign(1.0,dh1)*max(0.0,min(abs(dh1),b1*&
sign(1.0,dh1)*dh2,b1*sign(1.0,dh1)*dh0))
		 dbh2 = sign(1.0,dh2)*max(0.0,min(abs(dh2),b1*&
sign(1.0,dh2)*dh0,b1*sign(1.0,dh2)*dh1))

		 d3h =  dbh2-2*dbh1+dbh0
		 dsh1 = dh1-d3h/6

		 dbh1 = sign(1.0,dh1)*max(0.0,min(abs(dh1),b1*&
sign(1.0,dh1)*dh2,b1*sign(1.0,dh1)*dh3))
		 dbh2 = sign(1.0,dh2)*max(0.0,min(abs(dh2),b1*&
sign(1.0,dh2)*dh3,b1*sign(1.0,dh2)*dh1))
		 dbh3 = sign(1.0,dh3)*max(0.0,min(abs(dh3),b1*&
sign(1.0,dh3)*dh1,b1*sign(1.0,dh3)*dh2))

		 d3h =  dbh3-2*dbh2+dbh1
		 dsh2 = dh2-d3h/6

		 dbh2 = sign(1.0,dh2)*max(0.0,min(abs(dh2),b1*&
sign(1.0,dh2)*dh3,b1*sign(1.0,dh2)*dh4))
		 dbh3 = sign(1.0,dh3)*max(0.0,min(abs(dh3),b1*&
sign(1.0,dh3)*dh4,b1*sign(1.0,dh3)*dh2))
		 dbh4 = sign(1.0,dh4)*max(0.0,min(abs(dh4),b1*&
sign(1.0,dh4)*dh2,b1*sign(1.0,dh4)*dh3))

		 d3h =  dbh4-2*dbh3+dbh2
		 dsh3 = dh3-d3h/6

		 if (ilim_c.eq.0) then

		    phir1 = 1.
		    phir2 = 1.
		    phir1i = 1.
		    phir2i = 1.

!.... Extended variables

	        Hv_b(i,j) = Hv(i,j)   + (phir1*dsh1+phir1i*2*dsh2)/6.0
	        Hv_t(i,j) = Hv(i,j+1) - (2*phir2*dsh2+phir2i*dsh3)/6.0

		 elseif (ilim_c.eq.1) then

	        dshb1 = sign(1.0,dsh1)*max(0.0,min(abs(dsh1),b0*&
dsh2*sign(1.0,dsh1)))
            dshb2 = sign(1.0,dsh2)*max(0.0,min(abs(dsh2),b0*&
dsh1*sign(1.0,dsh2)))
	        dshb3 = sign(1.0,dsh2)*max(0.0,min(abs(dsh2),b0*&
dsh3*sign(1.0,dsh2)))	     
		dshb4 = sign(1.0,dsh3)*max(0.0,min(abs(dsh3),b0*&
dsh2*sign(1.0,dsh3)))

	        Hv_b(i,j) = Hv(i,j)   + (dshb1+2*dshb2)/6.0
            Hv_t(i,j) = Hv(i,j+1) - (2*dshb3+dshb4)/6.0

		 endif

   enddo
enddo


   do j = sy,ey
	  do i = sx,ex
	     Hcx(i,j) = 0.5*(H_r(i-1,j)+H_l(i,j))
	     Hcy(i,j) = 0.5*(H_t(i,j-1)+H_b(i,j))
      enddo
   enddo

!..... call approximate Riemann solvers
call Riemanns4 (H_l,H_r,Hu_l,Hu_r,Hv_l,Hv_r,H_b,H_t,Hu_b,&
Hu_t,Hv_b,Hv_t,flux,fluy,1,sx,ex,sy,ey,ix,iy,iriemann)


if(numerical_scheme.eq.0)then
! HO FD for leading order terms
!
    call interp_xei (H,H_i,1,ix)  ! fourth order interpolation
    call interp_xei (u,ui,1,ix)  ! fourth order interpolation
    call interp_xei (v,vi,1,ix)  ! fourth order interpolation
    call interp_yei (H,H_j,1,iy)  ! fourth order interpolation
    call interp_yei (u,uj,1,iy)  ! fourth order interpolation
    call interp_yei (v,vj,1,iy)  ! fourth order interpolation
	
	do j=sy-1,ey
	   do i=sx-1,ex
	     if(id(i,j).ne.0)then

			flux(i,j,1)=H_i(i,j)*ui(i,j)
			flux(i,j,2)=H_i(i,j)*ui(i,j)*ui(i,j) ! + 0.5*grav*Hs_i**2.
			flux(i,j,3)=H_i(i,j)*ui(i,j)*vi(i,j)

			fluy(i,j,1)=H_j(i,j)*vj(i,j)
			fluy(i,j,2)=H_j(i,j)*vj(i,j)*uj(i,j)
			fluy(i,j,3)=H_j(i,j)*vj(i,j)*vj(i,j) ! + 0.5*grav*Hs_i**2.

		 endif
		enddo
	enddo

endif



    call interp_xei (z,zi,1,ix)  ! fourth order interpolation
	call interp_yei (z,zj,1,iy) ! fourth order interpolation
	do j=sy,ey
	   do i=sx,ex
			zx_4th(i,j)  = (zi(i,j)-zi(i-1,j))/dx
			zy_4th(i,j)  = (zj(i,j)-zj(i,j-1))/dy	
		enddo
	enddo


    do i=1,nx-1  ! special for USGS MHDemo - to reduce internal source vorticity blowup1
! 		xc=3.*(x(i)-4500.)/500.
 		conv_reduc(i)=1. !0.5*(1.+tanh(xc))
! 		do j=1,ny-1
!			flux(i,j,1)=(H(i+1,j)+H(i,j))*(u(i+1,j)+u(i,j))/4.*(1.-conv_reduc(i))+flux(i,j,1)*conv_reduc(i)
! 			fluy(i,j,1)=(H(i,j+1)+H(i,j))*(v(i,j+1)+v(i,j))/4.*(1.-conv_reduc(i))+fluy(i,j,1)*conv_reduc(i)
!   			flux(i,j,2)=flux(i,j,2)*conv_reduc(i)
!   			flux(i,j,3)=flux(i,j,3)*conv_reduc(i)
!			fluy(i,j,2)=fluy(i,j,2)*conv_reduc(i)
! 			fluy(i,j,3)=fluy(i,j,3)*conv_reduc(i)
!	   enddo
	enddo

! correct near shoreline flux
	do j=sy-1,ey
	   do i=sx-1,ex
		if(id(i,j).eq.9)then
			flux(i,j,1)=0.
			flux(i,j,2)=0.
			flux(i,j,3)=0.
			fluy(i,j,1)=0.
			fluy(i,j,2)=0.
			fluy(i,j,3)=0.
			zx_4th(i,j)=0.
			zy_4th(i,j)=0.
		elseif(d(i,j).le.cutoff_mat(i,j).and.id(i,j).eq.0)then
			zx_4th(i,j)= (z(i+1,j)-z(i-1,j))/dx/2.	
			zy_4th(i,j)= (z(i,j+1)-z(i,j-1))/dy/2.	
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

		endif				

		if(id(i,j).eq.4.or.id(i,j).eq.14.or.id(i,j).eq.24)then
			
			call shorefluy(i,j,2,zy_4th,fluy,z,u,v,d,dy)  ! calcs flux at i,j for id=4

			call shorefluy(i,j-1,3,zy_4th,fluy,z,u,v,d,dy)  ! calcs upwinded flux at i,j

		elseif(id(i,j).eq.3.or.id(i,j).eq.13.or.id(i,j).eq.23)then

			call shorefluy(i,j,1,zy_4th,fluy,z,u,v,d,dy)  ! calcs flux at i,j for id=3

			call shorefluy(i,j,3,zy_4th,fluy,z,u,v,d,dy)  ! calcs upwinded flux at i,j

		endif				

	  enddo
	enddo

	! channels of single grid point width
	do j=sy,ey
	   do i=sx,ex
		if(level(i,j).eq.0)then

			if(level(i,j-1).eq.1.and.level(i,j+1)&
.eq.1.and.level(i+1,j).eq.1.and.level(i-1,j).eq.0)then  ! should be id=2, but with v=0
				
				call shoreflux(i,j,2,zx_4th,flux,z,u,v,d,dx)  ! calcs flux at i,j for id=2

				call shoreflux(i-1,j,3,zx_4th,flux,z,u,v,d,dx)  ! calcs upwinded flux at i,j

				if(level(i-2,j).ne.0) zx_4th(i-1,j)=zx_4th(i,j)
				flux(i-1,j,3)=0.
				flux(i,j,3)=0.


			elseif(level(i,j-1).eq.&
1.and.level(i,j+1).eq.1.and.level(i+1,j).eq.0.and.level(i-1,j).eq.1)then   ! should be id=1, but with v=0
				
				call shoreflux(i,j,1,zx_4th,flux,z,u,v,d,dx)  ! calcs flux at i-1,j for id=1

				call shoreflux(i,j,3,zx_4th,flux,z,u,v,d,dx)  ! calcs upwinded flux at i,j
				
				if(level(i+2,j).ne.0) zx_4th(i+1,j)=zx_4th(i,j)
				flux(i-1,j,3)=0.
				flux(i,j,3)=0.


			elseif(level(i,j-1).eq.1.and.&
level(i,j+1).eq.0.and.level(i+1,j).eq.1.and.level(i-1,j).eq.1)then   ! should be id=3, but with u=0
				
				call shorefluy(i,j,1,zy_4th,fluy,z,u,v,d,dy)  ! calcs flux at i,j-1 for id=3

				call shorefluy(i,j,3,zy_4th,fluy,z,u,v,d,dy)  ! calcs upwinded flux at i,j

				if(level(i,j+2).ne.0) zy_4th(i,j+1)=zy_4th(i,j)
				fluy(i,j-1,2)=0.
				fluy(i,j,2)=0.

						
			elseif(level(i,j-1).eq.0.and.level(i,j+1)&
.eq.1.and.level(i+1,j).eq.1.and.level(i-1,j).eq.1)then   ! should be id=4, but with u=0
				
				call shorefluy(i,j,2,zy_4th,fluy,z,u,v,d,dy)  ! calcs flux at i,j for id=4

				call shorefluy(i,j-1,3,zy_4th,fluy,z,u,v,d,dy)  ! calcs upwinded flux at i,j
		
				if(level(i,j-2).ne.0) zy_4th(i,j-1)=zy_4th(i,j)
				fluy(i,j-1,2)=0.
				fluy(i,j,2)=0.


			endif
		endif
	  enddo
	enddo

	
	! lone points
	do j=sy,ey
	   do i=sx,ex
		if(level(i,j).eq.0)then  
			if(level(i-1,j).eq.1.and.level(i+1,j).eq.1)then  ! one wet alone
			   
				zx_4th(i-1,j)= 0.	
				zx_4th(i,j)= 0.
				zx_4th(i+1,j)= 0.	

				flux(i-1,j,1)=0.
				flux(i-1,j,2)=0.
				flux(i-1,j,3)=0.

				flux(i,j,1)=0.
				flux(i,j,2)=0.
				flux(i,j,3)=0.	

			elseif(level(i-1,j).eq.1.and.level(i+2,j).eq.1)then  ! two wets alone
				if(id(i,j).gt.10.and.id(i+1,j).gt.10)then
					zx_4th(i,j)= 0.
					zx_4th(i+1,j)= 0.

					flux(i-1,j,1)=0.
					flux(i-1,j,2)=0.
					flux(i-1,j,3)=0.

					flux(i,j,1)=0.
					flux(i,j,2)=0.
					flux(i,j,3)=0.
					
					flux(i+1,j,1)=0.
					flux(i+1,j,2)=0.
					flux(i+1,j,3)=0.
				else
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
					flux(i,j,3)=0.
				endif
			elseif(level(i-1,j).eq.1.and.level(i+3,j).eq.1)then  ! three wets alone, giving me problems for some reason.., usually in T, freeze em
				if(id(i,j).gt.10.and.id(i+2,j).gt.10)then
					zx_4th(i,j)= 0.
					zx_4th(i+2,j)= 0.

					flux(i-1,j,1)=0.
					flux(i-1,j,2)=0.
					flux(i-1,j,3)=0.

					if(id(i+1,j).ne.0)then
						zx_4th(i+1,j)= 0.
						flux(i,j,1)=0.
						flux(i,j,2)=0.
						flux(i,j,3)=0.
					
						flux(i+1,j,1)=0.
						flux(i+1,j,2)=0.
						flux(i+1,j,3)=0.
					else
						zx_4th(i+1,j)= (z(i+2,j)-z(i,j))/dx/2.
					endif

					flux(i+2,j,1)=0.
					flux(i+2,j,2)=0.
					flux(i+2,j,3)=0.

				endif
			endif		

			if(level(i,j-1).eq.1.and.level(i,j+1).eq.1)then  ! one wet alone
				zy_4th(i,j-1)= 0.	
				zy_4th(i,j)= 0.
				zy_4th(i,j+1)= 0.	

				fluy(i,j-1,1)=0.
				fluy(i,j-1,2)=0.
				fluy(i,j-1,3)=0.

				fluy(i,j,1)=0.
				fluy(i,j,2)=0.
				fluy(i,j,3)=0.

			elseif(level(i,j-1).eq.1.and.level(i,j+2).eq.1)then! two wets alone
				if(id(i,j).gt.10.and.id(i,j+1).gt.10)then
					zy_4th(i,j)=0.
					zy_4th(i,j+1)= 0.

					fluy(i,j-1,1)=0.
					fluy(i,j-1,2)=0.
					fluy(i,j-1,3)=0.

					fluy(i,j,1)=0.
					fluy(i,j,2)=0.
					fluy(i,j,3)=0.
										
					fluy(i,j+1,1)=0.
					fluy(i,j+1,2)=0.
					fluy(i,j+1,3)=0.
				else
					zy_4th(i,j-1)= 0.	
					if(z(i,j).gt.z(i,j-1))then
						zy_4th(i,j)= (z(i,j+1)-z(i,j-1))/dy/2.
					else
						zy_4th(i,j)= (z(i,j+1)-z(i,j))/dy
					endif
					if(z(i,j+1).gt.z(i,j+2))then
						zy_4th(i,j+1)= (z(i,j+2)-z(i,j))/dy/2.
					else
						zy_4th(i,j+1)= (z(i,j+1)-z(i,j))/dy
					endif
					zy_4th(i,j+2)= 0.
				
					call shorefluy(i,j,3,zy_4th,fluy,z,u,v,d,dy)  ! calcs upwinded flux at i,j
					fluy(i,j,2)=0.
				endif
			elseif(level(i,j-1).eq.1.and.level(i,j+3).eq.1)then ! three wets alone, giving me problems for some reason.., usually in T, freeze em
				if(id(i,j).gt.10.and.id(i,j+2).gt.10)then
					zy_4th(i,j)=0.
					zy_4th(i,j+2)= 0.

					fluy(i,j-1,1)=0.
					fluy(i,j-1,2)=0.
					fluy(i,j-1,3)=0.

					if(id(i,j+1).ne.0)then
						zy_4th(i,j+1)= 0.
						fluy(i,j,1)=0.
						fluy(i,j,2)=0.
						fluy(i,j,3)=0.
										
						fluy(i,j+1,1)=0.
						fluy(i,j+1,2)=0.
						fluy(i,j+1,3)=0.
					else
						zy_4th(i,j+1)=(z(i,j+2)-z(i,j))/dy/2.
					endif

					fluy(i,j+2,1)=0.
					fluy(i,j+2,2)=0.
					fluy(i,j+2,3)=0.
				endif							
			endif

		endif
	  enddo
	enddo

! special boundary corrections
	do j=sy,ey
	   do i=sx,ex
	    if(level(i,j).eq.1)then
			if(level(i-1,j).eq.0.and.level(i+1,j)&
.eq.0.and.level(i,j-1).eq.0.and.level(i,j+1).eq.0)then  ! one dry alone	
				zx_4th(i,j)= (z(i+1,j)-z(i-1,j))/dx/2.
				zy_4th(i,j)= (z(i,j+1)-z(i,j-1))/dy/2.
			endif
		endif
	  enddo
	enddo

call interp_xei (z,zi,interp2,ix)
call interp_yei (z,zj,interp2,iy)
call interp_xei (d,hi,interp2,ix)
call interp_yei (d,hj,interp2,iy)
call interp_xoi (u,ui,interp2,ix) 
call interp_yoi (v,vj,interp2,iy)
call interp_xei (v,vi,interp2,ix)
call interp_yei (u,uj,interp2,iy)

!za  = alpha*d
!zai = alpha*hi
!zaj = alpha*hj

za = -d+(1+alpha)*H
zai = -hi+(1+alpha)*(zi+hi)
zaj = -hj+(1+alpha)*(zj+hj)

hu = d*u
hv = d*v

hui = hi*ui
hvi = hi*vi
huj = hj*uj
hvj = hj*vj

!....

if (icor .eq.1) then
do j=2,endy
   do i=2,endx
      ux(i,j)  = (ui(i,j)-ui(i-1,j))/dx
      vy(i,j)  = (vj(i,j)-vj(i,j-1))/dy
      uy(i,j)  = (uj(i,j)-uj(i,j-1))/dy
      vx(i,j)  = (vi(i,j)-vi(i-1,j))/dx
   enddo
enddo

call interpolations(1,limiter2,id,z,phil,phir)
zx = (phir-phil)/dx
call interpolations(4,limiter2,id,z,phil,phir)
zy = (phir-phil)/dy

!call interpolations(2,2,id,u,phil,phir)
!ux = (phir-phil)/dx
call interpolations(3,limiter2,id,ux,phil,phir)
uxx = (phir-phil)/dx

call interpolations(2,limiter2,id,hu,phil,phir)
hux = (phir-phil)/dx
call interpolations(3,limiter2,id,hux,phil,phir)
huxx = (phir-phil)/dx

!call interpolations(6,2,id,v,phil,phir)
!vy = (phir-phil)/dy
call interpolations(5,limiter2,id,vy,phil,phir)
vyy = (phir-phil)/dy

call interpolations(6,limiter2,id,hv,phil,phir)
hvy = (phir-phil)/dy
call interpolations(5,limiter2,id,hvy,phil,phir)
hvyy = (phir-phil)/dy

!call interpolations(5,2,id,u,phil,phir)
!uy = (phir-phil)/dy
call interpolations(6,limiter2,id,uy,phil,phir)
uyy = (phir-phil)/dy

!call interpolations(3,2,id,v,phil,phir)
!vx = (phir-phil)/dx
call interpolations(2,limiter2,id,vx,phil,phir)
vxx = (phir-phil)/dx

call interpolations(2,limiter2,id,uy,phil,phir)
uxy = (phir-phil)/dx
call interpolations(2,limiter2,id,d*uy,phil,phir)
huxy = (phir-phil)/dx
call interpolations(6,limiter2,id,vx,phil,phir)
vxy = (phir-phil)/dy
call interpolations(6,limiter2,id,d*vx,phil,phir)
hvxy = (phir-phil)/dy

else

do j=2,endy
   do i=2,endx
      zx(i,j)  = (zi(i,j)-zi(i-1,j))/dx
      zy(i,j)  = (zj(i,j)-zj(i,j-1))/dy
      ux(i,j)  = (ui(i,j)-ui(i-1,j))/dx
      vy(i,j)  = (vj(i,j)-vj(i,j-1))/dy
      uy(i,j)  = (uj(i,j)-uj(i,j-1))/dy
      vx(i,j)  = (vi(i,j)-vi(i-1,j))/dx
   enddo
enddo

call interp_xei (ux,tmpi1,interp2,ix)
call interp_yei (vy,tmpj1,interp2,iy)

do j=2,endy
   do i=2,endx
      uxx(i,j)  = (tmpi1(i,j)-tmpi1(i-1,j))/dx
      vyy(i,j)  = (tmpj1(i,j)-tmpj1(i,j-1))/dy
      hux(i,j)  = (hui(i,j)-hui(i-1,j))/dx
      hvy(i,j)  = (hvj(i,j)-hvj(i,j-1))/dy
   enddo
enddo

call interp_xei (hux,tmpi1,interp2,ix)
call interp_yei (hvy,tmpj1,interp2,iy)
call interp_xoi (vx,tmpi2,interp2,ix) 
call interp_yoi (uy,tmpj2,interp2,iy)

do j=2,endy
   do i=2,endx
      huxx(i,j)  = (tmpi1(i,j)-tmpi1(i-1,j))/dx
      hvyy(i,j)  = (tmpj1(i,j)-tmpj1(i,j-1))/dy
      vxx(i,j)  = (tmpi2(i,j)-tmpi2(i-1,j))/dx
      uyy(i,j)  = (tmpj2(i,j)-tmpj2(i,j-1))/dy
   enddo
enddo

call interp_xei (vy,tmpi1,interp2,ix) 
call interp_yei (ux,tmpj1,interp2,iy)
call interp_xei (hvy,tmpi2,interp2,ix) 
call interp_yei (hux,tmpj2,interp2,iy)

do j=2,endy
   do i=2,endx
      vxy(i,j) = (tmpi1(i,j)-tmpi1(i-1,j))/dx
      uxy(i,j) = (tmpj1(i,j)-tmpj1(i,j-1))/dy
      hvxy(i,j) = (tmpi2(i,j)-tmpi2(i-1,j))/dx
      huxy(i,j) = (tmpj2(i,j)-tmpj2(i,j-1))/dy
   enddo
enddo

endif

	if(bc_1.eq.1)then
		i=overlap+1
		do j=1,endy
			hvxy(i,j)=0.
			vxy(i,j)=0.
			huxy(i,j)=0.
			uxy(i,j)=0.
			uxx(i,j)=0.
			huxx(i,j)=0.
			vx(i,j)=0.
			vxx(i,j)=0.
			uy(i,j)=0.
			uyy(i,j)=0.
			zx(i,j)=0.
	 	enddo
	endif

	if(bc_2.eq.1)then
		i=endx-overlap
		do j=1,endy
			hvxy(i,j)=0.
			vxy(i,j)=0.
			huxy(i,j)=0.
			uxy(i,j)=0.
			uxx(i,j)=0.
			huxx(i,j)=0.
			vx(i,j)=0.
			vxx(i,j)=0.
			uy(i,j)=0.
			uyy(i,j)=0.
			zx(i,j)=0.
		enddo
	endif

	if(bc_3.eq.1)then
		j=overlap+1
		do i=1,endx
			hvxy(i,j)=0.
			vxy(i,j)=0.
			huxy(i,j)=0.
			uxy(i,j)=0.
			vyy(i,j)=0.
			hvyy(i,j)=0.
			vx(i,j)=0.
			vxx(i,j)=0.
			uy(i,j)=0.
			uyy(i,j)=0.
			zy(i,j)=0.
		enddo
	endif

	if(bc_4.eq.1)then
		j=endy-overlap
		do i=1,endx
			hvxy(i,j)=0.
			vxy(i,j)=0.
			huxy(i,j)=0.
			uxy(i,j)=0.
			vyy(i,j)=0.
			hvyy(i,j)=0.
			vx(i,j)=0.
			vxx(i,j)=0.
			uy(i,j)=0.
			uyy(i,j)=0.
			zy(i,j)=0.
		enddo
	endif

! add for external source layer subroutine
	dvdyy(:,:,1)=vyy
	dudxx(:,:,1)=uxx
	dvdxy(:,:,1)=vxy
	dudxy(:,:,1)=uxy

! add for velocity profile calculations
do j=1,endy
	do i=1,endx
		S_grp(i,j,1)=ux(i,j)+vy(i,j)
		T_grp(i,j,1)=hux(i,j)+hvy(i,j) !+dhdt(i,j,n_loc)
		dSdx(i,j,1)=uxx(i,j)+vxy(i,j)
        dTdx(i,j,1)=huxx(i,j)+hvxy(i,j) !+dhdxt(i,j,n_loc)
		dSdy(i,j,1)=uxy(i,j)+vyy(i,j)
        dTdy(i,j,1)=huxy(i,j)+hvyy(i,j) !+dhdyt(i,j,n_loc)
	enddo
enddo


!.....nonlinear terms at x-dir interface
do j=2,endy-2
   do i=2,endx-2
!	     uix(i,j) = (15*(u(i+1,j)-u(i,j))-(u(i+2,j)-u(i-1,j)))/12./dx
!	     huix(i,j) = (15*(hu(i+1,j)-hu(i,j))-(hu(i+2,j)-hu(i-1,j)))/12./dx
!		 uixx(i,j) = (u(i+2,j)+u(i-1,j)-u(i+1,j)-u(i,j))/dx/dx*0.5
!		 huixx(i,j) = (hu(i+2,j)+hu(i-1,j)-hu(i+1,j)-hu(i,j))/dx/dx*0.5
!	     vix(i,j) = (15*(v(i+1,j)-v(i,j))-(v(i+2,j)-v(i-1,j)))/12./dx
!	     hvix(i,j) = (15*(hv(i+1,j)-hv(i,j))-(hv(i+2,j)-hv(i-1,j)))/12./dx
	      
	     uix(i,j) = (u(i+1,j)-u(i,j))/dx
	     huix(i,j) = (hu(i+1,j)-hu(i,j))/dx
		 uixx(i,j) = (u(i+2,j)+u(i-1,j)-u(i+1,j)-u(i,j))/dx/dx*0.5
		 huixx(i,j) = (hu(i+2,j)+hu(i-1,j)-hu(i+1,j)-hu(i,j))/dx/dx*0.5
	     vix(i,j) = (v(i+1,j)-v(i,j))/dx
	     hvix(i,j) = (hv(i+1,j)-hv(i,j))/dx
   enddo
enddo

!1
call interp_yoi (vi,cor,3,iy2)

do j=3,endy-3
   do i=3,endx-3
      viy(i,j) = (cor(i,j)-cor(i,j-1))/dy 
   enddo
enddo

!2
call interp_yoi (hvi,cor,3,iy2)

do j=3,endy-3
   do i=3,endx-3
      hviy(i,j) = (cor(i,j)-cor(i,j-1))/dy   
   enddo
enddo

!3 
do j=2,endy-2
   do i=2,endx-2
	     cor(i,j) = (vi(i,j+1)-vi(i,j))/dy
   enddo
enddo

do j=3,endy-3
   do i=3,endx-3
      viyy(i,j) = (cor(i,j)-cor(i,j-1))/dy
   enddo
enddo

!4
do j=2,endy-2
   do i=2,endx-2
	     cor(i,j) = (hvi(i,j+1)-hvi(i,j))/dy
   enddo
enddo

do j=3,endy-3
   do i=3,endx-3
      hviyy(i,j) = (cor(i,j)-cor(i,j-1))/dy
   enddo
enddo

!5
call interp_yei (uix,cor,interp2,iy2)

do j=3,endy-3
   do i=3,endx-3
      uixy(i,j) = (cor(i,j)-cor(i,j-1))/dy  
   enddo
enddo

!6
call interp_yei (huix,cor,interp2,iy2)

do j=3,endy-3
   do i=3,endx-3
      huixy(i,j) = (cor(i,j)-cor(i,j-1))/dy  
   enddo
enddo

!7
call interp_yoi (vix,cor,interp2,iy2)

do j=3,endy-3
   do i=3,endx-3 
	  vixy(i,j) = (cor(i,j)-cor(i,j-1))/dy
   enddo
enddo

!8
call interp_yoi (hvix,cor,interp2,iy2)

do j=3,endy-3
   do i=3,endx-3
      hvixy(i,j) = (cor(i,j)-cor(i,j-1))/dy   
   enddo
enddo

!.....nonlinear terms at y-dir interface

do j=2,endy-2
   do i=2,endx-2
!	     vjy(i,j) = (15*(v(i,j+1)-v(i,j))-(v(i,j+2)-v(i,j-1)))/12./dy
!	     hvjy(i,j) = (15*(hv(i,j+1)-hv(i,j))-(hv(i,j+2)-hv(i,j-1)))/12./dy
!		 vjyy(i,j) = (v(i,j+2)+v(i,j-1)-v(i,j+1)-v(i,j))/dy/dy*0.5
!		 hvjyy(i,j) = (hv(i,j+2)+hv(i,j-1)-hv(i,j+1)-hv(i,j))/dy/dy*0.5
!	     ujy(i,j) = (15*(u(i,j+1)-u(i,j))-(u(i,j+2)-u(i,j-1)))/12./dy
!	     hujy(i,j) = (15*(hu(i,j+1)-hu(i,j))-(hu(i,j+2)-hu(i,j-1)))/12./dy		
	     	  
	     vjy(i,j) = (v(i,j+1)-v(i,j))/dy
	     hvjy(i,j) = (hv(i,j+1)-hv(i,j))/dy
		 vjyy(i,j) = (v(i,j+2)+v(i,j-1)-v(i,j+1)-v(i,j))/dy/dy*0.5
		 hvjyy(i,j) = (hv(i,j+2)+hv(i,j-1)-hv(i,j+1)-hv(i,j))/dy/dy*0.5
	     ujy(i,j) = (u(i,j+1)-u(i,j))/dy
	     hujy(i,j) = (hu(i,j+1)-hu(i,j))/dy
   enddo
enddo

!1
call interp_xoi (uj,cor,3,ix2)

do j=3,endy-3
   do i=3,endx-3
      ujx(i,j) = (cor(i,j)-cor(i-1,j))/dx
   enddo
enddo

!2
call interp_xoi (huj,cor,3,ix2)

do j=3,endy-3
   do i=3,endx-3
      hujx(i,j) = (cor(i,j)-cor(i-1,j))/dx
   enddo
enddo

!3
call interp_xoi (ujy,cor,interp2,ix2)

do j=3,endy-3
   do i=3,endx-3
      ujxy(i,j) = (cor(i,j)-cor(i-1,j))/dx
   enddo
enddo

!4
call interp_xoi (hujy,cor,interp2,ix2)

do j=3,endy-3
   do i=3,endx-3
      hujxy(i,j) = (cor(i,j)-cor(i-1,j))/dx
   enddo
enddo

!5
do j=3,endy-3
   do i=3,endx-3
	     cor(i,j) = (uj(i+1,j)-uj(i,j))/dx
   enddo
enddo

do j=4,endy-4
   do i=4,endx-4
     ujxx(i,j) = (cor(i,j)-cor(i-1,j))/dx
   enddo
enddo

!6
do j=3,endy-3
   do i=3,endx-3
	     cor(i,j) = (huj(i+1,j)-huj(i,j))/dx
   enddo
enddo

do j=4,endy-4
   do i=4,endx-4
     hujxx(i,j) = (cor(i,j)-cor(i-1,j))/dx
   enddo
enddo

!7
call interp_xei (vjy,cor,interp2,ix2)

do j=3,endy-3
   do i=3,endx-3
      vjxy(i,j) = (cor(i,j)-cor(i-1,j))/dx
   enddo
enddo

!8
call interp_xei (hvjy,cor,interp2,ix2)

do j=3,endy-3
   do i=3,endx-3
      hvjxy(i,j) = (cor(i,j)-cor(i-1,j))/dx
   enddo
enddo

do j=4,endy-4
   do i=4,endx-4
	  hvyzx(i,j) = (hviy(i,j)*zi(i,j)-hviy(i-1,j)*zi(i-1,j))/dx
	  huxzy(i,j) = (hujx(i,j)*zj(i,j)-hujx(i,j-1)*zj(i,j-1))/dy
   enddo
enddo


!...............................................................................

if(ibotfric.eq.1 .or. idissip.eq.1 )then
do j=2,endy-2
   do i=2,endx-2
		 Hs_i=max(cutoff_mat(i,j),H(i,j))


		 if(ish_mat(i,j).eq.0)then
			ua(i,j) = u(i,j)  - ((z(i,j)**2-z(i,j)&
*d(i,j)+d(i,j)**2)/6.0-0.5*za(i,j)**2)*(dSdx(i,j,1)) &
                 - (0.5*(z(i,j)-d(i,j))-za(i,j))*(dTdx(i,j,1)) ! &
                 ! + psix(i,j,n_loc)*(za(i,j)**2*0.5-za(i,j)*z(i,j)+(2*z(i,j)**2-2*z(i,j)*d(i,j)-d(i,j)**2)/6.0) ! depth-averaged velocity 
			
			va(i,j) = v(i,j)  - ((z(i,j)**2-z(i,j)&
*d(i,j)+d(i,j)**2)/6.0-0.5*za(i,j)**2)*(dSdy(i,j,1)) &
		       - (0.5*(z(i,j)-d(i,j))-za(i,j))*(dTdy(i,j,1)) ! &
                 ! + psiy(i,j,n_loc)*(za(i,j)**2*0.5-za(i,j)*z(i,j)+(2*z(i,j)**2-2*z(i,j)*d(i,j)-d(i,j)**2)/6.0) ! depth-averaged velocity 
 
         else
			ua(i,j) = u(i,j) 
			va(i,j) = v(i,j) 
		 endif
                    
         if(bf_type.eq.0)then	
			Rh = max(100.0,sqrt(va(i,j)**2+ua(i,j)**2)*Hs_i/dissipcoef)
			Cf_c = -1.8*log10(6.9/Rh+(ks/Hs_i/3.7)**1.11)  ! Moody approx
			Cf(i,j) = 0.25/Cf_c**2
			if (Rh.le.411.58) Cf(i,j) = 0.25*64./Rh
		 elseif(bf_type.eq.1)then
			if(d(i,j).le.0.)then
				ks_i=0.085 ! n for rocky bottom (breakwater)
			else
				ks_i=0.085*(1.-conv_reduc(i))+ks*conv_reduc(i)
			endif
!			ks_i=ks		 
		    Cf(i,j)=ks_i**2.*9.81/Hs_i**0.3333  ! Mannings
		 elseif(bf_type.eq.2)then		 
			Cf(i,j)=ks  ! constant friction factor
         endif

		 if(id(i,j).ne.0.or.Hs_i.lt.cutoff_mat_d(i,j))then
			if(ua(i,j)*hx(i,j).ge.0.or.va(i,j)*hy(i,j).ge.0)then
				Cf(i,j) = max(Cf(i,j),0.04)  ! Puleo & Holland downrush friction factor
			else
				Cf(i,j) = max(Cf(i,j),0.01)  ! Puleo & Holland uprush friction factor
			endif
		 endif
		          
         Rhx(i,j) = Cf(i,j)*ua(i,j)*sqrt(va(i,j)**2+ua(i,j)**2)
         Rhy(i,j) = Cf(i,j)*va(i,j)*sqrt(va(i,j)**2+ua(i,j)**2)

      enddo
   enddo

endif   


!--------- add dissipation --------------
if (idissip.eq.1) then
	do j=2,endy-2
	   do i=2,endx-2
		    Hs_i=max(cutoff_mat(i,j),H(i,j))
            vtv(i,j) = Ch*Hs_i*sqrt(Cf(i,j)*(ua(i,j)**2+va(i,j)**2)) + dissipcoef  ! vertical eddy visc

			cm_i=cm
			if(id(i,j).ne.0.or.Hs_i.lt.cutoff_mat_d(i,j)) cm_i=1.0
			vt(i,j) = cm_i*dx*dy*sqrt(4.*ux(i,j)**2 + 4.*vy(i,j)**2 + 4.*(ux(i,j)*vy(i,j)) + & 
                       (uy(i,j)+vx(i,j))**2 ) + dissipcoef  ! horizontal eddy visc
      enddo
   enddo


   call interp_xei (vt,vti,interp2,ix)
   call interp_yei (vt,vtj,interp2,iy)

endif
if (ihvor.eq.1) then
	do j=2,endy-2
	   do i=2,endx-2
		    if(H(i,j).le.Ch_length)then
				psix(i,j,n_loc) = 0.
				psiy(i,j,n_loc) = 0.
			else
				psix(i,j,n_loc) = Rhx(i,j)/(vtv(i,j)*H(i,j))
				psiy(i,j,n_loc) = Rhy(i,j)/(vtv(i,j)*H(i,j))
			endif
      enddo
   enddo
else
	do j=2,endy-2
	   do i=2,endx-2
				psix(i,j,n_loc) = Rhx(i,j)/(vtv(i,j)*H(i,j))
				psiy(i,j,n_loc) = Rhy(i,j)/(vtv(i,j)*H(i,j))
      enddo
   enddo   
endif

! Internal source (type 1) wave generation
	do j = sy,ey
		do i = sx,ex
			  E(i,j)=0.
			  F(i,j)=0.
			  G(i,j)=0.
		enddo
	enddo

	if(int_src.eq.1.or.plus_tide.eq.1)then  !^^ Internal source type 1
	 if(n_loc.eq.4.and.count_min.eq.1)then !&&
	  do j = sy,ey
	    do i = sx,ex

		  f_int_src(i,j)=0.
		  f_int_src_v(i,j)=0.
		  
          if(bl_hor_wall(i,j).le.50)then !%%
          
!%%%%%%%%%%%%%%%%%  INTERNAL SOURCE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          tmp=0.

		  tmp=min(1.,max(0.,d(i,j)/depth))**2.

          do theta=1,num_theta
			do freq=1,num_freq
      
                 if(nn_loc.le.4)then !@@
                  cosA(freq,theta)=cos(inc_ang_spec(freq,theta))
                  sinA(freq,theta)=sin(inc_ang_spec(freq,theta))
                  wA(freq,theta)=2.*3.1415/per_spec(freq,theta) 
                  kA(freq,theta)=wA(freq,theta)/tanh(depth*2.*3.1415/L_spec(freq,theta))  ! this IS NOT WAVENUMBER, used by velocity source          				
                 endif !@@
            
				 if(is_oreint.eq.1)then
					x_c=x(i)-end_x_t/2.
					y_c=y(j)-x0-L_spec(freq,theta)*real(theta/num_theta)
				 else
					x_c=x(i)-x0-L_spec(freq,theta)*real(theta/num_theta) 
					y_c=y(j)-end_y_t/2.-2.*L_spec(freq,theta)*real(freq/num_freq)
				 endif
				
				 loc=abs(((x_c-L_spec(freq,theta)*& 
                         sinA(freq,theta))*cosA(freq,theta)+& 
                       (y_c+L_spec(freq,theta)*& 
                         sinA(freq,theta))*sinA(freq,theta)))               
			     if(loc/L_spec(freq,theta).le.1.0)then !$$			   
				   is_coef=0.5*D_src(freq,theta)*exp(-beta_src(freq,theta)*loc**2)	
				   is_sin=sin(-wA(freq,theta)*(t(nn_loc)+shift_spec(freq,theta))) 		
				   src_amp=is_coef*is_sin*kA(freq,theta)			   
				   f_int_src(i,j)=f_int_src(i,j)+src_amp*cosA(freq,theta)                    
				   f_int_src_v(i,j)=f_int_src_v(i,j)+src_amp*sinA(freq,theta)  
			     endif !$$
			   enddo
            enddo
            
            if(plus_tide.eq.1)then
					loc=x_c
				    is_coef=0.5*D_tide*exp(-beta_tide*loc**2)	
				    is_sin=eta_in(nn_loc)
				    src_amp=is_coef*is_sin				
					f_int_src(i,j)=f_int_src(i,j)+src_amp*sqrt(9.81/depth)  
			endif          
		endif !%%
		
       enddo
     enddo
     
     endif !&&
     
     
	 do j = sy-1,ey
		do i = sx-1,ex
				   F(i,j) = F(i,j)+f_int_src(i,j)*H(i,j)
				   G(i,j) = G(i,j)+f_int_src_v(i,j)*H(i,j)
        enddo
     enddo
	endif !^^


!-------------start calculating E, F, and G
	
  !----------calculate E------------
   do j=sy,ey
      do i=sx-1,ex
         temp2 = 1./6.*(zi(i,j)**3+hi(i,j)**3)-0.5*zai(i,j)**2*(zi(i,j)+hi(i,j))
	     temp3 = 0.5*(zi(i,j)**2-hi(i,j)**2)-zai(i,j)*(zi(i,j)+hi(i,j))
	     E1(i,j) = temp2*(uixx(i,j)+vixy(i,j))+temp3*(huixx(i,j)+hvixy(i,j))
      enddo
   enddo

   do j=sy-1,ey
      do i=sx,ex
         temp2 = 1./6.*(zj(i,j)**3+hj(i,j)**3)-0.5*zaj(i,j)**2*(zj(i,j)+hj(i,j))
	     temp3 = 0.5*(zj(i,j)**2-hj(i,j)**2)-zaj(i,j)*(zj(i,j)+hj(i,j))
	     E2(i,j) = temp2*(ujxy(i,j)+vjyy(i,j))+temp3*(hujxy(i,j)+hvjyy(i,j))
      enddo
   enddo

   do j=sy,ey
      do i=sx,ex
	     E(i,j) = E(i,j)+(E1(i,j)-E1(i-1,j))/dx + (E2(i,j)-E2(i,j-1))/dy &
-(flux(i,j,1)-flux(i-1,j,1))/dx-(fluy(i,j,1)-fluy(i,j-1,1))/dy  !- dhdt(i,j,n_loc)
      enddo
   enddo

   if (ihvor.eq.1) then
      call interp_xoi (psix(:,:,n_loc),tmpi1,interp2,ix)
      call interp_yoi (psiy(:,:,n_loc),tmpj1,interp2,iy)
	  tmpi2 = tmpi1*(zi+hi)*(zai*zai*0.5-zai*zi+(2*zi*zi-2*zi*hi-hi*hi)/6.0)
	  tmpj2 = tmpj1*(zj+hj)*(zaj*zaj*0.5-zaj*zj+(2*zj*zj-2*zj*hj-hj*hj)/6.0)
	  do i = sx,ex
	     do j = sy,ey
		    E(i,j) = E(i,j) - (tmpi2(i,j)-tmpi2(i-1,j))/dx &
- (tmpj2(i,j)-tmpj2(i,j-1))/dy
		 enddo
      enddo
   endif


   !----------caclulate F,F1------------

   call interp_xei (E,Ei,interp2,ix)

   do j=sy,ey
      do i=sx-1,ex
         EzSTi(i,j) = Ei(i,j)*(zi(i,j)*(uix(i,j)+viy(i,j))+(huix(i,j)+hviy(i,j)))
	     TzS2i(i,j) = (zi(i,j)*(uix(i,j)+viy(i,j)))**2 &
                    + 2.0*zi(i,j)*(uix(i,j)+viy(i,j))*(huix(i,j)+hviy(i,j)) &
                    + (huix(i,j)+hviy(i,j))**2
      enddo
   enddo
   do j=sy,ey
      do i=sx-1,ex
	     uSxvSyi(i,j) = 0.5*(zai(i,j)**2-zi(i,j)**2)*&
(ui(i,j)*(uixx(i,j)+vixy(i,j))+vi(i,j)*(uixy(i,j)+viyy(i,j))) 
	     uTxvTyi(i,j) = (zai(i,j)-zi(i,j))*(ui(i,j)*&
(huixx(i,j)+hvixy(i,j))+vi(i,j)*(huixy(i,j)+hviyy(i,j)))
      enddo
   enddo	   		   

   do j=sy,ey
      do i=sx,ex
		 temp1 = -(flux(i,j,2)-flux(i-1,j,2))/dx&
-(fluy(i,j,2)-fluy(i,j-1,2))/dy + u(i,j)*((E1(i,j)-E1(i-1,j))/dx + (E2(i,j)-E2(i,j-1))/dy)
	     temp3 = -(EzSTi(i,j)-EzSTi(i-1,j))/dx*H(i,j) &
                    + E(i,j)*(0.5*(za(i,j)**2-z(i,j)**2)*&
dSdx(i,j,1)+(za(i,j)-z(i,j))*dTdx(i,j,1) &
                              -zx(i,j)*(z(i,j)*S_grp(i,j,1)+T_grp(i,j,1)))
         temp4 = -(uSxvSyi(i,j)-uSxvSyi(i-1,j))/dx*H(i,j)
	     temp5 = -(uTxvTyi(i,j)-uTxvTyi(i-1,j))/dx*H(i,j)
	     temp6 = -0.5*(TzS2i(i,j)-TzS2i(i-1,j))/dx*H(i,j)
         temp7 = -v(i,j)*(zai(i,j)-zai(i-1,j))/dx*(huxy(i,j)+hvyy(i,j)+za(i,j)*(uxy(i,j)+vyy(i,j))) &
			 +v(i,j)*(zaj(i,j)-zaj(i,j-1))/dy*&
(huxx(i,j)+hvxy(i,j)+za(i,j)*(uxx(i,j)+vxy(i,j))) &
				 -((vi(i,j)-vi(i-1,j))/dx-(uj(i,j)-uj(i,j-1))&
/dy) & 
				 *((za(i,j)-0.5*(z(i,j)-d(i,j)))*(huxy(i,j)+hvyy(i,j)) & 
				  +(0.5*za(i,j)**2-(z(i,j)**2-z(i,j)*&
d(i,j)+d(i,j)**2)/6.0)*(uxy(i,j)+vyy(i,j)))

		 F(i,j) = F(i,j) + (temp1+temp3+temp4+temp5+temp6+temp7*&
H(i,j))-grav*H(i,j)*zx_4th(i,j) !+grav*Hcx(i,j)*(d(i+1,j)-d(i-1,j))/dx*0.5
	     
		 F1(i,j) = 0.5*(z(i,j)**2-za(i,j)**2)*vxy(i,j)-za(i,j)&
*hvxy(i,j)+zx(i,j)*z(i,j)*vy(i,j)+hvyzx(i,j)
		 F1(i,j) = F1(i,j)*H(i,j)
	  enddo
   enddo


   !----------caculate G,G1-------------

   call interp_yei (E,Ej,interp2,iy)


   do j=sy-1,ey
     do i=sx,ex
         EzSTj(i,j) = Ej(i,j)*(zj(i,j)*(ujx(i,j)+vjy(i,j))+(hujx(i,j)+hvjy(i,j)))
	     TzS2j(i,j) = (zj(i,j)*(ujx(i,j)+vjy(i,j)))**2 &
                    + 2.0*zj(i,j)*(ujx(i,j)+vjy(i,j))*(hujx(i,j)+hvjy(i,j)) &
                    + (hujx(i,j)+hvjy(i,j))**2
      enddo
   enddo	
   
   do j=sy-1,ey   		   
     do i=sx,ex
	     uSxvSyj(i,j) = 0.5*(zaj(i,j)**2-zj(i,j)**2)*(uj(i,j)*&
(ujxx(i,j)+vjxy(i,j))+vj(i,j)*(ujxy(i,j)+vjyy(i,j))) 
	     uTxvTyj(i,j) = (zaj(i,j)-zj(i,j))*(uj(i,j)*&
(hujxx(i,j)+hvjxy(i,j))+vj(i,j)*(hujxy(i,j)+hvjyy(i,j)))
      enddo
   enddo	   		   

   do j=sy,ey
      do i=sx,ex
         temp1 = -(flux(i,j,3)-flux(i-1,j,3))/dx-&
(fluy(i,j,3)-fluy(i,j-1,3))/dy + v(i,j)*((E1(i,j)-E1(i-1,j))/dx + &
(E2(i,j)-E2(i,j-1))/dy)
	     temp3 = -(EzSTj(i,j)-EzSTj(i,j-1))/dy*H(i,j) &
                    + E(i,j)*(0.5*(za(i,j)**2-z(i,j)**2)*&
dSdy(i,j,1)+(za(i,j)-z(i,j))*dTdy(i,j,1) &
                              -zy(i,j)*(z(i,j)*S_grp(i,j,1)+T_grp(i,j,1)))
         temp4 = -(uSxvSyj(i,j)-uSxvSyj(i,j-1))/dy*H(i,j)
	     temp5 = -(uTxvTyj(i,j)-uTxvTyj(i,j-1))/dy*H(i,j)
	     temp6 = -0.5*(TzS2j(i,j)-TzS2j(i,j-1))/dy*H(i,j)
         temp7 =  u(i,j)*(zai(i,j)-zai(i-1,j))/dx*(huxy(i,j)+hvyy(i,j)+za(i,j)*(uxy(i,j)+vyy(i,j))) &
				 -u(i,j)*(zaj(i,j)-zaj(i,j-1))/dy*&
(huxx(i,j)+hvxy(i,j)+za(i,j)*(uxx(i,j)+vxy(i,j))) &
				 +((vi(i,j)-vi(i-1,j))/dx-(uj(i,j)-uj(i,j-1))/&
dy) & 
				 *((za(i,j)-0.5*(z(i,j)-d(i,j)))&
*(huxx(i,j)+hvxy(i,j))+(0.5*za(i,j)**2-(z(i,j)**2&
-z(i,j)*d(i,j)+d(i,j)**2)/6.0)*(uxx(i,j)+vxy(i,j)))
	  
		 G(i,j) = G(i,j) + (temp1+temp3+temp4+temp5+temp6+temp7*&
H(i,j))-grav*H(i,j)*zy_4th(i,j) !+ grav*Hcy(i,j)*(d(i,j+1)-d(i,j-1))/dy*0.5

	  	 G1(i,j) = 0.5*(z(i,j)**2-za(i,j)**2)*&
uxy(i,j)-za(i,j)*huxy(i,j)+zy(i,j)*z(i,j)*ux(i,j)+huxzy(i,j)
		 G1(i,j) = G1(i,j)*H(i,j)

		 dhzudx(i,j,1)=(flux(i,j,1)-flux(i-1,j,1))/dx
		 dhzudy(i,j,1)=(flux(i,j,1)-flux(i,j-1,1))/dy
		 dhzvdx(i,j,1)=(fluy(i,j,1)-fluy(i-1,j,1))/dx
		 dhzvdy(i,j,1)=(fluy(i,j,1)-fluy(i,j-1,1))/dy
      enddo
   enddo


! Horizontal vorticity correction terms
if(ihvor.eq.1)then

   call interp_xoi (psix(:,:,n_loc),tmpi1,interp2,ix) 
   call interp_xei (psiy(:,:,n_loc),tmpi2,interp2,ix)

   call interp_yei (psix(:,:,n_loc),tmpj1,interp2,iy)
   call interp_yoi (psiy(:,:,n_loc),tmpj2,interp2,iy)

   do j = sy,ey
      do i = sx,ex
	     tmp1(i,j) = (ui(i,j)*tmpi1(i,j)*zi(i,j)-ui(i-1,j)*tmpi1(i-1,j)*zi(i-1,j))/dx &
		           + (vi(i,j)*tmpi2(i,j)*zi(i,j)-vi(i-1,j)*tmpi2(i-1,j)*zi(i-1,j))/dx
	     tmp2(i,j) = (uj(i,j)*tmpj1(i,j)*zj(i,j)-uj(i,j-1)*tmpj1(i,j-1)*zj(i,j-1))/dy &
		           + (vj(i,j)*tmpj2(i,j)*zj(i,j)-vj(i,j-1)*tmpj2(i,j-1)*zj(i,j-1))/dy
      enddo
   enddo

   F = F - tmp1*(z-d)*0.5*H
   G = G - tmp2*(z-d)*0.5*H

   do j = sy,ey
      do i = sx,ex
	     tmp1(i,j) = (ui(i,j)*tmpi1(i,j)-ui(i-1,j)*tmpi1(i-1,j))/dx &
		           + (vi(i,j)*tmpi2(i,j)-vi(i-1,j)*tmpi2(i-1,j))/dx
	     tmp2(i,j) = (uj(i,j)*tmpj1(i,j)-uj(i,j-1)*tmpj1(i,j-1))/dy &
		           + (vj(i,j)*tmpj2(i,j)-vj(i,j-1)*tmpj2(i,j-1))/dy
      enddo
   enddo

   F = F + tmp1*(z*z-z*d+d*d)/6.0*H
   G = G + tmp2*(z*z-z*d+d*d)/6.0*H

   do j = sy,ey
      do i = sx,ex
	     tmp1(i,j) = (ui(i,j)*tmpi1(i,j)*(0.5*zai(i,j)*zai(i,j)-zi(i,j)*zai(i,j)) &
    -ui(i-1,j)*tmpi1(i-1,j)*(0.5*zai(i-1,j)*zai(i-1,j)-zi(i-1,j)*zai(i-1,j)))/dx &
		           + (vi(i,j)*tmpi2(i,j)*(0.5*zai(i,j)*zai(i,j)-zi(i,j)*zai(i,j)) &
 -vi(i-1,j)*tmpi2(i-1,j)*(0.5*zai(i-1,j)*zai(i-1,j)-zi(i-1,j)*zai(i-1,j)))/dx
	     tmp2(i,j) = (uj(i,j)*tmpj1(i,j)*(0.5*zaj(i,j)*zaj(i,j)-zj(i,j)*zaj(i,j)) &
		             -uj(i,j-1)*tmpj1(i,j-1)*(0.5*zaj(i,j-1)*zaj(i,j-1)-zj(i,j-1)*zaj(i,j-1)))/dy &
		           + (vj(i,j)*tmpj2(i,j)*(0.5*zaj(i,j)*zaj(i,j)-zj(i,j)*zaj(i,j)) &
	     -vj(i,j-1)*tmpj2(i,j-1)*(0.5*zaj(i,j-1)*zaj(i,j-1)-zj(i,j-1)*zaj(i,j-1)))/dy
      enddo
   enddo

   F = F - tmp1*H
   G = G - tmp2*H

!... w*W

   tmp1 = psix(:,:,n_loc)*((z*z+z*d-2*d*d)*(ux+vy)/6.0+0.5*H*(hux+hvy))
   tmp2 = psiy(:,:,n_loc)*((z*z+z*d-2*d*d)*(ux+vy)/6.0+0.5*H*(hux+hvy))

   F = F + tmp1*H
   G = G + tmp2*H

!... hori. in KXI term

   tmp1 = -(vx-uy)*psiy(:,:,n_loc)*(0.5*za*za-za*z+(2*z*z-2*z*d-d*d)/6.0)
   tmp2 =  (vx-uy)*psix(:,:,n_loc)*(0.5*za*za-za*z+(2*z*z-2*z*d-d*d)/6.0)

   F = F - tmp1*H
   G = G - tmp2*H

   do j = sy,ey
      do i = sx,ex
	     tmp1(i,j) = (tmpi2(i,j)-tmpi2(i-1,j))/dx*&
(z(i,j)*z(i,j)-z(i,j)*d(i,j)+d(i,j)*d(i,j))/6.0
         tmp2(i,j) = (tmpj1(i,j)-tmpj1(i,j-1))/dy*&
(z(i,j)*z(i,j)-z(i,j)*d(i,j)+d(i,j)*d(i,j))/6.0
      enddo
   enddo

   F = F - v*tmp1*H + v*tmp2*H
   G = G + u*tmp1*H - u*tmp2*H

   do j = sy,ey
      do i = sx,ex
	     tmp1(i,j) = (tmpi2(i,j)*zi(i,j)-tmpi2(i-1,j)*zi(i-1,j))/dx*(z(i,j)-d(i,j))*0.5
         tmp2(i,j) = (tmpj1(i,j)*zj(i,j)-tmpj1(i,j-1)*zj(i,j-1))/dy*(z(i,j)-d(i,j))*0.5
      enddo
   enddo

   F = F + v*tmp1*H - v*tmp2*H
   G = G - u*tmp1*H + u*tmp2*H

   do j = sy,ey
      do i = sx,ex
	     tmp1(i,j) = (tmpi2(i,j)*(0.5*zai(i,j)*zai(i,j)-zai(i,j)*zi(i,j)) &
		           -  tmpi2(i-1,j)*(0.5*zai(i-1,j)*zai(i-1,j)-zai(i-1,j)*zi(i-1,j)))/dx
	     tmp2(i,j) = (tmpj1(i,j)*(0.5*zaj(i,j)*zaj(i,j)-zaj(i,j)*zj(i,j)) &
		           -  tmpj1(i,j-1)*(0.5*zaj(i,j-1)*zaj(i,j-1)-zaj(i,j-1)*zj(i,j-1)))/dy
      enddo
   enddo

   F = F + v*tmp1*H - v*tmp2*H
   G = G - u*tmp1*H + u*tmp2*H

endif

! shoreline reduction to NLSW equations
	do j=sy,ey
		do i=sx,ex
			if(ish_mat(i,j).eq.1)then
				E(i,j)= -(flux(i,j,1)-flux(i-1,j,1))/dx&
-(fluy(i,j,1)-fluy(i,j-1,1))/dy !- dhdt(i,j,n_loc)
				F(i,j)= -(flux(i,j,2)-flux(i-1,j,2))/dx&
-(fluy(i,j,2)-fluy(i,j-1,2))/dy -grav*H(i,j)*zx_4th(i,j) !+grav*Hcx(i,j)*(d(i+1,j)-d(i-1,j))/dx*0.5
				F1(i,j)=0.
				G(i,j)= -(flux(i,j,3)-flux(i-1,j,3))/dx&
-(fluy(i,j,3)-fluy(i,j-1,3))/dy -grav*H(i,j)*zy_4th(i,j) !+grav*Hcy(i,j)*(d(i,j+1)-d(i,j-1))/dy*0.5
				G1(i,j)=0.
			endif
		enddo
	enddo

! add bottom stress	
	if(ibotfric.eq.1)then
		do j=sy,ey
			do i=sx,ex
		      F(i,j)=F(i,j)-Rhx(i,j)
		      G(i,j)=G(i,j)-Rhy(i,j)
			enddo
		enddo
	endif

! Dissipation terms
if (idissip.eq.1) then
   
   call interp_yoi (vx,tmpj1,interp2,iy)
   call interp_xoi (uy,tmpi1,interp2,ix) 
   
   
   do j=sy,ey      
      do i=sx,ex
         temp1 = (vti(i,j)*uix(i,j)-vti(i-1,j)*uix(i-1,j))/dx
         temp2 = (vtj(i,j)*ujy(i,j)-vtj(i,j-1)*ujy(i,j-1))/dy
         temp3 = vtv(i,j)*(uix(i,j)+viy(i,j) - uix(i-1,j)-viy(i-1,j))/dx
         temp4 = (vtj(i,j)*tmpj1(i,j)-vtj(i,j-1)*tmpj1(i,j-1))/dy
         F(i,j) = F(i,j) + (2*temp1+temp2-1*temp3+temp4)*H(i,j)
      enddo
   enddo

   do j=sy,ey
      do i=sx,ex
         temp1 = (vti(i,j)*vix(i,j)-vti(i-1,j)*vix(i-1,j))/dx
         temp2 = (vtj(i,j)*vjy(i,j)-vtj(i,j-1)*vjy(i,j-1))/dy
         temp3 = vtv(i,j)*(ujx(i,j)+vjy(i,j) - ujx(i,j-1)-vjy(i,j-1))/dy
         temp4 = (vti(i,j)*tmpi1(i,j)-vti(i-1,j)*tmpi1(i-1,j))/dx
         G(i,j) = G(i,j) + (2*temp1+temp2-1*temp3+temp4)*H(i,j)
      enddo
   enddo

endif  ! for dissipation terms


! Backscatter model
	if(backscatter.eq.1)then
		if(count_min.eq.1) call create_rand_matrix
		do j=sy,ey
			do i=sx,ex
			  BS = CB*sqrt(ua(i,j)**2+va(i,j)**2)*&
sqrt(dissipcoef*sqrt(Cf(i,j))/dt)
!		      Hs_i=max(cutoff_mat(i,j),H(i,j))
!			  BS_breaker = CB*nut(i,j,n_loc)/Hs_i*sqrt(dissipcoef/dt)  ! need to calibrate leading coef, can't assume it to be the same as for bottom turbulence
		      F(i,j)=F(i,j)+BS*rand(i,j,1) !+ BS_breaker*rand(i,j,1)
		      G(i,j)=G(i,j)+BS*rand(i,j,2) !+ BS_breaker*rand(i,j,2)
       enddo
     enddo	
	endif	

	     
return
end
