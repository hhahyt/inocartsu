C*******************************************************************************
C........    Groups the physical variables into the 
C........... conservative forms given in Wei & Kirby
      subroutine FV_var_grp(n_loc,nn_loc)
      use mainvar_module, only:zeta,h,u,v,E,F,F1,G,G1,t,H_total,
     -  i,j,endx,endy,overlap,UU,VV,cur_level,nx,ny,disp_prop,cutoff,
     -  dx,dy,dt,dim,dims,v_tmp,bl_hor_wall,level,xrank,xcomm,yrank,
     -  ycomm,cutoff_mat,breaker_type,psix,psiy,ihvor,z_alp,wvmk_loc_i,
     -  wave_breaking,x,y,swash_d
	use FV_var_module

      integer n_loc,nn_loc,sxT,exT,syT,eyT
	  real Fp_v, Fc_v, Gp_v, Gc_v, za, zc, coef1,coef2,term1,term2,term3,conv_u,pot_u,shore_cut
	cur_level=1
      nn_loc=nn_loc

	if(n_loc.le.4)then 
		call FV_ix_iy_eval(n_loc,nn_loc) !update near boundary booleans
	endif

	z_FV=zeta(:,:,n_loc,1)
	d_FV=h(:,:,n_loc)
	u_FV=u(:,:,n_loc,1)
	v_FV=v(:,:,n_loc,1)
	E_FV=E(:,:,n_loc,1)
	F_FV=F(:,:,n_loc,1)
	F1_FV=F1(:,:,n_loc,1)
	G_FV=G(:,:,n_loc,1)
	G1_FV=G1(:,:,n_loc,1)
	
c	if(n_loc.eq.4.or.nn_loc.le.4)then
	  if(dim.eq.1)then
		call FV_efg_calc_1D(z_FV,d_FV,u_FV,E_FV,F_FV,F1_FV,n_loc,nn_loc)
	  else
		call FV_efg_calc(z_FV,d_FV,u_FV,v_FV,E_FV,F_FV,F1_FV,
     -		G_FV,G1_FV,n_loc,nn_loc)
	  endif
c	endif
c	zeta(:,:,n_loc,1)=z_FV
c	h(:,:,n_loc)=d_FV
c	u(:,:,n_loc,1)=u_FV
c	v(:,:,n_loc,1)=v_FV
	E(:,:,n_loc,1)=E_FV
	F(:,:,n_loc,1)=F_FV
	F1(:,:,n_loc,1)=F1_FV
	G(:,:,n_loc,1)=G_FV
	G1(:,:,n_loc,1)=G1_FV


C%%%%%%%%%%%%%%%%%%%  SPONGE LAYER  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	call add_sponge(n_loc,nn_loc)
c	call add_bottom_friction(n_loc,nn_loc)  ! bottom friction included in FV_efg_calc

	if(wave_breaking.eq.1)then
	 if(breaker_type.eq.1)then
		call add_breaking(n_loc,nn_loc)
	 elseif(breaker_type.eq.0)then
		call add_breaking_ktrans(n_loc,nn_loc)
	 endif
	endif
	
      if(n_loc.ge.3)then  ! good for all # of levels

       do j=1+overlap,endy-overlap
        do i=wvmk_loc_i(j),endx-overlap

         if(ish_mat(i,j).ge.0)then
          if(n_loc.eq.3)then
C*************** Predictor equations for zeta, u, and v **************

            UU(i,j,n_loc+1,cur_level)=
     -            UU(i,j,n_loc,cur_level)+dt/12.*
     -           (23.*(F(i,j,n_loc,cur_level))-
     -           16.*(F(i,j,n_loc-1,cur_level))+
     -           5.*(F(i,j,n_loc-2,cur_level)))+
     -            (2.*F1(i,j,n_loc,cur_level)-
     -            3.*F1(i,j,n_loc-1,cur_level)+
     -              F1(i,j,n_loc-2,cur_level))

			if(ihvor.eq.1)then
			 za=z_alp(i,j,3,cur_level)
			 zc=zeta(i,j,3,cur_level)
			 H_total=zc+d_FV(i,j)
			 
			 coef1=H_total*(zc**2.-zc*d_FV(i,j)
     -					+d_FV(i,j)**2.)/6.
			 coef2=-H_total*(zc-d_FV(i,j)-2.*za)/2.   
			 
			 term1=z_alp(i,j,3,cur_level)**2./2.
     -     -zeta(i,j,3,cur_level)*z_alp(i,j,3,cur_level)
     
			 term2=z_alp(i,j,2,cur_level)**2./2.
     -     -zeta(i,j,2,cur_level)*z_alp(i,j,2,cur_level)
     
			 term3=z_alp(i,j,1,cur_level)**2./2.
     -     -zeta(i,j,1,cur_level)*z_alp(i,j,1,cur_level)
			 
			 Fp_v=coef1*  
     -              (2.*psix(i,j,3)  
     -              -3.*psix(i,j,2)  
     -              +1.*psix(i,j,1)) 
     -           +coef2*  
     -              (2.*psix(i,j,3)*zeta(i,j,3,cur_level)  
     -              -3.*psix(i,j,2)*zeta(i,j,2,cur_level)  
     -              -1.*psix(i,j,1)*zeta(i,j,1,cur_level)) 
     -          -H_total* 
     -              (2.*psix(i,j,3)*term1  
     -              -3.*psix(i,j,2)*term2  
     -              -1.*psix(i,j,1)*term3)
						 
			 UU(i,j,n_loc+1,cur_level)=
     -         UU(i,j,n_loc+1,cur_level) +Fp_v     
			endif
			
				
            if(dim.eq.2)then      ! only calc for 2D problem
              VV(i,j,n_loc+1,cur_level)=
     -              VV(i,j,n_loc,cur_level)+dt/12.*
     -              (23.*(G(i,j,n_loc,cur_level))-
     -              16.*(G(i,j,n_loc-1,cur_level))+
     -              5.*(G(i,j,n_loc-2,cur_level)))+
     -              (2.*G1(i,j,n_loc,cur_level)-
     -              3.*G1(i,j,n_loc-1,cur_level)+
     -                G1(i,j,n_loc-2,cur_level))
     
			  if(ihvor.eq.1)then
			 
				Gp_v=coef1*  
     -              (2.*psiy(i,j,3)  
     -              -3.*psiy(i,j,2)  
     -              +1.*psiy(i,j,1)) 
     -           +coef2*  
     -              (2.*psiy(i,j,3)*zeta(i,j,3,cur_level)  
     -              -3.*psiy(i,j,2)*zeta(i,j,2,cur_level)  
     -              -1.*psiy(i,j,1)*zeta(i,j,1,cur_level)) 
     -          -H_total* 
     -              (2.*psiy(i,j,3)*term1  
     -              -3.*psiy(i,j,2)*term2  
     -              -1.*psiy(i,j,1)*term3)
						 
				VV(i,j,n_loc+1,cur_level)=
     -				VV(i,j,n_loc+1,cur_level) +Gp_v     
     
				endif     
            endif    
		    
			
			if(level(i,j).ne.0)then
				UU(i,j,n_loc+1,cur_level)=0.
				VV(i,j,n_loc+1,cur_level)=0.
				u(i,j,n_loc+1,cur_level)=0.
				v(i,j,n_loc+1,cur_level)=0.
			endif     

          elseif(n_loc.eq.4)then
C*************** Corrector equations for zeta, u, and v ***************

            UU(i,j,n_loc,cur_level)=
     -            UU(i,j,n_loc-1,cur_level)+dt/24.*
     -            (9.*(F(i,j,n_loc,cur_level))+
     -            19.*(F(i,j,n_loc-1,cur_level))-
     -            5.*(F(i,j,n_loc-2,cur_level))+
     -            (F(i,j,n_loc-3,cur_level)))+
     -            (F1(i,j,n_loc,cur_level)-
     -            F1(i,j,n_loc-1,cur_level))

			if(ihvor.eq.1)then
			 za=z_alp(i,j,4,cur_level)
			 zc=zeta(i,j,4,cur_level)
			 H_total=zc+d_FV(i,j)
			 
			 coef1=H_total*(zc**2.-zc*d_FV(i,j)
     -					+d_FV(i,j)**2.)/6.
			 coef2=-H_total*(zc-d_FV(i,j)-2.*za)/2.   
			 
			 term1=z_alp(i,j,4,cur_level)**2./2.
     -     -zeta(i,j,4,cur_level)*z_alp(i,j,4,cur_level)
     
			 term2=z_alp(i,j,3,cur_level)**2./2.
     -     -zeta(i,j,3,cur_level)*z_alp(i,j,3,cur_level)
			 
			 Fc_v=coef1*  
     -              (psix(i,j,4)  
     -              -psix(i,j,3)) 
     -           +coef2*  
     -              (psix(i,j,4)*zeta(i,j,4,cur_level)  
     -              -psix(i,j,3)*zeta(i,j,3,cur_level)) 
     -          -H_total* 
     -              (psix(i,j,4)*term1  
     -              -psix(i,j,3)*term2)
						 
			 UU(i,j,n_loc,cur_level)=
     -         UU(i,j,n_loc,cur_level) +Fc_v     
			endif

            if(dim.eq.2)then
              VV(i,j,n_loc,cur_level)=
     -            VV(i,j,n_loc-1,cur_level)+dt/24.*
     -              (9.*(G(i,j,n_loc,cur_level))+
     -              19.*(G(i,j,n_loc-1,cur_level))-                
     -              5.*(G(i,j,n_loc-2,cur_level))+
     -              (G(i,j,n_loc-3,cur_level)))+
     -              (G1(i,j,n_loc,cur_level)-
     -              G1(i,j,n_loc-1,cur_level))
     
				if(ihvor.eq.1.and.ish_mat(i,j).eq.0)then
			 
				Gc_v=coef1*  
     -              (psiy(i,j,4)  
     -              -psiy(i,j,3)) 
     -           +coef2*  
     -              (psiy(i,j,4)*zeta(i,j,4,cur_level)  
     -              -psiy(i,j,3)*zeta(i,j,3,cur_level)) 
     -          -H_total* 
     -              (psiy(i,j,4)*term1  
     -              -psiy(i,j,3)*term2)
						 
				VV(i,j,n_loc,cur_level)=
     -				VV(i,j,n_loc,cur_level) +Gc_v     
				endif     
            endif

			if(level(i,j).ne.0)then
				UU(i,j,n_loc,cur_level)=0.
				VV(i,j,n_loc,cur_level)=0.
				u(i,j,n_loc,cur_level)=0.
				v(i,j,n_loc,cur_level)=0.
			endif 

          endif
               
         endif

        enddo
       enddo


ccccc Special for seaside
c       if(nn_loc.le.3)then
c       do j=1+overlap,endy-overlap
c        do i=wvmk_loc_i(j),endx-overlap
c			if(x(i).ge.32.5)then
c     			zeta(i,j,n_loc,cur_level)=-h(i,j,n_loc)+swash_d
c			endif
c		enddo
c	   enddo
c	   endif
ccccccccccc

		call FV_LDU(d_FV,zeta(:,:,4,1),4,alpha,dx,dy,sx,ex,sy,ey)

		sxT=1+overlap
		exT=endx-overlap
		syT=1+overlap
		eyT=endy-overlap

	  if(dims(1).eq.1)then
		    call tridiag(LOWERx(sxT:exT,syT:eyT),
     -					DIAGx(sxT:exT,syT:eyT),
     -					UPPERx(sxT:exT,syT:eyT),
     -					RHSx(sxT:exT,syT:eyT),
     -			exT-sxT+1,eyT-syT+1,
     -			u(sxT:exT,syT:eyT,4,cur_level),
     -			v(sxT:exT,syT:eyT,4,cur_level),1)
	  else
		
		    call tridagp(LOWERx(sxT:exT,syT:eyT),
     -					DIAGx(sxT:exT,syT:eyT),
     -					UPPERx(sxT:exT,syT:eyT),
     -					RHSx(sxT:exT,syT:eyT),
     -			exT-sxT+1,eyT-syT+1,xrank,xcomm,
     -			dims(1),u(sxT:exT,syT:eyT,4,cur_level))

	  endif 
        
        if(dim.eq.2)then
	      if(dims(2).eq.1)then
		        call tridiag(LOWERy(sxT:exT,syT:eyT),
     -					DIAGy(sxT:exT,syT:eyT),
     -					UPPERy(sxT:exT,syT:eyT),
     -					RHSy(sxT:exT,syT:eyT),
     -			exT-sxT+1,eyT-syT+1,
     -			u(sxT:exT,syT:eyT,4,cur_level),
     -			v(sxT:exT,syT:eyT,4,cur_level),2)
	      else

		        call tridagp(transpose(LOWERy(sxT:exT,syT:eyT)),
     -					transpose(DIAGy(sxT:exT,syT:eyT)),
     -					transpose(UPPERy(sxT:exT,syT:eyT)),
     -					transpose(RHSy(sxT:exT,syT:eyT)),
     -			eyT-syT+1,exT-sxT+1,yrank,ycomm,
     -			dims(2),v_tmp(syT:eyT,sxT:exT))

				v(sxT:exT,syT:eyT,4,cur_level)=
     -					transpose(v_tmp(syT:eyT,sxT:exT))

	      endif
        endif


	
	  shore_cut=0.90	
       do j=1+overlap,endy-overlap
        do i=wvmk_loc_i(j),endx-overlap

         if(ish_mat(i,j).ge.0)then
          if(n_loc.eq.3)then
C*************** Predictor equations for zeta, u, and v **************

            zeta(i,j,n_loc+1,cur_level) = 
     -              zeta(i,j,n_loc,cur_level)+
     -              dt/12.*
     -              (23*E(i,j,n_loc,cur_level)-
     -              16*E(i,j,n_loc-1,cur_level)+
     -              5*E(i,j,n_loc-2,cur_level))

		 if(zeta(i,j,n_loc+1,cur_level)+h(i,j,n_loc).
     -           le.shore_cut*cutoff_mat(i,j))then 
				zeta(i,j,n_loc+1,cur_level)=max(zeta(i,j,n_loc+1,cur_level),
     -            -h(i,j,n_loc)+cutoff_mat(i,j)*shore_cut)
				
				if(nn_loc.ge.4)
     -            E(i,j,n_loc,cur_level)= 1./23.*(12./dt*
     -            (zeta(i,j,n_loc+1,cur_level)-
     -            zeta(i,j,n_loc,cur_level))+
     -            16*E(i,j,n_loc-1,cur_level)-
     -            5*E(i,j,n_loc-2,cur_level))

             endif

          elseif(n_loc.eq.4)then
C*************** Corrector equations for zeta, u, and v ***************

           zeta(i,j,n_loc,cur_level)=
     -            zeta(i,j,n_loc-1,cur_level)+dt/24.*
     -            (9.*E(i,j,n_loc,cur_level)+
     -            19.*E(i,j,n_loc-1,cur_level)-
     -            5.*E(i,j,n_loc-2,cur_level)+
     -            E(i,j,n_loc-3,cur_level))

		 if(zeta(i,j,n_loc,cur_level)+h(i,j,n_loc).
     -           le.shore_cut*cutoff_mat(i,j))then 
				zeta(i,j,n_loc,cur_level)=max(zeta(i,j,n_loc,cur_level),
     -            -h(i,j,n_loc)+cutoff_mat(i,j)*shore_cut)

				  E(i,j,n_loc,cur_level)= 1./9.*(24./dt*
     -              (zeta(i,j,n_loc,cur_level)
     -            -zeta(i,j,n_loc-1,cur_level))-
     -            19.*E(i,j,n_loc-1,cur_level)+
     -            5.*E(i,j,n_loc-2,cur_level)-
     -            E(i,j,n_loc-3,cur_level))

			endif 

		  endif
		 endif
		enddo
	   enddo

	 endif

      return
      end      ! end subroutine
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


