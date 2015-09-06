C*******************************************************************************
C........    Groups the physical variables into the 
C........... conservative forms given in Wei & Kirby
      subroutine add_breaking_ktrans(n_loc,nn_loc)

      use mainvar_module
      integer n_loc,nn_loc,nc_loc,nc
	real Cd_brk,sigma_brk,Hloc(endx,endy),
     -              cloc(endx,endy),rat_loc(endx,endy),
     -              nu_src,nu_sink,nu_RHS,
     -              tmp_nut(endx,endy),dzdy_c,trigger

      cur_level=1
      
      if(nn_loc.le.4)then
       do nc=1,4 
        do j=js1,je1
         do i=1,endx
         	B_mat(i,j,nc)=0.	
			B_RHS(i,j,nc)=0.	
			k_RHS(i,j,nc)=0.
			nut(i,j,nc)=0.	
		  enddo
		 enddo
		enddo
	  endif
	  
	  
      do j=js1,je1
        do i=1,endx
            Hloc(i,j)=max(h(i,j,n_loc)+zeta(i,j,n_loc,1),
     -              cutoff_mat_d(i,j))
			cloc(i,j)=sqrt(9.81*Hloc(i,j))	
			tmp_nut(i,j)=0.  
        enddo
      enddo

	  conv_v=0.
	  dzdy_c=0.

      do j=js2,je2
        do i=2,endx-1	  
            if(u(i,j,n_loc,1).gt.0.)then
				conv_u=(cloc(i-1,j)+cloc(i,j))/2.*
     -              (B_mat(i,j,n_loc)-B_mat(i-1,j,n_loc))/(dx)
			else
				conv_u=-(cloc(i+1,j)+cloc(i,j))/2.*
     -              (B_mat(i+1,j,n_loc)-B_mat(i,j,n_loc))/(dx)
			endif

		    if(dim.eq.2)then
              if(v(i,j,n_loc,1).gt.0.)then
				conv_v=(cloc(i,j-1)+cloc(i,j))/2.*
     -              (B_mat(i,j,n_loc)-B_mat(i,j-1,n_loc))/(dy)
			  else
				conv_v=-(cloc(i,j+1)+cloc(i,j))/2.*
     -              (B_mat(i,j+1,n_loc)-B_mat(i,j,n_loc))/(dy)
			  endif
			endif
     
		    rat_loc(i,j)=0.
		    if(E(i,j,n_loc,1).gt.0.)then
		      rat_loc(i,j)=E(i,j,n_loc,1)/cloc(i,j)  
		    endif
			
			trigger=0.65-0.50*0.5*(1+tanh(8.*3.14*(B_mat(i,j,n_loc)-0.125)))
			B=0.5*(1+tanh(10.*(rat_loc(i,j)-trigger)))
			     
			B_RHS(i,j,n_loc)= -conv_u-conv_v
     -			   + 0.25*(3.*B-1.*B_mat(i,j,n_loc))/(Hloc(i,j)/cloc(i,j))	 	   
         enddo
      enddo	  
      
      do j=js2,je2
        do i=2,endx-1

          if(n_loc.eq.3)then
C*************** Predictor equations for B **************

            B_mat(i,j,n_loc+1) = B_mat(i,j,n_loc)+
     -         dt*(B_RHS(i,j,n_loc))  
            B_mat(i,j,n_loc+1)=max(B_mat(i,j,n_loc+1),1.e-16)  
          elseif(n_loc.eq.4)then
C*************** Corrector equations for B ***************

            B_mat(i,j,n_loc)=B_mat(i,j,n_loc-1)+
     -         dt/2.*(B_RHS(i,j,n_loc)+B_RHS(i,j,n_loc-1))   
            B_mat(i,j,n_loc)=max(B_mat(i,j,n_loc),1.e-16)	
          else
            B_mat(i,j,n_loc)=1.e-16
          endif    
          
          if(bl_hor_wall(i,j).eq.99) B_mat(i,j,4)=1.e-16
          
		  if(E(i,j,n_loc,1).gt.0.)then            
             nu_src=B_mat(i,j,n_loc)*Hloc(i,j)
     -           *abs(E(i,j,n_loc,1))   
		  else
			nu_src=0.
		  endif
     
          nu_sink=nut(i,j,n_loc)
          
          k_RHS(i,j,n_loc)=2.*(3.*nu_src-1.*nu_sink)/
     -              (Hloc(i,j)/cloc(i,j))
         enddo
      enddo	  
      
      
      
      do j=js2,je2
        do i=2,endx-1

          if(n_loc.eq.3)then
C*************** Predictor equations for nu **************
            tmp_nut(i,j) = nut(i,j,n_loc)+
     -         dt*(k_RHS(i,j,n_loc))  
            tmp_nut(i,j)=max(tmp_nut(i,j),1.e-16)  
          elseif(n_loc.eq.4)then
C*************** Corrector equations for nu ***************

            tmp_nut(i,j)=nut(i,j,n_loc-1)+
     -         dt/2.*(k_RHS(i,j,n_loc)+k_RHS(i,j,n_loc-1))   
            tmp_nut(i,j)=max(tmp_nut(i,j),1.e-16)	
          else
            tmp_nut(i,j)=1.e-16
          endif                                        

         enddo
      enddo	  
      
      
      if(n_loc.eq.3)then
         nc_loc=n_loc+1
      else
         nc_loc=n_loc
      endif  
      
      if(dim.eq.1)then        
       do j=js3,je3
        do i=3,endx-2

			nut(i,j,nc_loc)=tmp_nut(i,j)

         enddo
       enddo	  
      else     
       do j=js3,je3
        do i=3,endx-2
		
			nut(i,j,nc_loc)=tmp_nut(i,j)

         enddo
       enddo    
      endif
      
       do j=js3,je3
        do i=3,endx-2
          
          stress_xx(1,i,j)=nut(i,j,n_loc)*
     -           (Hloc(i+1,j)*u(i+1,j,n_loc,1) - 
     -            Hloc(i-1,j)*u(i-1,j,n_loc,1))/(2.*dx)
     
          stress_xx(2,i,j)=nut(i,j,n_loc)*
     -           (Hloc(i,j+1)*u(i,j+1,n_loc,1) - 
     -            Hloc(i,j-1)*u(i,j-1,n_loc,1))/(2.*dy)
     -           +nut(i,j,n_loc)*
     -           (Hloc(i+1,j)*v(i+1,j,n_loc,1) - 
     -            Hloc(i-1,j)*v(i-1,j,n_loc,1))/(2.*dx)
     
          stress_xx(3,i,j)=nut(i,j,n_loc)*
     -           (Hloc(i,j+1)*v(i,j+1,n_loc,1) - 
     -            Hloc(i,j-1)*v(i,j-1,n_loc,1))/(2.*dy)
               
       enddo
      enddo

      do j=1+overlap,endy-overlap
       do i=1+overlap,endx-overlap
C%%%%%%%%%%%%%%%%%%%%%%%  WAVE BREAKING AND EDDY VISC DISSIPATION %%%%%%%%%%%%%%%%%%%%%%
        WB_x(i,j)=0.
        WB_y(i,j)=0.
        deriv_1x=0.
        deriv_1y=0.
        deriv_2x=0.
        deriv_2y=0.

        if(bl_hor_wall(i,j).le.50)then

	    deriv_1x=(stress_xx(1,i+1,j)-
     -            stress_xx(1,i-1,j))/(2.*dx)
	    if(dim.eq.2)then
	        deriv_2y=(stress_xx(2,i+1,j)-
     -              stress_xx(2,i-1,j))/(2.*dx)

	        deriv_1y=(stress_xx(3,i,j+1)-
     -              stress_xx(3,i,j-1))/(2.*dy)

	        deriv_2x=(stress_xx(2,i,j+1)-
     -              stress_xx(2,i,j-1))/(2.*dy)

            WB_y(i,j)=(deriv_1y+ 0.5*deriv_2y)
          endif   

          WB_x(i,j)=(deriv_1x+ 0.5*deriv_2x)

          param_G=WB_y(i,j)
          param_F=WB_x(i,j)

		F(i,j,n_loc,cur_level)=F(i,j,n_loc,cur_level)+param_F
		G(i,j,n_loc,cur_level)=G(i,j,n_loc,cur_level)+param_G

          endif
        enddo
      enddo



      return 

      end      ! end subroutine
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


