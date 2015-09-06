subroutine FV_LDU(d,z,n_loc,alpha,dx,dy,sx,ex,sy,ey)

use mainvar_module, only:nx,ny,bl_hor_wall,bl_x_wall,&
bl_y_wall,overlap,endx,endy,v,u,VV,UU,disp_prop,&
cutoff,z_alp,cutoff_mat,loc_sum,bc_1,bc_2,bc_3,bc_4,wvmk_loc_i
use FV_var_module, only:LOWERx,DIAGx,UPPERx,RHSx,LOWERy,DIAGy,&
UPPERy,RHSy,id,ix,iy,icor,limiter2,interp2,zx,zy,ish_mat
implicit none

integer :: run,sum,ii,jj,n_loc,sum_0
integer :: sx,ex,sy,ey,cur_level
integer :: i,j,i_s,i_e,j_s,j_e,bc_x_wall,bc_y_wall
real :: d(nx,ny),z(nx,ny)	
real :: dx,dy,alpha,H_total,z_loc,zx_loc,zy_loc
real :: za,zi(nx,ny),zj(nx,ny)
real :: phil(nx,ny),phir(nx,ny)

cur_level=1

if(icor.eq.1)then
	call interpolations(1,limiter2,id,z,phil,phir)
	zx = (phir-phil)/dx
	call interpolations(4,limiter2,id,z,phil,phir)
	zy = (phir-phil)/dy
else
   call interp_xei (z,zi,interp2,ix)
   call interp_yei (z,zj,interp2,iy)
   do j=sy,ey
      do i=sx,ex
         zx(i,j)  = (zi(i,j)-zi(i-1,j))/dx
         zy(i,j)  = (zj(i,j)-zj(i,j-1))/dy
      enddo
   enddo
endif


	sum_0=0
      do j=3,endy-2
        do i=3,endx-2 
		 H_total=d(i,j)+z(i,j)
         za = -d(i,j)+(1+alpha)*H_total	  
 !        za = alpha*d(i,j)
		 z_alp(i,j,n_loc,1) = za
		 
		 z_loc=z(i,j)
		 zx_loc=zx(i,j)
		 zy_loc=zy(i,j)
		 
		 bc_x_wall=0
		 bc_y_wall=0
		 
		 if(bl_x_wall(i,j).eq.1.and.bc_1.eq.1) bc_x_wall=1
		 if(bl_x_wall(i,j).eq.2.and.bc_2.eq.1) bc_x_wall=1
		 if(bl_y_wall(i,j).eq.1.and.bc_3.eq.1) bc_y_wall=1
		 if(bl_y_wall(i,j).eq.2.and.bc_4.eq.1) bc_y_wall=1
		 
		 RHSy(i,j)=VV(i,j,n_loc,1)/H_total
         if(ish_mat(i,j).eq.1.or.bc_y_wall.eq.1.or.i.lt.wvmk_loc_i(j))then
	        LOWERy(i,j) = 0.0
		    DIAGy(i,j)  = 1.0
  		    UPPERy(i,j) = 0.0
			if(id(i,j).eq.9.or.bc_y_wall.eq.1.or.i.lt.wvmk_loc_i(j)) RHSy(i,j)=v(i,j,n_loc,1)		 
		 else
	        LOWERy(i,j) = (za**2-z_loc**2)/2./dy**2 + (za-z_loc)*&
d(i,j-1)/dy**2 + zy_loc*(z_loc+d(i,j-1))/dy/2.
		    DIAGy(i,j)  = 1-(za**2-z_loc**2)/dy**2 - &
2*(za-z_loc)*d(i,j)/dy**2
  		    UPPERy(i,j) = (za**2-z_loc**2)/2./dy**2 + &
(za-z_loc)*d(i,j+1)/dy**2 - zy_loc*(z_loc+d(i,j+1))/dy/2.
		 endif

             
		 RHSx(i,j)=UU(i,j,n_loc,1)/H_total 
         if(ish_mat(i,j).eq.1.or.bc_x_wall.eq.1.or.i.lt.wvmk_loc_i(j))then
	        LOWERx(i,j) = 0.0
		    DIAGx(i,j)  = 1.0
  		    UPPERx(i,j) = 0.0
			if(id(i,j).eq.9.or.bc_x_wall.eq.1.or.&
i.lt.wvmk_loc_i(j)) RHSx(i,j)=u(i,j,n_loc,1)		 
		 else
	        LOWERx(i,j) = (za**2-z_loc**2)/2./dx**2 + &
(za-z_loc)*d(i-1,j)/dx**2 + zx_loc*(z_loc+d(i-1,j))/dx/2.
		    DIAGx(i,j)  = 1-(za**2-z_loc**2)/dx**2 &
- 2*(za-z_loc)*d(i,j)/dx**2
  		    UPPERx(i,j) = (za**2-z_loc**2)/2./dx**2 + &
(za-z_loc)*d(i+1,j)/dx**2 - zx_loc*(z_loc+d(i+1,j))/dx/2.
		 endif 
                   
      enddo
   enddo
return
end
