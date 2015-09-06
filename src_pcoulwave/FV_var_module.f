      module FV_var_module
      implicit none

      integer sx,ex,sy,ey,idissip,ibotfric,ilim,
     -		limiter2,limiter3,iriemann,icor,interp2

	integer, ALLOCATABLE :: ix(:,:),iy(:,:),ix2(:,:),iy2(:,:),
     -		id(:,:),ish_mat(:,:)

      real b0,b1,cm,ks,alpha,BS
      
      real, ALLOCATABLE :: z_FV(:,:),d_FV(:,:),u_FV(:,:),rand(:,:,:),
     - v_FV(:,:),E_FV(:,:),F_FV(:,:),F1_FV(:,:),G_FV(:,:),
     - G1_FV(:,:),LOWERx(:,:),DIAGx(:,:),UPPERx(:,:),RHSx(:,:), 
     - LOWERy(:,:),DIAGy(:,:),UPPERy(:,:),RHSy(:,:),zx(:,:),zy(:,:)      
      
      end module FV_var_module


