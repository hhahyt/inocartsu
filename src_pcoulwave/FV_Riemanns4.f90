subroutine Riemanns4 (H_l,H_r,Hu_l,Hu_r,Hv_l,Hv_r,&
H_b,H_t,Hu_b,Hu_t,Hv_b,Hv_t,flux,fluy,dir,sx,ex,sy,ey,ix,iy,iriemann)

use mainvar_module, only:nx,ny,cutoff_mat
implicit none

integer dir,i,j,sx,sy,ex,ey,ix(nx,ny),iy(nx,ny),iriemann

real :: pi,sl,sm,sr,dl,dr,ul,ur,grav,vl,vr,H_l(nx,ny),&
H_r(nx,ny),Hu_l(nx,ny),Hu_r(nx,ny),Hv_l(nx,ny),Hv_r(nx,ny),H_b(nx,ny),H_t(nx,ny)
real :: Hu_b(nx,ny),Hu_t(nx,ny)
real :: Hv_b(nx,ny),Hv_t(nx,ny)
real :: cdl(6),cdr(6),fsl(6),fsr(6),csl(6),csr(6),fdl(6),fdr(6)
real :: flux(nx,ny,3),fluy(nx,ny,3)



pi = 3.141592
grav = 9.81

if (iriemann.eq.1) then

!..... X-direction

do j = sy,ey
   do i = sx-1,ex
         cdl(1) = H_l(i,j)
         cdr(1) = H_r(i,j)
         cdl(2) = Hu_l(i,j)
         cdr(2) = Hu_r(i,j)
         cdl(3) = Hv_l(i,j)
         cdr(3) = Hv_r(i,j)

!   ROTA-FW

      if (cdl(1).gt.cutoff_mat(i,j)) then
	     fdl(1) = cdl(2)
         fdl(2) = cdl(2)*cdl(2)/cdl(1) !+ 0.5*grav*cdl(1)*cdl(1)
         fdl(3) = cdl(2)*cdl(3)/cdl(1)

         dl = cdl(1)
         ul = cdl(2)/dl
         vl = cdl(3)/dl

      else
	     fdl(1) = 0
         fdl(2) = 0
         fdl(3) = 0

         dl = 0
         ul = 0
         vl = 0

      endif

!   ROTA-FW

      if (cdr(1).gt.cutoff_mat(i,j)) then
         fdr(1) = cdr(2)
         fdr(2) = cdr(2)*cdr(2)/cdr(1) !+ 0.5*grav*cdr(1)*cdr(1)
         fdr(3) = cdr(2)*cdr(3)/cdr(1)

         dr = cdr(1)
         ur = cdr(2)/dr
	     vr = cdr(3)/dr

      else
         fdr(1) = 0
         fdr(2) = 0
         fdr(3) = 0

         dr = 0
         ur = 0
	     vr = 0

      endif

      call estime(sl,sm,sr,dl,ul,dr,ur)

      csl(1) = dl*(sl - ul)/(sl - sm)
      csl(2) = csl(1)*sm
      csl(3) = csl(1)*vl

      csr(1) = dr*(sr - ur)/(sr - sm)
      csr(2) = csr(1)*sm
      csr(3) = csr(1)*vr

      fsl(1) = fdl(1) + sl*(csl(1) - cdl(1))
      fsr(1) = fdr(1) + sr*(csr(1) - cdr(1))
      fsl(2) = fdl(2) + sl*(csl(2) - cdl(2))
      fsr(2) = fdr(2) + sr*(csr(2) - cdr(2))

      fsl(3) = fdl(3) + sl*(csl(3) - cdl(3))
      fsr(3) = fdr(3) + sr*(csr(3) - cdr(3))


      if (sl.ge.0.d0) then 
         flux(i,j,1) = fdl(1)
         flux(i,j,2) = fdl(2)
         flux(i,j,3) = fdl(3)
	  elseif (sl.le.0.d0 .and. sm.ge.0.d0) then
         flux(i,j,1) = fsl(1)
         flux(i,j,2) = fsl(2)
         flux(i,j,3) = fsl(3)
      elseif (sm.le.0.d0 .and. sr.ge.0.d0) then
         flux(i,j,1) = fsr(1)
         flux(i,j,2) = fsr(2)
         flux(i,j,3) = fsr(3)
      elseif (sr.le.0.d0) then
         flux(i,j,1) = fdr(1)
         flux(i,j,2) = fdr(2)
         flux(i,j,3) = fdr(3)
      endif

!   ROTA-BK
!      if (sm .ge. 0.d0) then
!         flux(i,j,3) = flux(i,j,1)*vl !*sleng(i)
!      else
!         flux(i,j,3) = flux(i,j,1)*vr !*sleng(i)
!      endif

   enddo
enddo


!..... Y-direction

do i = sx,ex
   do j = sy-1,ey
         cdl(4) = H_b(i,j)
         cdr(4) = H_t(i,j)
         cdl(5) = Hu_b(i,j)
         cdr(5) = Hu_t(i,j)
         cdl(6) = Hv_b(i,j)
         cdr(6) = Hv_t(i,j)

!   ROTA-FW

      if (cdl(4).gt.cutoff_mat(i,j)) then
	  fdl(4) = cdl(6)
      fdl(5) = cdl(6)*cdl(5)/cdl(4)
      fdl(6) = cdl(6)*cdl(6)/cdl(4) !+ 0.5*grav*cdl(4)*cdl(4)

      dl = cdl(4)
      ul = cdl(5)/dl
      vl = cdl(6)/dl
      
	  else
	  fdl(4) = 0
      fdl(5) = 0
      fdl(6) = 0

      dl = 0
      ul = 0
      vl = 0
	  endif

!   ROTA-FW

      if (cdr(4).gt.cutoff_mat(i,j)) then
      fdr(4) = cdr(6)
      fdr(5) = cdr(6)*cdr(5)/cdr(4)
      fdr(6) = cdr(6)*cdr(6)/cdr(4) !+ 0.5*grav*cdr(4)*cdr(4)

      dr = cdr(4)
      ur = cdr(5)/dr
	  vr = cdr(6)/dr
      
	  else
      fdr(4) = 0
      fdr(5) = 0
      fdr(6) = 0

      dr = 0
      ur = 0
	  vr = 0
	  endif


      call estime(sl,sm,sr,dl,vl,dr,vr)

      csl(4) = dl*(sl - vl)/(sl - sm)
      csl(6) = csl(4)*sm
      csr(4) = dr*(sr - vr)/(sr - sm)
      csr(6) = csr(4)*sm

      csl(5) = csl(4)*ul
      csr(5) = csr(4)*ur

      fsl(4) = fdl(4) + sl*(csl(4) - cdl(4))
      fsr(4) = fdr(4) + sr*(csr(4) - cdr(4))
      fsl(6) = fdl(6) + sl*(csl(6) - cdl(6))
      fsr(6) = fdr(6) + sr*(csr(6) - cdr(6))

      fsl(5) = fdl(5) + sl*(csl(5) - cdl(5))
      fsr(5) = fdr(5) + sr*(csr(5) - cdr(5))

      if (sl.ge.0.d0) then 
         fluy(i,j,1) = fdl(4)
         fluy(i,j,3) = fdl(6)
         fluy(i,j,2) = fdl(5)
	  elseif (sl.le.0.d0 .and. sm.ge.0.d0) then
         fluy(i,j,1) = fsl(4)
         fluy(i,j,3) = fsl(6)
         fluy(i,j,2) = fsl(5)
      elseif (sm.le.0.d0 .and. sr.ge.0.d0) then
         fluy(i,j,1) = fsr(4)
         fluy(i,j,3) = fsr(6)
         fluy(i,j,2) = fsr(5)
      elseif (sr.le.0.d0) then
         fluy(i,j,1) = fdr(4)
         fluy(i,j,3) = fdr(6)
         fluy(i,j,2) = fdr(5)
      endif

!   ROTA-BK

!      if (sm .ge. 0.d0) then
!         fluy(i,j,2) = fluy(i,j,1)*ul !*sleng(i)
!      else
!         fluy(i,j,2) = fluy(i,j,1)*ur !*sleng(i)
!      endif
   enddo
enddo

elseif (iriemann.eq.2) then

do j = sy,ey
   do i = sx-1,ex
         cdl(1) = H_l(i,j)
         cdr(1) = H_r(i,j)
         cdl(2) = Hu_l(i,j)
         cdr(2) = Hu_r(i,j)
         cdl(3) = Hv_l(i,j)
         cdr(3) = Hv_r(i,j)

      if (cdl(1).gt.cutoff_mat(i,j)) then
	     fdl(1) = cdl(2)
         fdl(2) = cdl(2)*cdl(2)/cdl(1) !+ 0.5*grav*cdl(1)*cdl(1)
         fdl(3) = cdl(2)*cdl(3)/cdl(1)

         dl = cdl(1)
         ul = cdl(2)/dl
         vl = cdl(3)/dl

      else
	     fdl(1) = 0
         fdl(2) = 0
         fdl(3) = 0

         dl = 0
         ul = 0
         vl = 0

      endif

      if (cdr(1).gt.cutoff_mat(i,j)) then
         fdr(1) = cdr(2)
         fdr(2) = cdr(2)*cdr(2)/cdr(1) !+ 0.5*grav*cdr(1)*cdr(1)
         fdr(3) = cdr(2)*cdr(3)/cdr(1)

         dr = cdr(1)
         ur = cdr(2)/dr
	     vr = cdr(3)/dr

      else
         fdr(1) = 0
         fdr(2) = 0
         fdr(3) = 0

         dr = 0
         ur = 0
	     vr = 0

      endif

      call estime(sl,sm,sr,dl,ul,dr,ur)

      if (sl.ge.0.d0) then 
         flux(i,j,1) = fdl(1)
         flux(i,j,2) = fdl(2)
         flux(i,j,3) = fdl(3)
	  elseif (sl.le.0.d0 .and. sr.ge.0.d0) then
         flux(i,j,1) = (sr*fdl(1)-sl*fdr(1)+sr*sl*(cdr(1)-cdl(1)))/(sr-sl)
         flux(i,j,2) = (sr*fdl(2)-sl*fdr(2)+sr*sl*(cdr(2)-cdl(2)))/(sr-sl)
         flux(i,j,3) = (sr*fdl(3)-sl*fdr(3)+sr*sl*(cdr(3)-cdl(3)))/(sr-sl)
      elseif (sr.le.0.d0) then
         flux(i,j,1) = fdr(1)
         flux(i,j,2) = fdr(2)
         flux(i,j,3) = fdr(3)
      endif
   enddo
enddo

!..... Y-direction

do i = sx,ex
   do j = sy-1,ey
         cdl(4) = H_b(i,j)
         cdr(4) = H_t(i,j)
         cdl(5) = Hu_b(i,j)
         cdr(5) = Hu_t(i,j)
         cdl(6) = Hv_b(i,j)
         cdr(6) = Hv_t(i,j)

      if (cdl(4).gt.cutoff_mat(i,j)) then
	  fdl(4) = cdl(6)
      fdl(5) = cdl(6)*cdl(5)/cdl(4)
      fdl(6) = cdl(6)*cdl(6)/cdl(4) !+ 0.5*grav*cdl(4)*cdl(4)

      dl = cdl(4)
      ul = cdl(5)/dl
      vl = cdl(6)/dl
      
	  else
	  fdl(4) = 0
      fdl(5) = 0
      fdl(6) = 0

      dl = 0
      ul = 0
      vl = 0
	  endif

      if (cdr(4).gt.cutoff_mat(i,j)) then
      fdr(4) = cdr(6)
      fdr(5) = cdr(6)*cdr(5)/cdr(4)
      fdr(6) = cdr(6)*cdr(6)/cdr(4) !+ 0.5*grav*cdr(4)*cdr(4)

      dr = cdr(4)
      ur = cdr(5)/dr
	  vr = cdr(6)/dr
      
	  else
      fdr(4) = 0
      fdr(5) = 0
      fdr(6) = 0

      dr = 0
      ur = 0
	  vr = 0
	  endif

      call estime(sl,sm,sr,dl,vl,dr,vr)

      if (sl.ge.0.d0) then 
         fluy(i,j,1) = fdl(4)
         fluy(i,j,2) = fdl(5)
         fluy(i,j,3) = fdl(6)
	  elseif (sl.le.0.d0 .and. sr.ge.0.d0) then
         fluy(i,j,1) = (sr*fdl(4)-sl*fdr(4)+sr*sl*(cdr(4)-cdl(4)))/(sr-sl)
         fluy(i,j,2) = (sr*fdl(5)-sl*fdr(5)+sr*sl*(cdr(5)-cdl(5)))/(sr-sl)
         fluy(i,j,3) = (sr*fdl(6)-sl*fdr(6)+sr*sl*(cdr(6)-cdl(6)))/(sr-sl)
      elseif (sr.le.0.d0) then
         fluy(i,j,1) = fdr(4)
         fluy(i,j,2) = fdr(5)
         fluy(i,j,3) = fdr(6)
      endif

   enddo
enddo

elseif (iriemann.eq.3) then

!..... X-direction

do j = sy,ey
   do i = sx-1,ex
         cdl(1) = H_l(i,j)
         cdr(1) = H_r(i,j)
         cdl(2) = Hu_l(i,j)
         cdr(2) = Hu_r(i,j)
         cdl(3) = Hv_l(i,j)
         cdr(3) = Hv_r(i,j)

!   ROTA-FW

      if (cdl(1).gt.cutoff_mat(i,j)) then
	     fdl(1) = cdl(2)
         fdl(2) = cdl(2)*cdl(2)/cdl(1) !+ 0.5*grav*cdl(1)*cdl(1)
         fdl(3) = cdl(2)*cdl(3)/cdl(1)

         dl = cdl(1)
         ul = cdl(2)/dl
         vl = cdl(3)/dl

      else
	     fdl(1) = 0
         fdl(2) = 0
         fdl(3) = 0

         dl = 0
         ul = 0
         vl = 0

      endif

!   ROTA-FW

      if (cdr(1).gt.cutoff_mat(i,j)) then
         fdr(1) = cdr(2)
         fdr(2) = cdr(2)*cdr(2)/cdr(1) !+ 0.5*grav*cdr(1)*cdr(1)
         fdr(3) = cdr(2)*cdr(3)/cdr(1)

         dr = cdr(1)
         ur = cdr(2)/dr
	     vr = cdr(3)/dr

      else
         fdr(1) = 0
         fdr(2) = 0
         fdr(3) = 0

         dr = 0
         ur = 0
	     vr = 0

      endif

      call estime(sl,sm,sr,dl,ul,dr,ur)

      csl(1) = dl*(sl - ul)/(sl - sm)
      csl(2) = csl(1)*sm
      csl(3) = csl(1)*vl

      csr(1) = dr*(sr - ur)/(sr - sm)
      csr(2) = csr(1)*sm
      csr(3) = csr(1)*vr

      fsl(1) = fdl(1) + sl*(csl(1) - cdl(1))
      fsr(1) = fdr(1) + sr*(csr(1) - cdr(1))
      fsl(2) = fdl(2) + sl*(csl(2) - cdl(2))
      fsr(2) = fdr(2) + sr*(csr(2) - cdr(2))

      fsl(3) = fdl(3) + sl*(csl(3) - cdl(3))
      fsr(3) = fdr(3) + sr*(csr(3) - cdr(3))


      if (sl.ge.0.d0) then 
         flux(i,j,1) = fdl(1)
         flux(i,j,2) = fdl(2)
         flux(i,j,3) = fdl(3)
	  elseif (sl.le.0.d0 .and. sm.ge.0.d0) then
         flux(i,j,1) = fsl(1)
         flux(i,j,2) = fsl(2)
         flux(i,j,3) = fsl(3)
      elseif (sm.le.0.d0 .and. sr.ge.0.d0) then
         flux(i,j,1) = fsr(1)
         flux(i,j,2) = fsr(2)
         flux(i,j,3) = fsr(3)
      elseif (sr.le.0.d0) then
         flux(i,j,1) = fdr(1)
         flux(i,j,2) = fdr(2)
         flux(i,j,3) = fdr(3)
      endif

      if (sl.ge.0.d0) then 
         flux(i,j,3) = fdl(3)
	  elseif (sl.le.0.d0 .and. sr.ge.0.d0) then
         flux(i,j,3) = (sr*fdl(3)-sl*fdr(3)+sr*sl*(cdr(3)-cdl(3)))/(sr-sl)
      elseif (sr.le.0.d0) then
         flux(i,j,3) = fdr(3)
      endif

   enddo
enddo


!..... Y-direction

do i = sx,ex
   do j = sy-1,ey
         cdl(4) = H_b(i,j)
         cdr(4) = H_t(i,j)
         cdl(5) = Hu_b(i,j)
         cdr(5) = Hu_t(i,j)
         cdl(6) = Hv_b(i,j)
         cdr(6) = Hv_t(i,j)

!   ROTA-FW

      if (cdl(4).gt.cutoff_mat(i,j)) then
	  fdl(4) = cdl(6)
      fdl(5) = cdl(6)*cdl(5)/cdl(4)
      fdl(6) = cdl(6)*cdl(6)/cdl(4) !+ 0.5*grav*cdl(4)*cdl(4)

      dl = cdl(4)
      ul = cdl(5)/dl
      vl = cdl(6)/dl
      
	  else
	  fdl(4) = 0
      fdl(5) = 0
      fdl(6) = 0

      dl = 0
      ul = 0
      vl = 0
	  endif

!   ROTA-FW

      if (cdr(4).gt.cutoff_mat(i,j)) then
      fdr(4) = cdr(6)
      fdr(5) = cdr(6)*cdr(5)/cdr(4)
      fdr(6) = cdr(6)*cdr(6)/cdr(4) !+ 0.5*grav*cdr(4)*cdr(4)

      dr = cdr(4)
      ur = cdr(5)/dr
	  vr = cdr(6)/dr
      
	  else
      fdr(4) = 0
      fdr(5) = 0
      fdr(6) = 0

      dr = 0
      ur = 0
	  vr = 0
	  endif


      call estime(sl,sm,sr,dl,vl,dr,vr)

      csl(4) = dl*(sl - vl)/(sl - sm)
      csl(6) = csl(4)*sm
      csr(4) = dr*(sr - vr)/(sr - sm)
      csr(6) = csr(4)*sm

      csl(5) = csl(4)*ul
      csr(5) = csr(4)*ur

      fsl(4) = fdl(4) + sl*(csl(4) - cdl(4))
      fsr(4) = fdr(4) + sr*(csr(4) - cdr(4))
      fsl(6) = fdl(6) + sl*(csl(6) - cdl(6))
      fsr(6) = fdr(6) + sr*(csr(6) - cdr(6))

      fsl(5) = fdl(5) + sl*(csl(5) - cdl(5))
      fsr(5) = fdr(5) + sr*(csr(5) - cdr(5))

      if (sl.ge.0.d0) then 
         fluy(i,j,1) = fdl(4)
         fluy(i,j,3) = fdl(6)
         fluy(i,j,2) = fdl(5)
	  elseif (sl.le.0.d0 .and. sm.ge.0.d0) then
         fluy(i,j,1) = fsl(4)
         fluy(i,j,3) = fsl(6)
         fluy(i,j,2) = fsl(5)
      elseif (sm.le.0.d0 .and. sr.ge.0.d0) then
         fluy(i,j,1) = fsr(4)
         fluy(i,j,3) = fsr(6)
         fluy(i,j,2) = fsr(5)
      elseif (sr.le.0.d0) then
         fluy(i,j,1) = fdr(4)
         fluy(i,j,3) = fdr(6)
         fluy(i,j,2) = fdr(5)
      endif

      if (sl.ge.0.d0) then 
         fluy(i,j,2) = fdl(5)
	  elseif (sl.le.0.d0 .and. sr.ge.0.d0) then
         fluy(i,j,2) = (sr*fdl(5)-sl*fdr(5)+sr*sl*(cdr(5)-cdl(5)))/(sr-sl)
      elseif (sr.le.0.d0) then
         fluy(i,j,2) = fdr(5)
      endif
   enddo
enddo


elseif (iriemann.eq.4) then

!..... X-direction

do j = sy,ey
   do i = sx-1,ex
         cdl(1) = H_l(i,j)
         cdr(1) = H_r(i,j)
         cdl(2) = Hu_l(i,j)
         cdr(2) = Hu_r(i,j)
         cdl(3) = Hv_l(i,j)
         cdr(3) = Hv_r(i,j)

!   ROTA-FW

      if (cdl(1).gt.cutoff_mat(i,j)) then
	     fdl(1) = cdl(2)
         fdl(2) = cdl(2)*cdl(2)/cdl(1) !+ 0.5*grav*cdl(1)*cdl(1)
         fdl(3) = cdl(2)*cdl(3)/cdl(1)

         dl = cdl(1)
         ul = cdl(2)/dl
         vl = cdl(3)/dl

      else
	     fdl(1) = 0
         fdl(2) = 0
         fdl(3) = 0

         dl = 0
         ul = 0
         vl = 0

      endif

!   ROTA-FW

      if (cdr(1).gt.cutoff_mat(i,j)) then
         fdr(1) = cdr(2)
         fdr(2) = cdr(2)*cdr(2)/cdr(1) !+ 0.5*grav*cdr(1)*cdr(1)
         fdr(3) = cdr(2)*cdr(3)/cdr(1)

         dr = cdr(1)
         ur = cdr(2)/dr
	     vr = cdr(3)/dr

      else
         fdr(1) = 0
         fdr(2) = 0
         fdr(3) = 0

         dr = 0
         ur = 0
	     vr = 0

      endif

      call estime(sl,sm,sr,dl,ul,dr,ur)

      csl(1) = dl*(sl - ul)/(sl - sm)
      csl(2) = csl(1)*sm
      csl(3) = csl(1)*vl

      csr(1) = dr*(sr - ur)/(sr - sm)
      csr(2) = csr(1)*sm
      csr(3) = csr(1)*vr

      fsl(1) = fdl(1) + sl*(csl(1) - cdl(1))
      fsr(1) = fdr(1) + sr*(csr(1) - cdr(1))
      fsl(2) = fdl(2) + sl*(csl(2) - cdl(2))
      fsr(2) = fdr(2) + sr*(csr(2) - cdr(2))

      fsl(3) = fdl(3) + sl*(csl(3) - cdl(3))
      fsr(3) = fdr(3) + sr*(csr(3) - cdr(3))


      if (sl.ge.0.d0) then 
         flux(i,j,1) = fdl(1)
         flux(i,j,2) = fdl(2)
	  elseif (sl.le.0.d0 .and. sm.ge.0.d0) then
         flux(i,j,1) = fsl(1)
         flux(i,j,2) = fsl(2)
      elseif (sm.le.0.d0 .and. sr.ge.0.d0) then
         flux(i,j,1) = fsr(1)
         flux(i,j,2) = fsr(2)
      elseif (sr.le.0.d0) then
         flux(i,j,1) = fdr(1)
         flux(i,j,2) = fdr(2)
      endif


      if (sm .ge. 0.d0) then
         flux(i,j,3) = flux(i,j,1)*vl !*sleng(i)
      else
         flux(i,j,3) = flux(i,j,1)*vr !*sleng(i)
      endif

   enddo
enddo


!..... Y-direction

do i = sx,ex
   do j = sy-1,ey
         cdl(4) = H_b(i,j)
         cdr(4) = H_t(i,j)
         cdl(5) = Hu_b(i,j)
         cdr(5) = Hu_t(i,j)
         cdl(6) = Hv_b(i,j)
         cdr(6) = Hv_t(i,j)

!   ROTA-FW

      if (cdl(4).gt.cutoff_mat(i,j)) then
	  fdl(4) = cdl(6)
      fdl(5) = cdl(6)*cdl(5)/cdl(4)
      fdl(6) = cdl(6)*cdl(6)/cdl(4) !+ 0.5*grav*cdl(4)*cdl(4)

      dl = cdl(4)
      ul = cdl(5)/dl
      vl = cdl(6)/dl
      
	  else
	  fdl(4) = 0
      fdl(5) = 0
      fdl(6) = 0

      dl = 0
      ul = 0
      vl = 0
	  endif

!   ROTA-FW

      if (cdr(4).gt.cutoff_mat(i,j)) then
      fdr(4) = cdr(6)
      fdr(5) = cdr(6)*cdr(5)/cdr(4)
      fdr(6) = cdr(6)*cdr(6)/cdr(4) !+ 0.5*grav*cdr(4)*cdr(4)

      dr = cdr(4)
      ur = cdr(5)/dr
	  vr = cdr(6)/dr
      
	  else
      fdr(4) = 0
      fdr(5) = 0
      fdr(6) = 0

      dr = 0
      ur = 0
	  vr = 0
	  endif


      call estime(sl,sm,sr,dl,vl,dr,vr)

      csl(4) = dl*(sl - vl)/(sl - sm)
      csl(6) = csl(4)*sm
      csr(4) = dr*(sr - vr)/(sr - sm)
      csr(6) = csr(4)*sm

      csl(5) = csl(4)*ul
      csr(5) = csr(4)*ur

      fsl(4) = fdl(4) + sl*(csl(4) - cdl(4))
      fsr(4) = fdr(4) + sr*(csr(4) - cdr(4))
      fsl(6) = fdl(6) + sl*(csl(6) - cdl(6))
      fsr(6) = fdr(6) + sr*(csr(6) - cdr(6))

      fsl(5) = fdl(5) + sl*(csl(5) - cdl(5))
      fsr(5) = fdr(5) + sr*(csr(5) - cdr(5))

      if (sl.ge.0.d0) then 
         fluy(i,j,1) = fdl(4)
         fluy(i,j,3) = fdl(6)
	  elseif (sl.le.0.d0 .and. sm.ge.0.d0) then
         fluy(i,j,1) = fsl(4)
         fluy(i,j,3) = fsl(6)
      elseif (sm.le.0.d0 .and. sr.ge.0.d0) then
         fluy(i,j,1) = fsr(4)
         fluy(i,j,3) = fsr(6)
      elseif (sr.le.0.d0) then
         fluy(i,j,1) = fdr(4)
         fluy(i,j,3) = fdr(6)

      endif

!   ROTA-BK

      if (sm .ge. 0.d0) then
         fluy(i,j,2) = fluy(i,j,1)*ul !*sleng(i)
      else
         fluy(i,j,2) = fluy(i,j,1)*ur !*sleng(i)
      endif

   enddo
enddo

endif

end
