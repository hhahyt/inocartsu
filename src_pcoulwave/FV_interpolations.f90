subroutine interpolations(flag,limiter,id,u,ul,ur)

use mainvar_module, only:nx,ny,dim,overlap
implicit none

integer :: sx,ex,sy,ey,sx1,ex1,sy1,ey1,i,j,flag,limiter
integer :: id(nx,ny)
real :: omeg,dupw,dloc,up_ratio,delta,pil,pir
real :: u(nx,ny),ur(nx,ny),ul(nx,ny)


if(dim.eq.1)then
		sy=overlap+1
		ey=overlap+1
else
		sy=2
		ey=ny-2
endif
sx=2
ex=nx-2

if (flag.eq.1) then  ! eta_x

do i = sx,ex
   do j = sy,ey

	     omeg =  0.d0
		 dupw = u(i,j)-u(i-1,j)
		 dloc = u(i+1,j)-u(i,j)

	  delta = 0.5d0*(1.d0+omeg)*dupw + 0.5d0*(1.D0-omeg)*dloc
	  up_ratio = dupw/dloc

      IF(limiter.EQ.1) DELTA = 0.D0  ! 1st order
      IF(limiter.EQ.2) DELTA = DELTA ! Centred 2nd order
      IF(limiter.EQ.3) CALL Superb(up_RATIO,OMEG,DELTA) ! Superbee like limiter
      IF(limiter.EQ.4) CALL Minmod(up_RATIO,OMEG,DELTA) ! Minimod like limiter
      IF(limiter.EQ.5) CALL Minmax(DUPW, DLOC,DELTA)    ! Minimod like limiter 2

         PIL = u(i,j)-0.5D0*DELTA
         PIR = u(i,j)+0.5D0*DELTA

      ul(i,j)  = pil
      ur(i,j)  = pir

   enddo
enddo

elseif (flag.eq.2) then  ! u_x

do i = sx,ex
   do j = sy,ey

         DUPW = u(i,j) - u(i-1,j)
         DLOC = u(i+1,j) - u(i,j)


      DELTA = 0.5D0*(1.D0+OMEG)*DUPW+0.5D0*(1.D0-OMEG)*DLOC
      up_RATIO = DUPW/DLOC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      if (dloc.eq.0) ratio = 1

      IF(limiter.EQ.1) DELTA = 0.D0  ! 1st order
      IF(limiter.EQ.2) DELTA = DELTA ! Centred 2nd order
      IF(limiter.EQ.3) CALL Superb(up_RATIO,OMEG,DELTA) ! Superbee like limiter
      IF(limiter.EQ.4) CALL Minmod(up_RATIO,OMEG,DELTA) ! Minimod like limiter
      IF(limiter.EQ.5) CALL Minmax(DUPW, DLOC,DELTA)    ! Minimod like limiter 2

      PIL = u(i,j) - 0.5D0*DELTA
      PIR = u(i,j) + 0.5D0*DELTA

      ul(i,j)  = pil
      ur(i,j)  = pir

   enddo
enddo

elseif (flag.eq.3) then  ! v_x

do i = sx,ex
   do j = sy,ey

         DUPW = u(i,j) - u(i-1,j)
         DLOC = u(i+1,j) - u(i,j)

      DELTA = 0.5D0*(1.D0+OMEG)*DUPW+0.5D0*(1.D0-OMEG)*DLOC
      up_RATIO = DUPW/DLOC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      if (dloc.eq.0) ratio = 1

      IF(limiter.EQ.1) DELTA = 0.D0  ! 1st order
      IF(limiter.EQ.2) DELTA = DELTA ! Centred 2nd order
      IF(limiter.EQ.3) CALL Superb(up_RATIO,OMEG,DELTA) ! Superbee like limiter
      IF(limiter.EQ.4) CALL Minmod(up_RATIO,OMEG,DELTA) ! Minimod like limiter
      IF(limiter.EQ.5) CALL Minmax(DUPW, DLOC,DELTA)    ! Minimod like limiter 2

      PIL = u(i,j) - 0.5D0*DELTA
      PIR = u(i,j) + 0.5D0*DELTA

      ul(i,j)  = pil
      ur(i,j)  = pir

   enddo
enddo

elseif (flag.eq.4) then  ! eta_y

do i = sx,ex
   do j = sy,ey

	     omeg =  0.d0
		 dupw = u(i,j)-u(i,j-1)
		 dloc = u(i,j+1)-u(i,j)

	  delta = 0.5d0*(1.d0+omeg)*dupw + 0.5d0*(1.D0-omeg)*dloc
	  up_ratio = dupw/dloc

      IF(limiter.EQ.1) DELTA = 0.D0  ! 1st order
      IF(limiter.EQ.2) DELTA = DELTA ! Centred 2nd order
      IF(limiter.EQ.3) CALL Superb(up_RATIO,OMEG,DELTA) ! Superbee like limiter
      IF(limiter.EQ.4) CALL Minmod(up_RATIO,OMEG,DELTA) ! Minimod like limiter
      IF(limiter.EQ.5) CALL Minmax(DUPW, DLOC,DELTA)    ! Minimod like limiter 2

         PIL = u(i,j)-0.5D0*DELTA
         PIR = u(i,j)+0.5D0*DELTA

      ul(i,j)  = pil
      ur(i,j)  = pir

   enddo
enddo

elseif (flag.eq.5) then  ! u_y

do i = sx,ex
   do j = sy,ey

         DUPW = u(i,j) - u(i,j-1)
         DLOC = u(i,j+1) - u(i,j)

      DELTA = 0.5D0*(1.D0+OMEG)*DUPW+0.5D0*(1.D0-OMEG)*DLOC
      up_RATIO = DUPW/DLOC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      if (dloc.eq.0) ratio = 1

      IF(limiter.EQ.1) DELTA = 0.D0  ! 1st order
      IF(limiter.EQ.2) DELTA = DELTA ! Centred 2nd order
      IF(limiter.EQ.3) CALL Superb(up_RATIO,OMEG,DELTA) ! Superbee like limiter
      IF(limiter.EQ.4) CALL Minmod(up_RATIO,OMEG,DELTA) ! Minimod like limiter
      IF(limiter.EQ.5) CALL Minmax(DUPW, DLOC,DELTA)    ! Minimod like limiter 2

      PIL = u(i,j) - 0.5D0*DELTA
      PIR = u(i,j) + 0.5D0*DELTA

      ul(i,j)  = pil
      ur(i,j)  = pir

   enddo
enddo

elseif (flag.eq.6) then  ! v_y

do i = sx,ex
   do j = sy,ey
   
         DUPW = u(i,j) - u(i,j-1)
         DLOC = u(i,j+1) - u(i,j)

      DELTA = 0.5D0*(1.D0+OMEG)*DUPW+0.5D0*(1.D0-OMEG)*DLOC
      up_RATIO = DUPW/DLOC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      if (dloc.eq.0) ratio = 1

      IF(limiter.EQ.1) DELTA = 0.D0  ! 1st order
      IF(limiter.EQ.2) DELTA = DELTA ! Centred 2nd order
      IF(limiter.EQ.3) CALL Superb(up_RATIO,OMEG,DELTA) ! Superbee like limiter
      IF(limiter.EQ.4) CALL Minmod(up_RATIO,OMEG,DELTA) ! Minimod like limiter
      IF(limiter.EQ.5) CALL Minmax(DUPW, DLOC,DELTA)    ! Minimod like limiter 2

      PIL = u(i,j) - 0.5D0*DELTA
      PIR = u(i,j) + 0.5D0*DELTA

      ul(i,j)  = pil
      ur(i,j)  = pir

   enddo
enddo

endif 

return
end