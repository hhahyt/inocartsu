!----------------------------------------------------------------------*

!superbee type slope limiter

subroutine superb(r, omega, delta)



implicit none

real  delta, denor, omega, phi, phir, r

phi = 0.0
if(r.ge.0.0)phi = 2.0*r
if(r.ge.0.5)phi = 1.0

if(r.ge.1.0)then
   denor = 1.0 - omega + (1.0 + omega)*r
   phir  = 2.0/denor
   phi   = min(phir, r)
   phi   = min(phi, 2.0)
endif

delta = phi*delta
return
end

!----------------------------------------------------------------------*

subroutine minmod(r, omega, delta)

! minmod type slope limiter delta

implicit none

real  delta, denor, omega, phi, phir, r

phi = 0.0
if(r.ge.0.0)phi = r

if(r.ge.1.0)then
   denor = 2.0*(1.0 - omega + (1.0 + omega)*r)
   phir  = 4.0/denor
   phi   = min(1.0, phir)
endif

delta    = phi*delta
return
end

!----------------------------------------------------------------------*

subroutine minmax(dupW, dloc, delta)
!minmax type slope limiter delta.

implicit none

real  betal, betar, delta, dloc, dupW, signo

betal = 1.0
betar = 1.0
signo = 0.5*(sign(1.0,dupW) + sign(1.0,dloc))
delta = signo*(min(betal*abs(dupW),betar*abs(dloc)))

return
end
