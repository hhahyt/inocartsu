
subroutine estime(sl,sm,sr,dl,ul,dr,ur)

implicit none

real  dl,ul,cl,dr,ur,cr,grav,sl,sm,sr,um,dm,al,ar,ql,qr

grav = 9.81d0


dm = (0.5*(sqrt(grav*dl)+sqrt(grav*dr))+0.25*(ul-ur))**2/grav
um = 0.5*(ul+ur)+sqrt(grav*dl)-sqrt(grav*dr)

if (dl.gt.0.0) then
   if (dm.gt.dl) then
      ql = sqrt(0.5*(dm+dl)*dm/dl/dl)
   else
      ql = 1.0
   endif
else
      ql = 0.0
endif

if (dr.gt.0.0) then
if (dm.gt.dr) then
   qr = sqrt(0.5*(dm+dr)*dm/dr/dr)
else
   qr = 1.0
endif
else
   qr = 0.0
endif

if (dl.gt.0.d0.and.dr.le.0.d0) then ! right dry
   sl = ul - sqrt(grav*dl)*ql
   sr = ul + 2.d0*sqrt(grav*dl)*ql
	
elseif (dl.le.0.d0.and.dr.gt.0.d0) then !left dry
   sl = ur - 2.d0*sqrt(grav*dr)*qr
   sr = ur + sqrt(grav*dr)*qr
elseif (dl.le.0.d0.and.dr.le.0.d0) then
   sl = 0.d0
   sr = 0.d0

else
   sl = ul-sqrt(grav*dl)*ql
   sr = ur+sqrt(grav*dr)*qr

endif

!sm = 0.0

end

