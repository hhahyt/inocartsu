C*******************************************************************************
C........    Calculates free surface and velocity
C.....  based on position  (following expression given in Wei et al, (not JFM paper))
      subroutine solit(uc,vc,zc,c,co,t,x0,x,y,ho,alp,inc_ang,wave_type,
     -             wave_hgt,L,bf_ratio,bet,cur_level,
     -			 num_levels)
      integer wave_type,coef,cur_level,num_levels,i
      real uc,zc,c,co,t,x0,theta,loc,x,y,ho,alp,k,w,per,za,
     -   amp,g,inc_ang,vc,wave_hgt,L,bf_ratio,bet(num_levels),
     -   beta,depth,grav,delta,coef0,coef1,coef2,xp,yp,ra,
     -   coef3,c1,f3,f1,c3,error,c2,f2,a,a1,a2,b,ksi
	real L_1,L_2,k_1,k_2,xa,xb,ang,xc,yc,xo,yo,a_1,a_2

      theta=inc_ang*3.1415/180.
      loc=x*cos(theta)+y*sin(theta)
      coef=1
      g=9.81
 
      if(wave_hgt.lt.0)coef=-1

      if(wave_type.eq.1)then   ! solitary wave solution written by Khairil Irfan Sitanggang

		beta=alp
		depth=ho
		grav=g
		delta = wave_hgt/depth
		coef0 = beta+1./3
		coef1 = 2.*delta*(beta+1./3)
		coef2 = -(3.*beta+1./3+2.*beta*delta)
		coef3 = 2.*beta
! find the root of the 3rd-order polynomial 
		c1 = 40.
		f3 = 10.
		f1 = 10.
		do while(f1*f3>0.0)
		  c3 = c1
		  c1 = c3-0.1
		  f1 = coef0+coef1*c1+coef2*c1**2+coef3*c1**3
		  f3 = coef0+coef1*c3+coef2*c3**2+coef3*c3**3
		enddo
		error = 100.0
		do i=1,2000
		  c2 = 0.5*(c1+c3) 
		  f1=coef0+coef1*c1+coef2*c1**2+coef3*c1**3
		  f2=coef0+coef1*c2+coef2*c2**2+coef3*c2**3
		  f3=coef0+coef1*c3+coef2*c3**2+coef3*c3**3
		  if(f2*f1.lt.0.0)then
			c3=c2
		  elseif(f3*f2.lt.0.0)then
			c1=c2
		  endif
		  error = abs(c3-c1)
		enddo
		c2 = c2*grav*depth
		if((c2-grav*depth).lt.0.0) pause 'since c^2-g*h < 0, no solution'
		c = sqrt(c2)
		a=(c2-grav*depth)/c
		a1=depth*(c2-grav*depth)/(3.*((beta+1./3)*grav*depth-beta*c2))
		a2=-depth*(c2-grav*depth)**2/(2.*grav*depth*c2)*((beta+1./3)*
     &		grav*depth+2.*beta*c2)/((beta+1./3)*grav*depth-beta*c2)
		b=sqrt((c2-grav*depth)/(4.*((beta+1./3)*grav*depth**3-
     &		beta*depth**2.*c2)))

		ksi=loc-c2*t-x0

          if(abs(b*ksi).gt.10)then
                  zc=0
                  uc=0
          else
                  zc=a1/cosh(b*ksi)**2+a2/cosh(b*ksi)**4
                  uc=a/cosh(b*ksi)**2
          endif



CC	TEMP FOR GOM TSUNAMI -East Breaks 1D

c	xo=	1.908e5		
c	yo= 0.
c	L=25.e3
c	W=8.e3
c	a_1=160.
c	a_2=100.
c	L_1=30.e3
c	L_2=40.e3
c
c	k_1=3.1415/L_1
c	k_2=3.1415/L_2
c
c	xa=xo-L_1
c	xb=xo+L_2
c
c	ang=0.*3.1415/180.
c
c	xc=(x-xo)*cos(ang) !+(y-yo)*sin(ang)
c	yc=0. !sqrt(((x-xo)*cos(ang))**2+((y-yo)*sin(ang))**2)
c	zc=0.
c	uc=0.
c	if(xc.ge.-L_1.and.xc.le.0.)then
c	   zc=-a_1*sin(k_1*(L_1+xc)) !*exp(-((yc)/W)**2)
c	elseif(xc.ge.0.and.xc.le.L_2)then
c	   zc=a_2*sin(k_2*(xc)) !*exp(-((yc)/W)**2)
c	endif


CC	TEMP FOR GOM TSUNAMI -East Breaks 2D

c	xo=	2.2e5		
c	yo= 2.85e5
c	L=25.e3
c	W=20.e3
c	a_1=120.
c	a_2=100.
c	L_1=30.e3
c	L_2=40.e3
c
c	k_1=3.1415/L_1
c	k_2=3.1415/L_2
c
c	xa=xo-L_1
c	xb=xo+L_2
c
c	ang=-80.*3.1415/180.
c
c	xp=x-xo
c	yp=y-yo
c	ra=sqrt(xp**2.+yp**2.)
c	xc=(x-xo)*cos(ang)+(y-yo)*sin(ang)
c	yc=sqrt(ra**2.-xc**2.)
c	zc=0.
c	uc=0.
c	if(xc.ge.-L_1.and.xc.le.0.)then
c	   zc=-a_1*sin(k_1*(L_1+xc))*exp(-((yc)/W)**2.)
c	elseif(xc.ge.0.and.xc.le.L_2)then
c	   zc=a_2*sin(k_2*(xc))*exp(-((yc)/W)**2.)
c	endif


CC	TEMP FOR GOM TSUNAMI -Florida Escarp 1D

c	xo=	1.28e6		
c	yo= 0.
c	L=25.e3
c	W=8.e3
c	a_1=-100.
c	a_2=-150.
c	L_1=30.e3
c	L_2=20.e3
c
c	k_1=3.1415/L_1
c	k_2=3.1415/L_2
c
c	xa=xo-L_1
c	xb=xo+L_2
c
c	ang=0.*3.1415/180.
c
c	xc=(x-xo)*cos(ang) !+(y-yo)*sin(ang)
c	yc=0. !sqrt(((x-xo)*cos(ang))**2+((y-yo)*sin(ang))**2)
c	zc=0.
c	uc=0.
c	if(xc.ge.-L_1.and.xc.le.0.)then
c	   zc=-a_1*sin(k_1*(L_1+xc)) !*exp(-((yc)/W)**2)
c	elseif(xc.ge.0.and.xc.le.L_2)then
c	   zc=a_2*sin(k_2*(xc)) !*exp(-((yc)/W)**2)
c	endif
c

CC	TEMP FOR GOM TSUNAMI - Florida Escarp 2D

c	xo=	1.32e6		
c	yo= 0.8e5
c	L=25.e3
c	W=60.e3
c	a_1=-100.
c	a_2=-150.
c	L_1=30.e3
c	L_2=20.e3
c
c	k_1=3.1415/L_1
c	k_2=3.1415/L_2
c
c	xa=xo-L_1
c	xb=xo+L_2
c
c	ang=0.*3.1415/180.
c
c	xp=x-xo
c	yp=y-yo
c	ra=sqrt(xp**2.+yp**2.)
c	xc=(x-xo)*cos(ang)+(y-yo)*sin(ang)
c	yc=sqrt(ra**2.-xc**2.)
c	zc=0.
c	uc=0.
c	if(xc.ge.-L_1.and.xc.le.0.)then
c	   zc=-a_1*sin(k_1*(L_1+xc))*exp(-((yc)/W)**2.)
c	elseif(xc.ge.0.and.xc.le.L_2)then
c	   zc=a_2*sin(k_2*(xc))*exp(-((yc)/W)**2.)
c	endif

CC	TEMP FOR GOM TSUNAMI -Campeche Escarp 1D

c	xo=9.25e5		
c	yo= 0.
c	L=25.e3
c	W=8.e3
c	a_1=-100.
c	a_2=-150.
c	L_1=30.e3
c	L_2=20.e3
c
c	k_1=3.1415/L_1
c	k_2=3.1415/L_2
c
c	xa=xo-L_1
c	xb=xo+L_2
c
c	ang=0.*3.1415/180.
c
c	xc=(x-xo)*cos(ang) !+(y-yo)*sin(ang)
c	yc=0. !sqrt(((x-xo)*cos(ang))**2+((y-yo)*sin(ang))**2)
c	zc=0.
c	uc=0.
c	if(xc.ge.-L_1.and.xc.le.0.)then
c	   zc=-a_1*sin(k_1*(L_1+xc)) !*exp(-((yc)/W)**2)
c	elseif(xc.ge.0.and.xc.le.L_2)then
c	   zc=a_2*sin(k_2*(xc)) !*exp(-((yc)/W)**2)
c	endif




cCC TEMP
c		  tmp=200.
c            loc=(x-tmp*1000)/1000.*2*3.14/tmp
c            zc=5.*sin(loc)*exp(-(loc)**2)
c            uc=0
cCC END TEMP


CC	TEMP SRI LANKA TSUNAMI
cCC!!!!!!!
c	xo=	3.12e3		! Indian Ocean
c	yo= 1.5e3
c	L=	1.0e3
c	ao=	-8.
c
c	ang1=-5.*3.1415/180.
c	ang2=25.*3.1415/180.
cCC!!!!!!!!
c
c	xo=	750	   ! Cascadia
c	yo= 1.5e3
c	L=	1.0e3
c	ao=	-8.
c
c	ang1=15.*3.1415/180.
c	ang2=5.*3.1415/180.
c
c	ko=	0.01			
c	k_1=0.01			
c					
c	x_2=xo-120		
c	k_2=0.05			
c	a_2=4.5		
c
c	ya=yo-L/2.
c	yb=yo+L/2.
c
c	xc=x/1000.
c	yc=y/1000.
c	ang=ang2+(ang1-ang2)/(yb-ya)*(yc-ya)
c	loc=(xc-xo)*cos(ang)+(yc-yo)*sin(ang)
c	loc2=(xc-x_2)*cos(ang)+(yc-yo)*sin(ang) 
c
c	if(yc.ge.ya.and.yc.le.yb)then    
c		zc=ao*exp(-(loc**2*ko**2))*sin(k_1*loc)+a_2*exp(-(loc2**2*k_2**2))
c	else
c		zc=0.
c	endif

c	xo=	6200	   ! atlantic ocean
c	yo= 1490
c	L=	1.0e3
c	ao=	1000.
c
c	ko=	0.1	
c
c	xc=x/1000.
c	yc=y/1000.
c	loc=sqrt((xc-xo)**2+(yc-yo)**2)
c
c	zc=ao*exp(-(loc**2*ko**2))
c
c	uc=0.

CC	END TEMP		


CC	TEMP FOR GOM TSUNAMI !!  ITS TSUNAMI FEVER BABY

c	xo=	4384		
c	yo= 3182
c	L=		1000.
c	W=		1000.
c	a_1=	0.
c	a_2=	50.
c	L_1=	1000.
c	L_2=	1000.
c
c	k_1=3.1415/L_1
c	k_2=3.1415/L_2
c
c	xa=xo-L_1
c	xb=xo+L_2
c
c	ang=0.*3.1415/180.
c
c	xc=(x-xo)*cos(ang)+(y-yo)*sin(ang)
c	yc=sqrt(((x-xo)*cos(ang))**2+((y-yo)*sin(ang))**2)
c	zc=0.
c	uc=0.
c	if(xc.ge.-L_1.and.xc.le.0.)then
c	   zc=-a_1*sin(k_1*(L_1+xc))*exp(-((yc)/W)**2)
c	elseif(xc.ge.0.and.xc.le.L_2)then
c	   zc=a_2*sin(k_2*(xc))*exp(-((yc)/W)**2)
c	endif


CC	TEMP FOR STOREGGA TSUNAMI !!  
c	xo=	600.e3		
c	yo= 620.e3
c	L=		150.e3
c	W=		75.e3
c	a_1=	10.
c	a_2=	20.
c	L_1=	100.e3
c	L_2=	100.e3
c
c	k_1=3.1415/L_1
c	k_2=3.1415/L_2
c
c	xa=xo-L_1
c	xb=xo+L_2
c
c	ang=-45.*3.1415/180.
c
c	xc=(x-xo)*cos(ang)+(y-yo)*sin(ang)
c	yc=sqrt(((x-xo)*cos(ang))**2+((y-yo)*sin(ang))**2)
c	zc=0
c	if(xc.ge.-L_1.and.xc.le.0.)then
c	   zc=a_1*sin(k_1*(L_1+xc))*exp(-((yc)/W)**2)
c	elseif(xc.ge.0.and.xc.le.L_2)then
c	   zc=-a_2*sin(k_2*(xc))*exp(-((yc)/W)**2)
c	endif

c	uc=0.
CC	TEMP FOR GOM TSUNAMI - Miss Canyon 2D

c	xo=	8.3e5		
c	yo= 0.25e5
c	L=25.e3
c	W=30.e3
c	a_1=200.
c	a_2=200.
c	L_1=60.e3
c	L_2=100.e3
c
c	k_1=3.1415/L_1
c	k_2=3.1415/L_2
c
c	xa=xo-L_1
c	xb=xo+L_2
c
c	ang=-55.*3.1415/180.
c
c	xp=x-xo
c	yp=y-yo
c	ra=sqrt(xp**2.+yp**2.)
c	xc=(x-xo)*cos(ang)+(y-yo)*sin(ang)
c	yc=sqrt(ra**2.-xc**2.)
c	zc=0.
c	uc=0.
c	if(xc.ge.-L_1.and.xc.le.0.)then
c	   zc=-a_1*sin(k_1*(L_1+xc))*exp(-((yc)/W)**2.)
c	elseif(xc.ge.0.and.xc.le.L_2)then
c	   zc=a_2*sin(k_2*(xc))*exp(-((yc)/W)**2.)
c	endif


CC


      elseif(wave_type.eq.2)then
             c=2.5814
             if(abs(2*(loc-c*t-x0)).gt.10)then
                  zc=0
                  uc=0
            else

                  zc=0.2324/cosh(2.0119*(loc-c*t-x0))+
     -                  0.02719/cosh(5.769*(loc-c*t-x0))

                  uc=0.1485/cosh(1.7681*(loc-c*t-x0))+
     -                  0.680265/(cosh(0.8763*(loc-c*t-x0)))**(2.6432)
             endif
      elseif(wave_type.eq.3)then
            c=2.4875
            if(abs(2*(loc-c*t-x0)).gt.10)then
                  zc=0
                  uc=0
            else

                  zc=0.10726/(cosh(1.0539*(loc-c*t-x0)))**2.+
     -                   0.088/(cosh(3.2248*(loc-c*t-x0)))**(0.5809)

                  uc=0.75346/(cosh(1.006145*(loc-c*t-x0)))**2.-
     -                   0.0599/(cosh(.77964*(loc-c*t-x0)))**(6.542)

             endif
      elseif(wave_type.eq.5)then
		  k=2*3.1415/L
		  w=sqrt(g*k*tanh(k*ho))
		  per=2*3.1415/w
		  amp=wave_hgt/2.*tanh(t/per)**4.

            zc=amp*sin(k*loc-w*t)
            
		  za=bet(cur_level)*ho
		  
		  uc=amp*g*k/w*cosh(k*(za+ho))/cosh(k*ho)*sin(k*loc-w*t)
c                  uc=-zc*co/ho

      elseif(wave_type.eq.7)then
            if(x.lt.x0-L/2.or.x.gt.x0+L/2)then
                  zc=0
                  uc=0
            elseif(x.le.x0)then
                  k=2*3.1415/L
                  loc=x-x0-L/2.
                  zc=-wave_hgt*sin(k*loc)
                  uc=-zc*co/ho
            else
                  k=2*3.1415/L
                  loc=x-x0-L/2.
                  zc=-wave_hgt*bf_ratio*sin(k*loc)
                  uc=zc*co/ho
            endif
      endif      

      vc=uc*sin(theta)
      uc=uc*cos(theta)




      return
      end      ! end subroutine
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      subroutine cnoidal(ep,mk2,x,y,t,zc,uc,vc,yt,yt2,rp1,r2,
     -           CK,h,L,per,co,inc_ang,wave_type,za)
      integer wave_type
      real ep,mk2,x,t,zc,uc,CK,h,L,per,theta,loc,y,k,w,
     -     uu,emmc,sn,cn,dn,yt,yt2,rp1,r2,co,inc_ang,vc,za

      theta=inc_ang*3.1415/180.
      loc=x*cos(theta)+y*sin(theta)

      if(wave_type.eq.4)then
            uu=(2*CK*(loc/L-t/per))    
            emmc=1-(mk2)
            call sncndn(uu,emmc,sn,cn,dn)
            
            uc=co*(ep*(cn**2-yt)+ep**2*(rp1+r2*cn**2-(cn**2.)**2./4.))
            zc=h*(ep*(cn**2-yt)-ep**2*(3./4.*cn**2*(1-cn**2)+yt2))
      elseif(wave_type.eq.5)then
             k=2*3.1415/L
             w=2*3.1415/per

             zc=ep/2.*h*cos(k*loc-w*t)
             uc=zc*9.81*k/w*cosh(k*(za+h))/cosh(k*h)
      endif

      vc=uc*sin(theta)
      uc=uc*cos(theta)

      return
      end
c

c.....Subroutine for elliptic function
c
      SUBROUTINE sncndn(uu,emmc,sn,cn,dn)
      REAL cn,dn,emmc,sn,uu,CA
      PARAMETER (CA=.00003)
      INTEGER i,ii,l
      REAL a,b,c,d,emc,u,em(13),en(13)
      LOGICAL bo
      emc=emmc
      u=uu
      if(emc.ne.0.)then
        bo=(emc.lt.0.)
        if(bo)then
          d=1.-emc
          emc=-emc/d
          d=sqrt(d)
          u=d*u
        endif
        a=1.
        dn=1.
        do 11 i=1,13
          l=i
          em(i)=a
          emc=sqrt(emc)
          en(i)=emc
          c=0.5*(a+emc)
          if(abs(a-emc).le.CA*a)goto 1
          emc=a*emc
          a=c
11      continue
1        u=c*u
        sn=sin(u)
        cn=cos(u)
        if(sn.eq.0.)goto 2
        a=cn/sn
        c=a*c
        do 12 ii=l,1,-1
          b=em(ii)
          a=c*a
          c=dn*c
          dn=(en(ii)+a)/(b+a)
          a=c/b
12      continue
        a=1./sqrt(c**2+1.)
        if(sn.lt.0.)then
          sn=-a
        else
          sn=a
        endif
        cn=c*sn
2       if(bo)then
          a=dn
          dn=cn
          cn=a
          sn=sn/d
        endif
      else
        cn=1./cosh(u)
        dn=cn
        sn=tanh(u)
      endif

      return
      END

      SUBROUTINE cnoidal_rf(x,y,z,rf)
      REAL rf,x,y,z,ERRTOL,TINY,BIG,THIRD,C1,C2,C3,C4
      PARAMETER (ERRTOL=.0008,TINY=1.5e-38,BIG=3.E37,THIRD=1./3.,
     *C1=1./24.,C2=.1,C3=3./44.,C4=1./14.)
      REAL alamb,ave,delx,dely,delz,e2,e3,sqrtx,sqrty,sqrtz,xt,yt,zt
      if(min(x,y,z).lt.0..or.min(x+y,x+z,y+z).lt.TINY.or.max(x,y,
     *z).gt.BIG)pause 'invalid arguments in rf'
      xt=x
      yt=y
      zt=z
1     continue
        sqrtx=sqrt(xt)
        sqrty=sqrt(yt)
        sqrtz=sqrt(zt)
        alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
        xt=.25*(xt+alamb)
        yt=.25*(yt+alamb)
        zt=.25*(zt+alamb)
        ave=THIRD*(xt+yt+zt)
        delx=(ave-xt)/ave
        dely=(ave-yt)/ave
        delz=(ave-zt)/ave
      if(max(abs(delx),abs(dely),abs(delz)).gt.ERRTOL)goto 1
      e2=delx*dely-delz**2
      e3=delx*dely*delz
      rf=(1.+(C1*e2-C2-C3*e3)*e2+C4*e3)/sqrt(ave)
      return
      END
c....
c....
c....
      SUBROUTINE cnoidal_rd(x,y,z,rd)
      REAL rd,x,y,z,ERRTOL,TINY,BIG,C1,C2,C3,C4,C5,C6
      PARAMETER (ERRTOL=.05,TINY=1.e-25,BIG=4.5E21,C1=3./14.,C2=1./6.,
     *C3=9./22.,C4=3./26.,C5=.25*C3,C6=1.5*C4)
      REAL alamb,ave,delx,dely,delz,ea,eb,ec,ed,ee,fac,sqrtx,sqrty,
     *sqrtz,sum,xt,yt,zt
      if(min(x,y).lt.0..or.min(x+y,z).lt.TINY.or.max(x,y,
     *z).gt.BIG)pause 'invalid arguments in rd'
      xt=x
      yt=y
      zt=z
      sum=0.
      fac=1.
1     continue
        sqrtx=sqrt(xt)
        sqrty=sqrt(yt)
        sqrtz=sqrt(zt)
        alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
        sum=sum+fac/(sqrtz*(zt+alamb))
        fac=.25*fac
        xt=.25*(xt+alamb)
        yt=.25*(yt+alamb)
        zt=.25*(zt+alamb)
        ave=.2*(xt+yt+3.*zt)
        delx=(ave-xt)/ave
        dely=(ave-yt)/ave
        delz=(ave-zt)/ave
      if(max(abs(delx),abs(dely),abs(delz)).gt.ERRTOL)goto 1
      ea=delx*dely
      eb=delz*delz
      ec=ea-eb
      ed=ea-6.*eb
      ee=ed+ec+ec
      rd=3.*sum+fac*(1.+ed*(-C1+C5*ed-C6*delz*ee)+delz*(C2*ee+delz*(-C3*
     *ec+delz*C4*ea)))/(ave*sqrt(ave))
      return
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


C**************************************************************************
C..........  Calculates UU from initial velocity 
       subroutine u2UU(nx,ny,nt,n,endx,endy,u,UU,h,hxx,b1,b2,dx,
     -                      v,VV,hyy,dy,hx,hy,
     -                 cur_level,num_levels,zeta,z_alp,nonlin_ind,dim)

      use mainvar_module, only:numerical_scheme
      implicit none
      integer n,nx,ny,nt,endy,endx,i,j,nonlin_ind,
     -            cur_level,num_levels,dim
      real u(nx,ny,nt,num_levels),
     -          UU(nx,ny,nt,num_levels),
     -          h(nx,ny,nt),hxx(nx,ny),
     -          b1(nx,ny,num_levels),
     -          b2(nx,ny,num_levels),dx,dudxx,dvdyy,
     -          v(nx,ny,nt,num_levels),
     -          VV(nx,ny,nt,num_levels),
     -          hyy(nx,ny),dy,tmp,
     -          hx(nx,ny),hy(nx,ny),dudx,dvdy,
     -          z,dzdx,dzdy,z1,eta,zeta(nx,ny,nt,num_levels),
     -          z_alp(nx,ny,nt,num_levels),detadx,detady


      do j=1,endy
       do i=1,endx
          if(i.eq.1)then
            dudx=(-3.*u(i+4,j,n,cur_level)+
     -                16.*u(i+3,j,n,cur_level)-
     -                36.*u(i+2,j,n,cur_level)+
     -                48.*u(i+1,j,n,cur_level)-
     -                25.*u(i,j,n,cur_level))/(12.*dx)
               dudxx=0
     
          elseif(i.eq.2)then
            dudx=(u(i+3,j,n,cur_level)-
     -                6.*u(i+2,j,n,cur_level)+
     -                18.*u(i+1,j,n,cur_level)-
     -                10.*u(i,j,n,cur_level)-
     -                3.*u(i-1,j,n,cur_level))/(12.*dx)          

               dudxx=(u(i+1,j,n,cur_level)-
     -                2.*u(i,j,n,cur_level)+
     -                u(i-1,j,n,cur_level))/(dx**2)

          elseif(i.eq.endx)then
            dudx=-(-3.*u(i-4,j,n,cur_level)+
     -                16.*u(i-3,j,n,cur_level)-
     -                36.*u(i-2,j,n,cur_level)+
     -                48.*u(i-1,j,n,cur_level)-
     -                25.*u(i,j,n,cur_level))/(12.*dx)

               dudxx=0

          elseif(i.eq.endx-1)then
            dudx=-(u(i-3,j,n,cur_level)-
     -                6.*u(i-2,j,n,cur_level)+
     -                18.*u(i-1,j,n,cur_level)-
     -                10.*u(i,j,n,cur_level)-
     -                3.*u(i+1,j,n,cur_level))/(12.*dx)       

               dudxx=(u(i+1,j,n,cur_level)-
     -                2.*u(i,j,n,cur_level)+
     -                u(i-1,j,n,cur_level))/(dx**2)
          else
            dudx=(-u(i+2,j,n,cur_level)+
     -                8.*u(i+1,j,n,cur_level)-
     -                8.*u(i-1,j,n,cur_level)+
     -                u(i-2,j,n,cur_level))/(12.*dx)

               dudxx=(-(u(i+2,j,n,cur_level))+
     -                16.*(u(i+1,j,n,cur_level))-
     -                30.*(u(i,j,n,cur_level))+
     -                16.*(u(i-1,j,n,cur_level))-
     -                (u(i-2,j,n,cur_level)))/(12.*dx**2)
          endif

          if(dim.eq.2)then
            if(j.eq.1)then
                dvdy=(-3.*v(i,j+4,n,cur_level)+
     -                  16.*v(i,j+3,n,cur_level)-
     -                  36.*v(i,j+2,n,cur_level)+
     -                  48.*v(i,j+1,n,cur_level)-
     -                  25.*v(i,j,n,cur_level))/(12.*dy)

                dvdyy=0
            elseif(j.eq.2)then
                dvdy=(v(i,j+3,n,cur_level)-
     -                  6.*v(i,j+2,n,cur_level)+
     -                  18.*v(i,j+1,n,cur_level)-
     -                  10.*v(i,j,n,cur_level)-
     -                  3.*v(i,j-1,n,cur_level))/(12.*dy)      
     
                dvdyy=(v(i,j+1,n,cur_level)-
     -                  2.*v(i,j,n,cur_level)+
     -                  v(i,j-1,n,cur_level))/(dy**2)   
           
            elseif(j.eq.endy)then
                dvdy=-(-3.*v(i,j-4,n,cur_level)+
     -                  16.*v(i,j-3,n,cur_level)-
     -                  36.*v(i,j-2,n,cur_level)+
     -                  48.*v(i,j-1,n,cur_level)-
     -                  25.*v(i,j,n,cur_level))/(12.*dy)

                dvdyy=0

               elseif(j.eq.endy-1)then
                dvdy=-(v(i,j-3,n,cur_level)-
     -                  6.*v(i,j-2,n,cur_level)+
     -                  18.*v(i,j-1,n,cur_level)-
     -                  10.*v(i,j,n,cur_level)-
     -                  3.*v(i,j+1,n,cur_level))/(12.*dy)

                dvdyy=(v(i,j+1,n,cur_level)-
     -                  2.*v(i,j,n,cur_level)+
     -                  v(i,j-1,n,cur_level))/(dy**2)
            else
                dvdy=(-v(i,j+2,n,cur_level)+
     -                  8.*v(i,j+1,n,cur_level)-
     -                  8.*v(i,j-1,n,cur_level)+
     -                  v(i,j-2,n,cur_level))/(12.*dy)

                dvdyy=(-(v(i,j+2,n,cur_level))+
     -                  16.*(v(i,j+1,n,cur_level))-
     -                  30.*(v(i,j,n,cur_level))+
     -                  16.*(v(i,j-1,n,cur_level))-
     -                  (v(i,j-2,n,cur_level)))/(12.*dy**2)
            endif
          endif

          dzdy=0
          if(nonlin_ind.eq.2)then
            z=zeta(i,j,n,1)
            if(i.eq.1)then
                dzdx=0
            elseif(i.eq.2)then
                dzdx=(zeta(i+3,j,n,cur_level)-
     -                6.*zeta(i+2,j,n,cur_level)+
     -                18.*zeta(i+1,j,n,cur_level)-
     -                10.*zeta(i,j,n,cur_level)-
     -                3.*zeta(i-1,j,n,cur_level))/(12.*dx)          
            elseif(i.eq.endx)then
                dzdx=0
            elseif(i.eq.endx-1)then
                dzdx=-(zeta(i-3,j,n,cur_level)-
     -                6.*zeta(i-2,j,n,cur_level)+
     -                18.*zeta(i-1,j,n,cur_level)-
     -                10.*zeta(i,j,n,cur_level)-
     -                3.*zeta(i+1,j,n,cur_level))/(12.*dx)       
            else
                dzdx=(-zeta(i+2,j,n,cur_level)+
     -                8.*zeta(i+1,j,n,cur_level)-
     -                8.*zeta(i-1,j,n,cur_level)+
     -                zeta(i-2,j,n,cur_level))/(12.*dx)
            endif

            if(dim.eq.2)then
                if(j.eq.1)then
                  dzdy=0
                elseif(j.eq.2)then
                  dvdy=(zeta(i,j+3,n,cur_level)-
     -                  6.*zeta(i,j+2,n,cur_level)+
     -                  18.*zeta(i,j+1,n,cur_level)-
     -                  10.*zeta(i,j,n,cur_level)-
     -                  3.*zeta(i,j-1,n,cur_level))/(12.*dy)      
                elseif(j.eq.endy)then
                  dzdy=0
                  elseif(j.eq.endy-1)then
                  dzdy=-(zeta(i,j-3,n,cur_level)-
     -                  6.*zeta(i,j-2,n,cur_level)+
     -                  18.*zeta(i,j-1,n,cur_level)-
     -                  10.*zeta(i,j,n,cur_level)-
     -                  3.*zeta(i,j+1,n,cur_level))/(12.*dy)
                else
                  dzdy=(-zeta(i,j+2,n,cur_level)+
     -                  8.*zeta(i,j+1,n,cur_level)-
     -                  8.*zeta(i,j-1,n,cur_level)+
     -                  zeta(i,j-2,n,cur_level))/(12.*dy)
                endif
            endif
          else
            z=0
            dzdx=0
            dzdy=0
          endif

          z1=z_alp(i,j,n,1)

          if(num_levels.eq.1)then
           UU(i,j,n,cur_level)=u(i,j,n,cur_level)*
     -          (1+b2(i,j,cur_level)*
     -          h(i,j,n)*hxx(i,j))+
     -           dudxx*h(i,j,n)**2.*
     -          (b1(i,j,cur_level)+b2(i,j,cur_level))+
     -           b2(i,j,cur_level)*h(i,j,n)*
     -           dudx*hx(i,j)

           if(dim.eq.2)then
           VV(i,j,n,cur_level)=
     -           v(i,j,n,cur_level)*
     -          (1+hyy(i,j)*(z1-z)-dzdy*hy(i,j))+
     -          0.5*(z1**2.-z**2.+2.*z1*h(i,j,n)-
     -          2.*z*h(i,j,n))*dvdyy+
     -          (2.*hy(i,j)*(z1-z)-dzdy*z-
     -          dzdy*h(i,j,n))*dvdy
           endif
          else
            eta=zeta(i,j,n,2)
            if(i.eq.1)then
                detadx=0
            elseif(i.eq.2)then
                detadx=(zeta(i+3,j,n,2)-
     -                6.*zeta(i+2,j,n,2)+
     -                18.*zeta(i+1,j,n,2)-
     -                10.*zeta(i,j,n,2)-
     -                3.*zeta(i-1,j,n,2))/(12.*dx)          
            elseif(i.eq.endx)then
                detadx=0
            elseif(i.eq.endx-1)then
                detadx=-(zeta(i-3,j,n,2)-
     -                6.*zeta(i-2,j,n,2)+
     -                18.*zeta(i-1,j,n,2)-
     -                10.*zeta(i,j,n,2)-
     -                3.*zeta(i+1,j,n,2))/(12.*dx)       
            else
                detadx=(-zeta(i+2,j,n,2)+
     -                8.*zeta(i+1,j,n,2)-
     -                8.*zeta(i-1,j,n,2)+
     -                zeta(i-2,j,n,2))/(12.*dx)
            endif

            if(dim.eq.2)then
                if(j.eq.1)then
                  detady=0
                elseif(j.eq.2)then
                  detady=(zeta(i,j+3,n,2)-
     -                  6.*zeta(i,j+2,n,2)+
     -                  18.*zeta(i,j+1,n,2)-
     -                  10.*zeta(i,j,n,2)-
     -                  3.*zeta(i,j-1,n,2))/(12.*dy)      
                elseif(j.eq.endy)then
                  detady=0
                  elseif(j.eq.endy-1)then
                  detady=-(zeta(i,j-3,n,2)-
     -                  6.*zeta(i,j-2,n,2)+
     -                  18.*zeta(i,j-1,n,2)-
     -                  10.*zeta(i,j,n,2)-
     -                  3.*zeta(i,j+1,n,2))/(12.*dy)
                else
                  detady=(-zeta(i,j+2,n,2)+
     -                  8.*zeta(i,j+1,n,2)-
     -                  8.*zeta(i,j-1,n,2)+
     -                  zeta(i,j-2,n,2))/(12.*dy)
                endif
            endif

            UU(i,j,n,cur_level)=
     -            u(i,j,n,cur_level)+
     -            0.5*(z1**2.-2.*z1*eta-z**2+2.*z*eta)*dudxx+
     -            (eta*dzdx+z*detadx-z1*detadx-z*dzdx)*dudx

            if(dim.eq.2)then
              VV(i,j,n,cur_level)=
     -            v(i,j,n,cur_level)+
     -            0.5*(z1**2.-2.*z1*eta-z**2+2.*z*eta)*dvdyy+
     -            (eta*dzdy+z*detady-z1*detady-z*dzdy)*dvdy
            endif

          endif
       enddo
      enddo

	  do j=1,endy
		do i=1,endx

			tmp=zeta(i,j,n,1)+
     -                  h(i,j,n)

            UU(i,j,n,cur_level)=UU(i,j,n,cur_level)*tmp

            if(dim.eq.2)then
              VV(i,j,n,cur_level)=VV(i,j,n,cur_level)*tmp
            endif
		enddo
        enddo		

      return
      end      ! end subroutine
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

