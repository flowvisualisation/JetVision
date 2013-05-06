
loadct,0


; K

den=rho(0)
pres=pr(0)

vr=v1(0)
vz=v2(0)
vphi=v3(0)

br=b1(0)
bz=b2(0)
bphi=b3(0)

k= den*sqrt(vr^2+vz^2)/sqrt(br^2+bz^2)


rr=rebin(reform(x1,n1,1),n1,n2 )
zz=rebin(reform(x2,n2,1),n2,n1 )
zz=transpose(zz)



rad_sphere=sqrt(rr^2+zz^2)

lambda=vphi*rr-rr*bphi/k

omega= vphi/rr - k*bphi /den/rr



r=x1
z=x2
dr=dx1
dz=dx2

ke= 0.5*(vr^2  +vz^2 + vphi^2)
kep= 0.5*(vr^2  +vz^2 )
ket= 0.5*(vphi^2  )
enth=(5./2.)*pres/den
mag= -(omega)*rr*bphi/k
grav= -1/rad_sphere


!P.POSITION=0
!P.MULTI=[0,2,3]
!P.CHARSIZE=2.0

;window,20, xs=900, ys=600
i=2.0
field_line, br, bz, x1,x2,  i,0.05, rr,zz

xbeg=0
xend=40
;alf, b1,b2,b3,v1,v2,rho,pr,mf,ma,ms,mson
posalf, b1,b2,b3,v1,v2,rho,pr,mf,ma,ms,mson
;interplot,  ms ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
;plot, sqrt(rpl^2+zpl^2), qinter1,$
;	xrange=[xbeg,xend],$
;	;yrange=[0,0.02],$
;	title="ms"

interplot,  grav ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
qgrav=qinter1

interplot,  kep ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
qkep=qinter1

interplot,  ket ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
qket=qinter1

interplot,  ms ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
a=where(qinter1 gt 1)
slow=a(0)
;slow=0
xs=sqrt( rpl(slow)^2 +  zpl(slow)^2)
xs=zpl(slow)
	 ;plots, [xs,xs],[0,max(qinter1)]

interplot,  ma ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
qrad=sqrt(rpl^2+zpl^2)
qrad=zpl
;plot, qrad, qinter1,$
;	xrange=[xbeg,xend],$
;	;yrange=[0,0.02],$
;	xtitle="Altitude (Z)",$
;	title="ma"

print, size(qinter1)
qinter1=shift(qinter1,-20)
a=where(qinter1 gt 1)
alfven=a(0)+20
;oplot, sqrt(rpl^2+zpl^2), qinter1
print , alfven
print,rpl(alfven)
print,zpl(alfven)
xf=sqrt( rpl(alfven)^2 +  zpl(alfven)^2)
xf=zpl(alfven)
;	 plots, [xf,xf],[0,max(qinter1)]
;	 plots, [xs,xs],[0,max(qinter1)]
;
interplot,  mf ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
;plot, qrad, qinter1,$
;	xrange=[xbeg,xend],$
;	;yrange=[0,0.02],$
;	xtitle="Altitude (Z)",$
;	title="mf"

a=where(qinter1 gt 1)
print, a(0)
fast=a(0)
print,rpl(fast)
print,zpl(fast)
xa= zpl(fast)
;	 plots, [xs,xs],[0,max(qinter1)]
;	 plots, [xa,xa],[0,max(qinter1)]
;	 plots, [xf,xf],[0,max(qinter1)]

interplot,ke ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
;plot, qrad, qinter1,$
;	xrange=[xbeg,xend],$
;	xtitle="Altitude (Z)",$
;	title="u^2/2"
;	 plots, [xs,xs],[0,max(qinter1)]
;	 plots, [xa,xa],[0,max(qinter1)]
;	 plots, [xf,xf],[0,max(qinter1)]

qke=qinter1

interplot,enth ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
;plot, qrad, qinter1,$
;	xrange=[xbeg,xend],$
;	;yrange=[0,0.02],$
;	xtitle="Altitude (Z)",$
;	title="P/rho"
;	 plots, [xs,xs],[0,max(qinter1)]
;	 plots, [xa,xa],[0,max(qinter1)]
;	 plots, [xf,xf],[0,max(qinter1)]

qenth=qinter1


interplot,  mag ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
;plot, qrad, qinter1,$
;	xrange=[xbeg,xend],$
;	;yrange=[0,0.02],$
;	xtitle="Altitude (Z)",$
;	title="Mag"
;	 plots, [xs,xs],[0,max(qinter1)]
;	 plots, [xa,xa],[0,max(qinter1)]
;	 plots, [xf,xf],[0,max(qinter1)]
qmag=qinter1



;plot, qrad,qmag, yrange=[-0.75,1.5]
;oplot, sqrt(rpl^2+zpl^2),qenth
;oplot, sqrt(rpl^2+zpl^2),qkep
;oplot, sqrt(rpl^2+zpl^2),qke
;grav=-1/sqrt(rpl^2+zpl^2)
;oplot, sqrt(rpl^2+zpl^2),grav
;oplot, sqrt(rpl^2+zpl^2),qke+qenth+qmag+grav

;im=tvread(filename="Bernoulli1",/png,/nodialog)


!P.MULTI=0
!P.CHARSIZE=2.0
!P.CHARTHICK=2.0
!P.THICK=3.0


window,21, xs=900, ys=600

xbeg=1e-2
xend=80
ybeg=0
yend=1.4
dist=sqrt(rpl^2+zpl^2)
dist=zpl
plot, dist,qmag, $
	xrange=[xbeg,xend],$
	xstyle=1, $
	/xlog,$
	yrange=[ybeg,yend],$
	ystyle=1, $
	title="Bernoulli Invariant r="+string(i),$
	xtitle="Altitude (Z)"
oplot, dist ,qenth, linestyle=1
oplot, dist ,qkep, linestyle=2
oplot, dist ,qket, linestyle=3
oplot, dist,qke, linestyle=4
oplot, dist,-qgrav, linestyle=5

items=['Mag','Enth','Pol', 'Tor', 'Kin', '-Grav','Tot' ]
legend, items, linestyle=indgen(7), /right
qe=qke+qenth+qmag+qgrav
oplot, dist,qe, linestyle=6
	 plots, [xs,xs],[ybeg,yend]
	 plots, [xa,xa],[ybeg,yend]
	 plots, [xf,xf],[ybeg,yend]

pt=20
;xyouts, pt, max(qke) , "Kinetic"
;xyouts, pt, qmag(505), "Magnetic"
;xyouts, pt, max(b), "Bernoulli"
;xyouts, pt, max(qke) , "Kinetic"
;xyouts, pt, max(qkep) , "Poloidal Kinetic"
;xyouts, pt, max(qket) , "Toroidal Kinetic"
;xyouts, pt, max(qenth) , "Enthalpy"
;xyouts, 20, max(grav)*8.5 , "Gravitational"

im=tvread(filename="Bernoulli2",/png,/nodialog)
end
