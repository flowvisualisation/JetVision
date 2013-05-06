
loadct,0


; K

den=rho(0)
pres=pr(0)

vx=v1(0)
vy=v2(0)
vz=v3(0)

bx=b1(0)
by=b2(0)
bz=b3(0)

k=rho(0)*sqrt(v1(0)^2+v2(0)^2)/sqrt(b1(0)^2+b2(0)^2)

rr=rebin(reform(x1,n1,1),n1,n1 )
zz=transpose(rr)
rad_sphere=sqrt(rr^2+zz^2)
lambda=v3(0)*rr-rr*b3(0)/k

omega= v3(0)/rr - k*b3(0) /rho(0)/rr

s=alog(pr(0)/rho(0)^(5./3.))

e=(5./2.)*pr(0) - 1/sqrt(rr^2 +zz^2) - 0.5*omega^2 *rr^2 + 0.5*(v1(0)^2  +v2(0)^2)$
	+0.5* (v3(0)/rr- omega)^2 * rr^2


r=x1
z=x2
dr=dx1
dz=dx2

ke= 0.5*(v1(0)^2  +v2(0)^2 + v3(0)^2)
kep= 0.5*(v1(0)^2  +v2(0)^2 )
ket= 0.5*(v3(0)^2  )
enth=(5./2.)*pr(0)/rho(0)
mag= -(omega)*rr*b3(0)/k
grav= -1/rad_sphere

!P.POSITION=0
!P.MULTI=[0,2,3]
!P.CHARSIZE=2.0

i=5
field_line, bx, by, x1,x2,  i,0.05, rr,zz

xbeg=0
xend=40
alf, b1,b2,b3,v1,v2,rho,pr,mf,ma,ms,mson
;interplot,  ms ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
;plot, sqrt(rpl^2+zpl^2), qinter1,$
;	xrange=[xbeg,xend],$
;	;yrange=[0,0.02],$
;	xtitle="spherical radius, R",$
;	title="ms"

interplot,  grav ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
qgrav=grav

interplot,  kep ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
qkep=qinter1

interplot,  ket ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
qket=qinter1


interplot,  ma ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
plot, sqrt(rpl^2+zpl^2), qinter1,$
	xrange=[xbeg,xend],$
	;yrange=[0,0.02],$
	xtitle="spherical radius, R",$
	title="ma"


interplot,  mf ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
plot, sqrt(rpl^2+zpl^2), qinter1,$
	xrange=[xbeg,xend],$
	;yrange=[0,0.02],$
	xtitle="spherical radius, R",$
	title="mf"


interplot,ke ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
plot, sqrt(rpl^2+zpl^2), qinter1,$
	xrange=[xbeg,xend],$
	xtitle="spherical radius, R",$
	title="u^2/2"
	 plots, [xs,xs],[0,max(qinter1)]
	 plots, [xa,xa],[0,max(qinter1)]
	 plots, [xf,xf],[0,max(qinter1)]

qke=qinter1

interplot,enth ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
plot, sqrt(rpl^2+zpl^2), qinter1,$
	xrange=[xbeg,xend],$
	;yrange=[0,0.02],$
	xtitle="spherical radius, R",$
	title="P/rho"
	 plots, [xs,xs],[0,max(qinter1)]
	 plots, [xa,xa],[0,max(qinter1)]
	 plots, [xf,xf],[0,max(qinter1)]

qenth=qinter1


interplot,  mag ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
plot, sqrt(rpl^2+zpl^2), qinter1,$
	xrange=[xbeg,xend],$
	;yrange=[0,0.02],$
	xtitle="spherical radius, R",$
	title="Mag"
	 plots, [xs,xs],[0,max(qinter1)]
	 plots, [xa,xa],[0,max(qinter1)]
	 plots, [xf,xf],[0,max(qinter1)]
qmag=qinter1



plot, sqrt(rpl^2+zpl^2),qmag, yrange=[-0.75,1.5]
oplot, sqrt(rpl^2+zpl^2),qenth
oplot, sqrt(rpl^2+zpl^2),qkep
oplot, sqrt(rpl^2+zpl^2),qke
grav=-1/sqrt(rpl^2+zpl^2)
oplot, sqrt(rpl^2+zpl^2),grav
oplot, sqrt(rpl^2+zpl^2),qke+qenth+qmag+grav

im=tvread(filename="Bernouilli1",/png,/nodialog)


!P.MULTI=0
!P.CHARSIZE=2.0
!P.CHARTHICK=2.0
!P.THICK=3.0


window,22, xs=900, ys=600

dist=sqrt(rpl^2+zpl^2)
dist=zpl
plot, dist,qmag, yrange=[-0.5,0.5],ystyle=1, $
	title="Bernouilli Invariant r="+string(i),$
	xtitle="Altitude (Z)"
oplot, dist ,qenth, linestyle=1
oplot, dist ,qkep, linestyle=2
oplot, dist ,qket, linestyle=3
oplot,  dist,qke, linestyle=4
grav=-1/dist
oplot, dist,grav, linestyle=5

items=['Mag','Enth','Pol', 'Tor', 'Kin', 'Grav','Sum' ]
legend, items, linestyle=indgen(7), /right
b=qke+qenth+qmag+grav
oplot, dist,b, linestyle=6

pt=20
;xyouts, pt, max(qke) , "Kinetic"
;xyouts, pt, qmag(505), "Magnetic"
;xyouts, pt, max(b), "Bernoulli"
;xyouts, pt, max(qke) , "Kinetic"
;xyouts, pt, max(qkep) , "Poloidal Kinetic"
;xyouts, pt, max(qket) , "Toroidal Kinetic"
;xyouts, pt, max(qenth) , "Enthalpy"
;xyouts, 20, max(grav)*8.5 , "Gravitational"

im=tvread(filename="Bernouilli4",/png,/nodialog)
end
