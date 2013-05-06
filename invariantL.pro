
window,22, xs=900, ys=600
loadct,0

!P.POSITION=0
!P.MULTI=0
!P.CHARSIZE=2.0
!P.CHARSIZE=3.0
!P.CHARTHICK=2.0
!P.THICK=3.0


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
lambda=v3(0)*rr-rr*b3(0)/k
mag = -rr*b3(0)/k
matter =v3(0)*rr

omega= v3(0)/rr - k*b3(0) /rho(0)/rr

s=alog(pr(0)/rho(0)^(5./3.))

e=(5./2.)*pr(0) - 1/sqrt(rr^2 +zz^2) - 0.5*omega^2 *rr^2 + 0.5*(v1(0)^2  +v2(0)^2)$
	+0.5* (v3(0)/rr- omega)^2 * rr^2


r=x1
z=x2
dr=dx1
dz=dx2


i=2

field_line, bx, by, x1,x2,  i,0.05, rr,zz
interplot,lambda ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
plot, sqrt(rpl^2+zpl^2), qinter1/i^0.5,$
	xrange=[0,10],$
	xtitle="spherical radius, R",$
	title="!7K!3 (Cons. of ang. mom), r="+string(i)

interplot,mag ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
oplot,  sqrt(rpl^2+zpl^2), qinter1/i^0.5,$
	linestyle=1


interplot,matter ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
oplot,  sqrt(rpl^2+zpl^2), qinter1/i^0.5,$
	linestyle=2

items=['!7K!3','Mag','Matter !7x!3r!U2!N' ]
legend, items, linestyle=indgen(3),/right

im=tvread(filename='InvariantL2', /nodialog, /png)


end
