
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

window,20, xs=900, ys=600
i=2.4
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


interplot,  vr ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
qvel=qinter1

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
	 plots, [xs,xs],[0,max(qinter1)]

interplot,  ma ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
qrad=sqrt(rpl^2+zpl^2)
qrad=zpl
plot, qrad, qinter1,$
	xrange=[xbeg,xend],$
	;yrange=[0,0.02],$
	xtitle="Altitude (Z)",$
	title="ma"

print, size(qinter1)
qinter1=shift(qinter1,-20)
a=where(qinter1 gt 1)
alfven=a(0)+10
;oplot, sqrt(rpl^2+zpl^2), qinter1
print , alfven
print,rpl(alfven)
print,zpl(alfven)
xa=sqrt( rpl(alfven)^2 +  zpl(alfven)^2)
xa=zpl(alfven)
	 plots, [xa,xa],[0,max(qinter1)]
	 plots, [xs,xs],[0,max(qinter1)]

interplot,  mf ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
plot, qrad, qinter1,$
	xrange=[xbeg,xend],$
	;yrange=[0,0.02],$
	xtitle="Altitude (Z)",$
	title="mf"

a=where(qvel gt 0)
print, a(0)
disk=a(0)
print,rpl(disk)
print,zpl(disk)
xd= zpl(disk)

a=where(qinter1 gt 1)
print, a(0)
fast=a(0)
print,rpl(fast)
print,zpl(fast)
xf= zpl(fast)

	 plots, [xs,xs],[0,max(qinter1)]
	 plots, [xa,xa],[0,max(qinter1)]
	 plots, [xf,xf],[0,max(qinter1)]

interplot,ke ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
plot, qrad, qinter1,$
	xrange=[xbeg,xend],$
	xtitle="Altitude (Z)",$
	title="u^2/2"
	 plots, [xs,xs],[0,max(qinter1)]
	 plots, [xa,xa],[0,max(qinter1)]
	 plots, [xf,xf],[0,max(qinter1)]

qke=qinter1

interplot,enth ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
plot, qrad, qinter1,$
	xrange=[xbeg,xend],$
	;yrange=[0,0.02],$
	xtitle="Altitude (Z)",$
	title="P/rho"
	 plots, [xs,xs],[0,max(qinter1)]
	 plots, [xa,xa],[0,max(qinter1)]
	 plots, [xf,xf],[0,max(qinter1)]

qenth=qinter1


interplot,  mag ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
plot, qrad, qinter1,$
	xrange=[xbeg,xend],$
	;yrange=[0,0.02],$
	xtitle="Altitude (Z)",$
	title="Mag"
	 plots, [xs,xs],[0,max(qinter1)]
	 plots, [xa,xa],[0,max(qinter1)]
	 plots, [xf,xf],[0,max(qinter1)]
qmag=qinter1



plot, qrad,qmag, yrange=[-0.75,1.5]
oplot, sqrt(rpl^2+zpl^2),qenth
oplot, sqrt(rpl^2+zpl^2),qkep
oplot, sqrt(rpl^2+zpl^2),qke
grav=-1/sqrt(rpl^2+zpl^2)
oplot, sqrt(rpl^2+zpl^2),grav
oplot, sqrt(rpl^2+zpl^2),qke+qenth+qmag+grav

im=tvread(filename="Bernoulli1",/png,/nodialog)


!P.MULTI=0
csize=1.5
!P.CHARSIZE=csize
!P.CHARTHICK=2.0
!P.THICK=3.0



!p.font=0
set_plot,'ps'
device, /encapsulated, filename="bernoulli.eps"
device, /times


loadct,0
!P.COLOR=0
!P.BACKGROUND=255
xbeg=1e-1
xend=200
ybeg=0
yend=0.7
dist=sqrt(rpl^2+zpl^2)
dist=zpl
plot, dist,qmag, $
	xrange=[xbeg,xend],$
	psym=1,$
	xstyle=1, $
	/xlog,$
	yrange=[ybeg,yend],$
	ystyle=1, $
	ytitle="E, r="+string(i, format='(F5.1)')+', t='+string(time(nlast)/2/!PI,format='(F8.1)'),$
	xtitle="Altitude (Z)"
oplot, dist ,qkep, linestyle=0
oplot, dist ,qenth, linestyle=1
oplot, dist ,qket, linestyle=2
oplot, dist,qke, linestyle=3
oplot, dist,-qgrav, linestyle=4

items=['Mag' ]
lineitems=['Pol','Enth', 'Tor', 'Kin', '-Grav','E' ]
psym=[1]
linestyle=[0,1,2,3,4,5]


legend, [items,lineitems ], psym=[psym,0,0,0,0,0,0],line=[0,linestyle], /right, charsize=csize/2
qe=qke+qenth+qmag+qgrav
oplot, dist,qe, linestyle=5
	 plots, [xd,xd],[ybeg,yend]
	    xyouts, xd*1.1, 0.9*yend, 'S!DD!N'
	 plots, [xs,xs],[ybeg,yend]
	    xyouts, xs*1.1, 0.9*yend, 'S!DSM!N'

	 plots, [xa,xa],[ybeg,yend]
	    xyouts, xa*1.1, 0.9*yend, 'S!DA!N'
	 plots, [xf,xf],[ybeg,yend]
	    xyouts, xf*1.1, 0.9*yend, 'S!DFM!N'


device, /close
set_plot,'x'
!P.COLOR=0
!P.BACKGROUND=255


pt=20
;xyouts, pt, max(qke) , "Kinetic"
;xyouts, pt, qmag(505), "Magnetic"
;xyouts, pt, max(b), "Bernoulli"
;xyouts, pt, max(qke) , "Kinetic"
;xyouts, pt, max(qkep) , "Poloidal Kinetic"
;xyouts, pt, max(qket) , "Toroidal Kinetic"
;xyouts, pt, max(qenth) , "Enthalpy"
;xyouts, 20, max(grav)*8.5 , "Gravitational"

end
