

!P.POSITION=0
!P.MULTI=[0,3,2]
!P.CHARSIZE=2.0
!P.CHARTHICK=2.0
!P.THICK=3.0


; K

filenum=126
filenum=nlast
pload,out=filenum
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


s=alog10(pres/den^(5./3.))



ke= 0.5*(vr^2  +vz^2 + vphi^2)
enth=(5./2.)*pres/den
mag= -(omega)*rr*bphi/k
grav= -1/rad_sphere


e= ke + enth+ mag + grav

r=x1
z=x2
dr=dx1
dz=dx2



set_plot,'ps'
device, /encapsulated, filename="invariants.eps"
i=1.6

field_line, br, bz, x1,x2,  i,0.05, rr,zz
interplot,k ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
qrad=sqrt(rpl^2+zpl^2)
qrad=zpl
xbeg=0
xend=2
plot, qrad, qinter1, yrange=[0,0.02],$
xrange=[xbeg,xend],$
	xtitle="Altitude, (Z)",$
	ytitle="k (mass loading)"

for i=2,8 do begin
field_line, br, bz, x1,x2,  i,0.05, rr,zz
interplot,k ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
qrad=sqrt(rpl^2+zpl^2)
qrad=zpl
oplot, qrad, qinter1, linestyle=i-2
endfor


i=1.6

field_line, br, bz, x1,x2,  i,0.05, rr,zz
interplot,lambda ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
qrad=sqrt(rpl^2+zpl^2)
qrad=zpl
plot, qrad, qinter1,$
xrange=[xbeg,xend],$
	yrange=[0,10],$
	title="t="+string(time(filenum)/2/!PI , format='(F8.1)'),$
	xtitle="Altitude, (Z)",$
	ytitle="l (Specific angular mom.)"

for i=2,8 do begin
field_line, br, bz, x1,x2,  i,0.05, rr,zz
interplot,lambda ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
qrad=sqrt(rpl^2+zpl^2)
qrad=zpl
oplot, qrad, qinter1, linestyle=i-2
endfor

i=1.6


field_line, br, bz, x1,x2,  i,0.05, rr,zz
interplot,omega ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
qrad=sqrt(rpl^2+zpl^2)
qrad=zpl
plot, qrad, qinter1,$
xrange=[xbeg,xend],$
	;yrange=[0,0.05],$
	xtitle="Altitude, (Z)",$
	ytitle="!7X!3 (Mag. surface rotation rate)"

for i=2,8 do begin
field_line, br, bz, x1,x2,  i,0.05, rr,zz
interplot, omega ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
qrad=sqrt(rpl^2+zpl^2)
qrad=zpl
oplot, qrad, qinter1, linestyle=i-2
endfor

i=1.6

field_line, br, bz, x1,x2,  i,0.05, rr,zz
interplot,S ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
qrad=sqrt(rpl^2+zpl^2)
qrad=zpl
plot, qrad, qinter1,$
xrange=[xbeg,xend],$
	;yrange=[0,0.05],$
	xtitle="Altitude, (Z)",$
	ytitle="S (Entropy)"

for i=2,8 do begin
field_line, br, bz, x1,x2,  i,0.05, rr,zz
interplot,S ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
qrad=sqrt(rpl^2+zpl^2)
qrad=zpl
oplot, qrad, qinter1, linestyle=i-2
endfor


i=1.6
field_line, br, bz, x1,x2,  i,0.05, rr,zz
interplot, e ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
qrad=sqrt(rpl^2+zpl^2)
qrad=zpl
plot, qrad, qinter1,$
xrange=[xbeg,xend],$
	;yrange=[-2,0],$
	xtitle="Altitude, (Z)",$
	ytitle="E (Specific Energy)"

for i=2,8 do begin
field_line, br, bz, x1,x2,  i,0.05, rr,zz
interplot, e ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
qrad=sqrt(rpl^2+zpl^2)
qrad=zpl
oplot, qrad, qinter1, linestyle=i-2
endfor


device, /close
set_plot,'x'
!P.COLOR=0
!P.BACKGROUND=255
;im=tvread(filename='Invariants', /nodialog, /png)


end
