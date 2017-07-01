
window,22, xs=900, ys=600
loadct,0

!P.POSITION=0
!P.MULTI=[0,3,2]
!P.CHARSIZE=2.0
!P.POSITION=0
!P.CHARSIZE=3.0
!P.CHARTHICK=2.0
!P.THICK=3.0


; K

den=rho(0)
pres=pr(0)

vx=v1(0)
vy=v2(0)
vz=v3(0)

br=b1(0)
bz=b2(0)
bphi=b3(0)


rr=rebin(reform(x1,n1,1),n1,n1 )
zz=transpose(rr)


rbphi=rr*bphi

 drBphidr=pres
  for j=0,n2-1 do begin
	      drBphidr(*,j)= deriv(x1,rbphi(*,j))
  endfor

 dBzdr=pres
  for j=0,n2-1 do begin
	      dBzdr(*,j)= deriv(x1,bz(*,j))
  endfor

 dprdr=pres
  for j=0,n2-1 do begin
	      dprdr(*,j)= -deriv(x1,pres(*,j))
  endfor

 dprdz=pres
  for i=0,n1-1 do begin
	      dprdz(i,*)= -deriv(x2,pres(i,*))
  endfor

 dBrdz=br
  for i=0,n1-1 do begin
	      dBrdz(i,*)= deriv(x2,br(i,*))
  endfor


 dBphidz=br
  for i=0,n1-1 do begin
	      dBphidz(i,*)= deriv(x2,bphi(i,*))
  endfor


	alf, b1,b2,b3,v1,v2,rho,pr,mf,ma,ms,mson

jz  =  drBphidr/rr
jr  = -dBphidz
jphi=  -dBrdz + dBzdr

fr = jphi*bz- jz*bphi

fphi = -jr*bz+jz*br

fz = jr*bphi - jphi*br

r=x1
z=x2
dr=dx1
dz=dx2
i=1.5

xbeg=0
xend=6
ybeg=-1e-4
yend=1e-4
field_line, br, bz, x1,x2,  i,0.05, rr,zz
interplot,dprdz ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
plot, sqrt(rpl^2+zpl^2), qinter1,$
	;yrange=[-0.00001,2e-6],$
	yrange=[ybeg,yend],$
	xrange=[xbeg,xend],$
	xtitle="spherical radius, R",$
	title="-!MG!Xp!Dz!N Pressure Gradient"

for i=2,10 do begin
field_line, br, bz, x1,x2,  i,0.05, rr,zz
interplot, dprdz ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
oplot, sqrt(rpl^2+zpl^2), qinter1, linestyle=i-2
endfor


i=1.5

field_line, br, bz, x1,x2,  i,0.05, rr,zz
interplot,dprdr ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
plot, sqrt(rpl^2+zpl^2), qinter1,$
	;yrange=[-0.00001,0.0],$
	yrange=[ybeg,yend],$
	xrange=[xbeg,xend],$
	xtitle="spherical radius, R",$
	title="-!MG!Xp!Dr!N Pressure Gradient"

for i=2,10 do begin
field_line, br, bz, x1,x2,  i,0.05, rr,zz
interplot, dprdr ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
oplot, sqrt(rpl^2+zpl^2), qinter1, linestyle=i-2
endfor


i=1.5

field_line, br, bz, x1,x2,  i,0.05, rr,zz
interplot,fr ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
plot, sqrt(rpl^2+zpl^2), qinter1,$
	;yrange=[-0.00001,0.0],$
	yrange=[ybeg,yend],$
	xrange=[xbeg,xend],$
	xtitle="spherical radius, R",$
	title="Fr,Radial Lorentz Force"

for i=2,10 do begin
field_line, br, bz, x1,x2,  i,0.05, rr,zz
interplot, fr ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
oplot, sqrt(rpl^2+zpl^2), qinter1, linestyle=i-2
endfor

i=1.5

field_line, br, bz, x1,x2,  i,0.05, rr,zz
interplot,fz ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
plot, sqrt(rpl^2+zpl^2), qinter1,$
;	yrange=[0.0, 2e-5],$
	yrange=[ybeg,yend],$
	xrange=[xbeg,xend],$
	xtitle="spherical radius, R",$
	title="F!Dz!N Lorentz Force"

for i=2,10 do begin
field_line, br, bz, x1,x2,  i,0.05, rr,zz
interplot, fz ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
oplot, sqrt(rpl^2+zpl^2), qinter1, linestyle=i-2
endfor

i=1.5

field_line, br, bz, x1,x2,  i,0.05, rr,zz
interplot, fphi ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
plot, sqrt(rpl^2+zpl^2), qinter1,$
	;yrange=[-0.00001,0.0],$
	yrange=[ybeg,yend],$
	xrange=[xbeg,xend],$
	xtitle="spherical radius, R",$
	title="F!7u!3 Lorentz Force"

for i=2,10 do begin
field_line, br, bz, x1,x2,  i,0.05, rr,zz
interplot, fphi ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
oplot, sqrt(rpl^2+zpl^2), qinter1, linestyle=i-2
endfor

;plot, x2,alog10(pres(100,*))


i=1.5

field_line, br, bz, x1,x2,  i,0.05, rr,zz
interplot, ms ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
plot, sqrt(rpl^2+zpl^2), qinter1,$
	yrange=[0,2.0],$
	xrange=[xbeg,xend],$
	xtitle="spherical radius, R",$
	title="Slow Mag Point"

for i=2,10 do begin
field_line, br, bz, x1,x2,  i,0.05, rr,zz
interplot, ms ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
oplot, sqrt(rpl^2+zpl^2), qinter1, linestyle=i-2
endfor

im=tvread(filename='Forcesx', /nodialog, /png)

!P.POSITION=0
!P.MULTI=0

end
