
window,22, xs=900, ys=600
loadct,0

!P.POSITION=0
!P.MULTI=[0,2,2]
!P.CHARSIZE=2.0
!P.POSITION=0
!P.CHARSIZE=2.0
!P.CHARTHICK=2.0
!P.THICK=3.0


; K

den=rho(0)
pres=pr(0)

vx=v1(0)
vy=v2(0)
vphi=v3(0)

br=b1(0)
bz=b2(0)
bphi=b3(0)


rr=rebin(reform(x1,n1,1),n1,n1 )
zz=transpose(rr)


rbphi=rr*bphi

 drBphidr=pres
 dBphidr=pres
 phi=-1/sqrt(rr^2+zz^2)

  for j=0,n2-1 do begin
	      drBphidr(*,j)= deriv(x1,rbphi(*,j))
	      dBphidr(*,j)= deriv(x1,bphi(*,j))
  endfor

 dBzdr=pres
  for j=0,n2-1 do begin
	      dBzdr(*,j)= deriv(x1,bz(*,j))
  endfor

 dprdr=pres
  for j=0,n2-1 do begin
	      dprdr(*,j)= deriv(x1,pres(*,j))
  endfor

 dphidr=pres
  for j=0,n2-1 do begin
	      dphidr(*,j)= deriv(x1,phi(*,j))
  endfor

 dphidz=pres
  for i=0,n1-1 do begin
	      dphidz(i,*)= deriv(x2,phi(i,*))
  endfor

 dprdz=pres
  for i=0,n1-1 do begin
	      dprdz(i,*)= deriv(x2,pres(i,*))
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



bp=sqrt(br^2+bz^2)
fp = (-bz*dprdz -br*dprdr)/bp

fm =  (bz*fz +br*fr)/bp
feff =  den/bp*(br*vphi^2/rr  - bz*dphidz - br*dphidr )


ftot=feff + fm +fp
r=x1
z=x2
dr=dx1
dz=dx2
i=1.5

xbeg=0
xend=6
ybeg=-1e-4
yend=1e-4


fluff=fp
field_line, br, bz, x1,x2,  i,0.05, rr,zz
interplot,fluff ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
plot, sqrt(rpl^2+zpl^2), qinter1,$
	;yrange=[-0.00001,2e-6],$
	yrange=[ybeg,yend],$
	xrange=[xbeg,xend],$
	xtitle="spherical radius, R",$
	title="Fp"

for i=2,10 do begin
field_line, br, bz, x1,x2,  i,0.05, rr,zz
interplot, fluff ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
oplot, sqrt(rpl^2+zpl^2), qinter1, linestyle=i-2
endfor


i=1.5
fluff=fm
field_line, br, bz, x1,x2,  i,0.05, rr,zz
interplot,fluff ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
plot, sqrt(rpl^2+zpl^2), qinter1,$
	;yrange=[-0.00001,0.0],$
	yrange=[ybeg,yend],$
	xrange=[xbeg,xend],$
	xtitle="spherical radius, R",$
	title="Fm"

for i=2,10 do begin
field_line, br, bz, x1,x2,  i,0.05, rr,zz
interplot, fluff ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
oplot, sqrt(rpl^2+zpl^2), qinter1, linestyle=i-2
endfor


i=1.5

fluff=feff
field_line, br, bz, x1,x2,  i,0.05, rr,zz
interplot,fluff ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
plot, sqrt(rpl^2+zpl^2), qinter1,$
	;yrange=[-0.00001,0.0],$
	yrange=[ybeg,yend],$
	xrange=[xbeg,xend],$
	xtitle="spherical radius, R",$
	title="Feff"

for i=2,10 do begin
field_line, br, bz, x1,x2,  i,0.05, rr,zz
interplot, fluff ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
oplot, sqrt(rpl^2+zpl^2), qinter1, linestyle=i-2
endfor

i=1.5

fluff=ftot
field_line, br, bz, x1,x2,  i,0.05, rr,zz
interplot,fluff ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
plot, sqrt(rpl^2+zpl^2), qinter1,$
;	yrange=[0.0, 2e-5],$
	yrange=[ybeg,yend],$
	xrange=[xbeg,xend],$
	xtitle="spherical radius, R",$
	title="Ftotal"

for i=2,10 do begin
field_line, br, bz, x1,x2,  i,0.05, rr,zz
interplot, fluff ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
oplot, sqrt(rpl^2+zpl^2), qinter1, linestyle=i-2
endfor

im=tvread(filename='Forceproj', /nodialog, /png)

!P.POSITION=0
!P.MULTI=0

end
