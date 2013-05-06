
window,22, xs=900, ys=600
loadct,0

!P.POSITION=0
!P.MULTI=0
!P.CHARSIZE=2.0
!P.POSITION=0
!P.CHARSIZE=2.0
!P.CHARTHICK=2.0
!P.THICK=3.0
!P.COLOR=0
!P.BACKGROUND=255


; K

den=rho(0)
pres=pr(0)

vr=v1(0)
vz=v2(0)
vphi=v3(0)

br=b1(0)
bz=b2(0)
bphi=b3(0)


r=x1
z=x2
dr=dx1
dz=dx2

rr=rebin(reform(x1,n1,1),n1,n2 )
zz=rebin(reform(x2,n2,1),n2,n1 )
zz=transpose(zz)


rad=sqrt(rr^2+zz^2)
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


 dvphidr=pres
  for j=0,n2-1 do begin
         dvphidr(*,j)= deriv(x1,vphi(*,j))
  endfor

 dvphidz=pres
  for i=0,n1-1 do begin
         dvphidz(i,*)= deriv(x2,vphi(i,*))
  endfor

 dvzdr=pres
  for j=0,n2-1 do begin
         dvzdr(*,j)= deriv(x1,vz(*,j))
  endfor

 drvrdr=pres
  for j=0,n2-1 do begin
         drvrdr(*,j)= deriv(x1,rr(*,j)*vr(*,j))
  endfor

 dvrdr=pres
  for j=0,n2-1 do begin
         dvrdr(*,j)= deriv(x1,vr(*,j))
  endfor

 dvzdz=pres
  for i=0,n1-1 do begin
         dvzdz(i,*)= deriv(x2,vz(i,*))
  endfor

 dvrdz=pres
  for i=0,n1-1 do begin
         dvrdz(i,*)= deriv(x2,vr(i,*))
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
jphi=  dBrdz - dBzdr

fr = jphi*bz- jz*bphi

fphi = -jr*bz+jz*br

fz = jr*bphi - jphi*br


div_v=1./rr*drvrdr + dvzdz


vsound2=pres/den
vsoundmid=vsound2
  for j=0,n1-1 do begin
	  vsoundmid(*,j)=pres(*,0)/den(*,0)
	  endfor
	  eta=sqrt(rr^3)*(vsoundmid + (2.0/5.0)*(1/rad-1/rr))

eta(*,*) =eta(*,*)*(eta(*,*) gt 0.0)
kvisc=2.0/3.0*den*eta

tau_rr=kvisc*(2*dvrdr -2.0/3.0*div_v )
tau_zz=kvisc*(2*dvzdz -2.0/3.0*div_v )
tau_phiphi=kvisc*(2*vr/rr -2.0/3.0*div_v )
tau_rphi=kvisc*(dvphidr - vphi/rr)
tau_zphi=kvisc*(dvphidz)
tau_rz=kvisc*(dvzdr +dvrdz)






 drtau_rrdr=pres
 rtau_rr=rr*tau_rr
  for j=0,n2-1 do begin
         drtau_rrdr(*,j)= deriv(x1,rtau_rr(*,j))
  endfor

 drtau_rzdr=pres
 rtau_rz=rr*tau_rz
  for j=0,n2-1 do begin
         drtau_rzdr(*,j)= deriv(x1,rtau_rz(*,j))
  endfor

 drtau_rphidr=pres
 rtau_rphi=rr*tau_rphi
  for j=0,n2-1 do begin
         drtau_rphidr(*,j)= deriv(x1,rtau_rphi(*,j))
  endfor

 dtau_rzdz=pres
  for i=0,n1-1 do begin
         dtau_rzdz(i,*)= deriv(x2,tau_rz(i,*))
  endfor

 dtau_rphidz=pres
  for i=0,n1-1 do begin
         dtau_rphidz(i,*)= deriv(x2,tau_rphi(i,*))
  endfor

 dtau_zzdz=pres
  for i=0,n1-1 do begin
         dtau_zzdz(i,*)= deriv(x2,tau_zz(i,*))
  endfor

 dtau_zphidz=pres
  for i=0,n1-1 do begin
         dtau_zphidz(i,*)= deriv(x2,tau_zphi(i,*))
  endfor


  divTr= 1/rr*   drtau_rrdr +  dtau_rzdz - tau_phiphi/rr
  divTz= 1/rr*   drtau_rzdr +  dtau_zzdz 
  divTphi= 1/rr*   drtau_rphidr +  dtau_zphidz + tau_rphi/rr

;feff =  den*(vphi^2/rr  - dphidr )

;  
; Projection and Interpolation
;

bp=sqrt(br^2+bz^2)
fp = (-bz*dprdz -br*dprdr)/bp

fm =  (bz*fz +br*fr)/bp
feff =  den/bp*(br*vphi^2/rr  - bz*dphidz - br*dphidr )


fv =  (bz*divTz +br*divTr)/bp

ftot=feff + fm +fp +fv


i=4.6
xbeg=1e-1
xend=10
ybeg=-1e-5
yend=1e-5
field_line, br, bz, x1,x2,  i,0.0005, rr,zz
fluff=fp
interplot,fluff ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
qrad= sqrt(rpl^2+zpl^2)
qrad= zpl
plot, qrad, qinter1,$
title='r='+string(i),$
	/xlog,$
	xrange=[xbeg,xend],$
	yrange=[ybeg,yend]

fluff=fm
interplot,fluff ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
qrad= sqrt(rpl^2+zpl^2)
qrad= zpl
oplot, qrad, qinter1,linestyle=1

fluff=feff
interplot,fluff ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
qrad= sqrt(rpl^2+zpl^2)
qrad= zpl
oplot, qrad, qinter1,linestyle=2

fluff=fv
interplot,fluff ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
qrad= sqrt(rpl^2+zpl^2)
qrad= zpl
oplot,  qrad, qinter1,linestyle=3
oplot, qrad, qinter1,psym=4

fluff=ftot
interplot,fluff ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
qrad= sqrt(rpl^2+zpl^2)
qrad= zpl
oplot, qrad, qinter1,linestyle=4

items = ['Fp','Fm', 'feff', ' fv', 'sum']
lines=indgen(5)
legend,  items,linestyle=lines,/right


im=tvread(filename="vforceproj",/png,/nodialog)


end
