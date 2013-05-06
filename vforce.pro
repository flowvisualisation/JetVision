
window,22, xs=900, ys=600
loadct,0

!P.POSITION=0
!P.MULTI=[0,1,2]
!P.CHARSIZE=2.0
!P.POSITION=0
!P.CHARSIZE=2.0
!P.CHARTHICK=2.0
!P.THICK=3.0


; K

den=rho(0)
pres=pr(0)

vr=v1(0)
vz=v2(0)
vphi=v3(0)

br=b1(0)
bz=b2(0)
bphi=b3(0)



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


bp=sqrt(br^2+bz^2)
fp = (-bz*dprdz -br*dprdr)/bp

fm =  (bz*fz +br*fr)/bp
;feff =  den/bp*(br*vphi^2/rr  - bz*dphidz - br*dphidr )
feff =  den*(vphi^2/rr  - dphidr )


ftot=feff + fm +fp

nx=20

print, x1(nx)

xbeg=0
xend=1
ybeg=-2e-5
yend=1e-5
sum=fz-dprdz+divTz - den*dphidz
plot, x2,fz(nx,*),$
	linestyle=0,$
	title="Vertical forces, r="+string(x1(nx)) ,$
	xtitle="z/r!D0!N",$
	yrange=[ybeg,yend],$
	xrange=[xbeg,xend]
oplot, x2,-dprdz(nx,*),$
	linestyle=1
oplot, x2, divTz(nx,*),$
	linestyle=2
oplot, x2, sum(nx,*),$
	linestyle=3
items = ['Fz','dp/dz', 'div T z']
items = ['Fz','dp/dz', 'div T z','sum']
lines=indgen(4)
legend,  items,linestyle=lines,/right

sum=fr-dprdr+divTr +feff
plot,x2,fr(nx,*),$
	linestyle=0,$
	title="Radial forces",$
	xtitle="z/r!D0!N",$
	xrange=[xbeg,xend]
oplot,x2, -dprdr(nx,*),$
	linestyle=1
oplot,x2, feff(nx,*),$
	linestyle=2
oplot,x2,  divTr(nx,*),$
	linestyle=3
oplot, x2, sum(nx,*),$
	linestyle=4

items = ['Fr','dp/dr', 'feff', ' divTr', 'sum']
lines=indgen(5)
legend,  items,linestyle=lines,/right


im=tvread(filename="vforce",/png,/nodialog)


end
