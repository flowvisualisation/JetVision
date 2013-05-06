
;window,22, xs=1000, ys=900
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



rr=rebin(reform(x1,n1,1),n1,n2 )
zz=rebin(reform(x2,n2,1),n2,n1 )
zz=transpose(zz)


rad=sqrt(rr^2 +zz^2)

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


 dvphidr=pres
  for j=0,n2-1 do begin
         dvphidr(*,j)= deriv(x1,vphi(*,j))
  endfor

 dvphidz=pres
  for i=0,n1-1 do begin
         dvphidz(i,*)= deriv(x2,vphi(i,*))
  endfor

 dvrdr=pres
  for j=0,n2-1 do begin
         dvrdr(*,j)= deriv(x1,vr(*,j))
  endfor
 dvrdz=pres
  for i=0,n1-1 do begin
         dvrdz(i,*)= deriv(x2,vr(i,*))
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

jr  = -dBphidz
jphi=  dBrdz - dBzdr
jz  =  drBphidr/rr

fr 	= 	jphi*bz	- 	jz*bphi
fphi 	= -jr*bz		+	jz*br
fz 	= 	jr*bphi 	- 	jphi*br



vsound2=pres/den
vsoundmid=vsound2
  for j=0,n1-1 do begin
     vsoundmid(*,j)=pres(*,0)/den(*,0)
     endfor
     eta=sqrt(rr^3)*(vsoundmid + (2.0/5.0)*(1/rad-1/rr))
eta(*,*) =eta(*,*)*(eta(*,*) gt 0.0)
kvisc=2.0/3.0*den*eta





kvisc=2.0/3.0*den*eta
tau_rphi=kvisc*(dvphidr - vphi/rr)
tau_zphi=kvisc*(dvphidz)





 drtau_rphidr=pres
 rtau_rphi=rr*tau_rphi
  for j=0,n2-1 do begin
         drtau_rphidr(*,j)= deriv(x1,rtau_rphi(*,j))
  endfor

 dtau_zphidz=pres
  for i=0,n1-1 do begin
         dtau_zphidz(i,*)= deriv(x2,tau_zphi(i,*))
  endfor


  divTphi= 1/rr*   drtau_rphidr +  dtau_zphidz + tau_rphi/rr


	
	trphi= 1/rr*   drtau_rphidr + tau_rphi/rr
	tzphi=  dtau_zphidz

!P.POSITION=0
!P.MULTI=[0,1,4]
!P.CHARSIZE=4.0
!P.CHARTHICK=2.0
!P.THICK=3.0

window,22, xs=400, ys=900
xbeg=0
xend=1
ybeg=-0.01
yend=0
nx=60
plot, x2,x1(nx)*divTphi(nx,*),$
	title="Visc. Torque at r="+string(x1(nx)),$
	xtitle="z/r!D0!N",$
	xrange=[xbeg,xend],$
	yrange=[ybeg,yend],$
	xstyle=1
plot, x2,x1(nx)*fphi(nx,*),$
	title="Mag. Torque at r="+string(x1(nx)),$
	xtitle="z/r!D0!N",$
	xrange=[xbeg,xend],xstyle=1
plot, x2,x1(nx)*trphi(nx,*),$
	title="rphi at r="+string(x1(nx)),$
	xtitle="z/r!D0!N",$
	xrange=[xbeg,xend],xstyle=1
plot, x2,x1(nx)*tzphi(nx,*),$
	title="zphi at r="+string(x1(nx)),$
	xtitle="z/r!D0!N",$
	xrange=[xbeg,xend],xstyle=1

;tau_zz=dvphidz

;display,  divTphi, /vbar
im=tvread(filename="torques2",/png,/nodialog)

ybeg=-0.0025
yend=0.002
window,23, xs=900, ys=900
!P.MULTI=1
plot, x2,x1(nx)*divTphi(nx,*),$
	title="Visc. Torque at r="+string(x1(nx)),$
	xtitle="z/r!D0!N",$
	xrange=[xbeg,xend],$
	yrange=[ybeg,yend],$
	xstyle=1
oplot, x2,x1(nx)*fphi(nx,*), linestyle=1
oplot, x2,x1(nx)*trphi(nx,*), linestyle=2
oplot, x2,x1(nx)*tzphi(nx,*), linestyle=3
oplot, x2,x1(nx)*divTphi(nx,*), psym=4
items=['viscous','mag','trphi','tzphi']
legend,items,linestyle=indgen(4),/right

window,24, xs=900, ys=900
plot,x2,fphi(nx,*)/divTphi(nx,*), xrange=[0,2],psym=2




ybeg=-0.0002
yend=0.0002
xend=2
!P.CHARSIZE=1.0
!P.CHARTHICK=1.0
set_plot,'ps'
device,filename='torques_only.eps', /encapsulated, xsize=16, ysize=8
plot, x2,x1(nx)*divTphi(nx,*),$
   title="Visc. Torque at r="+string(x1(nx)),$
   xtitle="z/r!D0!N",$
   xrange=[xbeg,xend],$
   yrange=[ybeg,yend],$
   ystyle=1,$
   xstyle=1
oplot, x2,x1(nx)*fphi(nx,*), linestyle=1
oplot, x2,x1(nx)*trphi(nx,*), linestyle=2
oplot, x2,x1(nx)*tzphi(nx,*), linestyle=3
oplot, x2,x1(nx)*divTphi(nx,*), psym=4
items=['visc','mag','trphi','tzphi']
legend,items,linestyle=indgen(4), /right
device,/close


set_plot,'x'

!P.COLOR=0
!P.BACKGROUND=255
end
