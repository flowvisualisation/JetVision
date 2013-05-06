
loadct,0
!P.COLOR=0
!P.BACKGROUND=255

;window,22, xs=1000, ys=700
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
jz  =  drBphidr/rr
jphi=  dBrdz - dBzdr

fr 	= 	jphi*bz	- 	jz*bphi
fz 	= 	jr*bphi 	- 	jphi*br
fphi 	= -jr*bz		+	jz*br



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
!P.CHARSIZE=2.0
!P.CHARTHICK=1.0
!P.THICK=3.0

window,22, xs=400, ys=700
xbeg=0.005
xend=2
ybeg=-0.01
yend=0
nx=60
ybeg=-0.0002
yend=0.0002
plot, x2,x1(nx)*divTphi(nx,*),$
	ytitle="Torques,r="+string(x1(nx),format='(F3.1)' ),$
	xtitle="z/r!D0!N",$
	xrange=[xbeg,xend],$
	yrange=[ybeg,yend],$
	xstyle=1
oplot, x2,x1(nx)*fphi(nx,*), linestyle=1
oplot, x2,x1(nx)*trphi(nx,*), linestyle=2
oplot, x2,x1(nx)*tzphi(nx,*), linestyle=3

ybeg=-0.005
yend=0.005
plot, x2,jphi(nx,*),$
	yrange=[ybeg,yend],$
	ytitle="Jphi,r="+string(x1(nx),format='(F3.1)' ),$
	xrange=[xbeg,xend],xstyle=1
oplot, x2,jr(nx,*), linestyle=1
oplot, x2,jr(nx,*), psym=1

ybeg=-0.003
yend=0.001
plot, x2,b1(nx,*,0),$
	yrange=[ybeg,yend],$
	ytitle="B1,r="+string(x1(nx),format='(F3.1)' ),$
	xrange=[xbeg,xend],xstyle=1
oplot, x2,b2(nx,*,0), linestyle=1
oplot, x2,b3(nx,*,0), linestyle=2
plot, x2,kvisc(nx,*),$
	ytitle="visc,r="+string(x1(nx),format='(F3.1)' ),$
	xtitle="z/r!D0!N",$
	xrange=[xbeg,xend],xstyle=1

;tau_zz=dvphidz

;display,  divTphi, /vbar
im=tvread(filename="torques2",/png,/nodialog)

ybeg=-0.001
yend=0.001
window,23, xs=400, ys=600
!P.MULTI=1
plot, x2,x1(nx)*divTphi(nx,*),$
	ytitle="Torques at r="+string(x1(nx)),$
	xtitle="z/r!D0!N",$
	/xlog,$
	xrange=[xbeg,xend],$
	yrange=[ybeg,yend],$
	xstyle=1
oplot, x2,x1(nx)*fphi(nx,*), linestyle=1
oplot, x2,x1(nx)*trphi(nx,*), linestyle=2
oplot, x2,x1(nx)*tzphi(nx,*), linestyle=3
;oplot, x2,x1(nx)*divTphi(nx,*), psym=4
items=['visc','mag','trphi','tzphi']
legend,items,linestyle=indgen(4),/right

window,24, xs=400, ys=600
plot,x2,fphi(nx,*)/divTphi(nx,*), xrange=[0,2],psym=2



xd=where(kvisc(nx,*)  eq 0)
xd=x2(xd(0))

!P.CHARSIZE=1.0
!P.CHARTHICK=1.0
set_plot,'ps'
device,filename='torques.eps',XSIZE=16, YSIZE=40, /encapsulated

!P.POSITION=0
!P.MULTI=[0,1,5]
!P.CHARSIZE=2
csize=1
!P.CHARTHICK=1.0
!P.THICK=3.0

xbeg=0.005
xend=2
ybeg=-0.01
yend=0
nx=60
ybeg=-0.0002
yend=0.0002
plot, x2,x1(nx)*divTphi(nx,*),$
   ytitle="Torques,r="+string(x1(nx),format='(F3.1)' ),$
   xrange=[xbeg,xend],$
   yrange=[ybeg,yend],$
   xstyle=1
oplot, x2,x1(nx)*fphi(nx,*), linestyle=1
oplot, x2,x1(nx)*trphi(nx,*), linestyle=2
oplot, x2,x1(nx)*tzphi(nx,*), linestyle=3
items=['visc','mag','trphi','tzphi']
legend,items,linestyle=indgen(4),/right,charsize=csize
plots, [xd,xd],[ybeg,yend]

ybeg=-0.002
yend=0.007
plot, x2,jphi(nx,*),$
   yrange=[ybeg,yend],$
   ytitle="Jphi,r="+string(x1(nx),format='(F3.1)' ),$
   xrange=[xbeg,xend],xstyle=1
oplot, x2,jr(nx,*), linestyle=1
oplot, x2,jr(nx,*), psym=1
items=['Jr','Jphi']
legend,items,linestyle=indgen(2),/right,charsize=csize

plots, [xd,xd],[ybeg,yend]

ybeg=-0.003
yend=0.001
plot, x2,b1(nx,*,0),$
   yrange=[ybeg,yend],$
   ytitle="B1,r="+string(x1(nx),format='(F3.1)' ),$
   xrange=[xbeg,xend],xstyle=1
oplot, x2,b2(nx,*,0), linestyle=1
oplot, x2,b3(nx,*,0), linestyle=2
items=['Br','Bz','Bphi']
legend,items,linestyle=indgen(3),/bottom,charsize=csize

plots, [xd,xd],[ybeg,yend]

plot, x2,kvisc(nx,*),$
   ytitle="visc,r="+string(x1(nx),format='(F3.1)' ),$
   xrange=[xbeg,xend],xstyle=1

plots, [xd,xd],[0,0.0012]

ybeg=-0.003
yend=0.001
plot, x2,v1(nx,*,0),$
   ytitle="v1,r="+string(x1(nx),format='(F3.1)' ),$
   xtitle="z/r!D0!N",$
   xrange=[xbeg,xend],xstyle=1

plots, [xd,xd],[-0.05,0.15]


device,/close


set_plot,'x'

!P.COLOR=0
!P.BACKGROUND=255
end
