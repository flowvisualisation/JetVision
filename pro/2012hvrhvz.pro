
; Calculate the torques
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


omega=v3(*,0)/x1


omegamatrix=rebin(reform(omega,n1,1),n1,n2 )

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


msmag=fphi/den/omegamatrix^2

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

msvisc=divTphi/den/omegamatrix^2


hpeak=fltarr(n1)
hmax=fltarr(n1)
hvr=fltarr(n1)
hvz=fltarr(n1)
rmpeak=fltarr(n1)
jphipeak=fltarr(n1)
bzpeak=fltarr(n1)
vzpeak=fltarr(n1)
vrpeak=fltarr(n1)
brpeak=fltarr(n1)
etapeak=fltarr(n1)
msmagpeak=fltarr(n1)
msviscpeak=fltarr(n1)
lambdapeak=fltarr(n1)

cellheight=indgen(n1)
;h2=x2
;;;;;

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

xbeg=1.3
xend=10
ybeg=-0.002
yend=0.001
window,23, xs=400, ys=600
!P.MULTI=1
plot, x1,x1*divTphi(*,0),$
	ytitle="Torques",$
	xtitle="r/r!D0!N",$
;	/xlog,$
	xrange=[xbeg,xend],$
	yrange=[ybeg,yend],$
	xstyle=1
oplot, x1,x1*fphi(*,0), linestyle=1
oplot, x1,x1*trphi(*,0), linestyle=2
oplot, x1,x1*tzphi(*,0), linestyle=3
;oplot, x2,x1(nx)*divTphi(nx,*), psym=4
items=['visc','mag','trphi','tzphi']
legend,items,linestyle=indgen(4),/right


mag_torque= fphi
visc_torque= divTphi
lambda=mag_torque/visc_torque
;lambda_ave=0
;r_ave=0

;for i=14,100 do begin

;lambda_ave=lambda_ave+ lambda(i)+1/lambda(i)*x1(i)^2 * omega(i)*dx1(i)
;r_ave=r_ave+ x1(i)^2 *omega(i) *dx1(i)

;endfor

;print, 'Lambda+1/Lambda', lambda_ave/r_ave

for i=0,n1-1 do begin
      a=where(eta(i,*) eq 0)
      zvr=where(vr(i,*) gt 0.0) 
      zvz=where(vz(i,*) gt 0.0) 
 ;     print,a(0)
    ;  a=a-1
      a=a * (a gt 0)
      zvr=zvr * (zvr gt 0)
      zvz=zvz * (zvz gt 0)
      hpeak(i)=x2(a(0))
      hmax(i)=x2(a(0))
      hvr(i) = x2(zvr(0))
      hvz(i) = x2(zvz(0))
      cellheight[i]=zvr(0)
      qq=cellheight[i]
      jphipeak[i]=jphi[i, qq ]
      bzpeak[i]=bz[i, qq ]
      vzpeak[i]=vz[i, qq ]
      vrpeak[i]=vr[i, qq ]
      brpeak[i]=br[i, qq ]
      etapeak[i]=eta[i, qq ]
      msmagpeak[i]=msmag[i, qq ]
      msviscpeak[i]=msvisc[i, qq ]
      lambdapeak[i]=lambda[i, qq ]
   endfor


window, title="lambda ", xsize=1000, ysize=800

!x.range=[1,20]
!y.range=0
plot, x1, hvr, linestyle=0, title='Lambda, !7K!3 '
oplot, x1, hvz, linestyle=1
oplot, x1, hmax, linestyle=2
;items=['ur/ omega h','(JxB)!D!7 u!3!N/ !7q x!3!U2!N h ','DivT!D!7u!3!N / !7q x!3!U2!N h ']
items=['hvr','hvz','hmax']
legend,items,linestyle=indgen(3),/right

im=tvread(filename='2012hvrhvz',/png,/nodialog)






end
