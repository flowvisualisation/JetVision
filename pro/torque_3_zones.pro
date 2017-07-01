
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





ybeg=-0.0002
yend=0.0002
xend=2
!P.CHARSIZE=1.0
!P.CHARTHICK=1.0
!P.font=0
!P.multi=[0,1,3]
!x.range=[2,3]
set_plot,'ps'
device,filename='torques_3_zones.eps', /encapsulated, xsize=16, ysize=16
device, /times

nx=24
xh=x2/0.1/x1(nx)
a=x1(nx)*divTphi(nx,*)
b=x1(nx)*fphi(nx,*)
c=x1(nx)*trphi(nx,*)
d=x1(nx)*tzphi(nx,*)
!y.range=[min([[a],[b],[c],[d] ]),max([[a],[b],[c],[d] ])]
!y.range=[-1e-4,1e-4]
plot, xh,x1(nx)*divTphi(nx,*),$
   title="Torques at t=" $
   +string(time(nlast)/2/!PI, format='(F6.2)')+'!9t!X!DK0!N',$
   xtitle="Altitude/Thermal heightscale, z/h",$
   ytitle="r="+string(x1(nx), format='(F5.1)'),$
   ystyle=1,$
   xstyle=1
oplot, xh,x1(nx)*fphi(nx,*), linestyle=1
oplot, xh,x1(nx)*trphi(nx,*), linestyle=2
oplot, xh,x1(nx)*tzphi(nx,*), linestyle=3
items=['Visc','Mag','dT!Dr!9f!X!N/dr','dT!Dz!9f!X!N/dz']
legend,items,linestyle=indgen(4), /right

nx=100
xh=x2/0.1/x1(nx)
a=x1(nx)*divTphi(nx,*)
b=x1(nx)*fphi(nx,*)
c=x1(nx)*trphi(nx,*)
d=x1(nx)*tzphi(nx,*)
!y.range=[min([[a],[b],[c],[d] ]),max([[a],[b],[c],[d] ])]
!y.range=[-1e-5,1e-5]
plot, xh,x1(nx)*divTphi(nx,*),$
   xtitle="Altitude/Thermal heightscale, z/h",$
   ytitle="r="+string(x1(nx), format='(F5.1)'),$
   ystyle=1,$
   xstyle=1
oplot, xh,x1(nx)*fphi(nx,*), linestyle=1
oplot, xh,x1(nx)*trphi(nx,*), linestyle=2
oplot, xh,x1(nx)*tzphi(nx,*), linestyle=3
legend,items,linestyle=indgen(4), /right

nx=400
xh=x2/0.1/x1(nx)
a=x1(nx)*divTphi(nx,*)
b=x1(nx)*fphi(nx,*)
c=x1(nx)*trphi(nx,*)
d=x1(nx)*tzphi(nx,*)
!y.range=[min([[a],[b],[c],[d] ]),max([[a],[b],[c],[d] ])]
plot, xh,x1(nx)*divTphi(nx,*),$
   xtitle="Altitude/Thermal heightscale, z/h",$
   ytitle="r="+string(x1(nx), format='(F5.1)'),$
   ystyle=1,$
   xstyle=1
oplot, xh,x1(nx)*fphi(nx,*), linestyle=1
oplot, xh,x1(nx)*trphi(nx,*), linestyle=2
oplot, xh,x1(nx)*tzphi(nx,*), linestyle=3
legend,items,linestyle=indgen(4), /right

device,/close


set_plot,'x'

!P.COLOR=0
!P.BACKGROUND=255
end
