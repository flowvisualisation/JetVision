
nend=nlast
pload,out=nend


br=b1(0)
bz=b2(0)
bphi=b3(0)
den=rho(0)
pres=pr(0)
        
vr=v1(0)
vz=v2(0) 
vphi=v3(0)


mu=(b1(0)^2+b2(0)^2+b3(0)^2)/pr(0)

poynting=v3(*,*,0)*b3(*,*,0)*sqrt(b1(*,*,0)^2 + b2(*,*,0)^2) 
poynting = poynting- b3(*,*,0)^2*sqrt(v1(*,*,0)^2 +v2(*,*,0)^2)
poynting = -poynting



kinetic=(0.5*rho(*,*,0)*(v1(*,*,0)^2 $
+v2(*,*,0)^2+v3(*,*,0)^2    )*sqrt(v1(*,*,0)^2 +v2(*,*,0)^2))
thermal= (2.5)*(pr(*,*,0))* sqrt(v1(*,*,0)^2 +v2(*,*,0)^2)

michel=poynting/(kinetic +thermal )
michel=poynting/(kinetic)

minmichel=fltarr(n1)



r=x1
z=x2
dr=dx1
dz=dx2
rr=x1
zz=0.1*x1

rrr=rebin(reform(r,n1,1),n1,n1 )
zz=rebin(reform(x2,n2,1),n2,n1 )
zzzz=transpose(zz)

rad=sqrt(rrr^2 + zzzz^2)
;michel=zzzz

sound2=(pr(0)/rho(0))
sound=sqrt(pr(0)/rho(0))



den=rho(0)
pres=pr(0)

vsound2=pres/den
vsoundmid=vsound2
  for j=0,n1-1 do begin
     vsoundmid(*,j)=pres(*,0)/den(*,0)
     endfor
     eta=sqrt(rrr^3)*(vsoundmid + (2.0/5.0)*(1/rad-1/rrr))
eta(*,*) =eta(*,*)*(eta(*,*) gt 0.0)
kvisc=2.0/3.0*den*eta

h2=x2
for i=0,n1-1 do begin
      a=where(eta(i,*) eq 0)
 ;     print,a(0)
;      a=a-1
      a=a+1
      a=a * (a gt 0)
      h2(i)=x2(a(0))
   endfor


zz=smooth(h2,5)

;zz=2.*sound(*,0)*rr^1.5+0.05

;interplot,michel ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
;interplot,mu ,r,z,dr,dz,rr,zz,qinter2,rpl,zpl

;loadct,33
;display,alog10(rho(0)),x1=x1,x2=x2, xrange=[1,10], yrange=[0,5], ims=4
;oplot,rr,zz



!P.MULTI=0
pres=pr(0)
dens=rho(0)

r=x1
h1=sqrt(pres(*,0)/dens(*,0)) *r^1.5
h2=h1

vr=v1(0)

for i=0,n1-1 do begin

      a=where(vr(i,*)>0)
      print,a(0)
      h2(i)=x2(a(0))
   endfor


xbeg=0
xend=25

i=5
field_line, br, bz, x1,x2,  i,0.0005, rf1,zf1
i=13
field_line, br, bz, x1,x2,  i,0.0005, rf2,zf2
i=15
field_line, br, bz, x1,x2,  i,0.0005, rf3,zf3


for i=0,1 do begin
if ( i eq 0) then begin
!p.font=1
set_plot, 'ps'
device,filename='Y2012height.eps', /encapsulated, /color
device, /times
device, xsize=16, ysize=16
!p.charsize=1.0
endif else begin
set_plot,'x'
window, xsize=1000, ysize=600
endelse
!x.style=1
!P.THICK=1.0
loadct,33
tvlct,0,0,0,1
tvlct,255,255,255,0
!p.background=0
!p.color=1

plot,x1, h2,$
title="Heightscales, isothermal and where vr>0, eta=0",$
   linestyle=1,$
         xtitle="Radius, R",$
         ytitle="Altitude (Z)",$
   xrange=[xbeg,xend]
oplot,x1,h1,linestyle=0
oplot, rr,zz,linestyle=2 

items=['Isothermal','vr>0', 'eta=0']
legend,items, linestyle=indgen(3)



oplot, rf1,zf1, color=90
oplot, rf2,zf2, color=133
oplot, rf3,zf3, color=233



posalf, b1,b2,b3,v1,v2,rho,pr,mf,ma,ms,mson
contour, mf, x1,x2, levels=1, color=233, /overplot, c_thick=2.0
contour, ma, x1,x2, levels=1, color=133, /overplot, c_thick=2.0 



if ( i eq 0) then begin
device,/close

endif else begin

im=tvread(filename="Y2012height",/png,/nodialog)
endelse 
endfor

; resettign back to x windows
set_plot, 'x'
end

