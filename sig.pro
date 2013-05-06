
nend=nlast
pload,out=nend

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
      print,a(0)
;      a=a-1
      a=a+1
      a=a * (a gt 0)
      h2(i)=x2(a(0))
   endfor


zz=smooth(h2,5)

;zz=2.*sound(*,0)*rr^1.5+0.05

interplot,michel ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
interplot,mu ,r,z,dr,dz,rr,zz,qinter2,rpl,zpl

loadct,33
display,alog10(rho(0)),x1=x1,x2=x2, xrange=[1,10], yrange=[0,5], ims=4
oplot,rr,zz


loadct,0
!p.color=0
!p.background=255

!p.font=0
set_plot,'ps'
device, filename="sigmu.eps"
device, /encapsulated
device, /times
loadct,0
!P.POSITION=0
!P.CHARSIZE=1
!P.MULTI=[0,1,2]
plot, rpl,smooth(qinter1,20,/edge_truncate), $
xrange=[1.5,10],$
xtitle="Radius, R",$
xstyle=1,$
/xlog,$
ytitle="!9s!X!U+!N"

plot, rpl, smooth(qinter2,10,/edge_truncate)/2, $
xtitle="Radius, R",$
/xlog,$
ytitle="!9m!X!U+!N",$
xstyle=1,$
xrange=[1.5,10]
device, /close



set_plot, 'x'
end

