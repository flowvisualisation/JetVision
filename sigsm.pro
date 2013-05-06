


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
zzzz=transpose(rrr)

rad=sqrt(rrr^2 + zzzz^2)
;michel=zzzz

sound2=(pr(0)/rho(0))
sound=sqrt(pr(0)/rho(0))



den=rho(0)
pres=pr(0)
br=b1(0)
bz=b2(0)
temp=pres/den

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
      a=a * (a gt 0)
      h2(i)=x2(a(0))
   endfor


zz=h2

;zz=2.1*sound(*,0)*rr^1.5+0.5


posalf, b1,b2,b3,v1,v2,rho,pr,mf,ma,ms,mson


for i=0,n1-1 do begin
      a=where(ms (i,*) gt 1)
      print,a(0)
      a=a * (a gt 0)
      h2(i)=x2(a(0))
   endfor



;window,2


interplot,michel ,r,z,dr,dz,rr,h2,qinter1,rpl,zpl
qsig=qinter1
interplot,poynting ,r,z,dr,dz,rr,h2,qinter1,rpl,zpl
qpoy=qinter1
interplot,kinetic ,r,z,dr,dz,rr,h2,qinter1,rpl,zpl
qkin=qinter1
;plot, rr,qinter1

;display,alog10(rho(0)),x1=x1,x2=x2
;display, michel ,x1=x1,x2=x2, /vbar
;oplot,rr,zz


set_plot,'x'
!P.COLOR=0
!P.BACKGROUND=255
window,2,xs=600,ys=800
loadct,0
!P.POSITION=0
!P.CHARSIZE=2
!P.MULTI=[0,1,4]
plot, rpl,qsig, xrange=[1,10],xstyle=1,$
title="Poynting to kinetic flux ratio, !7r!3"
plot, rpl,qpoy, xrange=[1,10],xstyle=1,$
title="Poynting flux"
plot, rpl,qkin, xrange=[1,10],xstyle=1,$
title="Kinetic flux"
plot, rpl,zpl, xrange=[1,10],xstyle=1,$
title="Slow"

im=tvread(filename='sigsm',/nodialog,/png)



!P.POSITION=0
!P.MULTI=0
set_plot,'ps'
device, /encapsulated, filename="sigsm.eps"
plot,rpl, qsig,$
	;/xlog,$
   xrange=[1,10],xstyle=1,$
	xtitle="Radius, R",$
	ytitle="!7r!3!DSM!N"
device,/close
set_plot,'x'

end

