
nend=nlast
pload,out=nend


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
;      a=where(eta(i,*) eq 0)
;      print,a(0)
;      a=a * (a gt 0)
;      h2(i)=x2(a(0))
   endfor


zz=h2

zz=2.*sound(*,0)*rr^1.5+0.05




;window,2


interplot,michel ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
;plot, rr,qinter1

;display,alog10(rho(0)),x1=x1,x2=x2
;display, michel ,x1=x1,x2=x2, /vbar
;oplot,rr,zz


window,2,xs=600,ys=800
loadct,0
!P.POSITION=0
!P.CHARSIZE=2
!P.MULTI=[0,1,3]
plot, rpl,qinter1, xrange=[1,10],xstyle=1,$
title="Poynting to kinetic flux ratio, !7r!3"
oplot,  rpl,qinter1, psym=1
plot, rr,zz, xrange=[1,10],$
title="Disk Surface Height, z"
cut=80
plot, x2, michel(cut,*), $
title="Vertical cut of  !7r!3, at r="+string(x1(cut)),$
xrange=[1,10]


n=nend
ll=6
zero=''
nts=strcompress(string(n),/remove_all)
lnt=strlen(nts)
for j=1,ll-lnt do zero=zero+'0'
	        fname='poyntingtokinetic.'+zero+nts
			  im=tvread(filename=fname,/nodialog,/png)


;window,3
;plot,x2,kinetic(80,*)


end

