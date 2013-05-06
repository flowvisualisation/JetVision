
nend=nlast
pload,out=nend


poynting=v3(*,*,0)*b3(*,*,0)*sqrt(b1(*,*,0)^2 + b2(*,*,0)^2) 
poynting = poynting- b3(*,*,0)^2*sqrt(v1(*,*,0)^2 +v2(*,*,0)^2)
poynting = -poynting

;plot, x1, poynting(*,0), title="Poynting",$
;  /xlog, $
;  xrange=[1,10]


kinetic=(0.5*rho(*,*,0)*(v1(*,*,0)^2 $
+v2(*,*,0)^2+v3(*,*,0)^2    )*sqrt(v1(*,*,0)^2 +v2(*,*,0)^2))
thermal= (2.5)*(pr(*,*,0))* sqrt(v1(*,*,0)^2 +v2(*,*,0)^2)

michel=poynting/(kinetic +thermal )

minmichel=fltarr(n1)



r=x1
z=x2
dr=dx1
dz=dx2
rr=x1
zz=0.1*x1

rrr=rebin(reform(r,n1,1),n1,n1 )
zzzz=transpose(rrr)
;michel=zzzz

sound2=(pr(0)/rho(0))
sound=sqrt(pr(0)/rho(0))
;interplot,michel ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl

zz=2.1*sound(*,0)*rr^1.5+0.5
!P.MULTI=0
loadct,33
display, alog10(rho(0)), x1=x1,x2=x2,$
   title='log(!7q!3)'+ $
	   'Time (IDR)='+string(time(n)/2./!PI,format='(F8.3)'),$
	xrange=[1,40],yrange=[1,5],$
	ims=1, /hbar
oplot, rr,zz, color=200, thick=3, linestyle=2
;!P.MULTI=[0,1,3]


n=nend
ll=6
zero=''
nts=strcompress(string(n),/remove_all)
lnt=strlen(nts)
for j=1,ll-lnt do zero=zero+'0'
	        fname='weaklymagnetisedSAD.'+zero+nts
			  im=tvread(filename=fname,/nodialog,/png)
;plot, x1,qinter1



;window,2


interplot,michel ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
;plot, rr,qinter1

;display,alog10(rho(0)),x1=x1,x2=x2
;display, michel ,x1=x1,x2=x2, /vbar
;oplot,rr,zz

!P.color=0
!p.background=255
window,2,xs=600,ys=700
loadct,0
!P.POSITION=0
!P.CHARSIZE=2
!P.MULTI=[0,1,3]
plot, rr,qinter1, xrange=[1,10],xstyle=1,$
title="Poynting to kinetic flux ratio, !7r!3, at disk surface"
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

