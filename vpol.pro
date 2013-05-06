
loadct,0

!P.POSITION=0
!P.POSITION=0
!P.MULTI=[0,1,2]
!P.MULTI=0
!P.CHARSIZE=2.0


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
r=x1
z=x2
dr=dx1
dz=dx2



 posalf, b1,b2,b3,v1,v2,rho,pr,mf,ma,ms,mson




!P.FONT=0
set_plot,'ps'
device, filename="vpolvtor.eps" , /encapsulated
device, /time

!P.CHARSIZE=2.0
!P.CHARTHICK=2.0
!P.THICK=3.0


vpol=sqrt(vr^2+vz^2)

i=2.4
field_line, br, bz, x1,x2,  i,0.01, rr,zz
interplot, vpol ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
plot, zpl, qinter1*(i)^0.5 ,$
   ;yrange=[-2,0],$
   xrange=[0.1,200],$
   xstyle=1,$
	/xlog,$
   xtitle="Altitude, Z",$
   ytitle="V!DP!N, r=2.4r!D0!N, t="+string(time(nlast)/2/!PI, format='(F8.2)')



fluff=ms
interplot,fluff ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
a=where(qinter1 gt 1)
slow=a(0)
print, 'slow',slow
print, zpl(a(0))
xs=zpl(a(0))

fluff=ma
interplot,fluff ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
a=where(qinter1 gt 1)
alfven=a(0)
print ,'alfven', alfven
xa=zpl(a(0))


fluff=mf
interplot,fluff ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
a=where(qinter1 gt 1)
fast=a(0)
print,'fast',(fast)
print, zpl(a(0))
xf=zpl(a(0))


nzpl = n_elements(zpl)
oplot,  zpl-zpl+xs , findgen(nzpl)/nzpl*4.
oplot,  zpl-zpl+xa ,  findgen(nzpl)/nzpl*4.
oplot,  zpl-zpl+xf , findgen(nzpl)/nzpl*4.
ybeg=1

;   plots, [xs,xs],[ybeg,yend]
   xyouts, xs*1.1,  0.9*ybeg, 'S!DSM!N'
;   plots, [xa,xa],[ybeg,yend]
   xyouts, xa*1.1,  0.9*ybeg, 'S!DA!N'
;   plots, [xf,xf],[ybeg,0.5*yend]
   xyouts, xf*1.1,  0.9*ybeg, 'S!DFM!N'



device, /close
set_plot,'x'


end
