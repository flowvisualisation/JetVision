
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


vsound=sqrt(pres/den)



rr=rebin(reform(x1,n1,1),n1,n2 )
zz=rebin(reform(x2,n2,1),n2,n1 )
zz=transpose(zz)


rad=sqrt(rr^2 +zz^2)

vsound2=pres/den
vsoundmid=vsound2
  for j=0,n1-1 do begin
     vsoundmid(*,j)=pres(*,0)/den(*,0)
     endfor
     eta=sqrt(rr^3)*(vsoundmid + (2.0/5.0)*(1/rad-1/rr))
eta(*,*) =eta(*,*)*(eta(*,*) gt 0.0)





kvisc=2.0/3.0*eta
;kvisc=eta





nx=40
nx2=100
nx3=400
ybeg=-0.0006
yend=0.0002
xend=2
!P.CHARSIZE=1.0
!P.CHARTHICK=1.0
!P.position=0
!P.MULTI=[0,1,5]
!x.range=[0,4]
!x.style=1
!p.font=0
set_plot,'ps'
device,filename='rho_u_nu_3zones.eps', /encapsulated
device, xsize=16, ysize=16
device, /times
a=den(nx,*)/den(nx,0)
b=den(nx2,*)/den(nx2,0)
c=den(nx3,*)/den(nx3,0)
!y.range=[min([[a],[b],[c]]), max([[a],[b],[c]]) ]
xh=x2/0.1/x1(nx)
plot, xh,a,$
   title="!9r!x at r="+string(x1(nx)),$
   xtitle="Altitude/thermal heightscale, z/h",$
   xstyle=1
xh=x2/0.1/x1(nx2)
oplot, xh,b,linestyle=1 
xh=x2/0.1/x1(nx3)
oplot, xh,c,linestyle=2 
items=['r='+string(x1(nx),format='(F5.1)'),'r='+string(x1(nx2),format='(F5.1)'),'r='+string(x1(nx3),format='(F5.1)')]
legend,items,linestyle=indgen(3), /right

a=vphi(nx,0:10)/x1(nx)
b=vphi(nx2,0:200)/x1(nx2)
c=vphi(nx3,0:500)/x1(nx3)
!y.range=[min([[a],[b],[c]]), max([[a],[b],[c]]) ]
xh=x2/0.1/x1(nx)
plot, xh,a, linestyle=0,title='!9w!X',$
   xtitle="Altitude/thermal heightscale, z/h"
xh=x2/0.1/x1(nx2)
oplot, xh,b,linestyle=1 
xh=x2/0.1/x1(nx3)
oplot, xh,c,linestyle=2 


h1=vsound(nx,0)*x1(nx)^1.5
a=-vr(nx,0:2500)
h2=vsound(nx2,0)*x1(nx2)^1.5
b=-vr(nx2,0:2500)
h3=vsound(nx3,0)*x1(nx3)^1.5
c=-vr(nx3,0:2900)
!y.range=[min([[a],[b],[c]]), 0]
!y.range=0
xh=x2(0:2500)/h1
plot, xh,a, linestyle=0,title='u!Dr!N',$
        xrange=[0,10],$
   xtitle="Altitude/thermal heightscale, z/h"
xh=x2(0:2500)/h2
oplot, xh,b,linestyle=1 
xh=x2(0:2900)/h3
oplot, xh,c,linestyle=2 

h1=vsound(nx,0)*x1(nx)^1.5
a=-vr(nx,0:10)*h1/kvisc(nx,0:10)
h2=vsound(nx2,0)*x1(nx2)^1.5
b=-vr(nx2,0:200)*h2/kvisc(nx2,0:200)
h3=vsound(nx3,0)*x1(nx3)^1.5
c=-vr(nx3,0:500)*h3/kvisc(nx3,0:500)
!y.range=[min([[a],[b],[c]]), 0]
!y.range=0
xh=x2(0:10)/h1
plot, xh,a, linestyle=0,title='hR!De!N/r=u!Dr!N*h/!9n!X!Dv!N',$
   xtitle="Altitude/thermal heightscale, z/h"
xh=x2(0:200)/h2
oplot, xh,b,linestyle=1 
xh=x2(0:500)/h3
oplot, xh,c,linestyle=2 

a=kvisc(nx,*)/kvisc(nx,0)
b=kvisc(nx2,*)/kvisc(nx2,0)
c=kvisc(nx3,*)/kvisc(nx3,0)
!y.range=[min([[a],[b],[c]]), max([[a],[b],[c]]) ]
!y.range=0
xh=x2/0.1/x1(nx)
plot, xh,a, linestyle=0,title='!9n!X!Dv!N'
xh=x2/0.1/x1(nx2)
oplot, xh,b,linestyle=1 
xh=x2/0.1/x1(nx3)
oplot, xh,c,linestyle=2 



device,/close


set_plot,'x'

!P.COLOR=0
!P.BACKGROUND=255
end
