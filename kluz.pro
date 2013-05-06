
;pload,out=0

rad=60
print, 'Radius is ',x1(rad)

!P.CHARSIZE=2
!P.THICK=2
!x.thick=2.0
!y.thick=2.0
!P.CHARTHICK=2

;set_plot,'ps'
;device, filename="kluz.eps", /encapsulated

pload,out=0
;plot, x2/0.1/x2(rad), (pr(rad,*,0)/pr(rad,0,0)),$
;	xrange=[0,6],$
;;	;/xlog,$
;	/ylog,$
;	title="Pressure, t="+string(time(nlast)/2/!PI,format='(F8.2)') ,$
;	xtitle="Altitude/Heightscale, z/h",$
;	ytitle="Pressure/P(z=0)"

pload,out=nlast
;oplot, x2/0.1/x2(rad), (pr(rad,*,0)/pr(rad,0,0)),linestyle=1


rr=rebin(reform(x1,n1,1),n1,n2 )
zz=rebin(reform(x2,n2,1),n2,n1 )
zz=transpose(zz)

rad_sphere=sqrt(rr^2 + zz^2)

rad=50
cs2=pr(0)/rho(0)
r0=x1(rad)
cs02= cs2(rad,0) + 5./2.*(1/rad_sphere(rad,*) -1/r0 )
plot, x2, cs02, xrange=[0,2]
;oplot, x2/0.1/x2(rad), (exp(-x2^2/(1.3*0.1*x1(rad))^2)),linestyle=2

items=['P, t=0','P, t=106','Gaussian']
legend, items,linestyle=indgen(3),/right 

;device, /close

set_plot,'x'

end
