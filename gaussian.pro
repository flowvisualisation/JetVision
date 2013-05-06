
;pload,out=0

rad=60
print, 'Radius is ',x1(rad)

!P.CHARSIZE=2
!P.THICK=2
!x.thick=2.0
!y.thick=2.0
!P.CHARTHICK=2

set_plot,'ps'
device, filename="pressure.eps", /encapsulated

pload,out=0
plot, x2/0.1/x2(rad), (pr(rad,*,0)/pr(rad,0,0)),$
	xrange=[0,6],$
	;/xlog,$
	/ylog,$
	title="Pressure, t="+string(time(nlast)/2/!PI,format='(F8.2)') ,$
	xtitle="Altitude/Heightscale, z/h",$
	ytitle="Pressure/P(z=0)"

pload,out=nlast
oplot, x2/0.1/x2(rad), (pr(rad,*,0)/pr(rad,0,0)),linestyle=1
oplot, x2/0.1/x2(rad), (exp(-x2^2/(1.3*0.1*x1(rad))^2)),linestyle=2

items=['P, t=0','P, t=106','Gaussian']
legend, items,linestyle=indgen(3),/right 

device, /close

set_plot,'x'

end
