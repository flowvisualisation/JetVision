

!P.CHARSIZE=1.5
!P.THICK=2
!x.thick=2.0
!y.thick=2.0
!P.CHARTHICK=2



set_plot,'ps'
device, filename="denlaunchreg.eps", /encapsulated
rad1=30
rad2=70
rad3=180
rad=rad1
plot, x2/0.1/x1(rad), rho(rad,*,0)/ rho(rad,0,0), $
	/ylog,$
	xtitle="Altitude/Thermal Heightscale, z/h",$
	ytitle="Normalised Density",$
	xrange=[0,5]
rad=rad2
oplot, x2/0.1/x1(rad), rho(rad,*,0)/ rho(rad,0,0), linestyle=1
rad=rad3
oplot, x2/0.1/x1(rad), rho(rad,*,0)/ rho(rad,0,0), linestyle=2

items=[											$
'r='+string(x1(rad1),format='(F5.1)'), $
'r='+string(x1(rad2),format='(F5.1)'),	$
'r='+string(x1(rad3),format='(F5.1)') 	$
]
legend, items, linestyle=indgen(3),/right
device, /close


temp=pr(0)/rho(0)

set_plot,'ps'
device, filename="templaunchreg.eps", /encapsulated
rad=50
plot, x2/0.1/x1(rad), temp(rad,*)/ temp(rad,0), $
	/ylog,$
	xtitle="Altitude/Thermal Heightscale, z/h",$
	ytitle="Normalised Temperature, T/T(z=0)",$
	xrange=[0,5]
rad=150
oplot, x2/0.1/x1(rad), temp(rad,*)/ temp(rad,0), linestyle=1
items=['r='+string(x1(50),format='(F5.1)'),'r='+string(x1(150),format='(F5.1)')]
legend, items, linestyle=indgen(2),/right
device, /close


temp=pr(0)

set_plot,'ps'
device, filename="preslaunchreg.eps", /encapsulated
rad=50
plot, x2/0.1/x1(rad), temp(rad,*)/ temp(rad,0), $
	/ylog,$
	xtitle="Altitude/Thermal Heightscale, z/h",$
	ytitle="Normalised Pressure, T/T(z=0)",$
	xrange=[0,5]
rad=150
oplot, x2/0.1/x1(rad), temp(rad,*)/ temp(rad,0), linestyle=1
items=['r='+string(x1(50),format='(F5.1)'),'r='+string(x1(150),format='(F5.1)')]
legend, items, linestyle=indgen(2),/right
device, /close


set_plot,'ps'
device, filename="vradlaunchreg.eps", /encapsulated
rad=30
temp=v1(0)/sqrt(pr(rad,0,0)/rho(rad,0,0) )
plot, x2/0.1/x1(rad), temp(rad,*), $
	xtitle="Altitude/Thermal Heightscale, z/h",$
	ytitle="Normalised Radial Velocity, v!DR!N/v(z=0)",$
	xrange=[0,5]
rad=70
temp=v1(0)/sqrt(pr(rad,0,0)/rho(rad,0,0) )
oplot, x2/0.1/x1(rad), temp(rad,*), linestyle=1
rad=180
temp=v1(0)/sqrt(pr(rad,0,0)/rho(rad,0,0) )
oplot, x2/0.1/x1(rad), temp(rad,*), linestyle=2
items=['r='+string(x1(30),format='(F5.1)'),'r='+string(x1(70),format='(F5.1)'),'r='+string(x1(180),format='(F5.1)')]
legend, items, linestyle=indgen(3)
device, /close


temp=v2(0)
set_plot,'ps'
device, filename="vvertlaunchreg.eps", /encapsulated
rad=50
temp=v2(0)/sqrt(pr(rad,0,0)/rho(rad,0,0) )
plot, x2/0.1/x1(rad), temp(rad,*), $
	xtitle="Altitude/Thermal Heightscale, z/h",$
	ytitle="Normalised Vertical Velocity, v!DZ!N/v(z=0)",$
	xrange=[0,5]
rad=150
temp=v2(0)/sqrt(pr(rad,0,0)/rho(rad,0,0) )
oplot, x2/0.1/x1(rad), temp(rad,*), linestyle=1
items=['r='+string(x1(50),format='(F5.1)'),'r='+string(x1(150),format='(F5.1)')]
legend, items, linestyle=indgen(2)
device, /close


set_plot,'x'



end
