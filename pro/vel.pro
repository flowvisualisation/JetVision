
!P.MULTI=[0,1,5]
!P.CHARSIZE=2.5
!P.CHARTHICK=2.5


pload,out=nlast


nheight=2760

r=x1*0.1

plot, r, 94*v1(*,nheight,0),$
	xtitle="R, [AU]",$
	ytitle="vR, [km/s]",$
	title="vR, z="+string(x2(nheight))
plot, r, 94*v2(*,nheight,0),$
	xtitle="R, [AU]",$
	ytitle="vZ, [km/s]",$
	title="vZ"+string(x2(nheight))
plot, r, 94*v3(*,nheight,0),$
	xtitle="R, [AU]",$
	ytitle="v!7U!3, [km/s]",$
	title="v!7U!3"+string(x2(nheight))
plot, r, alog10(rho(*,nheight,0)),$
	xtitle="R, [AU]",$
	title="dens"+string(x2(nheight))
plot, r, alog10(pr(*,nheight,0)/rho(*,nheight,0)),$
	xtitle="R, [AU]",$
	title="T"+string(x2(nheight))

im=tvread(filename="velcuts", /png, /nodialog)


	end
