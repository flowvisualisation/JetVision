
br=b1(0)
bz=b2(0)

vr=v1(0)
vz=v2(0)




alpha=atan( bz/br)

beta=atan(vz/vr)

!P.MULTI=0

nx=40
xbeg=0
xend=1
ybeg=-90
yend=90
plot, x2,(alpha(nx,*)-beta(nx,*))*180.0/!PI,$
;psym=2,$
title='Angle between B and u, r='+string(x1(nx)),$
ytitle='Angle [degree]',$
xrange=[xbeg,xend],$
yrange=[ybeg,yend]

;plot, alpha-beta

im=tvread(filename='Angle',/png,/nodialog)

end
