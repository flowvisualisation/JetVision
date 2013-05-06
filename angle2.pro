
br=b1(0)
bz=b2(0)

vr=v1(0)
vz=v2(0)




alpha=abs(atan( bz/br))

beta=abs(atan(vz/vr))


diff=(alpha-beta)*180.0/!PI


nx=40

r=x1
z=x2
dr=dx1
dz=dx2
field_line, br, bz, x1,x2,  x1(nx), 1e-5, rr,zz
interplot,diff ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl




!P.MULTI=0
!P.CHARSIZE=1.5

xbeg=8e-2
xend=10
ybeg=10
yend=10
set_plot, 'x'
;plot, x2,(alpha(nx,*)-beta(nx,*))*180.0/!PI,$
plot, zpl, qinter1,$
;psym=2,$
title='Angle between B and u, r='+string(x1(nx)),$
ytitle='Angle [degree]',$
xrange=[xbeg,xend],$
/xlog,$
xstyle=1,$
yrange=[ybeg,yend]
oplot, zpl, qinter1,psym=4

;plot, alpha-beta

im=tvread(filename='Angle',/png,/nodialog)

set_plot,'ps'
device,filename='angle.eps', /encapsulated

plot, x2,(alpha(nx,*)-beta(nx,*))*180.0/!PI,$
;psym=2,$
title='Angle between B!DP!N and u!DP!N, r='+string(x1(nx)),$
ytitle='Angle [degree]',$

xrange=[xbeg,xend],$
yrange=[ybeg,yend]
device,/close
set_plot,'x'
end
