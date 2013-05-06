
d=alog10(rho(0:511,0:511,0))
xx=x1(0:511)
yy=x2(0:511)


loadct,0
loadct,33
!P.CHARSIZE=1.

pos=[0.15,0.25,0.95,0.95]
imin=min(d)
imax=max(d)

!p.position=0
!p.multi=0
!P.font=0
set_plot,'ps'
device, filename="current.eps", /encapsulated, /color
device, /times
DEVICE, BITS_PER_PIXEL=8, COLOR=1
DEVICE, XSIZE=4, YSIZE=4, /INCHES
r=scale_vector(d,0,255)
;display, d, x1=xx, x2=yy, /hbar
tvimage, r,position=pos
contour,r,xx,yy, /nodata, /noerase, position=pos,$
xtitle="Radius, R",$
ytitle="Altitude, Z"
;  title='log(!9r!X), B field, S!DA!N, S!DF!N.'+ $
 ;           'Time (IDR)='+string(time(nlast)/2./!PI,format='(F8.2)')

xf= x2(n2-1)
xf=5


!P.BACKGROUND=255

current, b3, x1, rbphi,n2

minc=min(rbphi)
maxc=max(rbphi)
num=32
lev=findgen(32)*(maxc-minc)/32+minc
for i=0,31 do begin
contour,rbphi,x1,x2,/overplot,levels=lev(i)
endfor



    alf, b1,b2,b3,v1,v2,rho,pr,mf,ma,ms,mson
	 contour, mf, x1,x2, levels=1, color=233, /overplot, c_thick=2.0, c_linestyle=1
	 contour, ma, x1,x2, levels=1, color=133, /overplot, c_thick=2.0, c_linestyle=2




pos=[0.15,0.06,0.95,0.10]
colorbar, position=pos, range=[imin,imax],$
   format = '(f8.1)'

device, /close
set_plot,'x'


end
