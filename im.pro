
d=alog10(rho(0:511,0:511,0))
xx=x1(0:511)
yy=x2(0:511)


loadct,0
loadct,33
;window, xs=500, ys=500
!P.CHARSIZE=1.

pos=[0.1,0.2,0.9,0.9]
imin=min(d)
imax=max(d)



set_plot,'ps'
device, filename="logden.eps", /encapsulated, /color
DEVICE, BITS_PER_PIXEL=8, COLOR=1
DEVICE, XSIZE=4, YSIZE=4, /INCHES
r=scale_vector(d,0,255)
;display, d, x1=xx, x2=yy, /hbar
tvimage, r,position=pos
contour,r,xx,yy, /nodata, /noerase, position=pos,$
  title='log(!7q!3), B field, S!DA!N, S!DF!N.'+ $
            'Time (IDR)='+string(time(nlast)/2./!PI,format='(F8.2)')

xf= x2(n2-1)
xf=5


for i=1.6, xf,1 do begin

field_line, b1(0), b2(0), x1,x2,  i,0.05, rr,zz
oplot, rr,zz,color=255
;field_line, b1(0), b2(0), x1,x2,  10,i, rr,zz
;oplot, rr,zz,color=255
endfor




    alf, b1,b2,b3,v1,v2,rho,pr,mf,ma,ms,mson
	 contour, mf, x1,x2, levels=1, color=233, /overplot, c_thick=2.0
	 contour, ma, x1,x2, levels=1, color=133, /overplot, c_thick=2.0




pos=[0.1,0.06,0.9,0.12]
colorbar, position=pos, range=[imin,imax],$
   format = '(f8.3)'

device, /close
set_plot,'x'


end
