
xbeg=0
xend=712
ybeg=0
yend=1736

contour, alog10(rho(xbeg:xend,ybeg:yend,0)), x1(xbeg:xend),x2(ybeg:yend), /fill, nlev=60, xstyle=1, ystyle=1

alf, b1,b2,b3,v1,v2,rho,pr,mf,ma,ms,mson


contour, mf, x1,x2, levels=1, color=233, /overplot, c_thick=2.0

contour, ma, x1,x2, levels=1, color=133, /overplot, c_thick=2.0

loadct,33
!P.CHARSIZE=0.5
!P.CHARTHICK=1.0
set_plot,'ps'
device,filename='SADjet.eps', /encapsulated, /color,$
	XSIZE=5000/1000.,$
	YSIZE=10000/1000.

d=alog10(rho(xbeg:xend,ybeg:yend,0))
contour, d, $
	x1(xbeg:xend),x2(ybeg:yend),$
	position=[0.1,0.2,0.9,0.9],$
	title="log(!7q!3), S!DA!N, S!DF!N., Time(IDR)="+string(time(nlast)/!PI/2.0),$
	/fill,$
	nlev=60,$
	xstyle=1,$
	ystyle=1

imin=min(d)
imax=max(d)
cbdiv=5
xb0=0.1
yb0=0.06
xb1=0.9
yb1=0.09
color=200
charsize=0.5

 colorbar,/horizontal,position=[xb0,yb0,xb1,yb1],$
            range=[imin,imax], format = '(F8.2)',$
            /right,  divisions = cbdiv,charsize=charsize


xf= x2(n2-1)
xf=4

!P.THICK=1.0
for i=1.6, xf,1 do begin

field_line, b1(0), b2(0), x1,x2,  i,0.05, rr,zz
oplot, rr,zz,color=255
;field_line, b1(0), b2(0), x1,x2,  10,i, rr,zz
;oplot, rr,zz,color=255
endfor



alf, b1,b2,b3,v1,v2,rho,pr,mf,ma,ms,mson


contour, mf, x1,x2, levels=1, color=233, /overplot, c_thick=2.0

contour, ma, x1,x2, levels=1, color=133, /overplot, c_thick=2.0


device,/close


set_plot,'x'



end
