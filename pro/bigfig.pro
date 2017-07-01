
pload,out=0
d=alog10(rho(0:511,0:511,0))
xx=x1(0:511)
yy=x2(0:511)


loadct,0
loadct,33
window, xs=500, ys=500
!P.CHARSIZE=1.




set_plot,'ps'
device, filename="bigfig.eps", /encapsulated, /color
DEVICE, BITS_PER_PIXEL=8, COLOR=1
DEVICE, XSIZE=8, YSIZE=8, /INCHES

pos=[0.1,0.6,0.45,0.95]
pload,out=0
d=alog10(rho(0:511,0:511,0))
imin=min(d)
imax=max(d)
r=scale_vector(d,0,255)
;display, d, x1=xx, x2=yy, /hbar
tvimage, r,position=pos
contour,r,xx,yy, /nodata, /noerase, position=pos, ytitle='(a)'; ,$
xf= x2(n2-1)
xf=5
for i=1.6, xf,1 do begin
;field_line, b1(0), b2(0), x1,x2,  i,0.05, rr,zz
;oplot, rr,zz,color=255
endfor
;
    posalf, b1,b2,b3,v1,v2,rho,pr,mf,ma,ms,mson
	 contour, mf, x1,x2, levels=1, color=233, /overplot, c_thick=2.0
	 contour, ma, x1,x2, levels=1, color=133, /overplot, c_thick=2.0
;
pos=[0.1,0.53,0.45,0.56]
colorbar, position=pos, range=[imin,imax],$
   format = '(f8.3)'



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pos=[0.6,0.6,0.95,0.95]
pload,out=67
d=alog10(rho(0:511,0:511,0))
imin=min(d)
imax=max(d)
r=scale_vector(d,0,255)
tvimage, r,position=pos
contour,r,xx,yy, /nodata, /noerase, position=pos, ytitle='(b)'; ,$
;  title='log(!7q!3), B field, S!DA!N, S!DF!N.'+ $
;            'Time (IDR)='+string(time(nlast)/2./!PI,format='(F8.2)')
xf= x2(n2-1)
xf=5
for i=1.6, xf,1 do begin
;field_line, b1(0), b2(0), x1,x2,  i,0.05, rr,zz
;oplot, rr,zz,color=255
endfor
;
    posalf, b1,b2,b3,v1,v2,rho,pr,mf,ma,ms,mson
	 contour, mf, x1,x2, levels=1, color=233, /overplot, c_thick=2.0
	 contour, ma, x1,x2, levels=1, color=133, /overplot, c_thick=2.0
;
pos=[0.6,0.53,0.95,0.56]
colorbar, position=pos, range=[imin,imax],$
   format = '(f8.3)'


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pos=[0.1,0.1,0.45,0.45]
pload,out=nlast
d=alog10(rho(0:511,0:511,0))
imin=min(d)
imax=max(d)
r=scale_vector(d,0,255)
tvimage, r,position=pos
contour,r,xx,yy, /nodata, /noerase, position=pos, ytitle='(c)'; ,$
xf= x2(n2-1)
xf=5
for i=1.6, xf,1 do begin
;field_line, b1(0), b2(0), x1,x2,  i,0.05, rr,zz
;oplot, rr,zz,color=255
endfor
;
    posalf, b1,b2,b3,v1,v2,rho,pr,mf,ma,ms,mson
	 contour, mf, x1,x2, levels=1, color=233, /overplot, c_thick=2.0
	 contour, ma, x1,x2, levels=1, color=133, /overplot, c_thick=2.0
;
pos=[0.1,0.03,0.45,0.06]
colorbar, position=pos, range=[imin,imax],$
   format = '(f8.3)'


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pos=[0.6,0.1,0.95,0.45]
pload,out=nlast
d=alog10(rho(0:511,0:511,0))
imin=min(d)
imax=max(d)
r=scale_vector(d,0,255)
tvimage, r,position=pos
contour,r,xx,yy, /nodata, /noerase, position=pos, ytitle='(d)'; ,$
xf= x2(n2-1)
xf=5
for i=1.6, xf,1 do begin
;field_line, b1(0), b2(0), x1,x2,  i,0.05, rr,zz
;oplot, rr,zz,color=255
endfor
;
    posalf, b1,b2,b3,v1,v2,rho,pr,mf,ma,ms,mson
	 contour, mf, x1,x2, levels=1, color=233, /overplot, c_thick=2.0
	 contour, ma, x1,x2, levels=1, color=133, /overplot, c_thick=2.0
;
pos=[0.6,0.03,0.95,0.06]
colorbar, position=pos, range=[imin,imax],$
   format = '(f8.3)'


device, /close
set_plot,'x'


end
