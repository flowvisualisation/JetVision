;display, alog10(rho(0)), x1=x1,x2=x2, /hbar
window,1
loadct,0
pload,out=nlast
!P.position=0
!P.multi=[0,2,1]
plot,X1, (b1(*,0,0)^2 +  b2(*,0,0)^2 )/pr(*,0,0), $
;/xlog,$
/ylog,$
xrange=[1,20],title='!7l!3=b2/pr',charsize=2
plot,X1, (b1(*,0,0)^2 +  b2(*,0,0)^2 ), $
/xlog,$
/ylog,$
xrange=[1,20],title='b2'

;r=tvread(filename="aslope",/png,/nodialog)
!P.multi=0


loadct,33
;pload,out=0
display, alog10( rho(0) ), x1=x1,x2=x2, /vbar
q1=b1(0)
q2=b2(0)
xa=x1
xb=x2
size=x1(n1-1)
for i=0.5,size-0.5,2 do begin

;field_line, b1(0), b2(0), x1,x2,  i,0.05, rr,zz
field_line, q1,q2, xa,xb,  size-i,0.05, rr,zz
oplot, rr,zz,color=100
;field_line, q1,q2, xa,xb,  size-i, size/2, rr,zz
;oplot, rr,zz,color=100
;field_line, b1(0), b2(0), x1,x2,  0.01,i, rr,zz
;field_line, q1,q2, xa,xb,  0.01,i, rr,zz
;oplot, rr,zz,color=100
endfor
;r=tvread(/png,/nodialog)
end
