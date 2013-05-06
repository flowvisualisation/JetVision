
!P.THICK=1.0

loadct,33
nbeg=99
nend=nlast
nstep=1
for n=nbeg,nend,nstep do begin
pload,out=n
display, alog10(rho(0)),x1=x1,x2=x2, ims=0.15,$
	/vbar, /gst, $
	;xrange=[0,40],$
	;yrange=[0,10],$
	title='log(!7q!3), B field, S!DA!N, S!DF!N.'+ $
	'Time (IDR)='+string(time(n)/2./!PI,format='(F8.3)'),$
	label1='r',$
	label2='z'
	


xf= x2(n2-1)
xf=10

for i=1.6, xf,1 do begin

field_line, b1(0), b2(0), x1,x2,  i,0.05, rr,zz
oplot, rr,zz,color=255
;field_line, b1(0), b2(0), x1,x2,  10,i, rr,zz
;oplot, rr,zz,color=255
endfor

 alf, b1,b2,b3,v1,v2,rho,pr,mf,ma,ms



contour, ma, x1,x2, levels=1, color=133, /overplot, c_thick=2.0

contour, mf, x1,x2, levels=1, color=233, /overplot, c_thick=2.0

ll=6
zero=''
nts=strcompress(string(n),/remove_all)
lnt=strlen(nts)
for j=1,ll-lnt do zero=zero+'0'
   fname='SADJED.'+zero+nts
	im=tvread(filename=fname,/nodialog,/png)

endfor

end
