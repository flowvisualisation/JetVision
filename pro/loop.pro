
loadct,33

nbeg=99
nend=nbeg
for n=nbeg,nend,1 do begin
pload,out=n
display, alog10(rho(0)),x1=x1,x2=x2,$
	ims=0.15,$
;	xrange=[0,40],$
;	yrange=[0,40],$
	/vbar, $
	/gst,$
	title='log(!7q!3) S!DF!N,S!DA!N  t='+string(time(n)/2/!PI,format='(F8.3)'),$
	label1='r',$
label2='z'
	

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



contour, ma, x1,x2, levels=1, color=133, /overplot, c_thick=2.0
;contour, ms, x1,x2, levels=1, color=66, /overplot, c_thick=2.0

;contour, mson, x1,x2, levels=1, color=22, /overplot, c_thick=2.0
contour, mf, x1,x2, levels=1, color=233, /overplot, c_thick=2.0


ll=6
zero=''
nts=strcompress(string(n),/remove_all)
lnt=strlen(nts)
for j=1,ll-lnt do zero=zero+'0'
   fname='Surface.SADJED.'+zero+nts
	im=tvread(filename=fname,/nodialog,/png)

endfor

end
