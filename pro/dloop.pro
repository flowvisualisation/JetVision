
for i=0,nlast do begin
   pload,out=i
   display, alog10(rho(0)), ims=2, /vbar,x1=x1,x2=x2,$
	title='T='+string(time(i))
	wait,0.1
endfor
end 

