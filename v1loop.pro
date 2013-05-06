
window,2
pload,out=0
 plot, x1,v1(*,0,0), xrange=[0.5,20], xstyle=1
for i=1,nlast,5 do begin
	pload,out=i
	oplot, x1,v1(*,0,0)
endfor
window,3
 plot, x1,v1(*,0,0), xrange=[0.5,20], xstyle=1

end
