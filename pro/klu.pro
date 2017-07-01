
loadct,33
display, alog10(rho(0)), x1=x1,x2=x2,title=string(time(nlast))


loadct,0
window,2

plot, x1, v1(*,0,0),xrange=[0,20],xstyle=1, title='v1'

;window,3
;plot, x1, (b1(*,0,0)^2 + b1(*,0,0)^2)/pr(*,0,0),xrange=[0,20],xstyle=1
end
