
loadct,33
tvlct,0,0,0,0
tvlct,255,255,255,1
!p.color=0
!p.background=1


rr=rebin(reform(x1,n1,1),n1,n2 )
zz=rebin(reform(x2,n2,1),n2,n1 )
zz=transpose(zz)

window, 21, xs=500, ys=800

contour, rbphi, rr,zz, xrange=[0,20], yrange=[0,80], nlev=20

field_line, br, bz, x1,x2,  10,0.05, rf,zf
oplot,rf,zf, color=33

field_line, br, bz, x1,x2,  5,0.05, rf,zf
oplot,rf,zf, color=33



oplot, x1, hmax, color=233

window, 22
plot, x2, jr(14*15,*), /xlog
end
