
set_plot, 'ps'
;window,22, xs=1000, ys=900
loadct,0


; K

den=rho(0)
pres=pr(0)

vr=v1(0)
vz=v2(0)
vphi=v3(0)

br=b1(0)
bz=b2(0)
bphi=b3(0)



rr=rebin(reform(x1,n1,1),n1,n2 )
zz=rebin(reform(x2,n2,1),n2,n1 )
zz=transpose(zz)

rad=sqrt(rr^2 +zz^2)

vsound2=pres/den
vsoundmid=vsound2
  for j=0,n1-1 do begin
     vsoundmid(*,j)=pres(*,0)/den(*,0)
     endfor
     eta=sqrt(rr^3)*(vsoundmid + (2.0/5.0)*(1/rad-1/rr))
eta(*,*) =eta(*,*)*(eta(*,*) gt 0.0)
kvisc=2.0/3.0*den*eta




va=sqrt(br^2+bz^2)/sqrt(den)
va2=va^2

vs2=((5./3.)* pres /den)
vs=sqrt(vs2)

btot2 = br^2 +  bz^2 + bphi^2

cfast2 = 0.5*( vs2 +btot2/den + sqrt( (vs2 +btot2/den)^2  -4.0*vs2*va2 )  )

cslow2 = 0.5*( vs2 +btot2/den - sqrt( (vs2 +btot2/den)^2  -4.0*vs2*va2 )  )
cfast=sqrt(cfast2)
cslow=sqrt(cslow2)

vr=vr*(vr gt 0)
vz=vz*(vz gt 0)
vpol=sqrt(vr^2 + vz^2)


mf=vpol/cfast
ma=vpol/va
ms=vpol/cslow
mson=vpol/vs


!P.CHARSIZE=1

 contour, mf, x1,x2, levels=1, xrange=[0,40], yrange=[0,40],$
 	xtitle="Radius, R",$
 	ytitle="Altitude, Z"
  contour, ma, x1,x2, levels=1, xrange=[0,40], yrange=[0,40], $
	 	/overplot,$
		c_linestyle=1
   contour, ms, x1,x2, levels=1, xrange=[0,40], yrange=[0,40],$ 
	 	/overplot,$
		c_linestyle=2
	 contour, eta, x1,x2, levels=0.00001, xrange=[0,40], yrange=[0,40], $
	 	/overplot,$
		c_linestyle=3

	items=['S!DFM!N','S!DA!N','S!DSM!N','disk']
	legend, items, linestyle=indgen(4)

!P.font=0
set_plot,'ps'
device, /times
device,/encapsulated, filename="surfaces.eps"

 contour, mf, x1,x2, levels=1, xrange=[0,40], yrange=[0,40],$
 	xtitle="Radius, R",$
 	ytitle="Altitude, Z"
  contour, ma, x1,x2, levels=1, xrange=[0,40], yrange=[0,40], $
	 	/overplot,$
		c_linestyle=1
   contour, ms, x1,x2, levels=1, xrange=[0,40], yrange=[0,40],$ 
	 	/overplot,$
		c_linestyle=2
	 contour, eta, x1,x2, levels=0.00001, xrange=[0,40], yrange=[0,40], $
	 	/overplot,$
		c_linestyle=3

field_line, br, bz, x1, x2, 1.0, 0.01, rr, zz
oplot, rr,zz
field_line, br, bz, x1, x2, 5.3, 0.01, rr, zz
oplot, rr,zz
field_line, br, bz, x1, x2, 12.9, 0.01, rr, zz
oplot, rr,zz
xyouts, 6,10, "Zone I"
xyouts, 22,10, "Zone II"
xyouts, 34,10, "Zone III"
	items=['S!DFM!N','S!DA!N','S!DSM!N','disk']
	legend, items, linestyle=indgen(4)
device,/close
set_plot,'x'

end
