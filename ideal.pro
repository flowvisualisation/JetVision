
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

h2=x2
for i=0,n1-1 do begin
      a=where(eta(i,*) eq 0)
      print,a(0)
		a=a * (a gt 0)
      h2(i)=x2(a(0))
   endfor

xbeg=0
xend=5
plot,x1,h2,linestyle=0, $
	xrange=[xbeg,xend],$
	title="Ideal, Slow, Sonic surfaces"




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


for i=0,n1-1 do begin
      a=where( ms(i,*) gt 1 )
      print,a(0)
		a=a * (a gt 0)
      h2(i)=x2(a(0))
   endfor

oplot,x1,h2,linestyle=1



for i=0,n1-1 do begin
      a=where( mson(i,*) gt 1 )
      print,a(0)
		a=a * (a gt 0)
      h2(i)=x2(a(0))
   endfor
oplot,x1,h2,linestyle=2

items=['z_ideal','z_sm','z_sonic']
legend, items, linestyle=indgen(3)

im=tvread(filename="ideal", /png, /nodialog)

end
