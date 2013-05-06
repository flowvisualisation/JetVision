

!P.MULTI=0
pres=pr(0)
dens=rho(0)

r=x1
h1=sqrt(pres(*,0)/dens(*,0)) *r^1.5
h2=h1

vr=v1(0)

for i=0,n2-1 do begin
		
		a=where(vr(i,*)>0)
		print,a(0)
		h2(i)=x2(a(0))
	endfor


xbeg=0
xend=20
plot,x1, h2,$
title="Heightscales, iosthermal and where vr>0",$
	linestyle=1,$
	xrange=[xbeg,xend]
oplot,x1,h1,linestyle=0

items=['Isothermal','vr>0']
legend,items, linestyle=indgen(2)

im=tvread(filename="heightscale",/png,/nodialog)
end
