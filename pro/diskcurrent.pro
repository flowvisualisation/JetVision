; Compute current 
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

rbphi=rr*bphi

vsound2=pres/den
vsoundmid=vsound2
  for j=0,n2-1 do begin
     vsoundmid(*,j)=pres(*,0)/den(*,0)
     endfor
     eta=sqrt(rr^3)*(vsoundmid + (2.0/5.0)*(1/rad-1/rr))
eta(*,*) =eta(*,*)*(eta(*,*) gt 0.0)
kvisc=2.0/3.0*den*eta
h2=x2

for i=0,n1-1 do begin
      a=where(eta(i,*) eq 0)
      ;print,a(0)
      a=a * (a gt 0)
      h2(i)=x2(a(0))
   endfor

h2=sqrt(pres(*,0)/den(*,0))*x1^(3./2.)*2.2


; Compute disk height


r=x1
z=x2
dr=dx1
dz=dx2
h=x1-x1+78.0

interplot,rbphi ,r,z,dr,dz,r,h2,qinter1,rpl,zpl


window, 1
!P.POSITION=0
plot, rpl,qinter1, $
	xrange=[1,20],$
	xstyle=1,$
	title="rbphi"

interplot,bphi/bz ,r,z,dr,dz,r,h2,qinter1,rpl,zpl
window, 2
plot, rpl,qinter1, $
	xrange=[1,20],$
	xstyle=1,$
	title="bphi/bz"
end


