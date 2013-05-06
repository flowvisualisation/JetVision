

pro posalf, b1,b2,b3,v1,v2,rho,pr,mf,ma,ms,mson
va=sqrt(b1(0)^2+b2(0)^2)/sqrt(rho(0))
va2=va^2

vs2=((5./3.)* pr(0) /rho(0))
vs=sqrt(vs2)

btot2 = b1(0)^2 +  b2(0)^2 + b3(0)^2 

cfast2 = 0.5*( vs2 +btot2/rho(0) + sqrt( (vs2 +btot2/rho(0))^2  -4.0*vs2*va2 )  )

cslow2 = 0.5*( vs2 +btot2/rho(0) - sqrt( (vs2 +btot2/rho(0))^2  -4.0*vs2*va2 )  )
cfast=sqrt(cfast2)
cslow=sqrt(cslow2)

vpol=sqrt(v1(0)^2 + v2(0)^2)

vr=v1(0)*(v1(0) gt 0)
vz=v2(0)*(v2(0) gt 0)
vpol=sqrt(vr^2+vz^2)
;vpol=sqrt(v1(0)^2 + v2(0)^2)


;vpol = (v1(0)*b1(0)+v2(0)*b2(0))/sqrt(b1(0)^2 +  b2(0)^2)

mf=vpol/cfast
ma=vpol/va
ms=vpol/cslow
mson=vpol/vs

print, min(abs(mf)),min(mf), max(mf)
print, min(abs(ma)),min(ma), max(ma)
print, min(abs(ms)),min(ms), max(ms)
end
