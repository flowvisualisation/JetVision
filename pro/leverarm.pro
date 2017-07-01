
;
; Plot the lever arm, using interplot to follow a field line until
; the Alfven surface is reached
; The alfven surface  is defined as the place where the poloidal velocity
; becomes greater than the alfven speed, but it must be the outgoing 
; poloidal velocity where vz> vr & vz > 0
;



poynting=v3(*,*,0)*b3(*,*,0)*sqrt(b1(*,*,0)^2 + b2(*,*,0)^2)
poynting = poynting- b3(*,*,0)^2*sqrt(v1(*,*,0)^2 +v2(*,*,0)^2)
poynting = -poynting

;plot, x1, poynting(*,0), title="Poynting",$
;  /xlog, $
;  xrange=[1,10]


kinetic=(0.5*rho(*,*,0)*(v1(*,*,0)^2 $
+v2(*,*,0)^2+v3(*,*,0)^2)*sqrt(v1(*,*,0)^2 +v2(*,*,0)^2))

michel=poynting/kinetic

r=x1
z=x2
dr=dx1
dz=dx2
rr=x1

field_line, b1(0), b2(0), x1,x2,  2,0.05, rr,zz


interplot,michel ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl




end
