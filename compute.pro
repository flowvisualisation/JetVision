pro compute
  COMMON GRID
  COMMON PLUTO_VAR
  COMMON TIME_OUT
  COMMON FNAMES



common prim,  br,bz,bphi, vr,vz,vphi, den, pres
common geom, rr,zz, rad, omega, omegamatrix
common current, jr,jz,jphi, fr,fz,fphi
common visc, trphi, tzphi, divTphi

; Calculate the torques
loadct,0
!P.COLOR=0
!P.BACKGROUND=255

;window,22, xs=1000, ys=700
loadct,0


pload,out=nlast
; K

den=rho(0)
pres=pr(0)

vr=v1(0)
vz=v2(0)
vphi=v3(0)

br=b1(0)
bz=b2(0)
bphi=b3(0)


omega=v3(*,0)/x1


omegamatrix=rebin(reform(omega,n1,1),n1,n2 )

rr=rebin(reform(x1,n1,1),n1,n2 )
zz=rebin(reform(x2,n2,1),n2,n1 )
zz=transpose(zz)


rad=sqrt(rr^2 +zz^2)

rbphi=rr*bphi

 drBphidr=pres
 dBphidr=pres
 phi=-1/sqrt(rr^2+zz^2)

  for j=0,n2-1 do begin
         drBphidr(*,j)= deriv(x1,rbphi(*,j))
         dBphidr(*,j)= deriv(x1,bphi(*,j))
  endfor

 dBzdr=pres
  for j=0,n2-1 do begin
         dBzdr(*,j)= deriv(x1,bz(*,j))
  endfor

 dprdr=pres
  for j=0,n2-1 do begin
         dprdr(*,j)= deriv(x1,pres(*,j))
  endfor

 dphidr=pres
  for j=0,n2-1 do begin
         dphidr(*,j)= deriv(x1,phi(*,j))
  endfor


 dvphidr=pres
  for j=0,n2-1 do begin
         dvphidr(*,j)= deriv(x1,vphi(*,j))
  endfor

 dvphidz=pres
  for i=0,n1-1 do begin
         dvphidz(i,*)= deriv(x2,vphi(i,*))
  endfor

 dvrdr=pres
  for j=0,n2-1 do begin
         dvrdr(*,j)= deriv(x1,vr(*,j))
  endfor
 dvrdz=pres
  for i=0,n1-1 do begin
         dvrdz(i,*)= deriv(x2,vr(i,*))
  endfor

 dphidz=pres
  for i=0,n1-1 do begin
         dphidz(i,*)= deriv(x2,phi(i,*))
  endfor

 dprdz=pres
  for i=0,n1-1 do begin
         dprdz(i,*)= deriv(x2,pres(i,*))
  endfor

 dBrdz=br
  for i=0,n1-1 do begin
         dBrdz(i,*)= deriv(x2,br(i,*))
  endfor

 dBphidz=br
  for i=0,n1-1 do begin
         dBphidz(i,*)= deriv(x2,bphi(i,*))
  endfor

   alf, b1,b2,b3,v1,v2,rho,pr,mf,ma,ms,mson

jr  = -dBphidz
jz  =  drBphidr/rr
jphi=  dBrdz - dBzdr

fr 	= 	jphi*bz	- 	jz*bphi
fz 	= 	jr*bphi 	- 	jphi*br
fphi 	= -jr*bz		+	jz*br


msmag=fphi/den/omegamatrix^2

vsound2=pres/den
vsoundmid=vsound2
  for j=0,n1-1 do begin
     vsoundmid(*,j)=pres(*,0)/den(*,0)
     endfor
     eta=sqrt(rr^3)*(vsoundmid + (2.0/5.0)*(1/rad-1/rr))
eta(*,*) =eta(*,*)*(eta(*,*) gt 0.0)
kvisc=2.0/3.0*den*eta





kvisc=2.0/3.0*den*eta
tau_rphi=kvisc*(dvphidr - vphi/rr)
tau_zphi=kvisc*(dvphidz)





 drtau_rphidr=pres
 rtau_rphi=rr*tau_rphi
  for j=0,n2-1 do begin
         drtau_rphidr(*,j)= deriv(x1,rtau_rphi(*,j))
  endfor

 dtau_zphidz=pres
  for i=0,n1-1 do begin
         dtau_zphidz(i,*)= deriv(x2,tau_zphi(i,*))
  endfor


  divTphi= 1/rr*   drtau_rphidr +  dtau_zphidz + tau_rphi/rr


	
	trphi= 1/rr*   drtau_rphidr + tau_rphi/rr
	tzphi=  dtau_zphidz
end
