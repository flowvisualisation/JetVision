
den=rho(0)
pres=pr(0)

vr=v1(0)
vz=v2(0)
vphi=v3(0)

br=b1(0)
bz=b2(0)
bphi=b3(0)


rr=rebin(reform(x1,n1,1),n1,n1 )
zz=transpose(rr)


rad=sqrt(rr^2 +zz^2)


omega=vphi/rr

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


 domegadr=pres
  for j=0,n2-1 do begin
         domegadr(*,j)= deriv(x1,omega(*,j))
  endfor

 domegadz=pres
  for i=0,n1-1 do begin
         domegadz(i,*)= deriv(x2,omega(i,*))
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

jz  =  drBphidr/rr
jr  = -dBphidz
jphi=  dBrdz - dBzdr


vsound2=pres/den
vsoundmid=vsound2
  for j=0,n1-1 do begin
	       vsoundmid(*,j)=pres(*,0)/den(*,0)
			      endfor
					     eta=sqrt(rr^3)*(vsoundmid + (2.0/5.0)*(1/rad-1/rr))
						  eta(*,*) =eta(*,*)*(eta(*,*) gt 0.0)
						  kvisc=2.0/3.0*den*eta
etajr=eta*jr




integrand=br*domegadr + bz*domegadz



nx=20
integral=fltarr(n2+2)
print, n2
integral(*)=0

dz=dx1
for j=2,n2-1 do begin
	;print, j, j-1
	integral(j)= integral(j) + integrand(nx,j)*dz(j)
	endfor

window,8,xs=600,ys=800
!P.POSITION=0
!P.CHARSIZE=3.
!P.MULTI=0


print,x1(nx) 

xbeg=0
xend=5

plot, x2, -x1(nx)*integral,color=100,$
;plot, x2, etajr(nx,*),$
		linestyle=0,$
		xtitle='z',$
		xrange=[xbeg,xend],$
	;	yrange=[-1e-3,1e-3],$
		title="integ ,r="+string(x1(nx))
oplot, x2, etajr(nx,*),$
		linestyle=2




end
