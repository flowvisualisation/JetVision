
loadct,0

!P.POSITION=0
!P.MULTI=0
!P.CHARSIZE=2.0
!P.POSITION=0
!P.CHARSIZE=2.0
!P.CHARTHICK=2.0
!P.THICK=3.0
!P.COLOR=0
!P.BACKGROUND=255


; K

den=rho(0)
pres=pr(0)

vr=v1(0)
vz=v2(0)
vphi=v3(0)

br=b1(0)
bz=b2(0)
bphi=b3(0)


r=x1
z=x2
dr=dx1
dz=dx2

rr=rebin(reform(x1,n1,1),n1,n2 )
zz=rebin(reform(x2,n2,1),n2,n1 )
zz=transpose(zz)


rad=sqrt(rr^2+zz^2)
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
i=2
field_line, br, bz, x1,x2,  i,0.0005, rr,zz
fluff=jphi
interplot,fluff ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
qrad= zpl
window,22, xs=900, ys=600
!x.range=[0,10]
plot, qrad, qinter1,$
title='r='+string(i)
i=4
field_line, br, bz, x1,x2,  i,0.0005, rr,zz
fluff=jphi
interplot,fluff ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
qrad= zpl
oplot, qrad, qinter1,linestyle=1

end
