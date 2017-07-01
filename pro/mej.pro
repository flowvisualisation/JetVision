

; mf



for n=0,600,100 do begin
	pload,out=n

alf, b1,b2,b3,v1,v2,rho,pr,mf,ma,ms

den=rho(0)
pres=pr(0)

vx=v1(0)
vy=v2(0)
vz=v3(0)

bx=b1(0)
by=b2(0)
bz=b3(0)

k=rho(0)*sqrt(v1(0)^2+v2(0)^2)/sqrt(b1(0)^2+b2(0)^2)

rr=rebin(reform(x1,n1,1),n1,n1 )
zz=transpose(rr)
lambda=v3(0)*rr-rr*b3(0)/k

omega= v3(0)/rr - k*b3(0) /rho(0)/rr

s=alog(pr(0)/rho(0)^(5./3.))

e=(5./2.)*pr(0) - 1/sqrt(rr^2 +zz^2) - 0.5*omega^2 *rr^2 + 0.5*(v1(0)^2 +v2(0)^2)$
   +0.5* (v3(0)/rr- omega)^2 * rr^2


 ;contour, mf, x1,x2, levels=1, color=133,  c_thick=2.0


mej=rho(*,n2-1,0)*v2(*,n2-1,0)
meja=rho(*,n2-1,0)*v2(*,n2-1,0)
mejp=rho(*,n2-1,0)*v2(*,n2-1,0)

tmej=0
tmeja=0
tmejp=0
for i=0,n1-1 do begin

	if (mf(i,n2-1) lt 1) then begin
		mej(i)=0
		endif

	tmej=tmej+ 2*!PI* mej(i)*x1(i)*dx1(i)

	if (ma(i,n2-1) lt 1) then begin
		meja(i)=0
		endif

	tmeja=tmeja+ 2*!PI* meja(i)*x1(i)*dx1(i)



	if (mf(i,n2-1) lt 1) then begin
		mejp(i)=0
		endif

	tmejp=tmejp+ 2*!PI* meja(i)*x1(i)*e(i,n2-1)*dx1(i)
		endfor

print, total(mej), total(meja)
print, tmej, tmeja, tmejp


endfor
;plot, mej
end
