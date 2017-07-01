
window,12,title="Invariants", xs=900, ys=600
loadct,0

!P.POSITION=0
!P.MULTI=[0,3,2]
!P.MULTI=0
!P.CHARSIZE=2.0
!P.CHARSIZE=1.0
!P.CHARTHICK=1.0
!P.THICK=1.0


; K

den=rho(0)
pres=pr(0)



vr=v1(0)
vz=v2(0)
vphi=v3(0)

br=b1(0)
bz=b2(0)
bphi=b3(0)

k= den*sqrt(vr^2+vz^2)/sqrt(br^2+bz^2)


rr=rebin(reform(x1,n1,1),n1,n2 )
zz=rebin(reform(x2,n2,1),n2,n1 )
zz=transpose(zz)

rad_sphere=sqrt(rr^2+zz^2)

lambda=vphi*rr-rr*bphi/k

omega= vphi/rr - k*bphi /den/rr


s=alog(pres/den^(5./3.))



ke= 0.5*(vr^2  +vz^2 + vphi^2)
enth=(5./2.)*pres/den
mag= -(omega)*rr*bphi/k
grav= -1/rad_sphere


e= ke + enth+ mag + grav

r=x1
z=x2
dr=dx1
dz=dx2


posalf, b1,b2,b3,v1,v2,rho,pr,mf,ma,ms,mson



; loop through grid finding each alf point
;

npts=12
karr=fltarr(npts)
lambdaarr=fltarr(npts)
flarr=fltarr(npts)

ibeg=1
iend=npts-2


for i=ibeg,iend do begin

flarr(i)=1.1+0.4*i
print, 'i', i, flarr(i)

field_line, br, bz, x1,x2,   flarr(i),0.03, rr,zz


fluff=ma
interplot,fluff ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
a=where(qinter1 gt 1)
if (a(0) < 0) then begin
	a(0)=0
	endif

print, a(0)

fluff=k
interplot,fluff ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
qk=qinter1
print, qk(a(0))

fluff=lambda
interplot,fluff ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
qlambda=qinter1
print, qlambda(a(0))

fluff=bz
interplot,fluff ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl

karr(i)=qk(a(0))/fluff(0)*sqrt(rpl(0))
lambdaarr(i)=qlambda(a(0))/rpl(0)^(0.5)

print, 'radius=', flarr(i), ', k=', karr(i), ', lambda=', lambdaarr(i)

	endfor

;plot, karr, lambdaarr;,$


plot, karr, lambdaarr, psym=4,$
		xtitle="!7j!3",$
		ytitle="!7k!3",$
		title="!7j!3-!7k!3"

 oplot, karr, 1.5*(1.0+karr^(-2./3.))



im=tvread(filename="klambda", /png, /nodialog)

!P.font=0
set_plot,'ps'
device,filename='klambda.eps', /encapsulated
device, /times
plot, karr, lambdaarr, psym=4,$
;      title="!9k!3-!9l!3"
		xtitle="!9k!3",$
      ytitle="!9l!3"

 oplot, karr, 1.5*(1.0+karr^(-2./3.))

device,/close
set_plot,'x'
!P.COLOR=0
!P.BACKGROUND=255


end
