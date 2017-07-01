
window,13,title="Fastness", xs=900, ys=600
loadct,0

!P.POSITION=0
!P.MULTI=[0,3,2]
!P.MULTI=0
!P.CHARSIZE=2.0
!P.CHARSIZE=2.0
!P.CHARTHICK=2.0
!P.THICK=3.0


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


vpol=sqrt(vr^2+vz^2)
fastness=omega*rr/vpol

r=x1
z=x2
dr=dx1
dz=dx2


posalf, b1,b2,b3,v1,v2,rho,pr,mf,ma,ms,mson



; loop through grid finding each alf point
;

npts=11
karr=fltarr(npts)
frad=fltarr(npts)
ibeg=2
iend=9
for i=ibeg,iend do begin


field_line, br, bz, x1,x2, 1.5+ i*0.1,0.03, rr,zz


fluff=ma
interplot,fluff ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
a=where(qinter1 gt 1)
if (a(0) < 0) then begin
	a(0)=0
	endif

print, a(0)
fluff=fastness
interplot,fluff ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
qk=qinter1
print, qk(a(0))
karr(i)=qk(a(0))
frad(i)=rpl(a(0))


	endfor

;plot, karr, lambdaarr;,$


plot, frad, karr,  psym=4,$
		xtitle="Radius, R",$
		ytitle="!7x!3",$
		title="Fastness parameter, !7x!3"


im=tvread(filename="fastness", /png, /nodialog)

set_plot,'ps'
device,filename='fastness.eps', /encapsulated
plot, frad, karr,  psym=4,$
		xtitle="Radius, R",$
		ytitle="!7x!3",$
		title="Fastness parameter, !7x!3"
device,/close
set_plot,'x'
!P.COLOR=0
!P.BACKGROUND=255


end
