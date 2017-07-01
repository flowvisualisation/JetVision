!P.POSITION=0
!P.MULTI=[0,1,3]

loadct,0
!P.CHARSIZE=2
window, 13, xs=800, ys=900

br=b1(0)
bz=b2(0)
bphi=b3(0)

pres=pr(0)

mu=alog10((br^2 +bz^2 +bphi^2 )/pres)
;mu=alog10((br^2 +bz^2+bphi^2 )/pres)

nx=20
plot, x2, mu(nx,*) ,$
	title='!7l!3, Magnetisation, r=' +string(x1(nx)) ,$
;yrange=[0,2],$
xrange=[0,2]
oplot,x2,mu(30,*)
oplot,x2,mu(40,*)


xmin=1.2


pres=pr(0)
dens=rho(0)

r=x1
h1=sqrt(pres(*,0)/dens(*,0)) *r^1.5
h2=h1

vr=v1(0)

for i=0,n2-1 do begin

;      a=where(vr(i,*)>0)
;      print,a(0)
;      h2(i)=x2(a(0))
   endfor

den=rho(0)
h2=sqrt(pres(*,0)/den(*,0))*x1^(3./2.)*2

plot,x1, h2, $
	xrange=[xmin,20],xstyle=1

zz=h2
zz=2.4*sqrt(pres(*,0)/dens(*,0)) *r^1.5

r=x1
z=x2
dr=dx1
dz=dx2
rr=x1
interplot,  mu ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl

plot, rpl, qinter1,$
	xrange=[xmin,20],xstyle=1




im=tvread(filename="mu", /png, /nodialog)
end
