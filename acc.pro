
;pload,out=0

acc_rate=fltarr(100)



for n=nlast,nlast do begin
	pload,out=n

inner_rad=17
r_int=inner_rad
r_ext=63


sound=sqrt(pr(0)/rho(0))
sound2=(pr(0)/rho(0))


rr=x1
zz=2.3*sound(*,0)*x1^1.5 ;+0.2



vr=v1(0)
h2=x1

for i=0,n1-1 do begin

      a=where(vr(i,*)>0)
      ;print,a(0)
		a(0)=a(0)*(a(0) gt 0)
      h2(i)=x2(a(0))
   endfor


zz=h2

;window,22
;plot,x1,zz

loadct,33
display, alog(rho(0)),x1=x1,x2=x2,  ims=0.25
oplot,x1,zz,thick=2,linestyle=2,color=100

macc_int=fltarr(r_ext+1)
inthscale=fltarr(r_ext+1)

dns = rho(0)
vrd = v1(0)
vvr = v2(0)
for r_int=15,r_ext do begin

i=0
while x2(i) lt zz(r_int) do begin
	i=i+1
	nend=i
endwhile

;print, x2(nend-1), nend, zz(r_int), x2(nend), x2(nend+1)
for i=0,nend do begin
	macc_int(r_int) = macc_int(r_int) - 4*!PI*dns( r_int, i)*vrd(r_int,i)*x1(r_int)*dx2(i)
endfor
inthscale(r_int )=nend
endfor

window,2
loadct,0
plot, x1,macc_int,xrange=[x1(15), x1(r_ext)],$
	title="macc_int"



p1=dns*vrd
p2=dns*vvr

interplot, p1 ,x1,x2,dx1,dx2,rr,zz,qinter1,rpl,zpl
interplot, p2 ,x1,x2,dx1,dx2,rr,zz,qinter2,rpl,zpl


window,3
!P.POSITION=0
!P.MULTI=[0,1,2]
plot, rpl, qinter1,xrange=[1,10], title="rho vr vs r"
plot, rpl,qinter2,xrange=[1,10], title="rho vz vs r"

r_int=inner_rad

slope=deriv(x1,zz)

mej=0
for i=r_int,r_ext do begin
	;theta=atan(slope(i) )
	mej = mej + 4*!PI*( qinter2(i) - qinter1(i)*slope(i)) *x1(i)*dx1(i)
endfor

print, r_int
print, "Mass accreting at inner radius, r=",x1(r_int), macc_int(r_int)
print, "Mass accreting at outer radius, r=",x1(r_ext),macc_int(r_ext)
print, "Difference , ",macc_int(r_ext)-macc_int(r_int)
print, "Mass ejected from jet", mej
print, "Ejection to Accretion ratio", mej/(macc_int(r_ext))
print, "Average Ejection Efficiency", mej/(macc_int(r_ext))/alog(x1(r_ext) /x1(r_int))

acc_rate(n/100)= macc_int(r_int)
!P.MULTI=0

;macc_ext=0

;nend2=where( x2(0:511) gt  zz(r_ext)  )
;nend=nend2(0)


endfor

window,17
plot,acc_rate
end
