
nend=nlast
pload,out=nend

mu=(b1(0)^2+b2(0)^2+b3(0)^2)/pr(0)




r=x1
z=x2
dr=dx1
dz=dx2
rr=x1
zz=0.1*x1

rrr=rebin(reform(r,n1,1),n1,n1 )
zz=rebin(reform(x2,n2,1),n2,n1 )
zzzz=transpose(zz)

rad=sqrt(rrr^2 + zzzz^2)

sound2=(pr(0)/rho(0))
sound=sqrt(pr(0)/rho(0))



den=rho(0)
pres=pr(0)

vsound2=pres/den
vsoundmid=vsound2
  for j=0,n1-1 do begin
     vsoundmid(*,j)=pres(*,0)/den(*,0)
     endfor
     eta=sqrt(rrr^3)*(vsoundmid + (2.0/5.0)*(1/rad-1/rrr))
eta(*,*) =eta(*,*)*(eta(*,*) gt 0.0)
kvisc=2.0/3.0*den*eta

h2=x2
for i=0,n1-1 do begin
      a=where(eta(i,*) eq 0)
      print,a(0)
;      a=a-1
      a=a+1
;      a=a * (a gt 0)
      h2(i)=x2(a(0))
   endfor


zz=smooth(h2,5)


interplot,mu ,r,z,dr,dz,rr,zz,qinter2,rpl,zpl

loadct,33
display,alog10(rho(0)),x1=x1,x2=x2, xrange=[1,10], yrange=[0,5], ims=4
oplot,rr,zz


loadct,0
!p.color=0
!p.background=255

!p.font=0
set_plot,'ps'
device, filename="muz0zh.eps"
device, /encapsulated
device, /times
loadct,0
!P.POSITION=0
!P.CHARSIZE=1
!P.MULTI=[0,1,2]
!x.range=[1.2,10]
!y.range=0
plot, x1, mu[*,0], $
/xlog,$
xtitle="Radius, R",$
ytitle="!9m!X!D0!N"

plot, rpl, smooth(qinter2,10,/edge_truncate)/2, $
/xlog,$
xtitle="Radius, R",$
ytitle="!9m!X!U+!N"
device, /close



!p.multi=0
!p.position=0
!x.range=0
!y.range=0
set_plot, 'x'
end

