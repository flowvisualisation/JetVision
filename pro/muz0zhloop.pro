
set_plot,'x'
bign=[0,299,599]
pload,out=0

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

qint=rr
mu0=rr

for iterate=0,2 do begin
nend=bign[iterate]
pload,out=nend

r=x1
z=x2
dr=dx1
dz=dx2
rr=x1
zz=0.1*x1
rrr=rebin(reform(r,n1,1),n1,n2 )
zz=rebin(reform(x2,n2,1),n2,n1 )
zzzz=transpose(zz)
rad=sqrt(rrr^2 + zzzz^2)
mu=(b1(0)^2+b2(0)^2+b3(0)^2)/pr(0)
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
      ;print,a(0)
      a=a+1
      h2(i)=x2(a(0))
   endfor
zz=smooth(h2,5)
interplot,mu,r,z,dr,dz,rr,zz,qinter2,rpl,zpl
help, qint, qinter2
qint=[[qint], [qinter2] ]
arg=reform(mu[*,0])
mu0=[[mu0], [arg] ]
loadct,33
display,alog10(rho(0)),x1=x1,x2=x2, xrange=[1,10], yrange=[0,5], ims=4
oplot,rr,zz
endfor

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
!x.style=1
!y.range=[min(mu0[0:100,1:3] ), max(mu0[*,1:3])]
!y.range=[0,0.04]
plot, x1, mu0[*,1], $
/xlog,$
xtitle="Radius, R",$
ytitle="!9m!X!D0!N"
oplot, x1, mu0[*,2], linestyle=1
oplot, x1, mu0[*,3], linestyle=2

!y.range=[0,5]
plot, rpl, smooth(qint[*,1],10,/edge_truncate)/2, $
/xlog,$
xtitle="Radius, R",$
ytitle="!9m!X!U+!N"
oplot, rpl, smooth(qint[*,2],10,/edge_truncate)/2, linestyle=1
oplot, rpl, smooth(qint[*,3],10,/edge_truncate)/2, linestyle=2

items=['Time='+string(time(0)), 'Time='+string(time(299)/2/!PI), 'Time='+string(time(nlast)/2/!PI)]
legend, items,linestyle=indgen(3)
device, /close



!p.multi=0
!p.position=0
!x.range=0
!y.range=0
set_plot, 'x'
end

