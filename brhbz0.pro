
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
vr=v1(0)
vz=v2(0)
vphi=v3(0)
br=b1(0)
bz=b2(0)
bphi=b3(0)

vsound2=pres/den
vsoundmid=vsound2
bzmid=vsound2
  for j=0,n2-1 do begin
     vsoundmid(*,j)=pres(*,0)/den(*,0)
     bzmid(*,j)=bz(*,0)
     endfor
     eta=sqrt(rrr^3)*(vsoundmid + (2.0/5.0)*(1/rad-1/rrr))
eta(*,*) =eta(*,*)*(eta(*,*) gt 0.0)
kvisc=2.0/3.0*eta

h2=vsound2/x1^1.5

zz=smooth(h2,5)

rey=vr*rrr/kvisc
irey=kvisc/vr/rrr
rey=rey(*,*)*(kvisc(*,*) gt 0.0)
rey_int=total(rey,2)

brhbz0=br/bzmid

interplot,brhbz0 ,r,z,dr,dz,rr,zz,brhbz_inter,rplbrhbz,zplbrhbz
;interplot,rey ,r,z,dr,dz,rr,zz,rey_inter,rplrey,zplrey

loadct,33
display,alog10(rho(0)),x1=x1,x2=x2, xrange=[1,10], yrange=[0,5], ims=4
oplot,rr,zz


loadct,0
!p.color=0
!p.background=255

!p.font=0
;set_plot,'ps'
;device, filename="muz0zh.eps"
;device, /encapsulated
;device, /times
loadct,0
!P.POSITION=0
!P.CHARSIZE=1
!P.MULTI=[0,1,2]
!x.range=[1.2,10]
!y.range=0
plot, x1, rey_int, $
xtitle="Radius, R",$
ytitle="!9m!X!D0!N"

plot, rpl, smooth(brhbz_inter,10,/edge_truncate)/2, $
/xlog,$
xtitle="Radius, R",$
ytitle="!9m!X!U+!N"
;device, /close



!p.multi=0
!p.position=0
!x.range=0
!y.range=0
set_plot, 'x'
end

