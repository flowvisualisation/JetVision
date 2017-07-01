

pload,out=nlast
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
dz=dr


rr=rebin(reform(x1,n1,1),n1,n2 )
zz=rebin(reform(x2,n2,1),n2,n1 )
zz=transpose(zz)


rad=sqrt(rr^2 +zz^2)

vsound2=pres/den
vsoundmid=vsound2
bzmid=vsound2
  for j=0,n2-1 do begin
     vsoundmid(*,j)=pres(*,0)/den(*,0)
     bzmid(*,j)=bz(*,0)
     endfor
     eta=sqrt(rr^3)*(vsoundmid + (2.0/5.0)*(1/rad-1/rr))
eta(*,*) =eta(*,*)*(eta(*,*) gt 0.0)
kvisc=eta



rm=-vr*rr/kvisc

brbz=br/bzmid

h2=0.1*x1
h2=1.5*sqrt(vsound2(*,0))*x1^(1.5)
interplot,brbz ,r,z,dr,dz,r,h2,brbz_inter,rplbrbz,zplbrbz
brbzuni=interpol(brbz_inter, rplbrbz,r)




z=x1
surf_index=intarr(n1)

rm_int=x1
rm_ave=x1
for i=0,n1-1 do begin

a=where((vr(i,*) gt 0) and (vz(i,*) gt 0))
z(i)=x2(a(0)+1)

surf_index(i)=a(0)-1

endfor

for i=0,n1-1 do begin
surf_index(i)=max(where(x2 lt h2(i)) )
;print, surf_index(i)
endfor



for i=0,n1-1 do begin
zeta_index=max([surf_index(i),1])
print, zeta_index
rmh2=interpol(rm(i,*), x2,h2(i))
abscissa=[x2(0:zeta_index), h2(i)]
ord=[reform(rm(i,0:zeta_index)),rmh2]
;rm_int(i)= total(  rm(i,0:zeta_index) *dx2[0:zeta_index])
rm_int(i)=int_tabulated( abscissa,ord )
rm_ave(i)= rm_int(i)/h2(i)
endfor


!p.color=0
!p.background=255
!x.range=[1.2,10]
!x.style=1
!p.multi=[0,1,3]
!p.charsize=2
plot , r,brbzuni*r/h2, title="Br!U+!N/Bz!D0!N * r/h"
plot , r,rm_ave, title="Rm averaged"
plot , r,-rm_ave+brbzuni*r/h2, title="sum"
oplot,r, deriv(alog(r), alog(bz(*,0))), linestyle=1
items=['sum', 'dlog(Bz)/dlog(r)']
legend, items, linestyle=indgen(2), charsize=1
im=tvread(filename='cur'+string(nlast, format='(I04)'), /png, /nodialog)

!p.font=0
set_plot, 'ps'
device, filename="brbz_rm_integrated_.eps"
device, /encapsulated
device, /times
!p.color=0
!p.background=255
!x.range=[1.2,10]
!x.style=1
!p.multi=[0,1,3]
!p.charsize=2
plot , r,brbzuni*r/h2, title="Br!U+!N/Bz!D0!N * r/h", xtitle="Radius, R"
plot , r,rm_ave, title="Rm averaged", xtitle="Radius, R"
plot , r,-rm_ave+brbzuni*r/h2, title="sum", xtitle="Radius, R"
oplot,r, deriv(alog(r), alog(bz(*,0))), linestyle=1
items=['sum', 'dlog(Bz)/dlog(r)']
legend, items, linestyle=indgen(2), charsize=1

device, /close
set_plot, 'x'

end
