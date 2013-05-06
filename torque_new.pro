
n=nlast
pload,out=n


pres=pr(0)
dens=rho(0)
vr=v1(0)
vz=v2(0)
vphi=v3(0)

br=b1(0)
bz=b2(0)
bphi=b3(0)

z=x1
h2=x1
surf_index=intarr(n1)

rr=rebin(reform(x1,n1,1),n1,n2 )
zz=rebin(reform(x2,n2,1),n2,n1 )
zz=transpose(zz)
rad=sqrt(rr^2+zz^2)


jphi=rr*bphi


snd=sqrt(pres/dens)

Energy=0.5*(vr^2+vz^2+vphi^2) + (5./2.)*snd^2  - 1./rad
KineticEnergy=0.5*(vr^2+vz^2+vphi^2)


vsoundmid=snd

  for j=0,n1-1 do begin
     vsoundmid(*,j)=pres(*,0)/dens(*,0)
     endfor
     eta=sqrt(rr^3)*(vsoundmid + (2.0/5.0)*(1/rad-1/rr))
eta(*,*) =eta(*,*)*(eta(*,*) gt 0.0)
kvisc=2.0/3.0*dens*eta*0.9


 domegadr=pres
  for j=0,n2-1 do begin
         domegadr(*,j)= deriv(x1,vphi(*,j)/x1)
  endfor

tau_rphi=kvisc*rr*domegadr


for i=0,n1-1 do begin

a=where((vr(i,*) gt 0) and (vz(i,*) gt 0))
z(i)=x2(a(0)+1)

surf_index(i)=a(0)-1

;h2(i)=snd(i,0)*x1(i)^(1.5)*2+0.3
endfor

h2 = smooth(z,10,/edge_truncate,/NAN)
h2 = z
for i=0,n1-1 do begin
surf_index(i)=min(where(x2 gt h2(i)) )-1
endfor

dh2=smooth(deriv(x1,h2),10,/edge_truncate)


mdot=x1
mdot2=x1
jv=x1
jphiint=x1

mout=x1
magtorqueout=x1
kintorqueout=x1
mechanicalflux=x1
kineticflux=x1

poyntingflux=x1
viscouspower=x1

jdot=x1
edot=x1

jv=x1


for i=0,n1-1 do begin
zeta_index=max([surf_index(i),0])


;mdot(i) =  4*!PI*x1(i)*int_tabulated(x2(0:zeta_index ),  dens(i,0:zeta_index )* vr(i,0:zeta_index) )
mdot2(i) =  4*!PI*x1(i)*total(  dens(i,0:zeta_index )* vr(i,0:zeta_index)*dx2[0:zeta_index] ) 
jdot(i)=  4*!PI*x1(i)*total(  dens(i,0:zeta_index )* vr(i,0:zeta_index)*dx2[0:zeta_index]*x1(i)*vphi(i,0:zeta_index) ) 
jv(i)=-4*!PI*x1(i)^2*total(  tau_rphi(i,0:zeta_index)*dx2[0:zeta_index] )
edot(i) =  4*!PI*x1(i)*total(  dens(i,0:zeta_index )* vr(i,0:zeta_index)*dx2[0:zeta_index] *Energy(i,0:zeta_index) ) 
viscouspower(i)= -4*!PI*x1(i)*total(  tau_rphi(i,0:zeta_index)*dx2[0:zeta_index] *vphi(i,0:zeta_index)  )
jphiint(i)= total(  jphi(i,0:zeta_index) *dx2[0:zeta_index])

magtorqueout(i)=-4*!PI*x1(i)^2 *interpol( bphi(i,*) *(bz(i,*)-br(i,*)*dh2(i) )  ,x2,h2(i) )
kintorqueout(i)=4*!PI * x1(i)^2 * interpol( dens(i,*) *vphi(i,*)* (vz(i,*)-vr(i,*)*dh2(i) )  ,x2,h2(i) )
poyntingflux(i)= -4*!PI*x1(i) *interpol( vphi(i,*)* bphi(i,*) *(bz(i,*)-br(i,*)*dh2(i) )  ,x2,h2(i) )

mout(i)=4*!PI * x1(i) * interpol( dens(i,*) *(vz(i,*)-vr(i,*)*dh2(i) )  ,x2,h2(i) )
mechanicalflux(i)=4*!PI * x1(i) * interpol( Energy(i,*)*dens(i,*) *(vz(i,*)-vr(i,*)*dh2(i) )  ,x2,h2(i) )
kineticflux(i)=4*!PI * x1(i) * interpol( KineticEnergy(i,*)*dens(i,*) *(vz(i,*)-vr(i,*)*dh2(i) )  ,x2,h2(i) )


endfor

rin = 1.3
rout = 5.

indin = min(where(x1 gt rin))
indout = max(where(x1 lt rout))
pacc = -edot(indout)+edot(indin)
pvisc = -viscouspower(indout)+viscouspower(indin)
pmagjet = total(poyntingflux(indin:indout)*dx1(indin:indout))
pmecjet = total(mechanicalflux(indin:indout)*dx1(indin:indout))
pkinjet = total(kineticflux(indin:indout)*dx1(indin:indout))

massoutflowrate=total(mout(indin:indout)*dx1(indin:indout))

print, 'pvisc', pvisc
print, 'pacc', pacc
print, 'pmagjet', pmagjet
print, 'pmecjet', pmecjet
print, 'pkinjet', pkinjet
print, 'massoutflowrate', massoutflowrate

print, 'pmagjet/pkinjet', pmagjet/pkinjet
print, 'lambda', pmagjet/pkinjet+1
;window,3
;plot, x1, -deriv(x1, smooth(mdot2,10,/edge_truncate)) , xrange=[1.3,5]
;oplot, x1, mout
window,4
plot, x1, -deriv(x1, smooth(jdot,10,/edge_truncate)) , xrange=[1.3,5]
oplot,x1, deriv(x1, smooth(jv,10,/edge_truncate)) , linestyle=2
oplot,x1, magtorqueout, linestyle=3
oplot,x1, kintorqueout, linestyle=4
oplot, x1, deriv(x1, smooth(jv,10,/edge_truncate)) +magtorqueout+kintorqueout, linestyle=1


window,5

biglambda=magtorqueout/deriv(x1,smooth(jv,10,/edge_truncate)) 
plot, x1, biglambda, xrange=[1.3,5]

window,6
plot,x1, -deriv(x1,smooth(edot,10,/edge_truncate)), xrange=[1.3,5]
oplot,x1, -deriv(x1,smooth(viscouspower,10,/edge_truncate)), linestyle=2
oplot,x1, poyntingflux, linestyle=3
oplot,x1, -mechanicalflux, linestyle=4
window,7
plot , x1, smooth(jphiint,5, /edge_truncate) , xrange=[1,10], xstyle=1, title='int rphi dz' 
im=tvread(filename='cur'+string(n, format='(I04)'), /png, /nodialog)
end
