
set_plot,'x'
simtime=nlast
simtime=399

pload,out=simtime

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
dz=dx2

psi=fltarr(n1,n2)
psidisk=fltarr(n1,n2)
unitydisk=fltarr(n1,n2)
unity=fltarr(n1,n2)
unity(*,*)=1
psiave=fltarr(n1)
bzmid=bz(*,0)
bzr=bzmid*r*dr
bflux=fltarr(n1)
for i=14,n1-1 do begin
bflux[i]=bflux[i-1]+bzr[i]
endfor

rr=rebin(reform(x1,n1,1),n1,n2 )
zz=rebin(reform(x2,n2,1),n2,n1 )
zz=transpose(zz)

drmatrix=rebin(reform(dr,n1,1),n1,n2 )
dzmatrix=rebin(reform(dz,n2,1),n2,n1 )
dzmatrix=transpose(dzmatrix)


rbzdrmatrix=bz*rr*drmatrix
for j=0,n2-1 do begin
for i=14,n1-1 do begin
psi[i,j]=psi[i-1,j]+rbzdrmatrix[i,j]
endfor
endfor


rad=sqrt(rr^2 +zz^2)

rbphi=rr*bphi

 drBphidr=pres
 dBphidr=pres
 phi=-1/sqrt(rr^2+zz^2)

  for j=0,n2-1 do begin
         drBphidr(*,j)= deriv(x1,rbphi(*,j))
         dBphidr(*,j)= deriv(x1,bphi(*,j))
  endfor

 dBzdr=pres
  for j=0,n2-1 do begin
         dBzdr(*,j)= deriv(x1,bz(*,j))
  endfor

 dBrdz=br
  for i=0,n1-1 do begin
         dBrdz(i,*)= deriv(x2,br(i,*))
  endfor
if (  0 ) then begin
 dBphidz=br
  for i=0,n1-1 do begin
         dBphidz(i,*)= deriv(x2,bphi(i,*))
  endfor
  endif

;  alf, b1,b2,b3,v1,v2,rho,pr,mf,ma,ms,mson

;jz  =  drBphidr/rr
;jr  = -dBphidz
jphi=  dBrdz - dBzdr

alpha=0.9

vsound2=pres/den
vsound2mid=fltarr(n1,n2)
bzmid=vsound2
  for j=0,n2-1 do begin
     vsound2mid(*,j)=pres(*,0)/den(*,0)
     bzmid(*,j)=bz(*,0)
  endfor
     eta=alpha*sqrt(rr^3)*(vsound2mid + (2.0/5.0)*(1/rad-1/rr))
eta(*,*) =eta(*,*)*(eta(*,*) gt 0.0)
kvisc=eta



rm=-vr*rr/kvisc









etajphi=eta*jphi/100.0
urbz=vr*bz
uzbr=vz*br
urbzdisk=fltarr(n1,n2)
uzbrdisk=fltarr(n1,n2)
uzdisk=fltarr(n1,n2)
urdisk=fltarr(n1,n2)
bzdisk=fltarr(n1,n2)
brdisk=fltarr(n1,n2)


index=where(vz lt 0, count)
;index=where(eta gt 0, count)
IF count NE 0 THEN urbzdisk[index] = urbz[index]
IF count NE 0 THEN uzbrdisk[index] = uzbr[index]
IF count NE 0 THEN unitydisk[index] = unity[index]
IF count NE 0 THEN bzdisk[index] = bz[index]
IF count NE 0 THEN brdisk[index] = br[index]
IF count NE 0 THEN urdisk[index] = vr[index]
IF count NE 0 THEN uzdisk[index] = vz[index]

;normetajphi=max(etajphi)
normetajphi=1.0


if (0) then begin
set_plot,'x'
!p.position=0
!p.multi=[0,1,4]
!p.color=0
!p.background=255
!p.charsize=2

window, 17,xs=600,ys=800
plot, x1, total(etajphi,2)/normetajphi, linestyle=0,$
title="etajphi integrated"+'Time='+string(time(simtime)/2/!PI),$
xtitle="Radius, (R)", /xlog
!y.range=0
plot, x1, total(urbzdisk,2)/normetajphi, linestyle=2, /xlog
!y.range=0
plot, x1, total(uzbrdisk,2)/normetajphi, linestyle=1 , /xlog
plot, x1, (total(uzbrdisk,2)-  total(urbzdisk,2))/normetajphi, linestyle=1


items=['etajphi','urbzdisk', 'uzbrdisk']
legend,items,linestyle=[0,2,1],/right

endif

urbzdiskdz= urbzdisk*dzmatrix
etajphidz= etajphi*dzmatrix
psidiskdz= psidisk*dzmatrix
unitydiskdz= unitydisk*dzmatrix

integraldz=total(unitydiskdz,2 )
psitotal=total(psidiskdz,2 )/integraldz

psi1=psi(*,0)
psishift=shift(psi1,1)

numerator=(psi1 - psishift)  / (2 * !PI * x1 )
denominator=   (- total(urbzdiskdz,2) - total(etajphidz,2) ) / integraldz
;denominator=   (- total(urbzdiskdz,2)  ) / integraldz

for usingps=0,1 do begin

fname='2012integrated'+'advectiondiffusiontimescales'

if ( usingps ) then begin
set_plot,'ps'
device,filename=fname+'.eps',/encapsulated
phistring='!9f!X'
!p.font=0
device, /times
xs=5.
ys=3
DEVICE, XSIZE=xs, YSIZE=ys, /INCHES
endif else begin
loadct,0
set_plot,'x'
!p.font=-1
!p.color=0
!p.background=255
window, title="integrated timescales"
endelse

!p.multi=0
!x.range=[1,100]
!y.range=[1e-2,1e6]
plot, x1, (abs(numerator/denominator))/2/!PI ,$
  /xlog, /ylog, title="integrated "+$
greek('psi')+"/(-u!Dr!Nb!Dz!N -"+greek('eta')+"j!D"+greek('phi')+"!N), time="+$
string(time(simtime)/2/!PI, format='(I3)')+greek('tau')+"!DK!N", $
xtitle="Radius (R)"


if ( usingps ) then begin
device,/close
set_plot,'x'
endif else begin
set_plot,'x'
im=tvread(filename=fname,/nodialog,/png)
endelse

endfor
set_plot,'x'
end
