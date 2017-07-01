
common prim,  br,bz,bphi, vr,vz,vphi, den, pres
common geom, rr,zz, rad, omega, omegamatrix
common current, jr,jz,jphi, fr,fz,fphi
common visc, trphi, tzphi, divTphi


simtime=199
simtime=nlast

pload,out=simtime

compute
;;;;;;;;;;;;;;;

rm=-vr*rr/kvisc

brbz=br/bzmid

;;; here we specify the interpolation surface.
;;; using variable h2


hpeak=0.1*x1
hmax=0.1*x1
cellheight=indgen(n1)
rmpeak=fltarr(n1)
jphipeak=fltarr(n1)
bzpeak=fltarr(n1)
vzpeak=fltarr(n1)
vrpeak=fltarr(n1)
brpeak=fltarr(n1)
etapeak=fltarr(n1)

;h2=x2
for i=0,n1-1 do begin
      a=where(eta(i,*) eq 0)
 ;     print,a(0)
    ;  a=a-1
      a=a * (a gt 0)
      hpeak(i)=x2(a(0)-3)
      hmax(i)=x2(a(0))
      cellheight[i]=a(0)-3
      qq=cellheight[i]
;      qq=1
     ; rmpeak[i]=rm[i, qq ]
      jphipeak[i]=jphi[i, qq ]
      bzpeak[i]=bz[i, qq ]
      vzpeak[i]=vz[i, qq ]
      vrpeak[i]=vr[i, qq ]
      brpeak[i]=br[i, qq ]
      etapeak[i]=eta[i, qq ]
   endfor


loadct,0
!p.color=1
!p.background=255

!p.multi=[0,1,2]
!p.multi=0


jphinorm=bz(*,0)/x1
vrbznorm=-vrpeak*bzpeak/etapeak/jphinorm 


!x.range=[0,20]
!y.range=[0,0.1]
;!y.range=0
plot, x1,-vrpeak*bzpeak/jphinorm , linestyle=0
oplot, x1,etapeak*jphipeak/jphinorm, linestyle=1
oplot, x1,vzpeak*brpeak/jphinorm, linestyle=2
items=["vrBz","etajphi","vzbr"]
legend,items,linestyle=indgen(3), /right

test=vzpeak*brpeak/bz[*,0]/x1


window, 3

etajphi=eta*jphi ;*2*!PI
urbz=vr*bz
uzbr=vz*br
urbzdisk=fltarr(n1,n2)
uzbrdisk=fltarr(n1,n2)


index=where(etajphi ne 0, count)
IF count NE 0 THEN urbzdisk[index] = urbz[index]
IF count NE 0 THEN uzbrdisk[index] = uzbr[index]

normetajphi=max(etajphi)

!x.range=[0,20]
!y.range=[1e-7,0.0004]
!y.range=[1e-4,2]
window, title="etajphi"
plot, x1, total(etajphi,2)/normetajphi, /ylog, linestyle=0,$
title='Time='+string(time(simtime)/2/!PI),$
xtitle="Radius, (R)"
oplot, x1, -total(urbzdisk,2)/normetajphi, linestyle=2
oplot, x1, -total(uzbrdisk,2)/normetajphi, linestyle=1

items=['etajphi','urbzdisk', 'uzbrdisk']
legend,items,linestyle=[0,2,1],/right


urbzdiskdz= urbzdisk*dzmatrix
etajphidz= etajphi*dzmatrix

numerator=bflux/2/!PI/x1
denominator= (- total(urbzdiskdz,2) - total(etajphidz,2) ) /hmax

window, 3

!x.range=[1,100]
!y.range=[1,1e6]
plot, x1, (abs(numerator/denominator)) ,  /xlog, /ylog

end
