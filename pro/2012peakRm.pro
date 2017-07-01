

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

;;; here we specify the interpolation surface.
;;; using variable h2


hpeak=0.1*x1
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
      cellheight[i]=a(0)-3
      qq=cellheight[i]
      qq=1
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

etajphi=eta*jphi
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
plot, x1, total(etajphi,2)/normetajphi, /ylog, linestyle=0,$
xtitle="Radius, (R)"
oplot, x1, -total(urbzdisk,2)/normetajphi, linestyle=2
oplot, x1, -total(uzbrdisk,2)/normetajphi, linestyle=1

items=['etajphi','urbzdisk', 'uzbrdisk']
legend,items,linestyle=[0,2,1],/right

end
