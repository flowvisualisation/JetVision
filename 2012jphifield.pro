
loadct,0
simtime=nlast
pload,out=simtime

!P.POSITION=0
!P.MULTI=0
!P.CHARSIZE=2.0
!P.POSITION=0
!P.CHARSIZE=2.0
!P.CHARTHICK=2.0
!P.THICK=3.0
!P.COLOR=0
!P.BACKGROUND=255


; K

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

rr=rebin(reform(x1,n1,1),n1,n2 )
zz=rebin(reform(x2,n2,1),n2,n1 )
help,rr
zz=transpose(zz)
help,zz


rad=sqrt(rr^2+zz^2)
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

;	alf, b1,b2,b3,v1,v2,rho,pr,mf,ma,ms,mson

;jz  =  drBphidr/rr
;jr  = -dBphidz
jphi=  dBrdz - dBzdr
i=5
field_line, br, bz, x1,x2,  i,0.0005, rf1,zf1
i=13
field_line, br, bz, x1,x2,  i,0.0005, rf2,zf2
i=15
field_line, br, bz, x1,x2,  i,0.0005, rf3,zf3

for usingps=0,1 do begin
fname='2012jphifield'
phistring='!7u!X'
if ( usingps ) then begin
set_plot,'ps'
device,filename=fname+'.eps',/encapsulated, /color
phistring='!9f!X'
!p.font=0
device, /times
xs=10.
ys=6
DEVICE, XSIZE=xs, YSIZE=ys, /INCHES
endif else begin
set_plot,'x'
!p.font=-1
!p.color=0
!p.background=255
window, title=fname
endelse

!x.range=[1,20]
!y.range=[0,20]
!x.style=1
!P.THICK=1.0
loadct,33
tvlct,0,0,0,1
tvlct,255,255,255,0
!p.background=0
!p.color=1
maxjphi=max(jphi)
minjphi= min(jphi)
num=20
lev=fltarr(num)
for qi=0,num-1 do  begin
lev[qi] = minjphi+qi*(maxjphi-minjphi)/num
endfor 
levels=0.0001*[-3,-2,-1,1,2,3]
labels=['J','J','J','J','J','J']
linestyles=[0]
contour,  jphi,rr,zz, levels=levels, color=1, $
         c_linestyle=linestyles,$
        C_LABELS = levels, $
        c_colors=[40,40,40,200,200,200],$
         xtitle="Radius, (R)",$
         ytitle="Altitude (Z)",$
         title='J!D'+greek('phi')+'!N, field in zones 1-3'+', time='+$
string(time(simtime)/2/!PI, format='(I3)')+greek('tau')+"!DK!N"


oplot, rf1,zf1, color=90, linestyle=4
oplot, rf2,zf2, color=133, linestyle=4
oplot, rf3,zf3, color=233, linestyle=4


;alf, b1,b2,b3,v1,v2,rho,pr,mf,ma,ms,mson
;
posalf, b1,b2,b3,v1,v2,rho,pr,mf,ma,ms,mson
contour, mf, x1,x2, levels=1, color=233, /overplot, c_thick=2.0, c_linestyle=2
contour, ma, x1,x2, levels=1, color=133, /overplot, c_thick=2.0, c_linestyle=2

;im=tvread(filename="jphiBzones123",/png,/nodialog)

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

