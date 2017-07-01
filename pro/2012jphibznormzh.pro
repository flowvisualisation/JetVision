simtime=nlast
pload,out=simtime
den=rho(0)
pres=pr(0)

vr=v1(0)
vz=v2(0)
vphi=v3(0)

br=b1(0)
bz=b2(0)
bphi=b3(0)


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

 dprdr=pres
  for j=0,n2-1 do begin
         dprdr(*,j)= deriv(x1,pres(*,j))
  endfor

 dphidr=pres
  for j=0,n2-1 do begin
         dphidr(*,j)= deriv(x1,phi(*,j))
  endfor

 dphidz=pres
  for i=0,n1-1 do begin
         dphidz(i,*)= deriv(x2,phi(i,*))
  endfor

 dprdz=pres
  for i=0,n1-1 do begin
         dprdz(i,*)= deriv(x2,pres(i,*))
  endfor

 dBrdz=br
  for i=0,n1-1 do begin
         dBrdz(i,*)= deriv(x2,br(i,*))
  endfor

 dBphidz=br
  for i=0,n1-1 do begin
         dBphidz(i,*)= deriv(x2,bphi(i,*))
  endfor

 ;  alf, b1,b2,b3,v1,v2,rho,pr,mf,ma,ms,mson

jz  =  drBphidr/rr
jr  = -dBphidz
jphi=  dBrdz - dBzdr




;window,8,xs=600,ys=800
!P.POSITION=0
csize=2.0
!P.CHARSIZE=csize
!P.MULTI=0
;!P.MULTI=[0,2,3]

!p.font=-1
usingps=1

phistring='!7u!X'
if ( usingps ) then begin
set_plot,'ps'
device,filename='2012jphib.eps',/encapsulated
phistring='!9f!X'
!p.font=0
device, /times

xs=10.
ys=6
DEVICE, XSIZE=xs, YSIZE=ys, /INCHES


endif

nx=50
print,x1(nx) 

xbeg=0
xend=5

ybeg=min(bphi(nx,*))
yend=max(bz(nx,*))


!p.color=1
!p.background=255


pos1=13*5
pos2=14*10
pos3=400

!y.range=[-5,40]
!p.thick=2
plot, x2/0.1/x1(pos1),jphi(pos1,*)/bz(pos1,0),$
		xtitle='z/h',$
		xrange=[xbeg,xend],$
		linestyle=0,$
		title=' J!D'+phistring+'!N, r='+string(x1(nx),format='(F4.1)') + $
', time='+$
string(time(simtime)/2/!PI, format='(I3)')+greek('tau')+"!DK!N"
oplot, x2/0.1/x1(pos2),jphi(pos2,*)/bz(pos2,0),$
		linestyle=1
oplot, x2/0.1/x1(pos3),jphi(pos3,*)/bz(pos3,0),$
		linestyle=2
items = [ $
		'J!D'+phistring+'!N, r='+string(x1(pos1),format='(F4.1)'),$
		'J!D'+phistring+'!N, r='+string(x1(pos2),format='(F4.1)'),$
		'J!D'+phistring+'!N, r='+string(x1(pos3),format='(F4.1)')]
lines=indgen(3)
legend,  items,linestyle=lines, charsize=csize

if ( usingps ) then begin
device, /close
endif
set_plot,'x'
end
