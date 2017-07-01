!P.COLOR=0
!P.BACKGROUND=255

loadct,0
den=rho(0)
pres=pr(0)

vr=v1(0)
vz=v2(0)
vphi=v3(0)

br=b1(0)
bz=b2(0)
bphi=b3(0)


rr=rebin(reform(x1,n1,1),n1,n2 )
rr=rebin(reform(x2,n2,1),n1,n2 )

rr=br
zz=br
  for j=0,n2-1 do begin
         rr(*,j)= x1
  endfor
  for i=0,n1-1 do begin
         zz(i,*)= x2
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

   alf, b1,b2,b3,v1,v2,rho,pr,mf,ma,ms,mson

jz  =  drBphidr/rr
jr  = -dBphidz
jphi=  dBrdz - dBzdr


  
  vsound2=pres/den
  vsoundmid=vsound2
    for j=0,n1-1 do begin
         vsoundmid(*,j)=pres(*,0)/den(*,0)
	      endfor
	           eta=sqrt(rr^3)*(vsoundmid + (2.0/5.0)*(1/rad-1/rr))
		   eta(*,*) =eta(*,*)*(eta(*,*) gt 0.0)
		   kvisc=2.0/3.0*den*eta



window,8,xs=600,ys=800
!P.POSITION=0
!P.CHARSIZE=1.
!P.MULTI=0
!P.MULTI=[0,1,2]


nx=200
print,x1(nx) 

xbeg=0
xend=1

ybeg=-0.005
yend=0.006
plot, x2,br(nx,*),color=100,$
		linestyle=0,$
		xtitle='z',$
		xrange=[xbeg,xend],$
		yrange=[ybeg,yend],$
		ystyle=1,$
		title="B fields"
oplot, x1,bz(nx,*),$
		linestyle=1
oplot, x1,bphi(nx,*),$
		linestyle=2

items = ['Br','Bz', 'Bphi']
lines=indgen(3)
legend,  items,linestyle=lines,/right

ybeg=-0.02
yend=0.05
plot, x2,jphi(nx,*),$
		xtitle='z ',$
		xrange=[xbeg,xend],$
		yrange=[ybeg,yend],$
		linestyle=0,$
		title="Jphi, r="+string(x1(nx))
oplot, x2,dBrdz(nx,*),linestyle=1
oplot, x2,-dBzdr(nx,*),linestyle=2
items = ['Jphi','dBrdz','-dBzdr']
lines=indgen(3)
legend,  items,linestyle=lines, /right
;im=tvread(filename="jphicomponents",/png,/nodialog)


window, 10

plot, x1, -x1*v1(*,0,0)/eta, $
	xrange=[1.2,10],$
	xtitle="Radius",$
	/xlog,$
	xstyle=1,$
	title="Magnetic Reynolds number"
plot, x1, (b2(*,0,0)),$
	/xlog, $
	xrange=[1.2,10],$
	xtitle="Radius",$
	xstyle=1,$
	title="Bz"
im=tvread(filename="RMandBz",/png,/nodialog)
end


