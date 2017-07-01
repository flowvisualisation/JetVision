
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

   alf, b1,b2,b3,v1,v2,rho,pr,mf,ma,ms,mson

jz  =  drBphidr/rr
jr  = -dBphidz
jphi=  dBrdz - dBzdr




window,8,xs=600,ys=800
!P.POSITION=0
!P.CHARSIZE=2.
!P.MULTI=0
!P.MULTI=[0,2,3]


nx=30
print,x1(nx) 

xbeg=0
xend=5

ybeg=-1e-2
yend=1e-2
plot, x2,br(nx,*),color=100,$
		linestyle=0,$
		xtitle='z',$
		xrange=[xbeg,xend],$
		yrange=[ybeg,yend],$
		title="B fields, r="+string(x1(nx),format='(F4.1)')
oplot, x1,bz(nx,*),$
		linestyle=1
oplot, x1,bphi(nx,*),$
		linestyle=2
items = ['Br','Bz', 'Bphi']
lines=indgen(3)
legend,  items,linestyle=lines,/right




plot, x2,jr(nx,*),$
		xtitle='z ',$
		xrange=[xbeg,xend],$
		linestyle=0,$
		title="Jr, Jphi, r="+string(x1(nx),format='(F4.1)')
oplot, x2,jr(nx,*),psym=4
oplot, x1,jphi(nx,*),$
		thick=2,$
		linestyle=1
items = ['Jr','Jphi']
lines=indgen(2)
legend,  items,linestyle=lines, /right

nx=70
print,x1(nx) 
ybeg=-1e-2
yend=1e-2
plot, x2,br(nx,*),color=100,$
		linestyle=0,$
		xtitle='z',$
		xrange=[xbeg,xend],$
		yrange=[ybeg,yend],$
		title="B fields, r="+string(x1(nx),format='(F4.1)')
oplot, x1,bz(nx,*),$
		linestyle=1
oplot, x1,bphi(nx,*),$
		linestyle=2
items = ['Br','Bz', 'Bphi']
lines=indgen(3)
legend,  items,linestyle=lines,/right



plot, x2,jr(nx,*),$
		xtitle='z ',$
		xrange=[xbeg,xend],$
		linestyle=0,$
		title="Jr, Jphi, r="+string(x1(nx),format='(F4.1)')
oplot, x2,jr(nx,*),psym=4
oplot, x1,jphi(nx,*),$
		thick=2,$
		linestyle=1
items = ['Jr','Jphi']
lines=indgen(2)
legend,  items,linestyle=lines, /right



nx=180
print,x1(nx) 
ybeg=-1e-2
yend=1e-2
plot, x2,br(nx,*),color=100,$
		linestyle=0,$
		xtitle='z',$
		xrange=[xbeg,xend],$
		yrange=[ybeg,yend],$
		title="B fields, r="+string(x1(nx),format='(F4.1)')
oplot, x1,bz(nx,*),$
		linestyle=1
oplot, x1,bphi(nx,*),$
		linestyle=2
items = ['Br','Bz', 'Bphi']
lines=indgen(3)
legend,  items,linestyle=lines,/right


plot, x2,jr(nx,*),$
		xtitle='z ',$
		xrange=[xbeg,xend],$
		linestyle=0,$
		title="Jr, Jphi, r="+string(x1(nx),format='(F4.1)')
oplot, x2,jr(nx,*),psym=4
items = ['Jr','Jphi']
lines=indgen(2)
legend,  items,linestyle=lines, /right
oplot, x1,jphi(nx,*),$
		thick=2,$
		linestyle=1

end
