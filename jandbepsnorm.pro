
pload,out=nlast
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




;window,8,xs=600,ys=800
!P.POSITION=0
csize=1.3
!P.CHARSIZE=csize
csize=0.6
!P.MULTI=0
!P.MULTI=[0,2,3]
!p.font=0
set_plot,'ps'
device, filename='jandb123.eps',/encapsulated
device, /times

xs=5.
ys=6*xs/5.
DEVICE, XSIZE=xs, YSIZE=ys, /INCHES

nx=30
print,x1(nx) 

xbeg=0
xend=5

ybeg=min(bphi(nx,*))
yend=max(bz(nx,*))
plot, x2/0.1/x1(nx),br(nx,*),$
		linestyle=0,$
		xtitle='z/h',$
		xrange=[xbeg,xend],$
		yrange=[ybeg,yend],$
		title="B fields, r="+string(x1(nx),format='(F4.1)')
oplot, x2/0.1/x1(nx),bz(nx,*),$
		linestyle=1
oplot, x2/0.1/x1(nx),bphi(nx,*),$
		linestyle=2
items = ['B!Dr!N','B!Dz!N', 'B!D!9f!X!N']
lines=indgen(3)
legend,  items,linestyle=lines,/right ,charsize=csize




plot, x2/0.1/x1(nx),jphi(nx,*),$
		xtitle='z/h',$
		xrange=[xbeg,xend],$
		linestyle=1,$
		title="J!Dr!N, J!D!9f!X!N, r="+string(x1(nx),format='(F4.1)')
;oplot, x2/0.1/x1(nx),jr(nx,*),psym=4
oplot, x2/0.1/x1(nx),jr(nx,*),$
		linestyle=0
items = ['Jr','J!9f!X']
items = ['J!Dr!N','J!D!9f!X!N']
lines=indgen(2)
legend,  items,linestyle=lines, /right ,charsize=csize

nx=100
print,x1(nx) 
ybeg=min(bphi(nx,*))
yend=max(br(nx,*))
plot, x2/0.1/x1(nx),br(nx,*),$
		linestyle=0,$
		xtitle='z/h',$
		xrange=[xbeg,xend],$
		yrange=[ybeg,yend],$
		title="B fields, r="+string(x1(nx),format='(F4.1)')
oplot, x2/0.1/x1(nx),bz(nx,*),$
		linestyle=1
oplot, x2/0.1/x1(nx),bphi(nx,*),$
		linestyle=2
items = ['Br','Bz', 'B!9f!X']
items = ['B!Dr!N','B!Dz!N', 'B!D!9f!X!N']
lines=indgen(3)
legend,  items,linestyle=lines,/right ,charsize=csize



plot, x2/0.1/x1(nx),jr(nx,*),$
		xtitle='z/h',$
		xrange=[xbeg,xend],$
		linestyle=0,$
		title="J!Dr!N, J!D!9f!X!N, r="+string(x1(nx),format='(F4.1)')
;oplot, x2/0.1/x1(nx),jr(nx,*),psym=4
oplot, x2/0.1/x1(nx),jphi(nx,*),$
		linestyle=1
items = ['Jr','J!9f!X']
items = ['J!Dr!N','J!D!9f!X!N']
lines=indgen(2)
legend,  items,linestyle=lines, /right ,charsize=csize



nx=400
print,x1(nx) 
ybeg=min(bphi(nx,*))
yend=max(br(nx,*))
plot, x2/0.1/x1(nx),br(nx,*),$
		linestyle=0,$
		xtitle='z/h',$
		xrange=[xbeg,xend],$
		yrange=[ybeg,yend],$
		title="B fields, r="+string(x1(nx),format='(F4.1)')
oplot, x2/0.1/x1(nx),bz(nx,*),$
		linestyle=1
oplot, x2/0.1/x1(nx),bphi(nx,*),$
		linestyle=2
items = ['Br','Bz', 'B!9f!X']
items = ['B!Dr!N','B!Dz!N', 'B!D!9f!X!N']
lines=indgen(3)
legend,  items,linestyle=lines,/bottom ,charsize=csize


plot, x2/0.1/x1(nx),jphi(nx,*),$
		xtitle='z/h',$
		xrange=[xbeg,xend],$
		linestyle=1,$
		title="J!Dr!N, J!D!9f!X!N, r="+string(x1(nx),format='(F4.1)')
;oplot, x2/0.1/x1(nx),jr(nx,*),psym=4
items = ['J!Dr!N','J!D!9f!X!N']
lines=indgen(2)
legend,  items,linestyle=lines, /left ,charsize=csize
oplot, x2/0.1/x1(nx),jr(nx,*),$
		linestyle=0

device, /close
set_plot,'x'
end
