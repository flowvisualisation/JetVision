

for qq=1,nlast,99 do begin
n=qq
pload,out=n
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




z=x1
h2=x1
surf_index=intarr(n1)

jphiint=x1
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



for i=0,n1-1 do begin
zeta_index=max([surf_index(i),0])


jphiint(i)= total(  jphi(i,0:zeta_index) *dx2[0:zeta_index])

endfor


!p.color=0
!p.background=255
plot , x1, smooth(jphiint,5, /edge_truncate) , xrange=[1,10], xstyle=1, title='int Jphi dz' 
im=tvread(filename='cur'+string(qq, format='(I04)'), /png, /nodialog)

endfor

end
