
;window,22, xs=1000, ys=900
loadct,0

!P.POSITION=0
!P.MULTI=0
!P.CHARSIZE=2.0

; K

den=rho(0)
pres=pr(0)

vr=v1(0)
vz=v2(0)
vphi=v3(0)

br=b1(0)
bz=b2(0)
bphi=b3(0)
mf=0
ma=0
ms=0
mson=0
alf, b1,b2,b3,v1,v2,rho,pr,mf,ma,ms,mson
window,2
contour, mf,levels=1





rr=rebin(reform(x1,n1,1),n1,n2 )
zz=rebin(reform(x2,n2,1),n2,n1 )
zz=transpose(zz)

rad_sphere=sqrt(rr^2 +zz^2)


phi= -1/rad_sphere


 dBrdz=br
  for i=0,n1-1 do begin
         dBrdz(i,*)= deriv(x2,br(i,*))
  endfor

 dBphidz=br
  for i=0,n1-1 do begin
         dBphidz(i,*)= deriv(x2,bphi(i,*))
  endfor


rbphi=rr*bphi
 drBphidr=pres
 dBphidr=pres


  for j=0,n2-1 do begin
         drBphidr(*,j)= deriv(x1,rbphi(*,j))
         dBphidr(*,j)= deriv(x1,bphi(*,j))
  endfor

 dBzdr=pres
  for j=0,n2-1 do begin
         dBzdr(*,j)= deriv(x1,bz(*,j))
  endfor


vsound2=pres/den
vsoundmid=vsound2
  for j=0,n1-1 do begin
     vsoundmid(*,j)=pres(*,0)/den(*,0)
     endfor
     eta=sqrt(rr^3)*(vsoundmid + (2.0/5.0)*(1/rad-1/rr))
eta(*,*) =eta(*,*)*(eta(*,*) gt 0.0)
kvisc=2.0/3.0*den*eta


jr  = -dBphidz
jz  =  drBphidr/rr
jphi=  dBrdz - dBzdr


er		=  vz*bphi 	-	vphi*bz  + jr*eta
ephi 	= 	vr*bz 	-	vz*br 	+ jphi*eta
ez		= -vr*bphi 	+ 	vphi *br + jz*eta

poyr	 = -	ez*bphi 	+ 	ephi*bz
poyz	 =	   er*bphi 	- 	ephi *br
poyphi =		-er*bz 	+	ez*br


r=x1
z=x2
dr=dx1
dz=dx2
i=1.5

power = ( vr^2 + vz^2 + vphi^2)/2 + (2.5)*pres/den + phi

acc_powerr = -den*vr*power 
acc_powerz =  den*vz*power 

!P.MULTI=[0,1,1]
window,23
plot, X1,acc_powerr(*,0)


!P.MULTI=0
;window,4
;contour, ms,levels=1
;window,3
;contour, ma,levels=1
;end
window,1


totpower=0
totpoynting=0
totkinetic=0

jheight=512

!P.MULTI=[0,1,4]
plot, mf(*,jheight )
ogg=fltarr(n1+1)
ogg(*)=0
for i=0,n1-1 do begin
	if ( mf(i, jheight) gt 1  ) then begin
		totpower 	 =  totpower	+ 2*!PI*x1(i)*dx1(i)* (acc_powerz(i,jheight) + poyz(i,jheight)  )
		totpoynting	 =  totpoynting+ 2*!PI*x1(i)*dx1(i)*  (poyz(i ,jheight) )
		totkinetic	 =  totkinetic	+ 2*!PI*x1(i)*dx1(i)* (acc_powerz(i,jheight)   )
		ogg(i) =  2*!PI*x1(i)*dx1(i)* (acc_powerz(i,jheight) + poyz(i,jheight)  ) 
		endif
	endfor
print, "totpower  , totpoynting , totkinetic"  
print, "jet ", totpower  , totpoynting , totkinetic  
pjet=totpower
plot, ogg

rad=20


print, "at radius, x1 ",x1(rad) 
totpower=0
totpoynting=0
totkinetic=0
ogg(*)=0
for j=0,n2/2-1 do begin
	if ( vr(rad,j) lt 0  ) then begin
		totpower 	 =  totpower	+ 2*!PI*x1(rad)*dx1(j)* (acc_powerr(rad,j) +$
		poyr(rad,j)  )
		totpoynting	 =  totpoynting+ 2*!PI*x1(rad)*dx1(j)*  (poyr(rad,j)  )
		totkinetic	 =  totkinetic	+ 2*!PI*x1(rad)*dx1(j)* (acc_powerr(rad,j)   )
		ogg(j)=2*!PI*x1(rad)*dx1(j)* (acc_powerr(rad,j) +$
		      poyr(rad,j)  )

		endif
	endfor
plot, x1,ogg
print, "disk ", totpower  , totpoynting , totkinetic  
paccin=totpower



totpower=0
totpoynting=0
totkinetic=0
ogg(*)=0
rad=30
print, "at radius, x1 ",x1(rad) 
for j=0,n2/2-1 do begin
	if ( vr(rad,j) lt 0  ) then begin
		totpower 	 =  totpower	+ 2*!PI*x1(rad)*dx1(j)* (acc_powerr(rad,j) +$
		poyr(rad,j)  )
		totpoynting	 =  totpoynting+ 2*!PI*x1(rad)*dx1(j)*  (poyr(rad,j)  )
		totkinetic	 =  totkinetic	+ 2*!PI*x1(rad)*dx1(j)* (acc_powerr(rad,j)   )
		ogg(j)=  2*!PI*x1(rad)*dx1(j)* (acc_powerr(rad,j) +$
		      poyr(rad,j)  )
		endif
	endfor
plot,x1, ogg
print, "disk", totpower  , totpoynting , totkinetic  
paccout=totpower
print, "PJet /Pacc", pjet/(paccout -paccin)



totpower		=fltarr(n1+1)
totpoynting	=fltarr(n1+1)
totkinetic	=fltarr(n1+1)
h	=fltarr(n1+1)
ogg(*)=0
rad=20
for rad=0,n1-1 do begin
for j=0,n2/2-1 do begin
	if ( vr(rad,j) lt 0  ) then begin
		totpower(rad) 	 = totpower(rad)	+ 2*!PI*x1(rad)*dx1(j)* (acc_powerr(rad,j) +$
		poyr(rad,j)  )
		totpoynting(rad)= totpoynting(rad)+ 2*!PI*x1(rad)*dx1(j)*  (poyr(rad,j)  )
		totkinetic(rad) = totkinetic(rad)	+ 2*!PI*x1(rad)*dx1(j)* (acc_powerr(rad,j)   )
		ogg(j)=  2*!PI*x1(rad)*dx1(j)* (acc_powerr(rad,j) +$
		      poyr(rad,j)  )
			h(rad)=x2(j)
		endif
	endfor
	endfor
;plot,x1, ogg
;print, "disk", totpower  , totpoynting , totkinetic  

window,16
!P.MULTI=[0,1,4]

plot, x1,totpower
plot, x1,totpoynting
plot, x1,totkinetic
plot, x1,h
end
