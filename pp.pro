 ;window,1
 !P.MULTI=[0,4,2]

 !P.CHARSIZE=3
 window,xs=1000,ys=800

loadct,0
!p.background=255
!p.color=0

nstep=10
 for n=0,nlast,nstep do begin
	 pload,out=n
 plot, x1,v1(*,0,0)/sqrt(5.0/3.0*pr(*,0,0)/rho(*,0,0)),$
 	title="Accretion Mach No."
;	window,2
 plot, x1,  (b1(*,0,0)^2+b2(*,0,0)^2 ) /pr(*,0,0) ,$
 	/xlog, /ylog,$
	xrange=[1,10],$
 	title="!7l!3=b!U2!N/p"
 plot, x1, pr(*,0,0) ,$
 	/xlog, /ylog,$
	xrange=[1,10],$
 	title="pressure"

 plot, x1,  (b1(*,0,0)^2+b2(*,0,0)^2 )  ,$
 	/xlog, /ylog,$
	xrange=[1,10],$
 	title="b!U2!N"

 plot, x1,  (rho(*,0,0)*v1(*,0,0))  ,$
 	title="Accreting Mass Flux"

 plot, x1,  v1(*,0,0)*b2(*,0,0 )  ,$
	xrange=[1,10],$
 	title="u!Dr!N B!Dz!N"

 dBrdz=b1(*,*,0)

 for i=0,n1-1 do begin
	 dBrdz(i,*)= deriv(x2,b1(i,*,0))
	 endfor
; plot,dBrdz(*,0),$
;	xrange=[1,10],$
; 	title="dBrdz"


 jr=b1(*,*,0)
 for i=0,n1-1 do begin
    jr(i,*)= deriv(x2,b3(i,*,0))
    endfor

torque=jr(*,*)*b2(*,*,0)
plot, X1, X1*torque(*,0),title='t='+string(time(n)/2./!PI)



rr=rebin(reform(x1,n1,1),n1,n2 )
zz=rebin(reform(x2,n2,1),n2,n1 )
zz=transpose(zz)

	
 eta=pr(*,*,0)/rho(*,*,0)
	r_sphere=sqrt(rr*rr+zz*zz)
	gmm=5./3.
 for j=0,n2-1 do begin
 eta(*,j)=pr(*,j,0)/rho(*,j,0)
	 endfor
 eta(*,*)= eta(*,*) + ((gmm-1)/gmm)* (1/r_sphere(*,*) - 1/rr(*,*))
 alpha=0.9
 eta(*,*)=alpha*sqrt(1/rr(*,*)^3 )*eta(*,*)*rho(*,*,0)
eta(*,*) =eta(*,*)*(eta(*,*) gt 0.0)



;plot, eta(*,0)
 plot, x1, eta(*,0)*(dBrdz(*,0)- deriv(x1,b2(*,0,0 )))  ,$
	xrange=[1,10],$
 	title="!7g!3J!D!7u!3!N"
	oplot,x1,- v1(*,0,0)*b2(*,0,0 )  




ll=6
zero=''
nts=strcompress(string(n),/remove_all)
lnt=strlen(nts)
for j=1,ll-lnt do zero=zero+'0'
        fname='accretion.plots.'+zero+nts
im=tvread(filename=fname,/nodialog,/png)


endfor
 !P.MULTI=0

;loadct,33
;display,(eta),x1=x1,x2=x2, /vbar,ims=3, xrange=[1,10],charsize=0.5
;contour,eta, /fill, nlevels=50
;window,2
;plot, eta(*,0)
 end
