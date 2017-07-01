step=100

 ;window,1
 !P.POSITION=0
 !P.MULTI=[0,4,2]

 !P.CHARSIZE=3

 pload,out=0
 bzinit = b2(*,*,0)
 pinit  = pr(*,*,0)
 dinit  = rho(*,*,0)
 v1init  = v1(*,*,0)
 window,xs=1000,ys=800

loadct,0
!p.background=255
!p.color=0

 for n=0,nlast,nstep do begin
	 pload,out=n

	 mach_acc = v1(*,0,0)/sqrt(5.0/3.0*pr(*,0,0)/rho(*,0,0))

dratio=rho(*,*,0)/dinit(*,*)
plot,x1,dratio(*,0),$
 	/xlog, $
	xrange=[1,10],$
 	title="!7q!3/!7q!3init"

sound=sqrt(5.0/3.0*pr(*,0,0)/rho(*,0,0))

 plot, x1,  mach_acc/0.9/(sound*(x1^(0.5))) ,$
 	/xlog, $
	xrange=[1,10],$
 	title="Acc. Mach No/(!7a!3 c!DS!N /(!7X!3!DK!N r ))"



;	window,2
 plot, x1,  (b1(*,0,0)^2+b2(*,0,0)^2 ) /pr(*,0,0) ,$
	/ylog,$
 	/xlog, $
	xrange=[1,10],$
 	title="!7l!3=b!U2!N/p"
 plot, x1, pr(*,0,0) ,$
 	/xlog, /ylog,$
	xrange=[1,10],$
 	title='pr,t(IDR)='+string(time(n)/2./!PI,format='(F8.1)')


v1ratio=v1(*,*,0)/v1init(*,*)
plot,x1,v1ratio(*,0),$
 	/xlog,$
	xrange=[1,10],$
 	title="v1/vinit"

 plot, x1,  (rho(*,0,0)*v1(*,0,0))  ,$
 	/xlog, $
	xrange=[1,10],$
 	title="  !7q!3 v"

; plot, x1,  v1(*,0,0)*b2(*,0,0 )  ,$
; 	title="u!Dr!N B!Dz!N"

 dBrdz=b1(*,*,0)

 for i=0,n1-1 do begin
	 dBrdz(i,*)= deriv(x2,b1(i,*,0))
	 endfor
; plot,dBrdz(*,0),$
; 	title="dBrdz"

 jr=b1(*,*,0)
 for i=0,n1-1 do begin
	 jr(i,*)= deriv(x2,b3(i,*,0))
	 endfor

;torque=jr(*,*)*b2(*,*,0)
;plot, X1, X1*torque(*,0)



pratio=pr(*,*,0)/pinit(*,*)
plot,x1,pratio(*,0),$
	title="P/P!Dinit!N",$
 	/xlog, $
	xrange=[1,10]

bratio=b2(*,*,0)/bzinit(*,*)
plot,x1,bratio(*,0)^2,$
	title="(B!Dz!N/B!Dz,init!N)!U2!N",$
 	/xlog, $
	xrange=[1,10]




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



xyouts, 0.01,0.01,$
	't(Inner disk rotations)='+string(time(n)/2./!PI,format='(F8.1)'),$
	/normal, charsize=3





ll=6
zero=''
nts=strcompress(string(n),/remove_all)
lnt=strlen(nts)
for j=1,ll-lnt do zero=zero+'0'
        fname='poynting.plots.'+zero+nts
im=tvread(filename=fname,/nodialog,/png)


endfor
 !P.MULTI=0

;loadct,33
;display,(eta),x1=x1,x2=x2, /vbar,ims=3, xrange=[1,10],charsize=0.5
;contour,eta, /fill, nlevels=50
;window,2
;plot, eta(*,0)
 end
