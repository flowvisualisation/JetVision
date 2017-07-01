nstep=1
nbeg=nlast
nend=nlast

 ;window,1,xs=900,ys=600
 !P.POSITION=0
 !P.MULTI=[0,3,2]

 !P.CHARSIZE=3
 !P.CHARTHICK=2.0
 !P.THICK=3.0



 pload,out=0
 bzinit = b2(*,*,0)
 pinit  = pr(*,*,0)
 dinit  = rho(*,*,0)
 v1init  = v1(*,*,0)
 muinit=bzinit^2/pinit
 window,xs=1000,ys=800

loadct,0
!p.background=255
!p.color=0


 for n=nbeg,nend,nstep do begin
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
 plot, x1,  (b1(*,0,0)^2+b2(*,0,0)^2 ) /pr(*,0,0)/muinit ,$
	/ylog,$
 	/xlog, $
	xrange=[1,10],$
 	title="!7l!3=b!U2!N/p"


v1ratio=v1(*,*,0)/v1init(*,*)
plot,x1,v1ratio(*,0),$
 	/xlog,$
	xrange=[1,10],$
 	title="v1/vinit"


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






xyouts, 0.01,0.01,$
	't(Inner disk rotations)='+string(time(n)/2./!PI,format='(F8.1)'),$
	/normal, charsize=3





ll=6
zero=''
nts=strcompress(string(n),/remove_all)
lnt=strlen(nts)
for j=1,ll-lnt do zero=zero+'0'
        fname='poynting.plots.'+zero+nts
;im=tvread(filename=fname,/nodialog,/png)


endfor
 !P.MULTI=0

;loadct,33
;display,(eta),x1=x1,x2=x2, /vbar,ims=3, xrange=[1,10],charsize=0.5
;contour,eta, /fill, nlevels=50
;window,2
;plot, eta(*,0)
 end
