nstep=1
nbeg=599
nend=nbeg

 ;window,1,xs=900,ys=600
 !P.POSITION=0
 !P.MULTI=0
; !P.MULTI=[0,3,2]

 !P.CHARSIZE=2
 !P.CHARTHICK=2.0
 !P.THICK=3.0



 pload,out=0
 bzinit = b2(*,*,0)
 pinit  = pr(*,*,0)
 dinit  = rho(*,*,0)
 v1init  = v1(*,*,0)
 muinit=bzinit^2/pinit
 window,xs=700,ys=700

loadct,0
!p.background=255
!p.color=0


xbeg=1.2
xend=10.0
 for n=nbeg,nend,nstep do begin
	 pload,out=n


plot,x1,b2(*,0,0),$
	xtitle="Radius",$
 	/xlog, $
	xstyle=1,xrange=[xbeg,xend],$
 	title="B!DZ!N"
set_plot,'ps'
device,filename='bz.eps', /encapsulated
plot,x1,b2(*,0,0),$
	xtitle="Radius",$
 	/xlog, $
	xstyle=1,xrange=[xbeg,xend],$
 	title="B!DZ!N"
device,/close
set_plot,'x'

	 mach_acc = v1(*,0,0)/sqrt(5.0/3.0*pr(*,0,0)/rho(*,0,0))


dratio=rho(*,*,0)/dinit(*,*)
plot,x1,dratio(*,0),$
	xtitle="Radius",$
 	/xlog, $
	xstyle=1,xrange=[xbeg,xend],$
 	title="!7q!3/!7q!3init"
set_plot,'ps'
device,filename='dratio.eps', /encapsulated
plot,x1,dratio(*,0),$
	xtitle="Radius",$
 	/xlog, $
	xstyle=1,xrange=[xbeg,xend],$
 	title="!7q!3/!7q!3init"
device,/close
set_plot,'x'

sound=sqrt(5.0/3.0*pr(*,0,0)/rho(*,0,0))

 plot, x1,  mach_acc/0.9/(sound*(x1^(0.5))) ,$
	xtitle="Radius",$
 	/xlog, $
	xstyle=1,xrange=[xbeg,xend],$
 	title="Acc. Mach No/(!7a!3 c!DS!N /(!7X!3!DK!N r ))"



!P.font=0
set_plot,'ps'
device,filename='mach_acc.eps', /encapsulated
device, xs=16,ys=8
device, /times
!P.CHARSIZE=1.0
 plot, x1,  mach_acc/0.9/(sound*(x1^(0.5))) ,$
	xtitle="Radius, R",$
 	/xlog, $
	xstyle=1,xrange=[xbeg,xend],$
 	ytitle="m!Ds!N/m!Dth.!N at t="+string(time(n)/2/!PI, format='(F8.1)')+'!9t!3!DK0!N'
 	;ytitle="M!Dacc!N/(!7a!3 c!DS!N /(!7X!3!DK!N r )) at t="+string(time(n)/2/!PI, format='(F8.1)')+'!7s!3!DK0!N'
device,/close

set_plot,'x'

 plot, x1,  (b1(*,0,0)^2+b2(*,0,0)^2 ) /pr(*,0,0)/muinit ,$
	xtitle="Radius",$
	/ylog,$
 	/xlog, $
	xstyle=1,xrange=[xbeg,xend],$
 	title="!7l!3=b!U2!N/p"


set_plot,'ps'
device,filename='mu.eps', /encapsulated
 plot, x1,  (b1(*,0,0)^2+b2(*,0,0)^2 ) /pr(*,0,0)/muinit ,$
	xtitle="Radius",$
	/ylog,$
 	/xlog, $
	xstyle=1,xrange=[xbeg,xend],$
 	title="!7l!3=b!U2!N/p"
device,/close
set_plot,'x'

v1ratio=v1(*,*,0)/v1init(*,*)
plot,x1,v1ratio(*,0),$
	xtitle="Radius",$
 	/xlog,$
	xstyle=1,xrange=[xbeg,xend],$
 	title="v1/vinit"


set_plot,'ps'
device,filename='v1ratio.eps', /encapsulated
plot,x1,v1ratio(*,0),$
	xtitle="Radius",$
 	/xlog,$
	xstyle=1,xrange=[xbeg,xend],$
 	title="v1/vinit"
device,/close
set_plot,'x'

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
	xtitle="Radius",$
	title="P/P!Dinit!N",$
 	/xlog, $
	xstyle=1,xrange=[xbeg,xend]


set_plot,'ps'
device,filename='pratio.eps', /encapsulated
plot,x1,pratio(*,0),$
	xtitle="Radius",$
	title="P/P!Dinit!N",$
 	/xlog, $
	xstyle=1,xrange=[xbeg,xend]
device,/close
set_plot,'x'

bratio=b2(*,*,0)/bzinit(*,*)
plot,x1,bratio(*,0)^2,$
	xtitle="Radius",$
	title="(B!Dz!N/B!Dz,init!N)!U2!N",$
 	/xlog, $
	xstyle=1,xrange=[xbeg,xend]




set_plot,'ps'
device,filename='bratio.eps', /encapsulated
plot,x1,bratio(*,0)^2,$
	xtitle="Radius",$
	title="(B!Dz!N/B!Dz,init!N)!U2!N",$
 	/xlog, $
	xstyle=1,xrange=[xbeg,xend]
device,/close
set_plot,'x'



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
;display,(eta),x1=x1,x2=x2, /vbar,ims=3, xstyle=1,xrange=[xbeg,xend],charsize=0.5
;contour,eta, /fill, nlevels=50
;window,2
;plot, eta(*,0)
 end
