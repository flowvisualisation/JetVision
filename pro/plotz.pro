
den=rho(0)
pres=pr(0)

vr=v1(0)
vz=v2(0)
vphi=v3(0)

br=b1(0)
bz=b2(0)
bphi=b3(0)


rr=rebin(reform(x1,n1,1),n1,n1 )
zz=transpose(rr)


rad=sqrt(rr^2 +zz^2)





window,8,xs=600,ys=800
!P.POSITION=0
!P.CHARSIZE=3.
!P.MULTI=0
!P.MULTI=[0,1,4]


nx=40
for nx=20,40,10 do begin
print,x1(nx) 

xbeg=0
xend=1
plot, x1, alog10(den(nx,*)),$
		title="log10(!7q!3) at r="+string(x1(nx)),$
		xrange=[xbeg,xend],$
		xtitle='z '
temp=pres/den
plot, x1,alog10(temp(nx,*)),$
		title="temp at r="+string(x1(nx)),$
		xrange=[xbeg,xend],$
		xtitle='z '

cs=sqrt(5.0/3.0*pres/den)
mach=vz/cs

plot, x1,mach(nx,*),$
		title="vz/cs at r="+string(x1(nx)),$
		xrange=[xbeg,xend],$
		xtitle='z '
omega=vphi/rr
ybeg=0.4
yend=0.6
plot, x1,omega(nx,*)*x1(nx)^1.5,$
		title="vphi/r at r="+string(x1(nx)),$
		xrange=[xbeg,xend],$
;		yrange=[ybeg,yend],$
		xtitle='z '

ll=6
zero=''
nts=strcompress(string(nx),/remove_all)
lnt=strlen(nts)
for j=1,ll-lnt do zero=zero+'0'
	        fname='zplot.'+zero+nts
			  im=tvread(filename=fname,/nodialog,/png)

endfor

end


