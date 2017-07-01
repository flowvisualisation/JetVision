
	vr=v1(0)
	pres=pr(0)
	den=rho(0)


	ms=vr/sqrt(pres/den)
	ms=vr;/sqrt(pres/den)


	

nx=20
xbeg=0
xend=2
plot, x2, ms(nx,*) , title=' ms, r=' +string(x1(nx)),$
xrange=[xbeg,xend]
oplot,x2,ms(30,*)
oplot,x2,ms(40,*)

im=tvread(filename="ms", /png, /nodialog)

end
