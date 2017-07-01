
nstep=5
loadct,0
!P.Charsize=2.0
!P.Multi=[0,3,2]
!P.background=255
!P.color=0
window,2,xs=1100,ys=750
pload,out=0
 xbeg=0.
 xend=40.
 i=0
 plot, x1,v1(*,0,0), xrange=[xbeg,xend], xstyle=1,title="v1, t="+string(time(i))
 plot, x1,v2(*,0,0), xrange=[xbeg,xend], xstyle=1,title="v2"
 plot, x1,v3(*,0,0), xrange=[xbeg,xend], xstyle=1,title="v3"
 plot, x1,pr(*,0,0), xrange=[xbeg,xend], xstyle=1,title="pr"
 plot, x1,rho(*,0,0), xrange=[xbeg,xend], xstyle=1,title="rho"
 plot, x1,rho(*,0,0), xrange=[xbeg,xend], xstyle=1,title="rho"
!P.Multi=[0,3,2]
for i=1,nlast,nstep do begin
	pload,out=i
 plot, x1,v1(*,0,0), xrange=[xbeg,xend], xstyle=1,title="v1, t="+string(time(i))
 plot, x1,v2(*,0,0), xrange=[xbeg,xend], xstyle=1,title="v2"
 plot, x1,v3(*,0,0)*sqrt(x1), xrange=[xbeg,xend], xstyle=1,title="v3"
 plot, x1,pr(*,0,0), xrange=[xbeg,xend], xstyle=1,title="pr"
 plot, x1,rho(*,0,0), xrange=[xbeg,xend], xstyle=1,title="rho"
 plot, x1,pr(*,0,0)/rho(*,0,0)^(5./3.), /ylog,xrange=[xbeg,xend], xstyle=1,title="entropy"
 print,i, ' of ', nlast
 wait,1
endfor
!P.Multi=0
;window,3
; plot, x1,v1(*,0,0), xrange=[0.5,xend], xstyle=1

end
