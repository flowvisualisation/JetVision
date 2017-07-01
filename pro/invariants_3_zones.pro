


; K

filenum=126
filenum=nlast
pload,out=filenum
den=rho(0)
pres=pr(0)



vr=v1(0)
vz=v2(0)
vphi=v3(0)

br=b1(0)
bz=b2(0)
bphi=b3(0)

k= den*sqrt(vr^2+vz^2)/sqrt(br^2+bz^2)


rr=rebin(reform(x1,n1,1),n1,n2 )
zz=rebin(reform(x2,n2,1),n2,n1 )
zz=transpose(zz)


rad_sphere=sqrt(rr^2+zz^2)

lambda=vphi*rr-rr*bphi/k

omega= vphi/rr - k*bphi /den/rr


s=alog10(pres/den^(5./3.))



ke= 0.5*(vr^2  +vz^2 + vphi^2)
enth=(5./2.)*pres/den
mag= -(omega)*rr*bphi/k
grav= -1/rad_sphere


e= ke + enth+ mag + grav

r=x1
z=x2
dr=dx1
dz=dx2



!p.font=0

!P.POSITION=0
!P.MULTI=0
!P.CHARSIZE=1.0
!P.CHARTHICK=1.0
!P.THICK=1.0
set_plot,'ps'
device, /times
device, /encapsulated, filename="invariants_3_zones.eps"
i=1.6

nx=40
nx2=100
nx3=400
!x.range=[0,5]
!y.range=[-0.06,0.66]
xh=x2/x1(nx)/0.1
plot,xh, omega(*,nx),$
        title="!9W!X, Magnetic surface rotation rate"$
         +", t="+string(time(nlast)/2/!PI, format='(F5.1)')+'!9t!X!DK0!N',$
        xtitle="Normalised radius, r/h",$
        ytitle="!9W!X"
        
xh=x2/x1(nx2)/0.1
oplot,xh, omega(*,nx2),linestyle=1
xh=x2/x1(nx3)/0.1
oplot,xh, omega(*,nx3),linestyle=2

items=['r='+string(x1(nx),format='(F5.1)'),'r='+string(x1(nx2),format='(F5.1)'),'r='+string(x1(nx3),format='(F5.1)')]
legend,items,linestyle=indgen(3), /right


device, /close
set_plot,'x'
!P.COLOR=0
!P.BACKGROUND=255
;im=tvread(filename='Invariants', /nodialog, /png)


end
