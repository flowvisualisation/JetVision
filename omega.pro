
den=rho(0)
pres=pr(0)



vr=v1(0)
vz=v2(0)
vphi=v3(0)

br=b1(0)
bz=b2(0)
bphi=b3(0)

k= den*sqrt(vr^2+vz^2)/sqrt(br^2+bz^2)

rr=rebin(reform(x1,n1,1),n1,n1 )
zz=transpose(rr)


 i=1.6
 
 r=x1
 z=x2
 dr=dx1
 dz=dx2

 field_line, br, bz, x1,x2,  i,0.05, r1,z1
 interplot,vphi/rr ,r,z,dr,dz,r1,z1,qinter1,rpl,zpl

!P.POSITION=0
 !P.MULTI=[0,1,3]
 !P.CHARSIZE=2
 plot, zpl, qinter1



omega=vphi/rr

nx=20
plot, x2, omega(nx,*),title="r="+string(x1(nx))
nx=40
plot, x2, omega(nx,*),title="r="+string(x1(nx))
 end

