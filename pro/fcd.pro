
;window,3,xs=900,ys=600
window,3
!P.Charsize=3.0
!P.thick=2.0
!x.thick=2.0
!y.thick=2.0
!P.charthick=2.0
loadct,0
!P.background=255

qq=250
;display, alog(rho(0)), x1=x1,x2=x2, /gstretch, imsize=0.25

field_line, b1(0), b2(0), x1,x2, 5,60,rr,zz
vpol=94*sqrt(v1(0)^2+v2(0)^2)
vphi=94*v3(0)
r=x1
z=x2
dr=dx1
dz=dx2
interplot,vpol ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
interplot,vphi ,r,z,dr,dz,rr,zz,qinter2,rpl,zpl
qinter2= 0.1*x1*qinter2

q1=35
plot, qinter1(q1:qq), qinter2(q1:qq), color=100, xtitle='Vpol [km/s]', ytitle='R V!7U!3 [ AU km/s]'
xyouts, qinter1(qq-60), qinter2(qq-20), '<0.1 AU', color=100


qq=400
field_line, b1(0), b2(0), x1,x2, 10,5,rr,zz
vpol=94*sqrt(v1(0)^2+v2(0)^2)
vphi=94*v3(0)
r=x1
z=x2
dr=dx1
dz=dx2
interplot,vpol ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
interplot,vphi ,r,z,dr,dz,rr,zz,qinter2,rpl,zpl
qinter2= 0.1*x1*qinter2

oplot, qinter1(0:qq), qinter2(0:qq), color=100
xyouts, qinter1(qq), qinter2(qq), '1 AU', color=100


field_line, b1(0), b2(0), x1,x2, 20,5,rr,zz
vpol=94*sqrt(v1(0)^2+v2(0)^2)
vphi=94*v3(0)
r=x1
z=x2
dr=dx1
dz=dx2
interplot,vpol ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
interplot,vphi ,r,z,dr,dz,rr,zz,qinter2,rpl,zpl
qinter2= 0.1*x1*qinter2


oplot, qinter1(0:qq), qinter2(0:qq), color=100
xyouts, qinter1(qq), qinter2(qq), '2 AU', color=100


field_line, b1(0), b2(0), x1,x2, 30,5,rr,zz
vpol=94*sqrt(v1(0)^2+v2(0)^2)
vphi=94*v3(0)
r=x1
z=x2
dr=dx1
dz=dx2
interplot,vpol ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
interplot,vphi ,r,z,dr,dz,rr,zz,qinter2,rpl,zpl
qinter2= 0.1*x1*qinter2

oplot, qinter1(0:qq), qinter2(0:qq), color=100
xyouts, qinter1(qq), qinter2(qq), '3 AU', color=100




field_line, b1(0), b2(0), x1,x2, 40,5,rr,zz
vpol=94*sqrt(v1(0)^2+v2(0)^2)
vphi=94*v3(0)
r=x1
z=x2
dr=dx1
dz=dx2
interplot,vpol ,r,z,dr,dz,rr,zz,qinter1,rpl,zpl
interplot,vphi ,r,z,dr,dz,rr,zz,qinter2,rpl,zpl
qinter2= 0.1*x1*qinter2

oplot, qinter1(0:qq), qinter2(0:qq), color=100
xyouts, qinter1(qq), qinter2(qq), '4 AU', color=100


;plot, qinter1, qinter2, color=100, xtitle='v_z', ytitle='r v_phi'

im=tvread(filename='fcd', /nodialog, /png)
end
