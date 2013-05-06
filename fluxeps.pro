
!P.COLOR=0
!P.CHARSIZE=1.
!P.THICK=2
!P.BACKGROUND=255
loadct,0
!p.font=0
set_plot,'ps'
device, /times
device, /encapsulated, filename="flux.eps", XSIZE=16,YSIZE=8

pload,out=0
bz=b2(*,0,0)
r=x1
dr=dx1
bzr=bz*r*dr
bflux=fltarr(n1+1)
for i=14,n1-1 do begin
bflux[i]=bflux[i-1]+bzr[i]
endfor
plot, x1,bflux,$
        /xlog,$
        xrange=[1,400], xstyle=1,$
        ;/ylog,$
        title="Magnetic Flux  !9y!X vs Radius, R, z=0",$
        xtitle="Radius, R",$
        ytitle="Magnetic Flux, !9y!X"


pload,out=299
bz=b2(*,0,0)
r=x1
dr=dx1
bzr=bz*r*dr
bflux=fltarr(n1+1)
for i=14,n1-1 do begin
bflux[i]=bflux[i-1]+bzr[i]
endfor
oplot, x1,bflux, linestyle=1

pload,out=nlast
bz=b2(*,0,0)
r=x1
dr=dx1
bzr=bz*r*dr
bflux=fltarr(n1+1)
for i=14,n1-1 do begin
bflux[i]=bflux[i-1]+bzr[i]
endfor
oplot, x1,bflux, linestyle=2

items=['Time='+string(time(0)), 'Time='+string(time(299)/2/!PI), 'Time='+string(time(nlast)/2/!PI)]
legend, items,linestyle=indgen(3)

device, /close
set_plot,'x'

end


