pro plotgen1d, fname, xaxis, yaxis, title, xtitle, usingps, xrange


;for usingps=0,1 do begin
phistring='!7u!X'
if ( usingps ) then begin
set_plot,'ps'
device,filename=fname+'.eps',/encapsulated, /color
phistring='!9f!X'
!p.font=0
device, /times
xs=5.
ys=3
DEVICE, XSIZE=xs, YSIZE=ys, /INCHES
!x.range=xrange
!x.style=1
endif else begin
set_plot,'x'
!p.font=-1
!p.color=0
!p.background=255
window, title=fname
endelse


plot, xaxis, yaxis , /xlog,  title=title, xtitle=xtitle


if ( usingps ) then begin
device,/close
set_plot,'x'
endif else begin
set_plot,'x'
im=tvread(filename=fname,/nodialog,/png)
endelse

;endfor
set_plot,'x'

end
