

nend=63
nend=2000
!P.COLOR=0
!P.BACKGROUND=255
!P.THICK=3
!P.CHARSIZE=1
!P.CHARTHICK=1
!Y.RANGE=0
;!Y.RANGE=[1e-9,1e-3]
!p.font=0
set_plot, 'ps'
device, Filename = 'bzlog.ps',XSIZE=16, YSIZE=8,/encapsulated
device, /times

loadct,0
pload,out=0
bz=b2(*,0,0)
 plot, x1, bz,$
        /xlog, /ylog, $
        title="Time ="+string(time(nlast)/2.0/!PI,format='(I6)')+'!9t!X!DK0!N',$
        xtitle="Log Radius",$
        ytitle="Log |B!Dz!N|",$
        xrange=[1.1,200], xstyle=1

 pload,out=299
bz=b2(*,0,0)
 oplot, x1, bz,linestyle=1
 pload,out=599
bz=b2(*,0,0)
 oplot, x1, bz,linestyle=2

ta=string(0,format='(I6)')
tb=string(time(299)/2/!PI,format='(I6)')
tc=string(time(599)/2/!PI,format='(I6)')
        items=['t='+ta,'t='+tb,'t='+tc]
 legend,items,linestyle=indgen(3), /right

device, /Close
set_plot, 'x'

end

