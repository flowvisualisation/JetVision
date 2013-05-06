

nend=63
nend=nlast
!P.COLOR=0
!P.BACKGROUND=255
!P.THICK=3
!P.CHARSIZE=1
!P.CHARTHICK=1
!Y.RANGE=[1e-5,1e-1]
!p.font=0
set_plot, 'ps'
device, Filename = 'mulog.eps',XSIZE=16, YSIZE=8,/encapsulated
device, /times

loadct,0
pload,out=0
 plot, x1, (b1(*,0,0)^2+b2(*,0,0)^2)/pr(*,0,0),$
        /xlog, /ylog, $
        title="Time ="+string(time(nend)/2.0/!PI,format='(I6)')+'!9t!X!DK!N',$
        xtitle="Log Radius",$
        ytitle="Log !9m!X",$
        xrange=[1.15,200], xstyle=1
 pload,out=nend/2
 oplot, x1, (b1(*,0,0)^2+b2(*,0,0)^2)/pr(*,0,0),linestyle=1
 pload,out=nend
 oplot, x1, (b1(*,0,0)^2+b2(*,0,0)^2)/pr(*,0,0),linestyle=2

ta=string(0,format='(I6)')
tb=string(time(nend/2)/2/!PI,format='(I6)')
tc=string(time(nend)/2/!PI,format='(I6)')
        items=['t='+ta,'t='+tb,'t='+tc]
 legend,items,linestyle=indgen(3), /right

device, /Close
set_plot, 'x'

end

