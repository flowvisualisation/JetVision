
!P.POSITION=0
;ycut=1536
;ycut=2700
window, xs=1000, ys=500
!P.multi=[0,3,2]
!P.Charsize=3.0
!P.thick=2.0

!X.STYLE=1



bphi=b3(0)
rbphi = bphi

for i=1,3071 do begin
	rbphi(*,i)=x1*bphi(*,i)
	endfor


xend=1023
!P.background=255
loadct,0
 plot,  x1(0:xend), (v2(*,ycut, 0)), color=100 , title='Vz, h='+string(x2(ycut)), xtitle='r'

;xyouts, 10,0.2, 'Vz',  color=100	
 plot,  x1(0:xend), (b3(*,ycut, 0)), color=100  , title='B!7U!3', xtitle='r'
 plot,  x1(0:xend), (rho(*,ycut, 0)), color=100 , title='!7q!3 ', xtitle='r'


 plot,  x1(0:xend), (v3(*,ycut, 0)), color=100  , title='V!7U!3', xtitle='r'
 plot,  x1(0:xend), (b2(*,ycut, 0)), color=100 , title='Bz', xtitle='r'
 plot,  x1(0:xend), (rbphi(*,ycut)), color=100 , title='I', xtitle='r'

;xyouts, 100,100, 'Vz', /device, color=100	
im=tvread(filename="JetProfiles"+string(ycut, format='(I04)'),/png,/nodialog)

!P.multi=0

!P.multi=[0,1,6]


!P.font=0
set_plot,'ps'
device, filename="JetProfiles"+string(ycut, format='(I04)')+".eps" 
device, /encapsulated
device, /times
device, xsize=4, ysize=8, /inches
!P.CHARSIZE=1.0

 plot,  x1(0:xend), (v2(*,ycut, 0)), color=100 , title='Vz, h='+string(x2(ycut)), xtitle='r'
 plot,  x1(0:xend), (b3(*,ycut, 0)), color=100  , title='B!9f!X', xtitle='r'
 plot,  x1(0:xend), (rho(*,ycut, 0)), color=100 , title='!9r!X ', xtitle='r'
 plot,  x1(0:xend), (v3(*,ycut, 0)), color=100  , title='V!9f!X', xtitle='r'
 plot,  x1(0:xend), (b2(*,ycut, 0)), color=100 , title='Bz', xtitle='r'
 plot,  x1(0:xend), (rbphi(*,ycut)), color=100 , title='I', xtitle='r'

device, /close

set_plot,'x'


!P.multi=0

end
