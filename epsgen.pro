pro epsgen, field, num
;field='bx'
;num=18000
strnum=string(num, format='(I05)')
fi=FINDFILE(field+'t*.data.gz')
nfi=size(fi,/N_ELEMENTS)

nend=nfi-1
;nend=20
;nbeg=nend


xres=800
yres=200

set_plot,'ps'

!p.font=0

device, /times
device, filename=field+strnum+".eps", /encapsulated, /color
DEVICE, BITS_PER_PIXEL=8, COLOR=1
DEVICE, XSIZE=8, YSIZE=4, /INCHES

   loadct,33   ; -- background in black and white --
	tvlct, 255,255,255, !D.table_size-1
	tvlct, 0,0,0, !D.table_size-2
;!p.background=255
!p.color=254


erase, 'ffffff'x

nmfc=field+'t'+strnum+".data.gz"
   @FIELDS_EVOL
   ndim=size(serie,/DIMENSIONS) 
   mini=min(serie)
   maxi=max(serie)
; data range to be plotted
   rmin=mini
   rmax=maxi
   x1=0
   x2=ndim(0)-1
   y1=0
   y2=ndim(1)-1
   z1=0
;   z2=ndim(2)-1
   z2=0

   stepx=x2-x1+2
   stepy=y2-y1+2
   stepz=z2-z1+2

   ydim=10000.0
   xdim=ydim*ndim(0)/ndim(1)

fxy=fltarr(x2-x1+1, y2-y1+1)

for k=0L,y2-y1 do begin
for j=0L,x2-x1 do begin
fxy(j,k)=serie(x1+j,y1+k,z2)
endfor
endfor


print, size(fxy)
fxy=fxy[4500:18000-4501,*]
bx=fxy

b=size(fxy)
pos=[0.,0.,1.0,1.0]
r=scale_vector(fxy,0,253)


imin=min(fxy)
imax=max(fxy)


!P.CHARSIZE=1.0
pos=[0.05,0.25,0.95,0.9]


print, b
xr=findgen(b(1))/b(1)*227./2-227./2/2
xz=findgen(b(2))/b(2)*5
tvimage,r,position=pos
!x.style=1
contour,r,xr,xz, /nodata, /noerase, position=pos,$
  title=field+ ', t ='+string(tm*5.3e-3,format='(I10)')+"!9w!X!Dp!N!U-1!N"



pos=[0.05,0.1,0.95,0.15]
colorbar, position=pos, range=[imin,imax],$
	format = '(f8.3)'






device, /close
;	   Set_Plot, thisDevice
end
