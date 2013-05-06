

SET_PLOT, 'x'
; Read an image file containing elevation data  
file = FILEPATH('worldelv.dat', $  
   SUBDIRECTORY = ['examples', 'data'])  
image = READ_BINARY(file, DATA_DIMS = [360, 360])  
  
; Store the original display device  
oldDevice = !D.NAME  
DEVICE, GET_DECOMPOSED=oldDecomposed  
  
; Write the image to the 24-bit Z buffer device, using a  
; color lookup table. Note that we must set the DECOMPOSED  
; keyword to zero to create the color image.  
SET_PLOT, 'Z'  
ERASE ,'ffffff'x  
ERASE ,'000000'x  

DEVICE, SET_PIXEL_DEPTH=24  
DEVICE, SET_RESOLUTION=[720,720]  
DEVICE, DECOMPOSED=0  
LOADCT, 33  
TVLCT, 255,255,255, !D.TABLE_SIZE - 1  
  
;  tvimage, dist(200), position=[0.1,0.1,0.9,0.9]
n=nlast
pload,out=n

r=alog10(rho(0:511,0:511,0))
xr=x1(0:511)
xz=x1(0:511)

imin=min(r)
imax=max(r)

r=scale_vector(r,0,255)

pos=[0.1,0.2,0.9,0.9]
tvimage, r, position=pos



!P.CHARSIZE=1.8
!P.COLOR=0
contour,r,xr,xz, /nodata, /noerase, position=pos,$
  title='log(!7q!3), B field, S!DA!N, S!DF!N.'+ $
            'Time (IDR)='+string(time(n)/2./!PI,format='(F8.2)')


pos=[0.1,0.05,0.9,0.12]
colorbar, position=pos, range=[imin,imax],$
   format = '(f8.3)'




for i=1.3,10.3,1 do begin
;field_line, b1(0), b2(0), x1,x2,  i,0.2, rr,zz
;oplot, rr,zz,color=0
endfor

for i=10.3,20,2 do begin
;field_line, b1(0), b2(0), x1,x2,  i,0.2, rr,zz
;oplot, rr,zz,color=0
endfor



;posalf, b1,b2,b3,v1,v2,rho,pr,mf,ma,ms, mson

;contour, mf, x1,x2, levels=1, color=233, /overplot, c_thick=2.0

;contour, ma, x1,x2, levels=1, color=133, /overplot, c_thick=2.0



; Read the image array back from the Z buffer device  
; and write it to a TIFF file.  
new_image = TVRD(/TRUE)  
new_file = GETENV('IDL_TMPDIR')+'world.tif'  



ll=6
zero=''
nts=strcompress(string(n),/remove_all)
lnt=strlen(nts)
for j=1,ll-lnt do begin
zero=zero+'0'
endfor

fname='Field.SAD.'+zero+nts
new_file=fname


WRITE_TIFF, new_file, new_image  
  
; Change back to the original device  
SET_PLOT, oldDevice  
DEVICE, DECOMPOSED=oldDecomposed  
  
; Read in the TIFF file and display it.  
tif_image = READ_TIFF(new_file)  
WINDOW, XSIZE=720, YSIZE=720, TITLE='Image read from TIFF file'  
TV, tif_image, /TRUE  

end
