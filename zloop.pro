nstep=1
nbeg=nlast
nend=nlast
loadct,33
nx=n1
ny=n2
;pload,out=nlast


xres=800
yres=900
   thisDevice = !D.Name
   Set_Plot, 'Z'
   Erase
   Device, Set_Resolution=[xres,yres], Z_Buffer=0
   tvlct,r,g,b,/get
   loadct,0   ; -- background in black and white --
   background=255
   old_background = !P.BACKGROUND
   !P.BACKGROUND  = background
   tvlct,r,g,b 
ERASE ,'ffffff'x

loadct,33

for n=nbeg,nend,nstep do begin
	pload,out=n

br=b1(0:511,0:511,0)
bz=b2(0:511,0:511,0)
xr=x1(0:511)
xz=x2(0:511)

pos=[0.,0.,1.0,1.0]
r=fltarr(200,200)
r=dist(200)
r=scale_vector(r,0,255)
r(*,*)=100
tvimage,r,position=pos


r=alog10(rho(0:511,0:511,0))
imin=min(r)
imax=max(r)

r=scale_vector(r,0,255)

!P.CHARSIZE=2.0
pos=[0.1,0.2,0.9,0.9]


tvimage,r,position=pos
contour,r,xr,xz, /nodata, /noerase, position=pos,$
  title='log(!7q!3), B field, S!DA!N, S!DF!N.'+ $
          'Time (IDR)='+string(time(n)/2./!PI,format='(F8.3)')



pos=[0.1,0.05,0.9,0.12]
colorbar, position=pos, range=[imin,imax],$
	format = '(f8.3)'


;current, b3, x1, rbphi,n2
;
;minc=min(rbphi)
;maxc=max(rbphi)
;num=32
;lev=findgen(32)*(maxc-minc)/32+minc
;for i=0,31 do begin
;	contour,rbphi,x1,x2,/overplot,levels=lev(i)
;	endfor



for i=1.3,10.3,1 do begin
;field_line, b1(0), b2(0), x1,x2,  i,0.2, rr,zz
;oplot, rr,zz,color=255
endfor

for i=10.3,20,2 do begin
;field_line, b1(0), b2(0), x1,x2,  i,0.2, rr,zz
;oplot, rr,zz,color=255
endfor

posalf, b1,b2,b3,v1,v2,rho,pr,mf,ma,ms, mson

contour, mf, x1,x2, levels=1, color=233, /overplot, c_thick=2.0

contour, ma, x1,x2, levels=1, color=133, /overplot, c_thick=2.0

ll=6
zero=''
nts=strcompress(string(n),/remove_all)
lnt=strlen(nts)
for j=1,ll-lnt do begin
zero=zero+'0'
endfor

fname='Field.SAD.'+zero+nts
print,fname
;im=tvread(filename=fname,/nodialog,/png)

TVLCT, r, g, b, /Get
  snapshot = TVRD()


   image24 = BytArr(3, xres, yres)
   image24[0,*,*] = r[snapshot]
   image24[1,*,*] = g[snapshot]
   image24[2,*,*] = b[snapshot]

	Write_PNG, fname+'.png', image24
endfor


   Device, Z_Buffer=1
	   Set_Plot, thisDevice
end
