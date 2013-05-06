pro interplot,quant,r,z,dr,dz,rr,zz,qinter,rpl,zpl,dpl

nr = n_elements(r)
nz = n_elements(z)

istart = max([min(where(rr ge r[0])),min(where(zz ge z[0]))])
iend   = min([max(where(rr le r[nr-1])),max(where(zz le z[nz-1]))])

rpl = rr[istart:iend]
zpl = zz[istart:iend]
dpl = rpl
npl = n_elements(rpl)
qinter = fltarr(npl)

for i = 0L,npl-1 do begin

 i0 = (where(abs(rpl[i]-r) lt 0.51*dr))[0]
 i0 = max([0,i0])
 i0 = min([nr-1,i0])
 scrh1 = interpol(quant[i0,*],z,zpl[i])

 j0 = (where(abs(zpl[i]-z) lt 0.51*dz))[0]
 j0 = max([0,j0])
 j0 = min([nz-1,j0])
 scrh2 = interpol(quant[*,j0],r,rpl[i])

 qinter[i] = -quant[i0,j0] + scrh1 + scrh2

 if (i eq 0 ) then begin
  dpl[i] = 0.0
 endif else begin
  distr = rpl[i]-rpl[i-1]
  distz = zpl[i]-zpl[i-1]
  dist = sqrt(distr^2+distz^2)
  dpl[i] = dpl[i-1]+dist
 endelse

endfor

return
end



