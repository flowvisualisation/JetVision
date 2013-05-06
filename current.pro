pro current, b3, x1, rbphi, n2

bphi=b3(0)
rbphi = bphi

for i=1,n2-1 do begin
rbphi(*,i)=x1*bphi(*,i)
endfor

end
