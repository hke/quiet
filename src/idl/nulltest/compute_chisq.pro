pro compute_chisq, modelfile, specfile, nbin, scale, nulltest

; Read model spectrum
fits2cl, cl_in, modelfile
lmax = n_elements(cl_in[*,0])-1
cl_fid = fltarr(lmax+1,6)
cl_fid[*,0] = cl_in[*,0]
cl_fid[*,1] = cl_in[*,3]
cl_fid[*,3] = cl_in[*,1]
for l = 0, lmax do begin
   cl_fid[l,*] = cl_fid[l,*] * l*(l+1)/2./!pi
end


; Read observed spectrum
c     = 'z'
bins  = fltarr(3)
cls   = fltarr(12)
close,58
openr,58,specfile
readf,58,c
chisq = 0.
n     = 0.
for i = 0L, nbin-1 do begin
   readf, 58, bins, cls
   for j = 0, 5 do begin
      if (cls[2*j+1] > 0.) then begin
;      if (cls[2*j+1] > 0. and j eq 5) then begin
         if (nulltest eq 1) then begin
            cl0 = 0.
         end else begin
            cl0 = mean(cl_fid[bins[1]:bins[2],j])
         end
         cls[2*j] = scale * cls[2*j]
         chisq = chisq + (cls[2*j]-cl0)^2 / (cls[2*j+1])^2
         print, cl0, cls[2*j], cls[2*j+1], (cls[2*j]-cl0)^2 / (cls[2*j+1])^2
         n     = n + 1.
      end
   end
end
close,58

pte = 1. - chisqr_pdf(chisq, n)
print, 'Chi-square = ', chisq, ', n = ', n, ', PTE = ', pte

end
