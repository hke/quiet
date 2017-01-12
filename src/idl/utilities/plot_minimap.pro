pro plot_minimap, infile, nside, order, npix, pol

map = fltarr(12L*nside^2,3)
openr,58,infile

a = 'a'
pix = 0L
val = 0.
pixmax = 0L
valmax = -1d30
readf,58,a
map[*,*] = -1.6375d30
for i = 0, npix-1 do begin
   readf,58, pix, val
   map[pix,1] = val
   if (val gt valmax) then begin
      pixmax = pix
      valmax = val
   end
end

for i = 0, npix-1 do begin
   readf,58, pix, val
   map[pix,2] = val
   if (val gt valmax) then begin
      pixmax = pix
      valmax = val
   end
end
close,58

pix2ang_nest, nside, pixmax, theta, phi
theta = 90.-theta*180./3.14
phi   = phi*180./3.14

gnomview, map[*,pol], /nest, rot=[phi,theta]

end
