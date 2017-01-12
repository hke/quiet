pro comm_gelman_rubin, infile, outfile, spec, burnin

cls_in = readfits(infile)
 
lmax     = size(cls_in(*,0,0,0),/n_elements)-1
numspec  = size(cls_in(0,*,0,0),/n_elements)
numchain = size(cls_in(0,0,*,0),/n_elements)
numsim   = size(cls_in(0,0,0,*),/n_elements)-1

cls = fltarr(lmax+1, numchain, numsim)

for l = 0, lmax do begin
    for i = 0, numchain-1 do begin
        for k = 0, numsim-1 do begin
            cls[l,i,k] = cls_in[l,spec,i,k+1L]
        end
    end
end

mean_chains = fltarr(numchain)

n = 1.*(numsim-burnin)
m = 1.*numchain

openw,38,outfile
for l = 2L, lmax do begin
    W = 0.
    for i = 0, numchain-1 do begin
        mean_chains[i] = mean(cls[l,i,burnin:numsim-1])
        W = W + total((cls[l,i,burnin:numsim-1]-mean_chains[i])^2)
    end
    W = W / (1.*m*(n-1.))

    totmean = mean(mean_chains)
    B = n/(m-1.) * total((mean_chains-totmean)^2)

    V = (1.-1./n)*W + 1./n * B

    R = V/W

    printf,38, l, R
end
close,38



end
