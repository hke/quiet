function compute_conf_region2, hist, percentile, minval, maxval

numbin = size(hist(0,*),/n_elements)
minval = hist(0,0)
maxval = hist(0,numbin-1)
hist_int = hist(1,*)

totnum = total(hist_int)
sumhist = max(hist_int,maxbin)

lower = maxbin
upper = maxbin
;plot, hist_int
while sumhist lt percentile*totnum do begin
    if (upper eq numbin-1) then begin
        lower = lower-1
        sumhist = sumhist + hist_int(lower)
    endif else if (lower eq 0) then begin
        upper = upper+1
        sumhist = sumhist + hist_int(upper)
    endif else if hist_int(upper) ge hist_int(lower) then begin
        upper = upper+1
        sumhist = sumhist + hist_int(upper)
    endif else begin
        lower = lower-1
        sumhist = sumhist + hist_int(lower)
    endelse
end

maxbin = minval + (maxbin+0.5) * (maxval-minval)/numbin
if (lower eq 0) then begin
;   lower = minval - 0.5 * (maxval-minval)/numbin
   lower = 0.
end else begin
   lower  = minval + (lower+0.5)  * (maxval-minval)/numbin
end
upper  = minval + (upper+0.5)  * (maxval-minval)/numbin

return, [maxbin,lower,upper]

end


function compute_upper_limit_from_dist, hist, percentile

numbin = size(hist(0,*),/n_elements)
hist_int = hist(1,*)

totnum = total(hist_int)

sumhist = 0.
upper   = 0
while sumhist lt percentile*totnum do begin
   sumhist = sumhist + hist[1,upper]
   upper   = upper + 1
end

return, hist[0,upper-1]

end

function compute_upper_limit, samples, percentile

sorted = samples[sort(samples)]

return, sorted[percentile*n_elements(sorted)]

end


pro comm_print_conf_regions_multi, infile1, infile2, outfile1, outfile2, outfile_tot, spec, numreg, numbin, numexclude, l_gauss_low, l_gauss_high, l_low, l_high, dl
;*****************************************************************************
;+
; NAME:
;       comm_print_conf_regions
;
; PURPOSE:
;       Output information from a Commander-based analysis
;
; CALLING SEQUENCE:
;       comm_analysis, filename
;
; INPUTS:
;       filename: FITS-file containing a set of power spectra
;;
; OUTPUTS:
;
; MODIFICATION HISTORY:
;       March 1, 2004    Hans Kristian Eriksen, JPL
;                        Version 1.0
;

; cls is a three-dimensional array, containing the output from the
; commander program in the form cls(0:lmax, 1:numchain, 1:numsim)
  cls_in = readfits(infile1)
  lmax1     = size(cls_in(*,0,0,0),/n_elements)-1
  numspec1  = size(cls_in(0,*,0,0),/n_elements)
  numchain1 = size(cls_in(0,0,*,0),/n_elements)
  numsim1   = size(cls_in(0,0,0,*),/n_elements)-1
  cls1 = fltarr(lmax1+1, numchain1, numsim1+1)
  
  for l = 0, lmax1 do begin
     for i = 0, numchain1-1 do begin
        for k = 0, numsim1-1 do begin
           cls1[l,i,k] = cls_in[l,spec,i,k+1]
        end
     end
  end
  
  cls_in = readfits(infile2)
  lmax2     = size(cls_in(*,0,0,0),/n_elements)-1
  numspec2  = size(cls_in(0,*,0,0),/n_elements)
  numchain2 = size(cls_in(0,0,*,0),/n_elements)
  numsim2   = size(cls_in(0,0,0,*),/n_elements)-1
  cls2 = fltarr(lmax2+1, numchain2, numsim2+1)
  
  for l = 0, lmax2 do begin
     for i = 0, numchain2-1 do begin
        for k = 0, numsim2-1 do begin
           cls2[l,i,k] = cls_in[l,spec,i,k+1]
        end
     end
  end
  absmax = 10000.
  
  cls_int1 = dblarr(lmax1+1,numchain1*(numsim1-numexclude))
  cls_int2 = dblarr(lmax2+1,numchain2*(numsim2-numexclude))
  
  numval_in_chain1 = numsim1 - numexclude
  numsamp1 = numval_in_chain1*numchain1
  numval_in_chain2 = numsim2 - numexclude
  numsamp2 = numval_in_chain2*numchain2
  
  j = 0L
  for l = 0L, lmax1 do begin
     for i = 0L, numchain1-1 do begin
        for j = 0L, numval_in_chain1-1L do begin
           cls_int1(l,i*numval_in_chain1+j) = cls1(l,i,1+numexclude+j)
        end
     end
  end

  j = 0L
  for l = 0L, lmax2 do begin
     for i = 0L, numchain2-1 do begin
        for j = 0L, numval_in_chain2-1L do begin
           cls_int2(l,i*numval_in_chain2+j) = cls2(l,i,1+numexclude+j)
        end
     end
  end
  
  
; Compute the confidence regions
  
  vals = fltarr(3)
  hist1 = dblarr(2,numbin)
  hist2 = dblarr(2,numbin)
  close, 38
  close, 39
  close, 40
  openw,38,outfile1
  openw,39,outfile2
  openw,40,outfile_tot
  for l = l_low, l_high, dl do begin
     
     for i = 0, numreg-1 do begin
        
        percentile    = gauss_pdf(i+1) - gauss_pdf(-i-1)
        percentile_90 = 0.9
        maxval = max([max(cls_int1[l,*]), max(cls_int2[l,*])])
        minval = 0.
        hist1[1,*] = histogram(cls_int1[l,*],max=maxval,min=minval,nbins=numbin)
        hist1[0,*] = minval + (findgen(numbin)+0.5) * (maxval-minval)/numbin
;        hist1[1,*] = smooth(hist1[1,*],5)
        hist1[1,*] = hist1[1,*] / total(hist1[1,*])
        
        hist2[1,*] = histogram(cls_int2[l,*],max=maxval,min=minval,nbins=numbin)
        hist2[0,*] = minval + (findgen(numbin)+0.5) * (maxval-minval)/numbin
;        hist2[1,*] = smooth(hist2[1,*],5)
        hist2[1,*] = hist2[1,*] / total(hist2[1,*])
        
        if (l lt l_gauss_low or l gt l_gauss_high) then begin
           vals    = compute_conf_region2(hist1, percentile, minval, maxval)
           vals_90 = compute_conf_region2(hist1, percentile, minval, maxval)
        end else begin
           mu1    = mean(cls_int1[l,*])
           sigma1 = stddev(cls_int1[l,*])
           vals[0] = mu1
           vals[1] = mu1-sigma1
           vals[2] = mu1+sigma1
        end
        
        printf,38,l-5, vals[0], vals[2]-vals[0], vals[0]-vals[1], compute_upper_limit(cls_int1[l,*], 0.9)
        
        if (l lt l_gauss_low or l gt l_gauss_high) then begin
           vals = compute_conf_region2(hist2, percentile, minval, maxval)
        end else begin
           mu2    = mean(cls_int2[l,*])
           sigma2 = stddev(cls_int2[l,*])
           vals[0] = mu2
           vals[1] = mu2-sigma2
           vals[2] = mu2+sigma2
        end
        printf,39,l+5, vals[0], vals[2]-vals[0], vals[0]-vals[1], compute_upper_limit(cls_int2[l,*], 0.9)
        

;           print, 'size a = ', n_elements(hist1[1,*])
;           print, 'size b = ', n_elements(hist2[1,*])
        plot, hist1[0,*],hist1[1,*]
        oplot, hist2[0,*],hist2[1,*]

        hist1[1,*] = smooth(hist1[1,*],10)
        hist2[1,*] = smooth(hist2[1,*],10)
;           return
        
        if (l lt l_gauss_low or l gt l_gauss_high) then begin
           hist1[1,*] = hist1[1,*] * hist2[1,*]
           vals = compute_conf_region2(hist1, percentile, minval, maxval)
           hist1[1,*] = hist1[1,*] / total(hist1[1,*])
           oplot, hist1[0,*],hist1[1,*], linestyle=2
        end else begin
           mu_tot = (1./sigma1^2 * mu1 + 1./sigma2^2 * mu2) / (1./sigma1^2 + 1./sigma2^2)
           sigma_tot = sqrt(1. / (1./sigma1^2 + 1./sigma2^2))
           vals[0] = mu_tot
           vals[1] = mu_tot-sigma_tot
           vals[2] = mu_tot+sigma_tot
        end
        
        printf,40,l, vals[0], vals[2]-vals[0], vals[0]-vals[1], compute_upper_limit_from_dist(hist1, 0.9)
        print, l, vals[0], vals[2]-vals[0], vals[0]-vals[1], compute_upper_limit_from_dist(hist1, 0.9)
        wait, 1
     end
  end
  close, 38
  close, 39
  close, 40
  
end


