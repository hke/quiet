pro nulltest2, infile1, infile2, outprefix, lmin, lmax, dl, burnin

  cls1 = readfits(infile1)
  nspec1    = n_elements(cls1[0,*,0,0])
  numsamp1  = n_elements(cls1[0,0,0,*])-1 
  numchain1 = n_elements(cls1[0,0,*,0])
  ntot1     = numsamp1 * numchain1
  cls1 = cls1[*,*,*,1+burnin:numsamp1]

  cls2 = readfits(infile2)
  nspec2    = n_elements(cls2[0,*,0,0])
  numsamp2  = n_elements(cls2[0,0,0,*])-1
  numchain2 = n_elements(cls2[0,0,*,0])
  ntot2     = numsamp2 * numchain2
  cls2 = cls2[*,*,*,1+burnin:numsamp2]
  
  numbin = (lmax-lmin)/dl+1

  ; Load WMAP spectrum
  wmapfile = 'wmap_lcdm_sz_lens_wmap5_cl_v3.dat'
  lmax_wmap = 1000
  wmap = fltarr(lmax_wmap+1,3)
  openr,58,wmapfile
  i = 0
  TT = 0.
  TE = 0.
  EE = 0.
  for l = 2, lmax do begin
     readf,58, i, TT, EE, TE
     wmap[l,0] = TT
     wmap[l,1] = TE
     wmap[l,2] = EE
  end
  close,58
  

  ; Plot histograms
  !P.font = 0
  
  mydevice = !D.name
  set_plot, 'PS'
  filename = outprefix+'_hist.eps'
  device, file=filename, $
          xsize = 5, ysize = 1.2*numbin, $ 
          xoffset = 0.6, yoffset = 1.5, /inch ;
  !P.Multi = [0, 2, numbin, 0, 1]
  n = sqrt(ntot1)
  for i = 0, 1 do begin
     spec = 3+2*i
     for j = 0, numbin-1 do begin
        l1 = lmin + dl * j - dl/2 +1
        l2 = lmin + dl * j + dl/2 
        print, l1, l2
        if (spec eq 3) then begin
           Cl = mean(wmap[l1:l2,2])
        end else begin
           Cl = 0.
        end
        h1 = histogram(cls1[l1,spec,*,*], nbins=n, locations=x1)
        h1 = h1 / (total(h1) * (x1(1)-x1(0)))
        h2 = histogram(cls2[l1,spec,*,*], nbins=n, locations=x2)
        h2 = h2 / (total(h2) * (x2(1)-x2(0)))
        if (spec eq 3) then begin
           xt = 'EE, l ='+strtrim(string(l1),2)+'-'+strtrim(string(l2),2)
        end else begin
           xt = 'BB, l ='+strtrim(string(l1),2)+'-'+strtrim(string(l2),2)
        end
        plot, x1, h1, psym=10, xrange=[0, max([1.3*Cl, max(x1)])], yrange=[0,1.1*max([max(h1),max(h2)])], $
              ytitle=xt, xtitle='Power spectrum, C!Ll!N, (uK!U2!N)', thick=2
        oplot, x2, h2, psym=10
        if (spec eq 3) then begin
           oplot, [Cl, Cl], [0., 10*max(h1), 10*max(h2)], linestyle=1
        end
     end
  end
  device, /close
  set_plot, mydevice



end
