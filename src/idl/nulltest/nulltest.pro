pro nulltest, infile, outprefix, lmin, lmax, dl

  cls = readfits(infile)
  
  nspec    = n_elements(cls[0,*,0,0])
  numsamp  = n_elements(cls[0,0,0,*])-1
  numchain = n_elements(cls[0,0,*,0])
  ntot     = numsamp * numchain
  cls = cls[*,*,*,1:numsamp]
  
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
  

  ; Plot chains
  !P.font = 0
  
  mydevice = !D.name
  set_plot, 'PS'
  filename = outprefix+'_chains.eps'
  device, file=filename, $
          xsize = 5, ysize = 1.6*numbin, $ 
          xoffset = 0.6, yoffset = 1.5, /inch ;
  !P.Multi = [0, 2, numbin, 0, 1]
  for i = 0, 1 do begin
     spec = 3+2*i
     for j = 0, numbin-1 do begin
        l = lmin + dl * j
        minrange = min(cls[l,spec,*,*])
        maxrange = max(cls[l,spec,*,*])
        plot, cls[l,spec,0,*], yrange=[minrange,maxrange], xrange=[0,numsamp-1]
        for k = 1, numchain-1 do oplot, cls[l,spec,k,*]
     end
  end
  device, /close
  set_plot, mydevice


  ; Plot histograms
  !P.font = 0
  
  mydevice = !D.name
  set_plot, 'PS'
  filename = outprefix+'_hist.eps'
  device, file=filename, $
          xsize = 5, ysize = 1.2*numbin, $ 
          xoffset = 0.6, yoffset = 1.5, /inch ;
  !P.Multi = [0, 2, numbin, 0, 1]
  n = sqrt(ntot)
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
        h = histogram(cls[l1,spec,*,*], nbins=n, locations=x)
        if (spec eq 3) then begin
           xt = 'EE, l ='+strtrim(string(l1),2)+'-'+strtrim(string(l2),2)
        end else begin
           xt = 'BB, l ='+strtrim(string(l1),2)+'-'+strtrim(string(l2),2)
        end
        plot, x, h, psym=10, xrange=[0, max([1.3*Cl, max(x)])], ytitle=xt, $
              xtitle='Power spectrum, C!Ll!N, (uK!U2!N)'
        if (spec eq 3) then begin
           oplot, [Cl, Cl], [0., 10*max(h)], linestyle=1
        end
     end
  end
  device, /close
  set_plot, mydevice



end
