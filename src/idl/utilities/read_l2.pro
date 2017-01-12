; Routine to read a temporary level 2 file (Oslo format)
; 
; Input:
;    filename -- name of level 2 file
;
; Output
;    modules  -- vector of module numbers
;    time     -- double precision array with timing (mjd)
;    pointing -- 3xN array with pointing in Healpix convention
;    demod    -- 4xN array with demod QUIET data
;
; Hans Kristian Eriksen
; April 15th, 2009
;
pro read_l2, filename, modules, time, pointing, demod

  nmod        = 0L
  num_samples = 0L

  openr,58,filename,/f77

  ; Read number of modules and samples
  readu,58, nmod
  readu,58, num_samples

  ; Allocate data arrays
  time         = dblarr(num_samples)
  pointing     = dblarr(3,num_samples,nmod)
  demod        = dblarr(4,num_samples,nmod)
  modules      = lonarr(nmod)
  p_buffer     = dblarr(3,num_samples)
  demod_buffer = dblarr(4,num_samples)


  ; Read time information
  readu,58,time

  j = 0L
  for i = 0L, nmod-1L do begin

     ; Read module number
     readu,58, j
     modules[i] = j
     
     ; Read pointing
     readu,58, p_buffer
     pointing[*,*,i] = p_buffer

     ; Read demod data
     readu,58, demod_buffer
     demod[*,*,i] = demod_buffer
  end
  close,58


end
