pro RIPS_Mercury_AEOS_v03, PART=part, NIGHT=night

; *********************************************************************************************************************
; *********************************************************************************************************************
; Quick routine to extract Na emission from spectral scans over Mercury's disk with AEOS/RIPS in order to create a 2-D
; image of Na above and near Mercury's disk.
; 
; The program is split into 4 "parts":
;   0 = break kinetic series into manageable datacube
;       input:   kinetic series FITS file(s)
;       output:  RIPS imaging (x,y,t) and spectral (wl,y,t) datacubes 
;   1 = find Mercury centroids by cross-corellation with previous image
;       input:   imaging datacube from Part 0
;       output:  co-aligned imaging datacube (x,y,t)
;   2 = isolate sodium emission in every frame
;       input:   spectral datacube from Part 0
;       output:  the calibrated spectral datacube (wl,y,t) for Na exosphere
;   3 = extract Na signal and place into 1D spectra
;       input:   spectral datacube from Part 2 (wl,y,t)
;       output:  Na brightness (wl,y) and linewidth (wl,y)
;   4 = put everything together and build exosphere image
;       input:   Na brightness and linewidth from Part 3 and imaging cube from Part 0
;       output:  none yet, a combined image called "test" for now...
;  99 = do all of 0-4
;
; Input variables
;   part           - the "part" of the program to execute
;   night          - the night of AEOS data to use; either '20' or '25' for June 20 or 25 (Default='25'); 
;                    contains numerous sub-variables: mercury_dir, dark_dir, flat_dir, sky_dir, Mercury_file, 
;                    sky_file, dark_file, flat_file, Na_D_rngs, ims, sps
;   thresh         - hotpixel removal threshold (sigma above background; Default=12)
;   width          - median smoothing width to get local background reference (for pixels above thresh; Default=5)
;
; 
; Created:  05 Jul 2018 - CS
; Edits:    09 Jul 2018 - "v02"; improved moon drift calibration in Part 2; added call to routine for 
;                         determining spatial and spectral platescale; include crude "lucky imaging" treatment in 
;                         Part 1 (LM)
;           15 Sep 2018 - "V03" edited to give a spatial aperture is probablistic based on the seeing,
;                         co-adding of frames should be done as a weighted mean of the probabilties (CS)                         
; *********************************************************************************************************************
; *********************************************************************************************************************
; =====================================================================================================================
; Define variables
; =====================================================================================================================
SetDefaultValue, part, 4                                                  ; Which part we doin?
SetDefaultValue, night, '25'                                              ; set to '20' or '25' for night of June 20 or June 25
SetDefaultValue, thresh, 12                                               ; hotpixel removal threshold (sigma above background)
SetDefaultValue, width, 5                                                 ; Median smoothing width to get local background reference (for pixels above thresh) 
SetDefaultValue, minfac, 0.8                                              ; maximum factor to apply to reference spectrum for Na emission extraction
SetDefaultValue, maxfac, 1.2                                              ; minimum factor to apply to reference spectrum for Na emission extraction
SetDefaultValue, dfac, 0.01                                               ; factor increment for minfac->maxfac
SetDefaultValue, ct, 20;27                                                   ; default color table
SetDefaultValue, max_spec, 2950.                                          ; partial hack --> the max value of the spectral cube after coadding (used for "prettier" movies)
SetDefaultValue, max_img, 1.6e6;2.7e6                                           ; partial hack --> the max value of the image cube after coadding (used for "prettier" movies)
SetDefaultValue, img_extraction, [310,470,105,265]                        ; coordinates for extracting a square Mercury image, [x1,x2,y1,y2]
SetDefaultValue, spec_extraction, [315,475,120,280]                       ; coordinates for extracting a square Mercury spectral image, [x1,x2,y1,y2]
SetDefaultValue, movie_scale, 4                                           ; scale the x-y dimensions of the extracted "Na image" and "Mercury disk" arrays for movies by this value
SetDefaultValue, do_realign, 0                                            ; 0=run part 1 as normal, 1=use previous cross-correlation as master template
SetDefaultValue, stddev_cutoff, 850.                                      ; only use Mercury frames with stddev > stddev_cutoff, these are defined as "good" 
SetDefaultValue, no_display, 1                                            ; no display = 1 skips the slow graphics, 0 shows everything
SetDefaultValue, root_dir, 'C:\Users\schmidt\Desktop\RIPS\2018 - AEOS with RIPS\'    ; directory where data live
SetDefaultValue, Spatial_sigma_width, 4                                   ; sigma width of the spatial dimemsion weighting where the spectral slice will be placed
debug = 0

if night eq '20' then begin  
  ; ===== These variables for the Mercury scans from 20 June HST  ===== 
  mercury_dir = root_dir+'June 20 HST\'    ; directory with Mercury FITS files
  dark_dir    = root_dir+'June 21 HST\'    ; directory with dark file
  flat_dir    = root_dir+'June 21 HST\'    ; directory with flat file
  sky_dir     = root_dir+'June 21 HST\'    ; directory with sky file
  Mercury_file = 'Mercury - 500 kinetic 1 sec scan - Na 3A + ND imaging - Na order (2).fits' ; Mercury kinetic series frames
  sky_file    = 'Sky flats - Na 3A + ND imaging - Na order.fits'          ; sky frame
  dark_file   = 'Dark - 200 sec - Na 3A + ND imaging - Na order.fits'     ; dark frame
  flat_file   = 'Tungsten lamp flat - 200 sec - Na 3A + ND imaging - Na order_Autosave_recovered.fits'  ; flat frame
  Na_D_rngs   = [266,297,528,538]                                         ; ranges of Na D Mercury emission (x pixels in spec domain for D1 and D2)
  ims         = [183,809,39,444]                                          ; variables for image analysis: x1,x2,y1,y2 (formerly "imaging_statsec")
  sps         = [99,990,502,978]                                          ; variables for spectral analysis: x1,x2,y1,y2 (formerly "spectra_statsec")
endif else if night eq '25' then begin             
  ; ===== These variables for the Mercury scans from 25 June HST  ===== 
  Mercury_dir = root_dir+'June 25 HST\'    ; directory with Mercury FITS files
  dark_dir    = root_dir+'June 25 HST\'    ; directory with dark file
  flat_dir    = root_dir+'June 25 HST\'    ; directory with flat file
  sky_dir     = root_dir+'June 23 HST\'    ; directory with sky/solar file
;  Mercury_file = 'RIPS_setup_215.fits'                                    ; Mercury kinetic series frames (210-215)
  Mercury_file = 'RIPS_setup_' + strcompress((indgen(6)+210), /remove_all) + '.fits'         ; Mercury kinetic series frames (210-215)
  sky_file    = 'RIPS_setup_174.fits'                                     ; sky/solar frame (164=sky Na, 174=lunar drift Na)
  dark_file   = 'RIPS_setup_203.fits'                                     ; dark frame
  flat_file   = 'RIPS_setup_205.fits'                                     ; flat frame
  
  Na_D_rngs   = [230,250,474,494]                                         ; ranges of Na D Mercury emission (x pixels in spec domain for D1 and D2)                       
  ims         = [183,809,39,414]                                          ; variables for image analysis: x1,x2,y1,y2 (formerly "imaging_statsec")
  sps         = [99,990,502,951]                                          ; variables for spectral analysis: x1,x2,y1,y2 (formerly "spectra_statsec")
  iref_factor = 1.                                                        ; factor to scale reference spectrum by
  integration = 1.                                                        ; HACK (in that we should just look up the integration time from the header... but I'm lazy)
endif else stop

outdir        =  Mercury_dir + 'Processed\'                               ; output directory
nfiles        =  n_elements(Mercury_file)                                 ; # of different Mercury kinetic frames

!p.multi = 0.
cgloadct, ct;, /reverse

; =====================================================================================================================
; Find display size, set windows accordingly
; =====================================================================================================================
DEVICE, GET_SCREEN_SIZE = ss ;& PRINT, ss                                  ; find screen size
;winpos_x      = [ss(0)/2, ss(0)/2, ss(0)/2+ss(0)/4, ss(0)/2+ss(0)/4, ss(0)] ; default x window locations
winpos_x      = [0, ss(0)/2, 0, ss(0)/2+ss(0)/4, ss(0)] ; default x window locations
winpos_y      = [0, ss(1)/2, 0, ss(1)/2, 0]                                ; default y window locations

; =====================================================================================================================
; Part 0 : break kinetic series into manageable datacube
; =====================================================================================================================
if part eq 0 or part eq 99 then begin 
    icube = MRDFITS(Mercury_dir + Mercury_file(0), 0, header, /unsigned, /silent )
    s = size(icube, /dimensions) 
    cube = fltarr(s(0),s(1),s(2)*nfiles)
    for ifile = 0, nfiles-1 do begin
      icube = MRDFITS(Mercury_dir + Mercury_file(ifile), 0, header, /unsigned, /silent )
      cube(*,*,s(2)*ifile:s(2)*ifile+99) = icube
    endfor
    
    imaging_cube = reform(cube[ims(0):ims(1),ims(2):ims(3), *])
    spectra_cube = reform(cube[sps(0):sps(1),sps(2):sps(3), *])
    MWRFITS, imaging_cube, outdir + 'imaging_cube.fits', header, /create, /silent
    MWRFITS, spectra_cube, outdir + 'spectra_cube.fits', header, /create, /silent
    Sky_cube = MRDFITS(sky_dir + sky_file, 0, header, /unsigned, /silent )
    imaging_sky_cube = reform(sky_cube[ims(0):ims(1),ims(2):ims(3), *])
    spectra_sky_cube = reform(sky_cube[sps(0):sps(1),sps(2):sps(3), *])
    MWRFITS, imaging_sky_cube, outdir + 'imaging_sky_cube.fits', header, /create, /silent
    MWRFITS, spectra_sky_cube, outdir + 'spectra_sky_cube.fits', header, /create, /silent
    beep
endif

; =====================================================================================================================
; Part 1 : find Mercury centroids by cross-corellation with previous image
; =====================================================================================================================
if part eq 1 or part eq 99 then begin ;Find the centroids by cross-correlation with the last image
    imaging_cube = MRDFITS(outdir + 'imaging_cube.fits', 0, header, /fscale, /silent )
    frame = reform(imaging_cube[*,*,0])
    s = size(imaging_cube, /dimensions) 
    aligned_imaging_cube = fltarr(s) 
;    window, 0, xs = s[0], ys = s[1]
    window, 0, xpos=winpos_x[0], ypos=winpos_y[0], xs=s[0], ys=s[1]

    Dark = MRDFITS(dark_dir + dark_file, 0, header, /fscale, /silent )
    dark = dark[ims(0):ims(1),ims(2):ims(3)]
    dark = sigma_filter( dark, 5, N_sigma=5)
    flat = MRDFITS(flat_dir + flat_file, 0, header, /fscale ) 
    flat = (flat[ims(0):ims(1),ims(2):ims(3)] - dark) / mean(flat[ims(0):ims(1),ims(2):ims(3)]) 
    flat[WHERE(flat lt .01, /NULL)] = !values.f_Nan ;reject unusual counts for centroid
    flat[WHERE(flat gt 2.0, /NULL)] = !values.f_Nan ;reject unusual counts for centroid
;    flat = smart_shift(flat, 1.25, .5,  /interp) ;hack
    flat = smart_shift(flat, .8, .5,  /interp) ;hack
    gooddata = where(Finite(flat), ngooddata, comp=baddata, ncomp=nbaddata)
    if nbaddata gt 0 then flat[baddata] = interpol(flat[gooddata], gooddata, baddata)

    if do_realign then begin
      restore, outdir+'shift_array.sav'          ; contains shift_array, std_devs, and aligned_imaging_cube
      shift_array_old = shift_array              ; keep old shift_array for comparison
      reference = total(aligned_imaging_cube,3)  ; create "master Mercury" template for cross-correlation
    endif else begin
      ireference = (reform(imaging_cube[*,*,s[2]-1]) - dark) / flat
      acre, ireference, reference, thresh, width
    endelse
    shift_array = intarr(s[2],2)
    std_devs = fltarr(s[2])
    for i = 0, s[2]-1 do begin
        iframe = (reform(imaging_cube[*,*,i]) - dark) / flat
        
        ;--------------hack to get interpol right--------
        iframe[337:343,*] = !values.F_Nan 
        tv, bytscl(iframe)
        stop
        FILL_MISSING, iframe, !values.F_Nan, 1
        tv, bytscl(iframe)
        ;------------------------------------------------
        
        stop
        acre, iframe, frame, thresh, width
        CORREL_OPTIMIZE, reference, frame, xoffset_optimum, yoffset_optimum, NUMPIX = 150
        shift_array[i,*] = [xoffset_optimum, yoffset_optimum]
    endfor
    if night eq '20' then begin
      shift_array[21:25,0] = -221                            ;hack
      shift_array = smooth(shift_array, [5,1], /EDGE_MIRROR) ;hack
    endif
    for i = 0, s[2]-1 do begin
        frame = (reform(imaging_cube[*,*,i]) - dark) / flat ;make this into a movie!
        aligned_imaging_cube[*,*,i] = shift(frame, [shift_array[i,*]])
        tv, bytscl(shift(frame, [shift_array[i,*]]), 0, 5000)
        img_frame = reform(frame(img_extraction(0)-50:img_extraction(1)+50,img_extraction(2)-50:img_extraction(3)+50))   ; extract just the image part of the frame
        top_rows = total( img_frame(*,n_elements(img_frame(0,*))-10:*),2 )
        slit_pos = (where(top_rows eq max(top_rows)))[0]                  ; find position of slit
        frame(slit_pos-5:slit_pos+5,*) = 0.0                              ; zero out slit position values so they don't interfere with determining image "sharpness"
        std_devs(i) = stddev(img_frame)                                   ; seems like anything < 1000 is not that good
;        print, i, n_elements(frame) - n_elements( where( frame lt (max(frame)*.8)/2d ) )
;        wait, 0.05
;        stop
    endfor    
    save, shift_array, aligned_imaging_cube, std_devs, filename = outdir + 'shift_array.sav'
    beep
    stop
;    x = reform(aligned_imaging_cube(img_extraction(0)-50:img_extraction(1)+50,img_extraction(2)-50:img_extraction(3)+50,i))
endif

; =====================================================================================================================
; Part 2 : Isolate the sodium emission in every frame
; =====================================================================================================================
if part eq 2 or part eq 99 then begin 
    spectra_cube = MRDFITS(outdir + 'spectra_cube.fits', 0, header, /fscale, /silent )
    frame = reform(spectra_cube[*,*,0])
    s = size(spectra_cube, /dimensions) 
    window, 0, xpos=winpos_x[0], ypos=winpos_y[0], xs=s[0], ys=s[1]
    tv, bytscl(Frame)
    
    ;For this data there is no dark or bias calibration frames
    ;Therefore, use the conventional mode date from the following night
    
    Dark = MRDFITS(dark_dir + dark_file, 0, header, /fscale, /silent )
    dark = dark[sps(0):sps(1),sps(2):sps(3)]
    dark = sigma_filter( dark, 5, N_sigma=5)
    acre, dark, dark, thresh, width                 ; try acre, still hot pixels in dark

;    xi = MRDFITS( outdir + 'imaging_sky_cube.fits', 0, header, /fscale, /silent )
;    xs = MRDFITS( outdir + 'spectra_sky_cube.fits', 0, header, /fscale, /silent )
;    stop
    if night eq '25' then begin
      spectra_sky_cube = MRDFITS(outdir + 'spectra_sky_cube.fits', 0, header, /fscale, /silent )
      ireference = spectra_sky_cube - dark
      acre, ireference, reference, thresh, width
      window, 1, xpos=winpos_x[1], ypos=winpos_y[1], xs=s[0], ys=s[1], title='IDL 1 - REFERENCE'
      tv, bytscl(reference)
      ref_spectrum = total(reference[*,120:220], 2)
;      amplifer_correction = mean(reference[*,s(1)-51:s(1)-1], dimension = 2)
;      amplifier_contribution = rebin(amplifer_correction, s[0], s[1])
;      dark = dark + amplifier_contribution                                ; correct for the fact that the dark correction from conventional mode is clearly insufficient for EMCCD mode

      ;get a spectral flat
      iflat = MRDFITS(flat_dir + flat_file, 0, header, /fscale ) 
      iflat = (iflat[sps(0):sps(1),sps(2):sps(3)] - dark)
      iflat = iflat / median(iflat) 
      acre, iflat, flat, thresh, width
      flat = smart_shift(flat, 0., -.65,  /interp) ;hack
      tv, bytscl(flat)

      reference = reference / flat
      tv, bytscl(reference) 
      ;ref_spectrum = total(reference[*,150:250], 2)
      ref_spectrum = total(reference, 2)
      
      ; Now find the xoffset in the reference spectrum (likely due to slightly different grating angles)
      y_Mercury = where( total(frame,1) ge 0.5*max(total(frame,1)) )    ; FWHM range of mercury
      frame_spectrum = total(frame(*,y_Mercury),2)
      frame_D1_x = ( where( frame_spectrum eq min(frame_spectrum) ) ) [0]
      ref_D1_x   = ( where( ref_spectrum   eq min(ref_spectrum)   ) ) [0]
      reference = shift(reference, frame_D1_x - ref_D1_x)
      ref_spectrum = shift(ref_spectrum, frame_D1_x - ref_D1_x)

      exosphere_spectra_cube = fltarr( size(spectra_cube, /dimensions) )
 
       ;now scale and subtract 
      Device, Window_State=theseWindows
      if theseWindows(2) ne 1 then window, 2, xpos=winpos_x[0]+s[0]+15, ypos=winpos_y[2], xs = s[0], ys = s[1], title='IDL 2 - FRAME'
      if theseWindows(3) ne 1 then window, 3, xpos=winpos_x[0]+s[0]+15, ypos=winpos_y[3], xs = s[0], ys = s[1], title='IDL 3 - SCALED REFERENCE'
      if theseWindows(4) ne 1 then window, 4, xpos=winpos_x[0]+s[0]*2+30, ypos=winpos_y[4], xs = s[0], ys = s[1], title='IDL 4 - EXOSPHERE'
      for i = 0, s[2]-1 do begin                                          ; LOOP over each frame in the series
        frame = reform(spectra_cube[*,*,i]) - dark                        ;dark subtract the spectrum
        frame = frame / flat                                              ;flatfield the spectrum
      
        ; Now find best scaling for ref_spectrum and subtract
;        y_Mercury = where( total(frame,1) ge 0.5*max(total(frame,1)) )    ; FWHM range of mercury
;        Mercury_spectrum = total(frame(*,y_Mercury),2)
        
        illum = frame / reference                                         ; scale the illumination against Mercury's disk
        for iD = 0, 2, 2 do illum[Na_D_rngs(iD+0):Na_D_rngs(iD+1),*] = !values.F_NaN  ; carefully avoid sodium emissions when scaling to the disk
      
;        illum_along_slit = mean(illum, dimension=1, /nan)
        illum_along_slit = MEDIAN(illum, dimension=1)
        illum_along_slit = MEDSMOOTH( illum_along_slit, 5 )               ;helps with hot pixels from bad flat-fielding where the slit has some dust
        scaled_reference = reference * rebin(transpose(illum_along_slit), s[0], s[1])
        Just_exosphere = Frame - scaled_reference * iref_factor
        
        ; Check that we've done a good job of minimizing spectrum between Na D lines
        between_D = smooth( just_exosphere(Na_D_rngs(1)+20:Na_D_rngs(2)-20,*), 6 )
        mean_vals = abs([ mean(between_D(*,min(Y_Mercury):max(Y_Mercury))), $ ; mean value across Mercury's disk
                          mean(between_D(*,0:min(Y_Mercury))), $              ; mean value below mercury's disk
                          mean(between_D(*,max(Y_Mercury):*)) ])              ; mean value above Mercury's disk
        
        if mean_vals(0) gt 2.*mean(mean_vals(1:2)) then begin             ; find the best "ref_factor" if we haven't
          nfacs = (maxfac - minfac)/dfac + 1                              ; # of factors to consider
          factors = cgscalevector(indgen(nfacs),minfac,maxfac)
          fac_diffs = fltarr(nfacs)                                       ; keep track of the differences in mean_vals for each factor
          for ifac = 0, nfacs-1 do begin
            Just_exosphere = frame - scaled_reference*factors(ifac)
            between_D = smooth( just_exosphere(Na_D_rngs(1)+20:Na_D_rngs(2)-20,*), 6 )
            mean_vals = abs([ mean(between_D(*,min(Y_Mercury):max(Y_Mercury))), $ ; mean value across Mercury's disk
                              mean(between_D(*,0:min(Y_Mercury))), $              ; mean value below mercury's disk
                              mean(between_D(*,max(Y_Mercury):*)) ])              ; mean value above Mercury's disk
            fac_diffs(ifac) = abs( mean_vals(0) - mean(mean_vals(1:2)) )            
          endfor
          best_fac = factors(( where( fac_diffs eq min(fac_diffs) ) )[0])
          print, i, best_fac
          just_exosphere = frame - scaled_reference * best_fac
        endif
;        stop
        
        
;        iminimize = 0
;        while mean_vals(0) gt 2.*mean(mean_vals(1:2)) do begin
;          illum = just_exosphere / reference                              ; scale the illumination against Mercury's disk
;          for iD = 0, 2, 2 do illum[Na_D_rngs(iD+0):Na_D_rngs(iD+1),*] = !values.F_NaN
;          illum_along_slit = MEDIAN(illum, dimension=1)
;          illum_along_slit = MEDSMOOTH( illum_along_slit, 5 ) 
;          scaled_reference = reference * rebin(transpose(illum_along_slit), s[0], s[1])
;          Just_exosphere = just_exosphere - scaled_reference * ref_factor 
;          between_D = smooth( just_exosphere(Na_D_rngs(1):Na_D_rngs(2),*), 5 )
;          mean_vals = [ mean(between_D(*,min(Y_Mercury):max(Y_Mercury))), $ ; mean value across Mercury's disk
;                        mean(between_D(*,0:min(Y_Mercury))), $              ; mean value below mercury's disk
;                        mean(between_D(*,max(Y_Mercury):*)) ]               ; mean value above Mercury's disk
;          iminimize++
;          tv, bytscl(just_exosphere,0,200)
;          wait, 0.25
;          if iminimize gt 10 then mean_vals(*) = 1.d
;          print, iminimize, mean(just_exosphere)
;          if iminimize gt 20 then stop
;        endwhile
;        stop
        
        
        wset, 2
        cgimage, bytscl(frame, 0, 200)
        wset, 3
        cgimage, bytscl(scaled_reference, 0, 200)
        wset, 4
        cgimage, bytscl(Just_exosphere, 0, 200)
        wait, 0.02
        exosphere_spectra_cube[ *, *, i] = Just_exosphere
        ;stop
      endfor ;i
      save, exosphere_spectra_cube, filename = outdir + 'exosphere_spectra_cube.sav'
    endif ;night=='25'


    if night eq '20' then begin
      ;For this data there is no reference spectrum without sodium
      ;The sky spectrum won't work, the lunar spectrum the following night has a wider slit
      ;Therefore, let's hack our way through using Mercury itself.
    
      reference = mean(spectra_cube, dimension = 3) - (dark)
      window, 1, xs = s[0], ys = s[1], xpos=winpos_x[1], ypos=winpos_y[1]
      tv, bytscl(reference)
      ref_spectrum = total(reference[*,120:220], 2)    
      amplifer_correction = mean(reference[*,426:476], dimension = 2)
      amplifier_contribution = rebin(amplifer_correction, s[0], s[1])
      dark = dark + amplifier_contribution                                ; correct for the fact that the dark correction from conventional mode is clearly insufficient for EMCCD mode
    
      ;get a spectral flat
      flat = MRDFITS(flat_dir + flat_file, 0, header, /fscale ) 
      flat = (flat[sps(0):sps(1),sps(2):sps(3)] - dark) 
      flat = flat / median(flat) 
    
      reference = (mean(spectra_cube, dimension = 3) - (dark)) / flat
      window, 1, xs = s[0], ys = s[1], xpos=winpos_x[1], ypos=winpos_y[1]
      tv, bytscl(reference, 0, 100) 
      ref_spectrum = total(reference[*,120:220], 2)    
      window, 2, xs = s[0], ys = s[1], xpos=winpos_x[2], ypos=winpos_y[2]
      cgplot, ref_spectrum, color = 'black', xr = [250,600]
      dummy = ref_spectrum
      dummy[286:297] = !values.f_Nan                                      ;remove Mercury D2 emission
      dummy[528:538] = !values.f_Nan                                      ;reject unusual counts for centroid
      gooddata = where(Finite(dummy), ngooddata, comp=baddata, ncomp=nbaddata)
      if nbaddata gt 0 then dummy[baddata] = interpol(dummy[gooddata], gooddata, baddata)
      cgplot, dummy, color = 'red', /overplot
      ref_spectrum = rebin(dummy, s[0], s[1])                             ; this is what a flat-fielded lunar drift scan reference spectrum might have looked like of they had been given the time to take one with this night  
    
      ;now scale and subtract 
      Device, Window_State=theseWindows
      if theseWindows(2) ne 1 then window, 2, xs = s[0], ys = s[1], xpos=winpos_x[2], ypos=winpos_y[2]
      if theseWindows(3) ne 1 then window, 3, xs = s[0], ys = s[1], xpos=winpos_x[3], ypos=winpos_y[3]
      if theseWindows(4) ne 1 then window, 4, xs = s[0], ys = s[1], xpos=winpos_x[4], ypos=winpos_y[4]
      stop
      exosphere_spectra_cube = fltarr( size(spectra_cube, /dimensions) )
      for i = 0, s[2]-1 do begin
        frame = reform(spectra_cube[*,*,i]) - dark                        ;dark subtract the spectrum
        frame = frame / flat                                              ;flatfield the spectrum
      
        illum = frame / ref_spectrum                                      ;scale the illumination against Mercury's disk
        for iD = 0, 2, 2 do illum[Na_D_rngs(iD+0):Na_D_rngs(iD+1),*] = !values.F_NaN  ; carefully avoid sodium emissions when scaling to the disk
;        illum[528:538, *] = !values.f_Nan                                 ;carefully avoid sodium emissions when scaling to the disk
;        illum[286:297, *] = !values.f_Nan                                 ;carefully avoid sodium emissions when scaling to the disk
      
        illum_along_slit = mean(illum, dimension=1, /nan)
        illum_along_slit = MEDSMOOTH( illum_along_slit, 5 )               ;helps with hot pixels from bad flat-fielding where the slit has some dust
        This_Frames_Fake_Mercury_Spectrum = ref_spectrum * rebin(transpose(illum_along_slit), s[0], s[1])
        Just_exosphere = Frame - This_Frames_Fake_Mercury_Spectrum
        wset, 2
        tv, bytscl(frame, 0, 200)
        wset, 3
        tv, bytscl(This_Frames_Fake_Mercury_Spectrum, 0, 200)
        wset, 4
        tv, bytscl(Just_exosphere, 0, 200)
        wait, 0.02
        exosphere_spectra_cube[ *, *, i] = Just_exosphere
        stop
      endfor
      save, exosphere_spectra_cube, filename = outdir + 'exosphere_spectra_cube.sav'
    endif ;night=='20'
    beep
;   stop
endif

; =====================================================================================================================
; Part 3 : extract and place the 1D spectra.
; =====================================================================================================================
if part eq 3 or part eq 99 then begin 
  restore, outdir + 'exosphere_spectra_cube.sav'
  s = size(exosphere_spectra_cube, /dimensions) 
  
  ;Dispersion calculation: at y of 265, lines occur at 390 and 632 in x 
  dispersion = (5895.92424 - 5889.95095) / (632. - 390.) ;Angstroms per pixel
  dispersion_Velocity = 299792.*dispersion / 5889.95095  ;Km/s per pixel at Na D2

   ;-------------------------------------------Extraction----------------------------------------------------------- 

    ;run MPFIT once to find the D2 line center
      img = total(exosphere_spectra_cube, 3)
      search = 16           ;pixels to search over
      expected_pixel = mean(Na_D_rngs(0:1))  ;Rough D2 pixel location
      D2_height = fltarr(n_elements(img[0,*]))
      D2_center = fltarr(n_elements(img[0,*]))
      D2_width  = fltarr(n_elements(img[0,*]))
      for i = 0, n_elements(img[0,*]) - 1 do begin
        result = mpfitpeak(findgen(search*2. + 1), img[expected_pixel-search:expected_pixel+search, i], a, STATUS = STATUS, /positive) 
        if STATUS ge 1 then D2_height[i] = A[0] else D2_height[i] = !values.F_nan
        if STATUS ge 1 then D2_center[i] = A[1] else D2_center[i] = !values.F_nan
        if STATUS ge 1 then D2_width[i] = A[2] else D2_width[i] = !values.F_nan
      endfor
      y = findgen(n_elements(img[0,*]))
      height = smooth(D2_height, 5, /edge_mirror) / s[2]
      height[0:25] = 1. & height[280:*] = 1.
      real = where(finite(D2_center), /NULL)
      COEFF = ROBUST_POLY_FIT(y[real], D2_center[real], 2)
      location = poly(findgen(n_elements(img[0,*])), coeff) + expected_pixel-search 
      real = where(finite(D2_width), /NULL)
      COEFF = ROBUST_POLY_FIT(y[real], D2_width[real], 2)
      width = poly(findgen(n_elements(img[0,*])), coeff) 
      dummy = img
      dummy[location, y] = max(dummy) 
      tv, bytscl(dummy, -100, 200)

    ;setup MPFIT Gaussian parameters
      guess = [80.,location[0],1.5,0.]
      A = guess
      parinfo = replicate( {value:0.D, fixed: 0b, limited: [0b,0b], limits: dblarr(2) }, 3)
      parinfo[0].limited = 1b                                 ;limit amplitude 
      parinfo[0].limits = [0.1, 1000. ]                       ;positive only
      parinfo[1].limited = 1b                                 ;limit sigma   
      ;parinfo[1].limits = [14., 18.]                          ;limit sigma width in pixels      
      parinfo[1].limits = [10., 22.]                          ;limit sigma width in pixels
      parinfo[2].limited = 1b                                 ;limit sigma   
      parinfo[2].limits = [1., 5.]                            ;limit sigma width in pixels
  
  if debug then window, 1, xs = 892, ys=450
  brightness = fltarr(s[2], s[1]) & err_brightness = fltarr(s[2], s[1])
  linewidth = fltarr(s[2], s[1]) & err_linewidth = fltarr(s[2], s[1])
  for n = 0, s[2]-1 do begin
      img = exosphere_spectra_cube[*,*,n]
      err_img = sqrt( abs(img) )     
      img = smooth(img,[1,5], /Nan)                   ;hack for now smoothing needed but not optimized
      err_img = smooth(err_img,[1,5], /Nan) / sqrt(5) ;hack for now smoothing needed but not optimized
      for i = 0, s[1]-1 do begin
               a[0] = height[i]
               a[1] = location[i] - expected_pixel + search
               a[2] = Width[i]
               if debug then begin
                  wset, 1
                  dummy = img
                  dummy[location[i]-search:location[i]+search, i] = max(dummy, /nan)
                  tv, bytscl(dummy)
                  wset, 0
                  cgplot, findgen(search*2. + 1), img[location[i]-search:location[i]+search, i]
               endif 
               result = mpfitpeak(findgen(search*2. + 1), img[location[i]-search:location[i]+search, i], A, PERROR = Err_A, $
                                  error = err_img[location[i]-search:location[i]+search, i], /positive, /NAN, STATUS = STATUS, PARINFO = PARINFO, NTERMS =3)
               if debug then begin
                  cgoplot, findgen(search*2. + 1), result, color = 'red' 
                  wait, 0.02
                  if i eq 200 then stop
               endif
               if STATUS ge 1 then begin
                  brightness[n,i] = A[0]*A[2]*SQRT(2*!DPI) 
                  err_brightness[n,i] = brightness[n,i] * sqrt( (err_A[0]/A[0])^2. + (err_A[2]/A[2])^2. )
                  linewidth[n,i] = dispersion_Velocity*sqrt((2.355*A[2])^2 - (2.45)^2) ;FWHM KM/S instrumental measured width is 2.45 pixels 
                  err_linewidth[n,i] = dispersion_Velocity*2.355*Err_A[2]
               endif else begin
                  brightness[n,i] = !values.F_nan
                  linewidth[n,i] = !values.F_nan
               endelse 
      endfor
      if not debug then print, 'Finished extracting:', n
  endfor  
  tv, bytscl(brightness, 0, 1000)
  save, brightness, linewidth, filename = outdir + 'brightness.sav'
  beep
endif  

; =====================================================================================================================
; Part 4 : put everything together and build images and movies of the exosphere.
; =====================================================================================================================
if part eq 4 or part eq 99 then begin  
   restore, outdir + 'brightness.sav'            ; contains brightness, linewidth
   restore, outdir + 'shift_array.sav'           ; contains shift_array, std_devs, and aligned_imaging_cube
   imaging_cube = MRDFITS(outdir + 'imaging_cube.fits', 0, header, /fscale, /silent )
   s = size(imaging_cube, /dimensions)
   
   ;brightness = congrid(brightness, s[2], s[1]) ;HACK NEED PROPER SCALING OF THE SPECTRUM TO IMAGING PLATESCALES
   brightness = congrid(brightness, s[2], s[1]) ;Jeff says the platescale ratio, based on lense focal lengths is   try that!!!
;stop
   home = [341,238] ;Mercury centroid of the last frame number 500, to which all the other frames are aligned
   
   big_cube = fltarr(s[0], 2*s[1], s[2])
   big_weighting_cube = fltarr(s[0], 2*s[1], s[2])
   img_cube = fltarr(s[0], s[1])
   spec_cube = fltarr(s[0], s[1])
   spec_dwell = fltarr(s[0], s[1])                                        ; keep track of the # of frames in each pixel
   spec_temp  = fltarr(s[0], s[1]) 
   big_frame_total = fltarr(s[0], 2*s[1])
   cgdisplay, s[0], s[1]*2, wid=0;, xpos=winpos_x[0], ypos=winpos_y[0]

   xs_spec = spec_extraction(1) - spec_extraction(0) + 1
   ys_spec = spec_extraction(3) - spec_extraction(2) + 1
   xs_img= img_extraction(1) - img_extraction(0) + 1
   ys_img = img_extraction(3) - img_extraction(2) + 1  
   cgdisplay, xs_spec*movie_scale, ys_spec*movie_scale, xpos=winpos_x[2], ypos=winpos_y[2], wid=1
   window, 2, xs = xs_img*movie_scale, ys = ys_img*movie_scale;, xpos=winpos_x[2], ypos=winpos_y[2]+ys_img*movie_scale+35
   window, 4, xs = s[0], ys = 2*s[1];, xpos=winpos_x[4], ypos=winpos_y[4]

   ; Some video related variables
   mpgFilename = outdir + 'Mercury Na image.mp4'
   mpgFilename2 = outdir + 'Mercury disk.mp4'
   mpgFilename3 = outdir + 'Mercury scan.mp4'
   video = IDLffVideoWrite(mpgFilename, format='mp4')
   video2 = IDLffVideoWrite(mpgFilename2, format='mp4')
   video3 = IDLffVideoWrite(mpgFilename3, format='mp4')
   framerate = 30
   framedims = [xs_spec*movie_scale, ys_spec*movie_scale]
   framedims2 = [xs_img*movie_scale, ys_img*movie_scale]
   framedims3 = [s[0], 2*s[1]]
   stream = video.AddVideoStream(framedims[0], framedims[1], framerate)
   stream2 = video2.AddVideoStream(framedims2[0], framedims2[1], framerate)
   stream3 = video3.AddVideoStream(framedims3[0], framedims3[1], framerate)

   good_frames = where(std_devs gt stddev_cutoff, complement=bad_frames)  ; these are the good frames we'll use below
   for ig = 0, n_elements(good_frames)-1 do begin                         ; loop over good frames
        i = good_frames[ig]
        big_frame = fltarr(s[0], 2*s[1])
        big_weighting_frame = fltarr(s[0], 2*s[1])
        empty_frame = make_array(s[0], s[1], /float, value = !values.F_NaN)
        empty_weighting_frame = make_array(s[0], s[1], /float, value = !values.F_NaN)
        
        ;assign this to 3 spatial pixels---unweighted
        empty_frame[home[0]+shift_array[i,0]-1:home[0]+shift_array[i,0]+1,*] = [brightness[i, *], brightness[i, *], brightness[i, *]] ;Make it bigger?

        ;assign this to a gaussian of spatial pixels---weighted
        weighting = GAUSSIAN_FUNCTION(Spatial_sigma_width, /NORMALIZE)
        empty_weighting_frame[home[0] + shift_array[i,0] - (N_elements(weighting)-1) / 2: home[0] + shift_array[i,0] + (N_elements(weighting)-1) / 2, *] $
                             = rebin(weighting, n_elements(weighting), s[1])

        ;big_frame[*, s[1]:*] = rotate(empty_frame, 7)                        ;place upper
        big_frame[*, s[1]:*] = rotate(empty_frame, 6)                         ;place upper
        ;big_weighting_frame[*, s[1]:*] = rotate(empty_weighting_frame, 7)    ;place upper weighting 
        big_weighting_frame[*, s[1]:*] = rotate(empty_weighting_frame, 6)     ;place upper weighting        
        big_frame[*, 0:s[1]-1] = reform(aligned_imaging_cube[*,*,i])          ;place lower
        big_weighting_frame[*, 0:s[1]-1] = 1.                                 ;not sure what to do with the lower weighting
        big_frame_total = big_frame_total + big_frame                         ;not sure what's going on here. for display? 

        wset, 0
        cgloadct, ct, /silent
        cgimage, bytscl(big_frame, 0, 4000)
        cgtext, 0.02, 0.96, 'Extracted Na emission', charsize=0.5*movie_Scale, col='black', /normal, charthick=1.25*movie_scale
        cgtext, 0.02, 0.96, 'Extracted Na emission', charsize=0.5*movie_Scale, col='white', /normal, charthick=0.5*movie_Scale
        cgtext, 0.02, 0.45, 'IMAGING CHANNEL', charsize=0.5*movie_Scale, col='black', /normal, charthick=1.25*movie_Scale
        cgtext, 0.02, 0.45, 'IMAGING CHANNEL', charsize=0.5*movie_Scale, col='white', /normal, charthick=0.5*movie_Scale
        cgtext, .82, .97, 'Frame # ' + strcompress(i+1,/remove_all), /data, col='black', /normal, charthick=1.5
        if not no_display then begin
          movie_png = cgsnapshot(filename = outdir + 'Mercury scan.png', /nodialog)
          image = read_png(outdir + 'Mercury scan.png')
          void = video3 -> Put(stream3, image)
        endif       
        big_cube[*,*,i] = big_frame
        big_weighting_cube[*,*,i] = big_weighting_frame
        img_cube = img_cube + reform(big_frame(*,0:s[1]))
        spec_frame = reform(big_frame(*,s[1]+1:*))
        spec_real = where(finite(spec_frame), complement=spec_bkg)
        ;spec_dwell[spec_real] = spec_dwell[spec_real] + integration       ; dwell time
        spec_dwell[spec_real] = (spec_dwell[spec_real] + integration)*empty_weighting_frame
        ;spec_frame[spec_bkg] = 0d
        spec_cube[spec_real] = spec_cube(spec_real) + spec_frame(spec_real)

        wset, 1
        spec_dwell_num = where(spec_dwell gt 0.)                          ; find pixels with spectral information included
        spec_frame = smooth(spec_frame, [0,8])                            ; HACK (just in the sense that maybe we don't want to smooth here)
        if ig lt n_elements(good_frames)-1 then $
          spec_temp[spec_dwell_num] = spec_cube[spec_dwell_num] / spec_dwell[spec_dwell_num] + spec_frame[spec_dwell_num] else $   
          spec_temp[spec_dwell_num] = spec_cube[spec_dwell_num] / spec_dwell[spec_dwell_num]  ; ie just don't plot slit contribution for last frame
        
        ;spec_movie = median(big_cube, dimension=3, /even)
        ;spec_movie = spec_movie[spec_extraction[0]:spec_extraction[1],spec_extraction[2]:spec_extraction[3]]
        ;spec_movie = bytscl(spec_movie, max(spec_movie)*0.05, max(spec_movie)*.8)
        
        spec_movie = bytscl(spec_temp[spec_extraction[0]:spec_extraction[1],spec_extraction[2]:spec_extraction[3]], max_spec*0.05, max_spec*.8)  
        
        spec_movie = congrid(spec_movie, xs_spec*movie_scale, ys_spec*movie_scale)
        cgimage, spec_movie
        cgcolorbar, /vertical, charsize=1.e-3;, oob_low='white'
        cgtext, .7, .05, '# frames = ' + strcompress(ig+1,/remove_all), /normal, charthick=0.5*movie_Scale, font=-1, charsize=0.5*movie_Scale
        cgtext, 0.02, 0.96, 'Mercury Na', /normal, charthick=0.75*movie_Scale, charsize=0.625*movie_Scale, font=-1
        if not no_display then begin
          movie_png = cgsnapshot(filename = outdir + 'Mercury Na.png', /nodialog)
          image = read_png(outdir + 'Mercury Na.png')
          void = video -> Put(stream, image)
        endif
        wset, 2
;        img_movie = bytscl(img_cube(img_extraction(0):img_extraction(1),img_extraction(2):img_extraction(3)), max_img*0.2, max_img)
        img_cube_cut = img_cube(img_extraction(0):img_extraction(1),img_extraction(2):img_extraction(3))
        img_movie = bytscl(img_cube_cut, 0.2*max(img_cube_cut), max(img_cube_cut))
        img_movie = congrid(img_movie, xs_img*movie_scale, ys_img*movie_scale)
        cgloadct, 0, /silent
        cgimage, img_movie
        cgcolorbar, /vertical, charsize=1.e-3
        cgtext, .7, .05, '# frames = ' + strcompress(ig+1,/remove_all), /data, /normal, charthick=0.375*movie_Scale, charsize=0.5*movie_scale
        cgtext, 0.02, 0.96, 'Solar continuum', /normal, charthick=0.75*movie_Scale, charsize=0.625*movie_Scale
        if not no_display then begin
          movie_png = cgsnapshot(filename = outdir + 'Mercury disk.png', /nodialog)
          image = read_png(outdir + 'Mercury disk.png')
          void = video2 -> Put(stream2, image)
        endif
        wset, 4
        cgimage, bytscl(big_frame_total / double(i),0,4000)
        ;if ig eq 3 then break
   endfor ;ig

   ; LOOP over bad frames too for comparison figure ???
   dobad = 1
   if dobad then begin
     img_cube_bad = fltarr(s[0], s[1])
     for ib = 0, n_elements(bad_frames)-1 do begin                         ; loop over bad frames
        i = bad_frames(ib)
        big_frame = fltarr(s[0], 2*s[1])
        empty_frame = make_array(s[0], s[1], /float, value = !values.F_NaN)
        ;empty_frame[home[0]+shift_array[i,0]-1:home[0]+shift_array[i,0]+1,*] = [brightness[i, *], brightness[i, *], brightness[i, *]] ;Make it bigger?
        empty_frame[home[0]+shift_array[i,0],*] = [brightness[i, *]]
        ;big_frame[*, s[1]:*] = rotate(empty_frame, 7) 
        big_frame[*, s[1]:*] = rotate(empty_frame, 6)
        big_frame[*, 0:s[1]-1] = reform(aligned_imaging_cube[*,*,i])
        big_frame_total = big_frame_total + big_frame
        big_cube[*,*,i] = big_frame
        img_cube_bad = img_cube_bad + reform(big_frame(*,0:s[1]))
      endfor ;ib     
    img_cube_cut_bad = img_cube_bad(img_extraction(0):img_extraction(1),img_extraction(2):img_extraction(3))
   endif ;dobad==1
     
   wset, 4
   test = median(big_cube, dimension = 3, /even) ; no weights
   test2 = total(big_cube*big_weighting_cube, 3, /NaN) / total(big_weighting_cube, 3, /NaN)
   test3 = test2
   test3[*, s[1]:*] = test3[*, s[1]:*] / spec_dwell ;to normalize the dwell time
   stop
   ;tv, bytscl(smooth(test,[3,3], /NaN), 0, 10000)
   cgimage, bytscl(test, 0, 4500)
   spec_dwell_num = where(spec_dwell gt 0.)                               ; find pixels with spectral information included
   spec_cube(spec_dwell_num) = spec_cube(spec_dwell_num) / spec_dwell(spec_dwell_num)  ; normalize by # of frames summed
   ;stop
   ; Now add some blank frames to the Na image movie before plotting a smoothed frame
   spec_raw = spec_cube(spec_extraction(0):spec_extraction(1),spec_extraction(2):spec_extraction(3))
   spec_raw = congrid(spec_raw, xs_spec*movie_scale, ys_spec*movie_scale) ; embiggen
   spec_movie = bytscl(spec_raw, max_spec*0.05, max_spec*.8)
   wset, 1
   loadct, ct, /silent
   cgimage, spec_movie
   cgcolorbar, /vertical, charsize=1.e-3
   cgtext, .7, .05, '# frames = ' + strcompress(ig,/remove_all), /normal, charthick=0.5*movie_Scale, font=-1, charsize=0.5*movie_Scale
   cgtext, 0.02, 0.96, 'Mercury Na', /normal, charthick=0.75*movie_Scale, charsize=0.625*movie_Scale, font=-1
   movie_png = cgsnapshot(filename = outdir + 'Mercury Na.png', /nodialog)
   image = read_png(outdir + 'Mercury Na.png')
   for i = 0, 30 do void = video -> Put(stream, image)
   smooth_width = 10
   cgimage, smooth(spec_raw,smooth_width*movie_scale,/edge_truncate)                ; HACK
   cgcolorbar, /vertical, charsize=1.e-3
   cgtext, .7, .05, '# frames = ' + strcompress(ig,/remove_all), /normal, charthick=0.5*movie_Scale, font=-1, charsize=0.5*movie_Scale
   cgtext, 0.02, 0.96, 'Mercury Na (smoothed=' + strcompress(smooth_width,/remove_all) + ')', /normal, charthick=0.75*movie_Scale, charsize=0.625*movie_Scale, font=-1
   movie_png = cgsnapshot(filename = outdir + 'Mercury Na.png', /nodialog)
   image = read_png(outdir + 'Mercury Na.png')
   for i = 0, 30 do void = video -> Put(stream, image)
      
   video -> Cleanup
   video2 -> Cleanup
   video3 -> Cleanup

   ; Overplot Na contours on top of Mercury disk
   loadct, 0, /silent
   cgimage, bytscl(img_movie), /axes
   spec_smoothed = bytscl(smooth(spec_raw,smooth_width*movie_scale,/edge_truncate))
   cgcontour, spec_smoothed, /onimage, color='black', levels=[50,100,150,175,200,225,250], xthick=3 
   cgcontour, spec_smoothed, /onimage, color='gold', levels=[50,100,150,175,200,225,250]
   al_legend, ['Disk','Na'], box=0, /bottom, /right, charsize=0.5*movie_scale, line=[0,0], linsize=0.2, thick=[6,2], color=['gray','gold'], textcolor='white'
   combined = cgsnapshot(filename = outdir + 'Mercury Na + disk.png', /nodialog)

   ; Generate the "bad images" Mercury disk for comparison
   if dobad then begin
     bad_img = bytscl(img_cube_cut_bad, 0.2*max(img_cube_cut_bad), max(img_cube_cut_bad))
     cgimage, bad_img
;    cgimage, img_cube_cut_bad
     cgcolorbar, /vertical, charsize=1.e-3
     cgtext, .7, .05, '# frames = ' + strcompress(ib,/remove_all), /data, /normal, charthick=0.375*movie_Scale, charsize=0.5*movie_Scale
     cgtext, 0.02, 0.96, 'Solar continuum', /normal, charthick=0.75*movie_Scale, charsize=0.625*movie_Scale
     bad_png = cgsnapshot(filename = outdir + 'Mercury disk - bad frames.png', /nodialog)

     ; And might as well compare the good vs bad contours while we're here...
     cgimage, img_movie, /axes
     cgcontour, img_movie, /onimage, color='gold'
     cgcontour, bad_img, /onimage, color='red'
     al_legend, ['good','bad'], box=0, /bottom, /right, charsize=2, line=[0,0], linsize=0.2, thick=[4,4], color=['orange','red'], textcolor=['orange','red']
     compare = cgsnapshot(filename = outdir + 'Mercury disk - good vs bad contours.png', /nodialog)
   endif ;dobad==1

   ; Create image showing the spatial coverage of the slit position
   wset, 2
   dwell_img = spec_dwell(spec_extraction(0):spec_extraction(1),spec_extraction(2):spec_extraction(3))
   cgloadct, 22, /silent
   cgimage, bytscl(dwell_img, 0, max(dwell_img))
   cgloadct, 22, ncolors=(1+fix(max(spec_dwell))), bottom=0
   cgcolorbar, /vertical, maxrange=max(spec_dwell), title='Total Na data coverage (sec)', /discrete, ncolors=fix(max(spec_dwell))
   combined = cgsnapshot(filename = outdir + 'Mercury slit coverage.png', /nodialog)

   beep
   stop
endif

; cgPS2Raster, 'test_2.ps', /JPEG
;   IDL> cgPS_Open, 'test_2.ps'
;   IDL> Graphic_Display
;   IDL> cgPS_Close, /PNG
  
end



