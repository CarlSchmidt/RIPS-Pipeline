pro RIPS_Mercury_Perkins_v02, PART=part, NIGHT=night

; *********************************************************************************************************************
; *********************************************************************************************************************
; Routine to extract Na emission from spectral scans over Mercury's disk with Perkins/RIPS in order to create a 2-D
; image of Na above and near Mercury's disk.
; 
; THIS REDUCTION PROGRAM IS SPLIT INTO 5 PARTS:
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
;       input:   spectral datacube from Part 2 (wl,y,t) with surface reflectance removed
;       output:  Na brightness (wl,y) and linewidth (wl,y)
;   4 = extract Surface reflectance and compare with a Hapke Model. Flux calibrate the data
;       input:   Na brightness and linewidth from Part 3 and imaging cube from Part 0
;       output:  Hapke Model, calibrated Na D images in Rayleighs and continuum in R/A
;  99 = do all of 0-4
;
; Input variables
;   part           - the "part" of the program to execute
;   night          - the night of RIPS data to use; only '15' is viable from the March 2018 run, but we've kept
;                    the option to specify night in case future runs are more successful!  There are numerous 
;                    sub-variables associated with each night: e.g., mercury_dir, dark_dir, flat_dir, sky_dir, 
;                    Mercury_file, sky_file, dark_file, flat_file, Na_D_rngs, ims, sps
;   width          - median smoothing width to get local background reference (for pixels above thresh; Default=5)
; 
; Created:  18 Oct 2018 - copied from RIPS_Mercury_AEOS_v03.pro to start (CS and LM)
; Edits:    25 Oct 2018 - rennamed to "v02", cleaned up with the plan to leave v02 as a minimal, yet still hard-coded
;                         version, and then to improve the automation in v03 (plus add additional improvements) (LM)  
;           26 Oct 2018 - re-worked for native IDL where simple 
;                         changed slit interpolation method used when finding the centroid in the imaging channel (CS)
;           3 July 2019 - re-write and simplify. Improved spatial alignment. (CS)

Start_time = SYSTIME(/SECONDS)

; Choose which drive we're pulling data from:
  Local = FILE_SEARCH('D:\DATA\Perkins\Perkins RIPS - March 2018\*')    ; Carl's Computer
  Cloud = FILE_SEARCH('Y:\obs_18\Perkins_RIPS_March\*')                 ; Luke's Cloud
  if local[0] ne '' then base_dir = 'D:\DATA\Perkins\Perkins RIPS - March 2018\'
  if Cloud[0] ne '' then base_dir = 'Y:\obs_18\Perkins_RIPS_March\'
                                                      
; =====================================================================================================================
; Define variables
; =====================================================================================================================
SetDefaultValue, part, 4.
SetDefaultValue, night, '15'                                              ; the night of the run; default '15' (only kept for possible future runs with multiple nights)
SetDefaultValue, thresh, 12                                               ; hotpixel removal threshold (sigma above background)
SetDefaultValue, width, 5                                                 ; Median smoothing width to get local background reference (for pixels above thresh) 
SetDefaultValue, minfac, 0.9                                              ; maximum factor to apply to reference spectrum for Na emission extraction
SetDefaultValue, maxfac, 1.3                                              ; minimum factor to apply to reference spectrum for Na emission extraction
SetDefaultValue, dfac, 0.01                                               ; factor increment for minfac->maxfac
SetDefaultValue, ct, 20                                                   ; default color table
SetDefaultValue, max_spec, 2500.                                          ; partial hack --> the max value of the spectral Na cube after coadding (used for "prettier" movies)
SetDefaultValue, max_img, 1.6e6                                           ; partial hack --> the max value of the image cube after coadding (used for "prettier" movies)
SetDefaultValue, img_extraction, [145,295,85,235]                         ; coordinates for extracting a square Mercury image, [x1,x2,y1,y2]
SetDefaultValue, spec_extraction, [145,295,77,227]                        ; coordinates for extracting a square Mercury spectral image, [x1,x2,y1,y2]
SetDefaultValue, movie_scale, 4                                           ; scale the x-y dimensions of the extracted "Na image" and "Mercury disk" arrays for movies by this value (used for bigger and therefore "prettier" movies)
SetDefaultValue, do_rotate, 0                                             ; 0=no rotation of image, 1=rotate Na image (inverted, =5 -- see below)
SetDefaultValue, do_median, 0                                             ; 0=use mean when constructing Na Mercury image, 1=use median instead
SetDefaultValue, do_smooth, 1                                             ; set to 1 to smooth the imaging channel before spatially coaligning
SetDefaultValue, smooth_width, 6                                          ; if do_smooth = 1, the pixel width of the smoothing kernel, HACK: CURRENTLY THE SAME FOR SPATIAL AND SPECTRAL. BE CAREFUL!
SetDefaultValue, stddev_cutoff, 1.                                        ; only use Mercury frames with stddev > [stddev_cutoff]_sigma, these are defined as "good" 
SetDefaultValue, Flux_Cal_Wavelength, 5893.                               ; Wavelength to extract and compare to Hapke code (Angstroms)

mercury_dir   = base_dir+'15\'                                            ; directory with Mercury FITS file(s)
dark_dir      = base_dir+'14\Carl_Keep\'                                  ; directory with dark file
flat_dir      = base_dir+'15\'                                            ; directory with flat file
sky_dir       = base_dir+'15\'                                            ; directory with sky file
Mercury_file  = strcompress(789+indgen(17), /remove_all) + '.fits'       ; SUFFIX of the Mercury kinetic series to use (SKIPPING 788 BECAUSE OF IT'S DIFFERENT ROTISERIZER ANGLE)
;Mercury_file  = strcompress(788, /remove_all) + '.fits'                   ; SUFFIX of the Mercury kinetic series to use (SKIPPING 788 BECAUSE OF IT'S DIFFERENT ROTISERIZER ANGLE)
sky_file      = 'RIPS1_Thu Mar 15 2018_01.12.12_785.fits'                 ; sky frame
dark_file     = 'Dark.fits'                                               ; dark frame
flat_file     = 'Flat_NaSpectra_624slitwidth_NaND1imaging'                ; flat frame withouth number and fits  e.g., filaname1.fits' 
n_flats       = 3                                                         ; Number of flat fields  
Na_D_rngs     = [235,255,475,495]                                         ; ranges of Na D Mercury emission (x pixels in spec domain for D1 and D2)
Na_D_rng      = 10                                                        ; number of pixels around Na D Mercury emission to avoid (in illum scaling for part 2)
do_realign    = 0                                                         ; 0=run part 1 as normal, 1=use previous cross-correlation as master template
;frame_ref     = 203                                                       ; this number identifies a frame # in imaging_cube used as the reference key spatially aligning all other frames, works best if object is off-slit
frame_ref     = 152
;frame_ref     = 52
offset_maxes  = [60,60]                                                   ; if the [x,y] spatial offset returned by image correlation is larger than these values, correlation fails and use the brightness centroid 
ims           = [283,759,39,374]                                          ; variables for image analysis: x1,x2,y1,y2 (formerly "imaging_statsec")
sps           = [99,899,502,941]                                          ; variables for spectral analysis: x1,x2,y1,y2 (formerly "spectra_statsec")
iref_factor   = 1.                                                        ; factor to scale reference spectrum by
count_factor  = 1.e4                                                      ; this is the scaling factor to apply to the white light Mercury image (i.e., the slit covers some regions of the disk, so we scale based on total slit coverage at each location)
outdir        =  Mercury_dir + 'Processed\'                               ; output directory
Spectral_platescale = 0.09669                                             ; SPECTRAL channel "/pix Jeff's Jun 5th 2018 email: Using the equatoial diameter of 44.48" on april 26, and the measured size on the spectral channel (460 pix) we get 0.09669 ["]/pix
imaging_platescale  = 0.114820                                            ; IMAGING channel "/pix from measureed Io - Europa seperation on 2018-04-27T07:29:33
Na_D1_WL      = 5895.92424                                                ; Rest wavelength from NIST
Na_D2_WL      = 5889.95095                                                ; Rest wavelength from NIST
Body          = 'Mercury'                                                 ; For SPICE purposes

; COULD OMIT THIS AND CALCULATE IT IN PART 3 AND 4!
; Dispersion calculation: at y of 265, lines occur at 390 and 632 in x 
  dispersion          = (Na_D1_WL - Na_D2_WL) / (632. - 390.)  ; Angstroms per pixel
  dispersion_Velocity = cspice_clight()*dispersion / Na_D2_WL  ; Km/s per pixel at Na D2
  
; Imaging_platescale:
  ;  Io - Europa Seperation on 2018-04-27T07:29:33 is...
  ;  pixel_sep = sqrt((794.-520.)^2+(422.-228.)^2)
  ;  angular_sep = 3600.*sqrt( ((227.76500 - 227.75408)*cos(-16.46847/!radeg))^2 + (16.46847-16.47068)^2 )
  ;  print, angular_sep / pixel_sep
  plate_scale_ratio  =  Spectral_platescale / imaging_platescale            ; Even if Jeff and Carl are both a few pixels off, this is accurate within ~ 3%
  
; ======================================LOAD SPICE========================================================================
  kernel_directory    = 'C:\SPICE\'
  CSPICE_KTOTAL, 'all', count
  PRINT, STRCOMPRESS('Cleaning ' + STRING(count) + ' lingering kernels out of memory . . .')
  i = 0
  
  WHILE i LT count DO BEGIN
    CSPICE_KDATA, 0, 'all', file, type, source, handle, found
    CSPICE_UNLOAD, file
    i = i + 1
  ENDWHILE

  CSPICE_FURNSH, STRCOMPRESS(FILE_SEARCH(kernel_directory, 'pck00009.tpc')) ; Body rotational states
  CSPICE_FURNSH, STRCOMPRESS(FILE_SEARCH(kernel_directory, 'naif0009.tls')) ; Leap seconds
  CSPICE_FURNSH, STRCOMPRESS(FILE_SEARCH(kernel_directory, 'de431.bsp'))    ; SPK (ephemeris kernel) for planets MOST CURRENT
  CSPICE_FURNSH, STRCOMPRESS(FILE_SEARCH(kernel_directory, 'sat317.bsp'))   ; SPK (ephemeris kernel) for satellites (most applicable data)
  CSPICE_FURNSH, STRCOMPRESS(FILE_SEARCH(kernel_directory, 'gm_de431.tpc')) ; Body Masses
  CSPICE_KTOTAL, 'all', count
; =================================SPICE KERNELS ALL LOADED===============================================================
; 
; Find Mercury's True Anomaly Angle
  filenames         = FILE_SEARCH(Mercury_dir + '*_'+Mercury_file)
  header = headfits(filenames[0])
  UTC = sxpar(header, 'DATE')
  cspice_str2et, UTC, et
  cspice_spkezr, Body, et, 'J2000', 'LT+S', 'Sun', state, ltime   
   cspice_bodvrd, 'MERCURY', 'GM', 1, mu
   mu = 1.3271246207371930E+11 ;Mercury hack from horizons!
   cspice_oscelt, state, et, mu, elts
   ecc = ELTS[1]        ; Eccentricity
   MA  = ELTS[5]        ; Mean Anomaly radians
   ; Now get true anomaly from Mean anomaly and Eccentricity, see https://en.wikipedia.org/wiki/True_anomaly 
   ; Roy, A.E. (1988). Orbital Motion (1 ed.). Bristol, UK; Philadelphia, PA: A. Hilger. ISBN 0852743602.
     TAA =  MA + (2.*ecc - 0.25*ecc^3)*sin(MA) + 1.25*(ecc^2)*sin(2.*MA)+ (13./12.)*(ecc^3)*sin(3.*MA) 
   True_Anomaly = TAA*!radeg ; True Anomaly in degrees 
   print, 'Mercury''S True Anomaly Angle:', True_Anomaly

; =====================================================================================================================
; Find display size, set window defaults, reset plotting variables just in case
; =====================================================================================================================
  cgloadct, ct
  CLEANPLOT                                                               ; Reset all plotting system variables (!P,!X,!Y,!Z) to their default values
  DEVICE, GET_SCREEN_SIZE = ss                                            ; find screen size
  winpos_x      = 0                                                       ; default initial x window location
  winpos_y      = 0                                                       ; default initial y window location

; =====================================================================================================================
; Part 0 : break kinetic series into two manageable sized datacubes
; =====================================================================================================================
if part eq 0 or part eq 99 then begin 
    header            = headfits(filenames[0])
    integration       = sxpar(header, 'EXPOSURE')                         ; integration time for individual frames
    n_frames_per_file = sxpar(header, 'NUMKIN')                           ; # of frames in the "kinetic series" for each "Mercury_file"
    nfiles            =  n_elements(Mercury_file)                         ; # of different Mercury kinetic frames
    s    = [1024, 1024, n_frames_per_file]
    cube = intarr(s[0], s[1], s[2]*nfiles)                                ; define a ***big** cube, where multiple kinetic series are stacked in the third dimension, INTEGER TYPE for now
    print, 'Found: ', N_elements(filenames), ' files'
    window, 0, xs = ims[1] - ims[0], ys = ims[3]-ims[2], title = 'Inspection: First frame, Imaging Channel' 
    for ifile = 0, nfiles-1 do begin                                      ; loop over kinetic series and add each to a stack
      icube = MRDFITS(filenames[ifile], 0, header, /unsigned, /silent )
      cube(*,*,s[2]*ifile:s[2]*ifile+n_frames_per_file-1) = icube
      integration       = sxpar(header, 'EXPOSURE')                       
      n_frames_per_file = sxpar(header, 'NUMKIN')   
      print, 'Loading: '+Mercury_file[ifile],' containing', strcompress(n_frames_per_file)+' frames of'+strcompress(integration)+' seconds each'    
      tv, bytscl(icube[ ims[0]:ims[1], ims[2]:ims[3], 0])
    endfor
    imaging_cube = float(reform(cube[ims[0]:ims[1],ims[2]:ims[3], *]))           ; extract the IMAGING portion of the kinetic series and combine
    spectra_cube = float(reform(cube[sps[0]:sps[1],sps[2]:sps[3], *]))           ; extract the SPECTRAL portion of the kinetic series and combine
    MWRFITS, imaging_cube, outdir + 'imaging_cube.fits', header, /create, /silent
    MWRFITS, spectra_cube, outdir + 'spectra_cube.fits', header, /create, /silent
    
  ; We'll want a "sky" cube as well  
    Sky_cube          = MRDFITS(sky_dir + sky_file, 0, header, /unsigned, /silent )
    imaging_sky_cube  = float(reform(sky_cube[ims[0]:ims[1],ims[2]:ims[3], *]))
    spectra_sky_cube  = float(reform(sky_cube[sps[0]:sps[1],sps[2]:sps[3], *]))
    MWRFITS, imaging_sky_cube, outdir + 'imaging_sky_cube.fits', header, /create, /silent
    MWRFITS, spectra_sky_cube, outdir + 'spectra_sky_cube.fits', header, /create, /silent
    beep
endif ; Note that above needs modification in case the integration times change between the fits files we're co-adding

; =====================================================================================================================
; Part 1 : find Mercury centroids by cross-correlation with a reference image
; =====================================================================================================================
if part eq 1 or part eq 99 then begin ;Find the centroids by cross-correlation with the last image
    if part ne 99 then imaging_cube = MRDFITS(outdir + 'imaging_cube.fits', 0, header, /fscale, /silent )  ; don't need to read in if we're continuing
    if part ne 99 then spectra_cube = MRDFITS(outdir + 'spectra_cube.fits', 0, header, /fscale, /silent )
    frame     = reform(imaging_cube[*,*,0])                               ; individual frame from kinetic series
    s         = size(imaging_cube, /dimensions)                           ; size of extracted IMAGING portion of frame
    ss        = size(spectra_cube, /dimensions)                           ; size of extracted SPECTRAL portion of frame

    ; Define variables based on frame sizes
      Imaging_shift_array = intarr(s[2],2)                                ; array of x, y shift values for aligning images
      sharpness_metric    = fltarr(s[2])                                  ; standard deviations of imaging frames (TEMPORARILY used to define image quality)
    
    ; Set up windows based on image sizes
      window, 0, xs = s[0], ys = s[1], xpos=winpos_x,             ypos=winpos_y,  title = 'Flat Field Inspection: Red = Behind Slit --> interpolate'
      window, 2, xs = s[0], ys = s[1], xpos=winpos_x+s[0]+20,     ypos=0,         title='Raw Reference Frame for spatial alignment'
      window, 3, xs = s[0], ys = s[1], xpos=winpos_x+s[0]+20,     ypos=s[1]+40,   title='Filtered/Smoothed Comparison'
      window, 4, xs = s[0], ys = s[1], xpos=winpos_x+2*(s[0]+20), ypos=0,         title='Aligned Raw Frame'
      window, 5, xs = s[0], ys = s[1], xpos=winpos_x+2*(s[0]+20), ypos=s[1]+40,   title='Aligned Filtered/Smoothed Comparison'

    ; Prepare Darks
      iDark = MRDFITS(dark_dir + dark_file, 0, header, /fscale, /silent )   ; read in dark
      dark  = idark[ims[0]:ims[1],ims[2]:ims[3]]
      dark  = sigma_filter( dark, width, N_sigma=thresh)
    
    ; Prepare Imaging Channel Flats 
      big_array = fltarr(1024, 1024, N_flats)
      for i = 0, N_flats-1 do big_array[*,*,i] = MRDFITS(flat_dir + flat_file + strcompress(i+1, /remove_all) + '.fits', 0, header, /fscale, /silent ) ; read in flats
      iflat = median(big_array, dimension = 3)  ; combine flat flat
      flat  = ( iflat[ims[0]:ims[1],ims[2]:ims[3]] - dark ) / max( iflat[ims[0]:ims[1],ims[2]:ims[3]] - dark ) ; Crop to just imaging portion
      flat  = flat + (1. - median(flat))                                    ; this normalizes such that the median of the flat is 1
      flat[WHERE(flat lt .01, /NULL)] = !values.f_Nan                       ; reject unusual counts for centroid
      flat[WHERE(flat gt 2.0, /NULL)] = !values.f_Nan                       ; reject unusual counts for centroid
    gooddata = where(Finite(flat), ngooddata, comp=baddata, ncomp=nbaddata)
    if nbaddata gt 0 then flat[baddata] = interpol(flat[gooddata], gooddata, baddata)
;    if do_realign then begin
;        restore, outdir + 'shift_array.sav'                               ; contains shift_array, sharpness_metric, and aligned_imaging_cube
;        shift_array_old = shift_array                                     ; keep old shift_array for comparison
;        reference = total(aligned_imaging_cube,3)                         ; create "master Mercury" template for cross-correlation
;    endif else begin
        ireference = (reform(imaging_cube[*,*,frame_ref]) - dark) / flat
        acre, ireference, reference, thresh, width                         ; clean any hot pixels
;    endelse

    ; Find indicies of all pixels under the slit in the imaging channel
      junk = min(total(flat, 2), slit_center)                         ; roughly find the slit center
      h    = histogram(flat, binsize = .05, REVERSE_INDICES=ri)       ; bin the flat into 0.05 bins
      slit_indicies = ri[ri[0]:ri[9]-1]                               ; find inidices of pixels in the lowest few bins of the flat's histogram 
      slit_indicies = array_indices(flat, slit_indicies)              ; convert into x & y indices
      keep = where(abs(slit_indicies[0,*] - slit_center) lt 10, /Null); keep only inidices within 10 pixels of slit center       
      slit_indicies = slit_indicies[*,keep]
      dummy = flat
      dummy[ slit_indicies[0,*], slit_indicies[1,*] ] = max(flat)     ; show the location of pixels fully under the slit.
      wset, 0      
      cgimage, hist_equal(dummy), title='flat'
      
    ; Interpolate pixels over the slit in Reference image that we're aligning to, best to choose a 'frame_ref' reference frame number with the slit fully off the disk
      reference[ slit_indicies[0,*], slit_indicies[1,*] ] = !values.F_Nan  ; define slit pixels as NaN
      fill_missing, reference, !values.F_Nan, 1                            ; interpolate over slit

    ; Log the correlations to the reference image (correlation post-alignment)
      Correl     = fltarr(s[2]) 
      Spectra_DN = fltarr(s[2]) 
      Imaging_DN = fltarr(s[2]) 

    print, 'Calculating Frame-By-Frame Spatial Shifts for Co-alignment...'
    print, '-----Frame-#---X-------Y-----Centroid-WRT-Frame-#', strcompress(frame_ref,/remove_all)
    for i = 0, s[2]-1 do begin
      iframe = (reform(imaging_cube[*,*,i]) - dark) / flat              ; "raw" imaging frame
      iframe[ slit_indicies[0,*], slit_indicies[1,*] ] = !values.F_Nan  ; define slit pixels as NaN
      fill_missing, iframe, !values.F_Nan, 1                            ; interpolate over slit
      frame = iframe
      if do_smooth then begin                                           ; set up bandpass filter to aid image alignment
        img_bandpass = float(bandpass_filter(smooth(frame,smooth_width,/edge_truncate), 0., 0.15, /butterworth))
        ref_bandpass = float(bandpass_filter(smooth(reference,smooth_width,/edge_truncate), 0., 0.15, /butterworth))
      endif else begin
        img_bandpass = float(bandpass_filter(frame, 0., 0.15, /butterworth))
        ref_bandpass = float(bandpass_filter(reference, 0., 0.15, /butterworth))
      endelse
            
      junk = max(total(ref_bandpass,2), ref_bandpass_centeroid_x) & junkx = max(total(img_bandpass,2), img_bandpass_centeroid_x)
      junk = max(total(ref_bandpass,1), ref_bandpass_centeroid_y) & junky = max(total(img_bandpass,1), img_bandpass_centeroid_y)
      CORREL_OPTIMIZE, ref_bandpass[ref_bandpass_centeroid_x-80:ref_bandpass_centeroid_x+80, ref_bandpass_centeroid_y-80:ref_bandpass_centeroid_y+80], $
                       img_bandpass[ref_bandpass_centeroid_x-80:ref_bandpass_centeroid_x+80, ref_bandpass_centeroid_y-80:ref_bandpass_centeroid_y+80], xoffset_optimum, yoffset_optimum, /NUMPIX
      Imaging_shift_array[i,*] = [xoffset_optimum, yoffset_optimum]             ; track required alignment offsets
      
      CORREL[i] = CORREL_IMAGES( Ref_bandpass[ref_bandpass_centeroid_x-80:ref_bandpass_centeroid_x+80, ref_bandpass_centeroid_y-80:ref_bandpass_centeroid_y+80], $
                                 img_bandpass[ref_bandpass_centeroid_x-80:ref_bandpass_centeroid_x+80, ref_bandpass_centeroid_y-80:ref_bandpass_centeroid_y+80], $
                                 xoffset_b = xoffset_optimum, yoffset_B = yoffset_optimum, XSHIFT =0, ySHIFT = 0) ;2-D Correlation at this alignment
      Spectra_DN[i] = total(spectra_cube[*,*,i])
      Imaging_DN[i] = total(img_bandpass[ref_bandpass_centeroid_x-80:ref_bandpass_centeroid_x+80, ref_bandpass_centeroid_y-80:ref_bandpass_centeroid_y+80])
        
      ; If correlation fails then use the brightness centroid of the filtered/smoothed frames
        if abs(xoffset_optimum) gt offset_maxes[0] or abs(yoffset_optimum) gt offset_maxes[1] then begin ; If Failed Correlation 
          Imaging_shift_array[i,*] = [ref_bandpass_centeroid_x - img_bandpass_centeroid_x, ref_bandpass_centeroid_y - img_bandpass_centeroid_y]
          print, i, Imaging_shift_array[i,0], Imaging_shift_array[i,1], '  <--- Correlation alignment failed, used brightness centroid'   
        endif else print, i, Imaging_shift_array[i,0], Imaging_shift_array[i,1]           ; update progress 

      ; Plot & Inspect progress to track the quality of the alignments
        wset, 2
          cgimage, reference, /keep_aspect
        wset, 4
          f = shift(frame, [Imaging_shift_array[i,*]])
          cgimage, f, /keep_aspect
          cgcontour, ref_bandpass/max(ref_bandpass), levels=[0.2,0.4,0.6,0.8], position = [0,0,1,1], /noerase
        wset, 3
          cgcontour, ref_bandpass/max(ref_bandpass), levels=[0.2,0.4,0.6,0.8], position = [0,0,1,1], /fill
          cgcontour, ref_bandpass/max(ref_bandpass), levels=[0.2,0.4,0.6,0.8], position = [0,0,1,1], /noerase
          cgcontour, img_bandpass/max(img_bandpass), levels=[0.2,0.4,0.6,0.8], c_color='red', position = [0,0,1,1], /noerase
        wset, 5
          cgcontour, ref_bandpass/max(ref_bandpass), levels=[0.2,0.4,0.6,0.8], position = [0,0,1,1], /fill
          cgcontour, ref_bandpass/max(ref_bandpass), title='aligned bandpass images', levels=[0.2,0.4,0.6,0.8], position = [0,0,1,1], /noerase
          f = shift(img_bandpass, [Imaging_shift_array[i,*]])
          cgcontour, f/max(f), /onimage, levels=[0.2,0.4,0.6,0.8], c_color='red', position = [0,0,1,1], /noerase
    endfor  ; end loop to determine Imaging_shift_array values to align images

    Spectra_shift_array = fix(round(float(Imaging_shift_array) / plate_scale_ratio))  ; Spectral and imaging channels have different platescales, hence different shifts
    Spectra_shift_array[*,1] = -Spectra_shift_array[*,1]                              ; X_i = X_s, but +Y_i in the imaging channel is -Y_s in the spectral channel
    aligned_imaging_cube = fltarr(s)                                                  ; Imaging cube after everything is aligned
    aligned_spectra_cube = fltarr(ss)                                                 ; Spectral cube after alignment in Y only (for now) (note that 

    ; Now loop over all frames to align both channels
      for i = 0, s[2]-1 do begin
        frame     = (reform(imaging_cube[*,*,i]) - dark) / flat                       ; "raw" imaging frame
        specframe = (reform(spectra_cube[*,*,i]))                                     ; "raw" spectral frame
        aligned_imaging_cube[*,*,i] = shift(frame, [Imaging_shift_array[i,*]])        ; Co-align imaging frame in X and Y
        aligned_spectra_cube[*,*,i] = shift(specframe, [0, Spectra_shift_array[i,1]]) ; Co-align Spectral frame **in Y only**. Y_Spec = -Y_img/PSR. Later on, we do the X_Spec = X_img/PSR.                     

        ; *****************************Sharpness criterion for lucky imaging************************************
          iframe[ slit_indicies[0,*], slit_indicies[1,*] ] = !values.F_Nan            ; mask slit pixels as NaN
          fill_missing, iframe, !values.F_Nan, 1                                      ; interpolate over slit
          frame = shift(iframe, [Imaging_shift_array[i,*]]) 
          subframe = frame[ref_bandpass_centeroid_x-80:ref_bandpass_centeroid_x+80-1, ref_bandpass_centeroid_y-80:ref_bandpass_centeroid_y+80-1] ; extract just the image part of the frame
          sharpness_metric[i] = total(abs(deriv(rebin(subframe, 40, 40))))  ; total the absolute value of the spatial derivative. higher --> sharper image (TEMP, we can do better!)
      endfor  
      sharpness_metric = sharpness_metric / mean(sharpness_metric) ;normalize this to 1 for convenience 
    save, Imaging_shift_array, Spectra_shift_array, correl, aligned_imaging_cube, aligned_spectra_cube, sharpness_metric, filename = outdir + 'shift_array.sav'
    beep
endif

;******************************testing sharpness metric***********************
if part eq 1.5 then begin
  restore, outdir + 'shift_array.sav'  
  ranked = SORT(correl)
  window, 0
  cgimage, total(aligned_imaging_cube[*,*, ranked[0:99]], 3), /keep_aspect
  window, 1
  cgimage, total(aligned_imaging_cube[*,*, ranked[1500:1599]], 3), /keep_aspect
  ;seems to work fairly well.
  stop
endif
;******************************testing sharpness metric***********************

; =====================================================================================================================
; Part 2 : Isolate the sodium emission in every specctral frame
; =====================================================================================================================
if part eq 2 or part eq 99 then begin 
    img_cube = MRDFITS(outdir + 'imaging_cube.fits', 0, header2, /fscale, /silent )
    if part ne 99 then restore, outdir + 'shift_array.sav'  ; contains Imaging_shift_array, Spectra_shift_array, Correl, sharpness_metric, aligned_imaging_cube and aligned_spectra_cube
    frame     = reform(aligned_spectra_cube[*,*,0])
    s         = size(aligned_spectra_cube, /dimensions) 
    si        = size(img_cube, /dimensions)
    window, 0, xpos=winpos_x,         ypos=winpos_y,           xs=s[0],  ys=s[1],  title='IDL 0 - SPECTRAL FRAME'
    window, 1, xpos=winpos_x,         ypos=winpos_y+s[1]+40,   xs=s[0],  ys=s[1],  title='IDL 1 - SCALED SODIUM FREE REFERENCE FOR REMOVING REFLECTANCE (Sky, Moon etc...)'
    window, 2, xpos=winpos_x+s[0]+20, ypos=winpos_y,           xs=s[0],  ys=s[1],  title='IDL 2 - EXOSPHERE'
    window, 3, xpos=winpos_x+s[0]+20, ypos=winpos_y+s[1]+40,   xs=si[0], ys=si[1], title='IDL 3 - IMAGING CHANNEL FRAME'

    Dark      = MRDFITS(dark_dir + dark_file, 0, header, /fscale, /silent )
    dark      = dark[sps[0]:sps[1],sps[2]:sps[3]]
    dark      = sigma_filter( dark, width, N_sigma=thresh)
    acre, dark, dark, thresh, width                                  ; clean up any remaining hot pixels

    if night eq '15' then begin
      print, 'Applying a sky spectrum as the solar reflectance spectrum...'
      spectra_sky_cube = MRDFITS(outdir + 'spectra_sky_cube.fits', /silent ) 
      ireference = spectra_sky_cube - dark
      acre, ireference, reference, thresh, width

      ; Make a spectral flat and flatten the spectral cube
        big_array = fltarr(1024, 1024, N_flats)
        for i = 0, N_flats-1 do big_array[*,*,i] = MRDFITS(flat_dir + flat_file + strcompress(i+1, /remove_all) + '.fits', 0, header, /fscale, /silent ) ; read in flats
        iflat = median(big_array, dimension = 3)  ; combine flat
        iflat     = (iflat[sps[0]:sps[1],sps[2]:sps[3]] - dark)
        iflat     = iflat / median(iflat) 
        acre, iflat, flat, thresh, width
        flat      = flat + (1. - median(flat))                       ; this normalizes again such that median(flat) = 1
        reference = reference / flat

      ; Now find the xoffset in the reference spectrum (likely due to slightly different grating angles)
        y_Mercury       = where( total(frame,1) ge 0.5*max(total(frame,1)) )    ; FWHM range of Mercury
        frame_spectrum  = total(frame(*,y_Mercury),2)
        ref_spectrum    = total(reference(*,y_Mercury), 2)
        frame_D1_x      = ( where( frame_spectrum(0:s[0]/2) eq min(frame_spectrum(0:s[0]/2))) ) [0]    ; hack (we only want the D1 line, so we just check the 1st half of spectrum)
        ref_D1_x        = ( where( ref_spectrum(0:s[0]/2)   eq min(ref_spectrum(0:s[0]/2))   ) ) [0]   ; hack (we only want the D1 line, so we just check the 1st half of spectrum)
        reference       = shift(reference, frame_D1_x - ref_D1_x)
        ref_spectrum    = shift(ref_spectrum, frame_D1_x - ref_D1_x)

      ; Find ranges of pixels surrounding Mercury Na D lines
        D1_x            = (where(frame_spectrum eq max(frame_spectrum)))[0]    ; x pixel for D1 line
        frame_temp      = frame_spectrum
        frame_temp(D1_x-Na_D_rng:D1_x+Na_D_rng) = mean(frame_spectrum)    ; get rid of Na emission peak to find next Na line
        D2_x            = (where(frame_temp eq max(frame_temp)))[0]       ; x pixel for D2 line

      ; Now scale and subtract 
        Surf_refl_spectra_cube = fltarr( size(aligned_spectra_cube, /dim)); variable to hold the raw Mercury spectra with exosphere + surface reflectance
        exosphere_spectra_cube = fltarr( size(aligned_spectra_cube, /dim)); variable to hold the residual Mercury spectra with exosphere only
        for i = 0, s[2]-1 do begin                                        ; LOOP over each frame in the series
          frame       = reform(aligned_spectra_cube[*,*,i]) - dark        ; dark subtract each spectrum
          frame       = frame / flat                                      ; flatfield each spectrum
      
          ; Now find best scaling for ref_spectrum and subtract (based on solar reflected emission between D lines        
            illum = frame / reference                                     ; scale the illumination against Mercury's disk
            for iD = 0, 2, 2 do illum[Na_D_rngs[iD+0]:Na_D_rngs[iD+1],*] = !values.F_NaN  ; carefully avoid sodium emissions when scaling to the disk
      
            illum_along_slit = MEDIAN(illum, dimension=1)
            illum_along_slit = MEDSMOOTH( illum_along_slit, 5 )           ; helps with hot pixels from bad flat-fielding where the slit has some dust
            scaled_reference = reference * rebin(transpose(illum_along_slit), s[0], s[1])
            Just_exosphere = Frame - scaled_reference * iref_factor

            ; Check that we've done a good job of minimizing the spectrum between the Na D lines
              between_D = smooth( just_exosphere[ Na_D_rngs[1]+2*Na_D_rng:Na_D_rngs[2]-2*Na_D_rng,*], smooth_width, /edge_truncate )
              mean_vals = abs([ mean(between_D[*,min(Y_Mercury):max(Y_Mercury)]), $ ; mean value across Mercury's disk
                                mean(between_D[*,0:min(Y_Mercury)]), $              ; mean value below mercury's disk
                                mean(between_D[*,max(Y_Mercury):*]) ])              ; mean value above Mercury's disk
                                
            ; find the best "ref_factor" if we haven't                            
              if mean_vals[0] gt 2.*mean(mean_vals[1:2]) then begin                 ; IFF on-disk is twice as bright as off-disk
                  nfacs = (maxfac - minfac)/dfac + 1                                ; # of factors to consider
                  factors = cgscalevector(indgen(nfacs),minfac,maxfac)
                  fac_diffs = fltarr(nfacs)                                         ; keep track of the differences in mean_vals for each factor
                  for ifac = 0, nfacs-1 do begin
                      Just_exosphere = frame - scaled_reference*factors(ifac)
                      between_D = smooth( just_exosphere(Na_D_rngs[1]+2*Na_D_rng:Na_D_rngs[2]-2*Na_D_rng,*), smooth_width )
                      mean_vals = abs([ mean(between_D(*,min(Y_Mercury):max(Y_Mercury))), $ ; mean value across Mercury's disk
                                        mean(between_D(*,0:min(Y_Mercury))), $              ; mean value below mercury's disk
                                        mean(between_D(*,max(Y_Mercury):*)) ])              ; mean value above Mercury's disk
                      fac_diffs(ifac) = abs( mean_vals[0] - mean(mean_vals(1:2)) )            
                  endfor ;ifac
                  best_fac = factors(( where( fac_diffs eq min(fac_diffs) ) )[0])
                  print, i, best_fac
                  just_exosphere = frame - scaled_reference * best_fac
              endif else begin
                  best_fac = 1.
                  just_exosphere = frame - scaled_reference * best_fac
              endelse
              surf_refl_spectra_cube[ *, *, i] = frame                        ; write each dark subtracted and flat divded spectra
              exosphere_spectra_cube[ *, *, i] = Just_exosphere               ; write each residual exosphere only spectra 
            
            ; Inspection 
              wset, 0
                cgimage, bytscl(frame, 0, 200), /axes                         ; plot the original spectral frame
                for iNa = 0, 3 do cgoplot, [1,1]*Na_D_rngs[iNa], !Y.Crange    ; show the Na exclusion zones
              wset, 1
                if do_smooth then cgimage, bytscl(scaled_reference, 0, 200), /axes, title='smoothed and scaled by ' + strtrim(best_fac) $
                             else cgimage, bytscl(scaled_reference, 0, 200), /axes, title='scaled by ' + strtrim(best_fac)
                for iNa = 0, 3 do cgoplot, [1,1]*Na_D_rngs[iNa], !Y.Crange
              wset, 2
                cgimage, bytscl(Just_exosphere, 0, 200), /axes                ; plot the Mercury Na emission after subtracting disk
                for iNa = 0, 3 do cgoplot, [1,1]*Na_D_rngs[iNa], !Y.Crange
              wset, 3
                img_frame = shift(reform(img_cube[*,*,i]), [Imaging_shift_array[i,*]]); show frame just for visual reference
                cgimage, bytscl( img_frame, 0, 3500 )
                cgtext, 1, 1, 'Frame ' + strtrim(1+i)

        endfor ;i (loop over frames)
        save, exosphere_spectra_cube, surf_refl_spectra_cube, filename = outdir + 'spectral_cubes.sav'
    endif ;night=='15'
    beep
endif

; =====================================================================================================================
; Part 3 : extract the 1D exosphere spectra.
; =====================================================================================================================
if part eq 3 or part eq 99 then begin 
  if part ne 99 then restore, outdir + 'spectral_cubes.sav'
  s = size(exosphere_spectra_cube, /dimensions) 
  
  window, 0, xpos=winpos_x,         ypos=winpos_y,           xs=s[0],  ys=s[1],  title='CO-ADDED RESIDUAL SPECTRUM (INSPECTION & LINE CENTERING)'
  window, 1, xpos=winpos_x,         ypos=winpos_y+s[1]+40,   xs=s[2],  ys=s[1],  title='EXTRACTED D2 BRIGHTNESS ALONG THE SLIT (Y) AS A TIME SERIES (X)'
  window, 2, xpos=winpos_x+s[2]+20, ypos=winpos_y+s[1]+40,   xs=s[2],  ys=s[1],  title='EXTRACTED D1 BRIGHTNESS ALONG THE SLIT (Y) AS A TIME SERIES (X)'

  ; Run MPFIT twice to find the D2 line center, and the D1 line center traces
    img            = total(exosphere_spectra_cube, 3)      ; Co-add all the spectra that have surface reflectance removed 
    dummy          = img                                   ; A duplicate dummy image for inspection purposes
    search         = 16                                    ; Pixels to search over when finding the emission line center

    D2_height = fltarr(n_elements(img[0,*])) & D1_height = fltarr(n_elements(img[0,*]))
    D2_center = fltarr(n_elements(img[0,*])) & D1_center = fltarr(n_elements(img[0,*]))
    D2_width  = fltarr(n_elements(img[0,*])) & D1_width  = fltarr(n_elements(img[0,*]))

    ; D2 Tracing 
      D2_expected_pixel = mean(Na_D_rngs[0:1])             ; Rough D2 pixel location
      for i = 0, n_elements(img[0,*]) - 1 do begin
        result = mpfitpeak(findgen(search*2. + 1), img[D2_expected_pixel-search:D2_expected_pixel+search, i], a, STATUS = STATUS, /positive) 
        if STATUS ge 1 then D2_height[i] = A[0] else D2_height[i] = !values.F_nan
        if STATUS ge 1 then D2_center[i] = A[1] else D2_center[i] = !values.F_nan
        if STATUS ge 1 then D2_width[i] = A[2] else D2_width[i] = !values.F_nan
      endfor
      y         = findgen(n_elements(img[0,*]))
      height    = smooth(D2_height, smooth_width, /edge_mirror) / s[2]
      height[0:25] = 1. & height[280:*] = 1.
      real      = where(finite(D2_center), /NULL)
      COEFF     = ROBUST_POLY_FIT(y[real], D2_center[real], 2)
      D2_trace  = poly(findgen(n_elements(img[0,*])), coeff) +D2_expected_pixel-search 
      real      = where(finite(D2_width), /NULL)
      COEFF     = ROBUST_POLY_FIT(y[real], D2_width[real], 2)
      width     = poly(findgen(n_elements(img[0,*])), coeff) 
      
    ; D1 Tracing
      D1_expected_pixel = mean(Na_D_rngs[2:3])              ; Rough D1 pixel location
      for i = 0, n_elements(img[0,*]) - 1 do begin
        result = mpfitpeak(findgen(search*2. + 1), img[D1_expected_pixel-search:D1_expected_pixel+search, i], a, STATUS = STATUS, /positive)
        if STATUS ge 1 then D1_height[i] = A[0] else D1_height[i] = !values.F_nan
        if STATUS ge 1 then D1_center[i] = A[1] else D1_center[i] = !values.F_nan
        if STATUS ge 1 then D1_width[i] = A[2] else D1_width[i] = !values.F_nan
      endfor
      y         = findgen(n_elements(img[0,*]))
      height    = smooth(D1_height, smooth_width, /edge_mirror) / s[2]
      height[0:25] = 1. & height[280:*] = 1.
      real      = where(finite(D1_center), /NULL)
      COEFF     = ROBUST_POLY_FIT(y[real], D1_center[real], 2)
      D1_trace  = poly(findgen(n_elements(img[0,*])), coeff) + D1_expected_pixel-search
      real      = where(finite(D1_width), /NULL)
      COEFF     = ROBUST_POLY_FIT(y[real], D1_width[real], 2)
      width     = poly(findgen(n_elements(img[0,*])), coeff)  
      
  ; Inspection
    dummy[round(D2_trace), y] = max(dummy) 
    dummy[round(D1_trace), y] = max(dummy) 
    wset, 0
    tv, bytscl(dummy, 0, 100000)

  ; Setup MPFIT Gaussian parameters
    parinfo = replicate( {value:0.D, fixed: 0b, limited: [0b,0b], limits: dblarr(2) }, 3)
    parinfo[0].limited = 1b                                 ; limit line amplitude 
    parinfo[0].limits = [0.1, 1000. ]                       ; positive amplitudes only
    parinfo[1].limited = 1b                                 ; limit center
    parinfo[1].limits = [14., 18.]                          ; line center in pixels
    parinfo[2].limited = 1b                                 ; limit line sigma   
    parinfo[2].limits = [1., 5.]                            ; limit sigma width in pixels      

  ; Define the arrays we'll write to:
    D2_brightness = fltarr(s[2], s[1]) & D2_err_brightness = fltarr(s[2], s[1])
    D2_linewidth  = fltarr(s[2], s[1]) & D2_err_linewidth = fltarr(s[2], s[1])
    D1_brightness = fltarr(s[2], s[1]) & D1_err_brightness = fltarr(s[2], s[1])
    D1_linewidth  = fltarr(s[2], s[1]) & D1_err_linewidth = fltarr(s[2], s[1])
  
  Print, 'Fitting exosphere D2 in every spatial bin of every frame (slow, ~200 frames/min)...'
  A = [80.,D2_trace[0],1.5,0.] ; Initial guess for fitting
  for n = 0, s[2]-1 do begin
    if do_smooth then img = smooth( exosphere_spectra_cube[*,*,n], [2, smooth_width] ) $
                 else img = exosphere_spectra_cube[*,*,n]
    err_img = sqrt( abs(img) )     
    for i = 0, s[1]-1 do begin
       a[0] = height[i]
       a[1] = D2_trace[i] - D2_expected_pixel + search
       a[2] = Width[i]
       
    ;               cgplot, findgen(search*2. + 1), img[D2_trace[i]-search:D2_trace[i]+search, i]
       result = mpfitpeak(findgen(search*2. + 1), img[D2_trace[i]-search:D2_trace[i]+search, i], A, PERROR = Err_A, $
                          error = err_img[D2_trace[i]-search:D2_trace[i]+search, i], /positive, /NAN, STATUS = STATUS, PARINFO = PARINFO, NTERMS =3)
                          ;               cgoplot, findgen(search*2. + 1), result, color = 'red'
                          ;               stop
       if STATUS ge 1 then begin
          D2_brightness[n,i] = A[0]*A[2]*SQRT(2*!DPI) 
          D2_err_brightness[n,i] = D2_brightness[n,i] * sqrt( (err_A[0]/A[0])^2. + (err_A[2]/A[2])^2. )
          D2_linewidth[n,i] = dispersion_Velocity*sqrt((2.355*A[2])^2 - (2.45)^2) ;FWHM KM/S instrumental measured width is 2.45 pixels 
          D2_err_linewidth[n,i] = dispersion_Velocity*2.355*Err_A[2]
       endif else begin
          D2_brightness[n,i] = !values.F_nan
          D2_linewidth[n,i] = !values.F_nan
       endelse 
    endfor
  endfor 
   
  Print, 'Fitting exosphere D1 in every spatial bin of every frame (slow, ~200 frames/min)...'
  A = [80.,D1_trace[0],1.5,0.] ; Initial guess for fitting
  for n = 0, s[2]-1 do begin
    if do_smooth then img = smooth( exosphere_spectra_cube[*,*,n], [2, smooth_width] ) $
                 else img = exosphere_spectra_cube[*,*,n]
    err_img = sqrt( abs(img) )
    for i = 0, s[1]-1 do begin
      a[0] = height[i]
      a[1] = D1_trace[i] - D1_expected_pixel + search
      a[2] = Width[i]

      ;               cgplot, findgen(search*2. + 1), img[D1_trace[i]-search:D1_trace[i]+search, i]
      result = mpfitpeak(findgen(search*2. + 1), img[D1_trace[i]-search:D1_trace[i]+search, i], A, PERROR = Err_A, $
        error = err_img[D1_trace[i]-search:D1_trace[i]+search, i], /positive, /NAN, STATUS = STATUS, PARINFO = PARINFO, NTERMS =3)
      ;               cgoplot, findgen(search*2. + 1), result, color = 'red'
      ;               stop
      if STATUS ge 1 then begin
        D1_brightness[n,i] = A[0]*A[2]*SQRT(2*!DPI)
        D1_err_brightness[n,i] = D1_brightness[n,i] * sqrt( (err_A[0]/A[0])^2. + (err_A[2]/A[2])^2. )
        D1_linewidth[n,i] = dispersion_Velocity*sqrt((2.355*A[2])^2 - (2.45)^2) ;FWHM KM/S instrumental measured width is 2.45 pixels
        D1_err_linewidth[n,i] = dispersion_Velocity*2.355*Err_A[2]
      endif else begin
        D1_brightness[n,i] = !values.F_nan
        D1_linewidth[n,i] = !values.F_nan
      endelse
    endfor
  endfor 
   
  wset, 1
    cgimage, D2_brightness, minvalue = 0, maxvalue = 1000, /axes, /keep_aspect, xtitle = 'frame number', ytitle = 'D2 Brightness along slit' 
  wset, 2
    cgimage, D1_brightness, minvalue = 0, maxvalue = 1000, /axes, /keep_aspect, xtitle = 'frame number', ytitle = 'D1 Brightness along slit'
    
  save, D2_brightness, D2_linewidth, D1_brightness, D1_linewidth, D2_Trace, D1_Trace, filename = outdir + 'brightness.sav'
  beep
endif  

; =====================================================================================================================
; Part 4 : Extract the 1D surface reflectance, and flux calibrate the data
; =====================================================================================================================
if part eq 4 or part eq 99 then begin
  if part ne 99 then restore, outdir + 'spectral_cubes.sav' ; We'll need the cube of (co-aligned?) surface reflectance spectra
  if part ne 99 then restore, outdir + 'brightness.sav'     ; We'll need the D1 & D2 Traces to figure out which pixel the calibration wavelength is located at.
  if part ne 99 then restore, outdir + 'shift_array.sav'    ; contains shift_arrays and aligned_cubes of both channels, sharpness_metric & correl 
  s  = size(surf_refl_spectra_cube , /dimensions)
  
  ; Reconstruct an image from the spectral channel at the calibration wavelength
    ; Find the pixel location of the calibration wavelength
      Cal_wavelength_trace    = fltarr(N_elements(D2_Trace))
      for i = 0, N_elements(D2_Trace)-1 do Cal_wavelength_trace[i] = interpol( [D2_Trace[i], D1_trace[i]], [Na_D2_Wl, Na_D1_Wl], Flux_Cal_Wavelength)
      Cal_wavelength_trace = round(Cal_wavelength_trace)
    
    ; As in part 3 for sodium, extact brightness slices. Use the surface reflectance cube, the exosphere cube already has the reflectance component removed
      DN_per_A  = fltarr(s[2], s[1])
      for n = 0, s[2]-1 do begin
        if do_smooth then img = smooth( surf_refl_spectra_cube[*,*,n], [2, smooth_width] ) $
                     else img = surf_refl_spectra_cube[*,*,n]
        for i = 0, s[1]-1 do begin              
          DN_per_A[n, i] = TOTAL(img[Cal_wavelength_trace[i]-(0.5/dispersion):Cal_wavelength_trace[i]+(0.5/dispersion), i]) 
        endfor
      endfor
      
  ; Reconstruct Calibration Image: Place this time series of 1D spectral slices along the slit into a (time averaged) 2D image
    flux_cube         = fltarr(s)   ; This datacube will hold the image of the slit for each frame, slit is extacted at the calibration wavelength over 1 Angstrom 
    flux_cube_weights = fltarr(s)   ; This datacube will hold the weights from the point spread function of the slit for each frame 
    Na_D1_cube        = fltarr(s)   ; This datacube will hold the image of the slit for each frame, slit is extacted at exosphere Na D1
    Na_D2_cube        = fltarr(s)   ; This datacube will hold the image of the slit for each frame, slit is extacted at exosphere Na D2
    
    ;    ; Get the point spread function to do!!!
    ;      big_array = fltarr(1024, 1024, 3)
    ;      for i = 0, 2 do  big_array[*,*,i] = MRDfits('D:\DATA\Perkins\Perkins RIPS - March 2018\15\Neon_NaSpectra_624slitwidth_NaND1imaging'+strcompress(i+1, /remove_all)+'.fits', 0)
    ;      LSF = median(big_array, dimension = 3)
    ;      LSF = reform(LSF[sps[0]:sps[1],sps[2]:sps[3], *])
    
    ; Use a flatfield to get the PSF, then adjust it for the different platescale in the spectral channel. HACK THIS NEGLECTS ANY (PROBABLY SMALL) FOCUS DIFFERENCES BETWEEN THE IMAGING AND SPECTRAL CHANNELS. 
      iDark = MRDFITS(dark_dir + dark_file, 0, header, /fscale, /silent )   ; read in dark
      dark  = idark[ims[0]:ims[1],ims[2]:ims[3]]
      dark  = sigma_filter( dark, width, N_sigma=thresh)
      big_array = fltarr(1024, 1024, N_flats)
      for i = 0, N_flats-1 do big_array[*,*,i] = MRDFITS(flat_dir + flat_file + strcompress(i+1, /remove_all) + '.fits', 0, header, /fscale, /silent ) ; read flat
      iflat = median(big_array, dimension = 3)  ; combine flat flat
      flat  = ( iflat[ims[0]:ims[1],ims[2]:ims[3]] - dark ) / max( iflat[ims[0]:ims[1],ims[2]:ims[3]] - dark ) ; Crop to just imaging portion
      flat  = flat + (1. - median(flat))                                    ; this normalizes such that the median of the flat is 1
      profile = total(flat[*, mean(ims[2:3])-10 : mean(ims[2:3])+10], 2)
      x = findgen(N_elements(profile))
      keep = where(x lt 220 or x gt 250)       
      coeffs = POLY_FIT( x[keep], profile[keep], 5 )
      PSF = GAUSSFIT( X, poly(x, coeffs)-profile, PSF_coeffs, nterms = 3 ) 
      window, 4, title = 'SPATIAL PSF ESTIMATE: Green fit will be used for weighting scheme when placing the slit'
        cgplot, profile
        cgplot, x[keep],  profile[keep], color = 'red', /overplot
        cgplot, poly(x, coeffs), color = 'blue', /overplot
        cgplot, poly(x, coeffs)-profile, /overplot
        cgplot, x, PSF, color = 'green', /overplot
      PSF_1D = GAUSSIAN_FUNCTION( PSF_Coeffs[2]/plate_scale_ratio, /NORMALIZE)  ; The spectral channel has a smaller plate scale than the imaging channel, correct for this here
      PSF_Size = N_elements(PSF_1D)
      spatial_weights = rebin(PSF_1D, psf_size, s[1])

    home = [241,175] ;Hack: Needed for now. Do it right!
       
    for i = 0, s[2]-1 do begin; loop over frames / time. Note that all spectral cubes have already been aligned in -Y/plate_scale_ratio. Only +X/plate_scale_ratio is needed  
      flux_cube[*,*,i]         = rebin(DN_per_A[i,*], s[0], s[1])
      Na_D2_cube[*,*,i]        = rebin(D2_Brightness[i,*], s[0], s[1])
      Na_D1_cube[*,*,i]        = rebin(D1_Brightness[i,*], s[0], s[1])
      weighting_frame          = fltarr(s[0], s[1]) 
      weighting_frame[home[0] + Spectra_shift_array[i,0] - PSF_Size/2.: home[0] + Spectra_shift_array[i,0] + PSF_Size/2.-1,*] = spatial_weights      
      flux_cube_weights[*,*,i] = weighting_frame
      cgimage, flux_cube_weights[*,*,i] * flux_cube[*,*,i]
      wait, 0.04
    endfor    
stop
    ; now take the geometric mean using the weighting
      Calib_img    = total(flux_cube*flux_cube_weights, 3, /NAN)  / total(flux_cube_weights, 3)
      Na_D2_img    = total(Na_D2_cube*flux_cube_weights, 3, /NAN) / total(flux_cube_weights, 3)
      Na_D1_img    = total(Na_D1_cube*flux_cube_weights, 3, /NAN) / total(flux_cube_weights, 3)
      some_valid_data       = where(finite(Calib_img), complement = no_brightness_information)
      Calib_img[no_brightness_information] = 0.
      Na_D2_img[no_brightness_information] = 0. 
      Na_D1_img[no_brightness_information] = 0.
   
    ; Gather top 10%
      Best_Fraction_of_Images = 0.50
      ranked = reverse(SORT(correl)) ; Rank highest to lowest correlations
      cut_index = fix(s[2]*Best_Fraction_of_Images)
   
    ; now take the geometric mean of Just the lucky ones
      Lucky_Calib_img    = total(flux_cube[*,*,ranked[0:cut_index]]*flux_cube_weights[*,*,ranked[0:cut_index]], 3, /NAN)  / total(flux_cube_weights[*,*,ranked[0:cut_index]], 3)
      Lucky_Na_D2_img    = total(Na_D2_cube[*,*,ranked[0:cut_index]]*flux_cube_weights[*,*,ranked[0:cut_index]], 3, /NAN) / total(flux_cube_weights[*,*,ranked[0:cut_index]], 3)
      Lucky_Na_D1_img    = total(Na_D1_cube[*,*,ranked[0:cut_index]]*flux_cube_weights[*,*,ranked[0:cut_index]], 3, /NAN) / total(flux_cube_weights[*,*,ranked[0:cut_index]], 3)
      some_valid_data    = where(finite(Lucky_Calib_img), complement = no_brightness_information)
      Lucky_Calib_img[no_brightness_information] = 0.
      Lucky_Na_D2_img[no_brightness_information] = 0.
      Lucky_Na_D1_img[no_brightness_information] = 0.    

  ; Generate a Hapke Model of Mercury at the observation time/geometry. Match its plate scale to the RIPS spectral channel
    header = headfits(outdir + 'imaging_cube.fits')  ; Grab a file header for the timestamp
    Surface_Flux_Calibration, Body = 'Mercury', Timestamp = sxpar(header, 'DATE'), Wavelength = Flux_Cal_Wavelength, Output_image = outdir + 'Hapke_For_Perkins.eps', /align_celestial_north, $
      Hapke_Platescale = Hapke_Platescale, MR_per_A = MR_per_A, Make_picture = Make_picture 
    Flux_Cal_Img_Size = size(MR_per_A, /dim)
    Flux_Cal_Spectral = congrid( MR_per_A, Flux_Cal_Img_Size[0] * Hapke_Platescale/Spectral_platescale, Flux_Cal_Img_Size[1] * Hapke_Platescale/Spectral_platescale, cubic = -0.5)
    effective_seeing = 2.7                                  ; "effective" seeing FWHM arcseconds
    effective_seeing = effective_seeing / 2.355             ; convert FWHM to sigma
    effective_seeing = effective_seeing/Spectral_platescale ; convert to Gaussian sigma in pixels
    visual_inspection = convol(rot(Flux_Cal_Spectral, 203, MISSING=0., /interp), GAUSSIAN_FUNCTION([effective_seeing,effective_seeing]), /EDGE_ZERO, /normalize)
    
    junk_s = max(Calib_img, loc_spec)
    spec_centroid = array_indices(Calib_img, loc_spec)
    junk_h = max(visual_inspection, loc_hapke)
    hapke_centroid = array_indices(visual_inspection, loc_hapke)
    
    sh = size(visual_inspection, /dimensions)              ; Hapke image dimensions
    si = size(aligned_imaging_cube, /dimensions)           ; Imaging channel image dimensions
    
    ; Need to ratio the blurred Hapke image and the data-generated calibration image over the same to determine sensitivity
      ;MR_per_DN = TOTAL(visual_inspection[hapke_centroid[0]-30:hapke_centroid[0]+30, $
      ;                                    hapke_centroid[1]-40:hapke_centroid[1]+40]) $  
      ;          / TOTAL(Calib_img[spec_centroid[0]-30:spec_centroid[0]+30, $
      ;                                     spec_centroid[1]-40:spec_centroid[1]+40]) ; this should work if seeing is correct, and yet...
                                           
    MR_per_DN = max(visual_inspection[hapke_centroid[0]-30:hapke_centroid[0]+30, $
                                           hapke_centroid[1]-40:hapke_centroid[1]+40]) $
              / max(Calib_img[spec_centroid[0]-30:spec_centroid[0]+30, $
                                           spec_centroid[1]-40:spec_centroid[1]+40]) ; ... and yet this looks better
     
    ; Assemble the imaging channel results 
      Imaging_Ch = total(aligned_imaging_cube, 3)                                                      & Lucky_Imaging_Ch = total(aligned_imaging_cube[*,*,ranked[0:cut_index]], 3)
      Imaging_Ch = rotate(Imaging_Ch, 7)                                                               & Lucky_Imaging_Ch = rotate(Lucky_Imaging_Ch, 7)
      Imaging_Ch = congrid(Imaging_Ch, si[0]/plate_scale_ratio, si[1]/plate_scale_ratio, cubic = -0.5) & Lucky_Imaging_Ch = congrid(Lucky_Imaging_Ch, si[0]/plate_scale_ratio, si[1]/plate_scale_ratio, cubic = -0.5)
      
    ; Calibrate the data in absolute flux, MegaRayleighs (line emission), MegaRayleighs/Angstrom (continuum) 
      Calib_img       = Calib_img * MR_per_DN
      Na_D2_img       = Na_D2_img * MR_per_DN
      Na_D1_img       = Na_D1_img * MR_per_DN
      Lucky_Calib_img = Lucky_Calib_img * MR_per_DN
      Lucky_Na_D2_img = Lucky_Na_D2_img * MR_per_DN
      Lucky_Na_D1_img = Lucky_Na_D1_img * MR_per_DN
                          
  ; Setup inspection windows
    window, 0, xpos=winpos_x,             ypos=winpos_y,         xs=s[0],  ys=s[1],  title='EXTRACTED DN ALONG IN A 1 ANGSTROM BAND (DN/A)'
    window, 1, xpos=winpos_x+s[0]+20,     ypos=winpos_y,         xs=s[0],  ys=s[1],  title='EXTRACTED NA D2 DN'
    window, 2, xpos=winpos_x+s[0]+20,     ypos=winpos_y+s[1]+40, xs=s[0],  ys=s[1],  title='EXTRACTED NA D1 DN'
    window, 3, xpos=winpos_x,             ypos=winpos_y+s[1]+40, xs=sh[0], ys=sh[1], title='HAPKE MODEL (MR/A)'
    window, 4, xpos=winpos_x+sh[0]+20,    ypos=winpos_y+s[1]+40, xs=si[0]/plate_scale_ratio, ys=si[1]/plate_scale_ratio, title='ALIGNED & COADDED IMAGING CHANNEL, MATCHED TO SPECTRAL CHANNEL PLATESCALE'
    window, 5, xpos=winpos_x,             ypos=winpos_y+2.*s[1]+40, xs=s[0],ys=s[1], title='LUCKY EXTRACTED DN ALONG IN A 1 ANGSTROM BAND (DN/A)'
    loadct, 22
    
    Continuum_colorbar_top = 40. ; MR/A
    wset, 0
      cgimage, Calib_img, minvalue = 0, maxvalue = Continuum_colorbar_top 
      cgColorbar, minrange = 0., maxrange = Continuum_colorbar_top, title = cgsymbol('Sigma')+'1'+cgsymbol('Angstrom') +' at '+ strcompress(Flux_Cal_Wavelength)+cgsymbol('Angstrom') + ' (MegaRayleighs /'+ cgsymbol('Angstrom')+')' 
    wset, 1 
      cgimage, Na_D2_img*1000., minvalue = 0, maxvalue = 2000.
      cgColorbar, minrange = 0, maxrange = 2000., title = 'Na D2 (KiloRayleighs)'
    wset, 2
      cgimage, Na_D1_img*1000., minvalue = 0, maxvalue = 1200.
      cgColorbar, minrange = 0, maxrange = 1200., title = 'Na D1 (KiloRayleighs)'
    wset, 3
      cgimage, visual_inspection, minvalue = 0, maxvalue = Continuum_colorbar_top ;Blurred_Hapke
    wset, 4
      cgimage, Imaging_Ch, /keep_aspect
      cgimage, Lucky_Imaging_Ch, /keep_aspect
    wset, 5
      cgimage, lucky_Calib_img, minvalue = 0, maxvalue = 40.
      cgColorbar, minrange = 0., maxrange = Continuum_colorbar_top, title = cgsymbol('Sigma')+'1'+cgsymbol('Angstrom') +' at '+ strcompress(Flux_Cal_Wavelength)+cgsymbol('Angstrom') + ' (MegaRayleighs /'+ cgsymbol('Angstrom')+')'   

    save, Lucky_Imaging_Ch, Imaging_Ch, Lucky_Calib_img, Calib_img, Lucky_Na_D2_img, Lucky_Na_D1_img, Na_D2_img, Na_D1_img, visual_inspection, Flux_Cal_Spectral, $
          filename = outdir + 'images_to_plot.sav' 
          stop
endif

; =====================================================================================================================
; Part 5 : Plot all images and and make movies of the exosphere.
; =====================================================================================================================
if part eq 5 or part eq 99 then begin
  if part ne 99 then restore, outdir + 'brightness.sav'            ; contains brightness, linewidth of both D2 & D1 lines
  
  cgPS_Open, filename = outdir+'Mercury_Perkins.eps', /ENCAPSULATED, xsize = 6.5, ysize = 4
  !P.font=1
  device, SET_FONT = 'Helvetica Bold', /TT_FONT
  !p.charsize = 1.2
  
  pos = cgLayout([2,1], OXMargin=[5.8,1], OYMargin=[5,1], XGap=1, YGap=10)
  thick = 4.
  
  cgimage, lucky_Calib_img, ytitle = 'Sodium D!L2!N (kR)', xtitle = 'Tangent Altitude (km)', /ynozero, color = 'blue', $ ;Hacked brightness to match Sprague et al 2012
    XStyle=9, Position=[POS[0,0], POS[1,0], POS[2,0], 0.52], Thick=thick, xthick = thick, ythick = thick, yr = [0.,1.5], xr = minmax(tangent_alt)
  cgPS_close
endif
stop



;      ; now we need to align the Hapke frame to the Spectral channel's image extracted dn along in a 1 Angstrom band
;        ; start with the brightness centroid as a guess
;          match_hapke_to_calib                         = FLTARR(s[0],s[1])
;          match_hapke_to_calib[0:sh[0]-1, 0:sh[1]-1]   = visual_inspection
;          junk_s = max(weighted_Calib_img, loc_spec)
;          spec_centroid = array_indices(weighted_Calib_img, loc_spec)
;          junk_h = max(match_hapke_to_calib, loc_hapk)
;          hapk_centroid = array_indices(match_hapke_to_calib, loc_hapk)
;          match_hapke_to_calib = shift(match_hapke_to_calib, spec_centroid[0] - hapk_centroid[0], spec_centroid[1] - hapk_centroid[1] )
;
;        ; fine tune the alignment
;          ;CORREL_OPTIMIZE, weighted_Calib_img, match_hapke_to_calib, xoffset_optimum, yoffset_optimum, /NUMPIX
;          ;match_hapke_to_calib = shift(match_hapke_to_calib, xoffset_optimum, yoffset_optimum )
;
;        window, 5, xpos=winpos_x+900, ypos=winpos_y+500,  xs=s[2],  ys=s[1],  title='hapke fake'
;        loadct, 22
;        cgimage, MR_per_DN
;        print, 'Stddev hapke match', stddev(MR_per_DN, /NAN)
;;        loadct, 0















;------------------------------------Old Code Below-----------------------------------




; =====================================================================================================================
; Part 5 : put everything together and build images and movies of the exosphere.
; =====================================================================================================================
if part eq 6 or part eq 99 then begin  
   if part ne 99 then restore, outdir + 'brightness.sav'            ; contains brightness, linewidth of both D2 & D1 lines
   if part ne 99 then restore, outdir + 'shift_array.sav'           ; contains shift_array, sharpness_metric, and aligned_imaging_cube
   if part ne 99 then imaging_cube = MRDFITS(outdir + 'imaging_cube.fits', 0, header, /fscale, /silent )
   s = size(imaging_cube, /dimensions)
   
   ; HACK! IF DIFFERENT FILES USE DIFFERENT INTEGRATION TIMES THIS WILL NOT BE ACCOUNTED FOR! 
     integration = sxpar(headfits(Mercury_dir + Mercury_file[0]), 'EXPOSURE') ; integration time for individual frames

   brightness = congrid(D2_brightness, s[2], s[1]) ;HACK NEED PROPER SCALING OF THE SPECTRUM TO IMAGING PLATESCALES
   ;brightness = congrid(brightness, s[2], round(S[2]*plate_scale_ratio))
;stop
   
   home = [341,238] ;Mercury centroid of the last frame number 500, to which all the other frames are aligned
   home = [236,165]
   home = [231,165]
   home = [241,175]
;   home = [131,102]  ; trying home values for a smaller cutout of the imaging and spectral frames to speed things up
   
;   big_cube = fltarr(s[0], 2*s[1], s[2])
   big_cube_spec = fltarr(s[0], s[1], s[2])
;   big_cube(*,*,*) = !Values.F_NAN
   big_cube_spec(*,*,*) = !Values.F_NAN
   img_cube = fltarr(s[0], s[1])
   img_count = fltarr(s[0], s[1])                ; array to track number of "real" pixels (ie not obscured by slit)
   spec_cube = fltarr(s[0], s[1])
   spec_dwell = fltarr(s[0], s[1])                                        ; keep track of the # of frames in each pixel
   spec_temp  = fltarr(s[0], s[1]) 
   big_frame_total = fltarr(s[0], 2*s[1])
   window, 0, xs = s[0], ys = 2*s[1], xpos=winpos_x, ypos=winpos_y
;   cgdisplay, s[0], s[1]*2, xpos=winpos_x, ypos=winpos_y, wid=0
;   window, 1, xs = s[0], ys = s[1], xpos=winpos_x[2], ypos=winpos_y[2]
   xs_spec = spec_extraction[1] - spec_extraction[0] + 1
   ys_spec = spec_extraction[3] - spec_extraction[2] + 1
   xs_img= img_extraction[1] - img_extraction[0] + 1
   ys_img = img_extraction[3] - img_extraction[2] + 1
   window, 1, xs = xs_img*movie_scale, ys = ys_img*movie_scale, xpos=winpos_x+s[0]+20, ypos=winpos_y   
;   cgdisplay, xs_spec*movie_scale, ys_spec*movie_scale, xpos=winpos_x+40, ypos=winpos_y, wid=1
;   window, 2, xs = s[0], ys = s[1], xpos=winpos_x[2], ypos=winpos_y[2]+s[1]+35
   window, 2, xs = xs_img*movie_scale, ys = ys_img*movie_scale, xpos=winpos_x+s[0]+20, ypos=winpos_y+ys_img*movie_scale+40
   window, 3, xs = s[0], ys = 2*s[1], xpos=winpos_x, ypos=winpos_y+2*s[1]+40

   ; Some video related variables
   mpgFilename = outdir + 'Mercury Na image.mp4'
   mpgFilename2 = outdir + 'Mercury disk.mp4'
   mpgFilename3 = outdir + 'Mercury scan.mp4'
   video = IDLffVideoWrite(mpgFilename, format='mp4')
   video2 = IDLffVideoWrite(mpgFilename2, format='mp4')
   video3 = IDLffVideoWrite(mpgFilename3, format='mp4')
   framerate = 30
   framedims = [xs_img*movie_scale, ys_img*movie_scale]
   framedims2 = [xs_img*movie_scale, ys_img*movie_scale]
   framedims3 = [s[0], 2*s[1]]
   stream  = video.AddVideoStream(framedims[0], framedims[1], framerate)
   stream2 = video2.AddVideoStream(framedims2[0], framedims2[1], framerate)
   stream3 = video3.AddVideoStream(framedims3[0], framedims3[1], framerate)

   good_frames = where(sharpness_metric gt (mean(sharpness_metric)-stddev_cutoff*stddev(sharpness_metric)), complement=bad_frames)  ; these are the good frames we'll use below (in terms of good seeing)
   ngood_shifts = n_elements(where(shift_array(*,0) ne -666))             ; these are the good frames (in terms of a reasonable shift_array value)

;   for ig = 253, n_elements(good_frames)-1 do begin                         ; loop over good frames
   for ig = 0, n_elements(good_frames)-1 do begin                         ; loop over good frames
    
        i = good_frames(ig)
        if shift_array[i,0] eq -666 then goto, skipbad
        
;   for i = 0, 20 do begin;s[2]-1 do begin
        big_frame = fltarr(s[0], 2*s[1])
        empty_frame = make_array(s[0], s[1], /float, value = !values.F_NaN)
        empty_frame[home[0]+shift_array[i,0]-1:home[0]+shift_array[i,0]+1,*] = [brightness[i, *], brightness[i, *], brightness[i, *]] ;Make it bigger?
        big_frame[*, s[1]:*] = rotate(empty_frame, 7) ;* 5. ;the factor of 5 is a HACK FOR NOW
        big_frame[*, 0:s[1]-1] = reform(aligned_imaging_cube[*,*,i])
        big_frame_total = big_frame_total + big_frame
;        stop
        wset, 0
        cgloadct, ct, /silent
        big_frame_scaled = big_frame / (total(big_frame,/nan)/7.427e7)
        cgimage, bytscl(big_Frame_scaled, 0, 10000)   ; 7e7 is hack for total of all frames
        cgtext, 0.02, 0.96, 'Extracted Na emission', charsize=0.5*movie_Scale, col='black', /normal, charthick=1.25*movie_scale
        cgtext, 0.02, 0.96, 'Extracted Na emission', charsize=0.5*movie_Scale, col='white', /normal, charthick=0.5*movie_Scale
        cgtext, 0.02, 0.45, 'IMAGING CHANNEL', charsize=0.5*movie_Scale, col='black', /normal, charthick=1.25*movie_Scale
        cgtext, 0.02, 0.45, 'IMAGING CHANNEL', charsize=0.5*movie_Scale, col='white', /normal, charthick=0.5*movie_Scale
        cgtext, .75, .97, 'Frame # ' + strtrim(i+1), /data, col='black', /normal, charthick=1.5
        movie_png = cgsnapshot(filename = outdir + 'Mercury scan.png', /nodialog)
        image = read_png(outdir + 'Mercury scan.png')
        void = video3 -> Put(stream3, image)
               
;        big_cube[*,*,i] = big_frame
        spec_cut_for_big_cube_spec = reform(big_frame(*,s[1]:*))
        big_cube_spec[*,*,i] = smooth( spec_cut_for_big_cube_spec, [0,8], /nan )
        img_cube = img_cube + reform(big_frame(*,0:s[1]))
        spec_frame = reform(big_frame(*,s[1]:*))
        spec_real = where(finite(spec_frame),complement=spec_bkg)         ; find elements with real values, otherwise assign to background
        spec_summed = total(finite(spec_frame),2)                         ; by summing rows we can find the columns corresponding to the slit location
        slit_columns = where(spec_summed eq max(spec_summed))             ; columns corresponding to slit location
        if n_elements(slit_columns) le 3 then begin                       ; sometimes we have a bad frame, in which case slit_columns thinks everything is beneath the slit
          img_count = img_count + 1.                                        ; increment "real" pixels
          img_count(slit_columns,*) = img_count(slit_columns,*) - 1.        ; but also account for the pixels covered by the slit
;        spec_dwell(spec_real) = spec_dwell(spec_real) + integration       ; dwell time
          spec_dwell(min(slit_columns):max(slit_columns),*) = spec_dwell(min(slit_columns):max(slit_columns),*) + integration
        endif
        spec_frame(spec_bkg) = 0d
        spec_frame = smooth(spec_frame, [0,8])                            ; HACK (just in the sense that maybe we don't want to smooth here)
        spec_cube(spec_real) = spec_cube(spec_real) + spec_frame(spec_real)

;wset, 1
;cgimage, spec_frame
;wset, 2
;cgimage, spec_dwell
;stop

        wset, 1
        spec_dwell_num = where(spec_dwell gt 0.)                          ; find pixels with spectral information included
        ;if i lt s[2]-1 then $
        if ig lt n_elements(good_frames)-1 then $
          spec_temp(spec_dwell_num) = spec_cube(spec_dwell_num) / spec_dwell(spec_dwell_num) + spec_frame(spec_dwell_num) else $   
          spec_temp(spec_dwell_num) = spec_cube(spec_dwell_num) / spec_dwell(spec_dwell_num)  ; ie just don't plot slit contribution for last frame
        spec_frame_slit = spec_frame(spec_extraction[0]:spec_extraction[1],spec_extraction[2]:spec_extraction[3])
        spec_frame_slit_value = where(spec_frame_slit gt 0., complement=spec_frame_slit_bkg)
        spec_frame_slit_2 = spec_frame_slit
        spec_frame_slit_2(spec_frame_slit_bkg) = 0.;!Values.F_NAN
        blah = big_cube_spec(spec_extraction[0]:spec_extraction[1],spec_extraction[2]:spec_extraction[3],*)
;        blah = big_cube_spec(img_extraction[0]:img_extraction[1],img_extraction[2]:img_extraction[3],*)
;        wslit = where( total(spec_frame_slit,2) eq max(total(spec_frame_slit,2)) )
        if do_median then spec_movie_temp = median(blah, dimension=3, /even) $
                     else spec_movie_temp = mean(blah, dimension=3, /nan)
        spec_movie = bytscl(spec_movie_temp, max_spec*0.05, max_spec*.8)
;        spec_movie = bytscl(spec_temp(spec_extraction[0]:spec_extraction[1],spec_extraction[2]:spec_extraction[3]), max_spec*0.05, max_spec*.8)
;        spec_movie = congrid(spec_movie, xs_spec*movie_scale, ys_spec*movie_scale)
;        if n_elements(wslit) eq 3 then spec_movie(min(wslit)>1:max(wslit)<(n_elements(spec_movie(*,0))-1),*) = 1e3 ; hack to make slit show up
        axis_format = {XTicklen:0.0001, YTickLen:0.0001, Xthick:4, Ythick:4}
        if do_rotate then cgimage, rotate(spec_movie,5), /axes, charsize=1e-3, position=[0.0,0.0,1.,1.], axkeywords=axis_format $
                     else cgimage, spec_movie, /axes, charsize=1e-3, position=[0.0,0.0,1.,1.], axkeywords=axis_format
        if do_rotate then cgimage, rotate(spec_frame_slit_2,5), transparent=20, missing_Value=0.0, ctindex=0 $
                     else cgimage, spec_frame_slit_2,  transparent=20, missing_Value=0.0, ctindex=0
        img_cube_cut = img_cube(img_extraction[0]:img_extraction[1],img_extraction[2]:img_extraction[3])
        ;cgcontour, normalize(img_cube_cut), /onimage, levels=[0.2], label=0 ;WTF is normalize?
        cgcontour, img_cube_cut / max(img_cube_cut, /NaN), /onimage, levels=[0.2], label=0 

;if ig gt 31 then stop        
        ; Next 2 lines draw slit
;        cgoplot, [1,1]*min(wslit), !Y.Crange, line=0, thick=2
;        cgoplot, [1,1]*max(wslit), !Y.Crange, line=0, thick=2

;        cgcolorbar, /vertical, charsize=1.e-3, position=[0.09,0.1,0.13,0.8];, oob_low='white'
        cgcolorbar, /vertical, charsize=1.e-3, position=[0.9,0.15,0.95,0.85]
        cgtext, .65, .05, '# frames = ' + strtrim(ig+1), /normal, charthick=0.5*movie_Scale, font=-1, charsize=0.5*movie_Scale
        cgtext, 0.02, 0.96, 'Mercury Na', /normal, charthick=0.75*movie_Scale, charsize=0.625*movie_Scale, font=-1
        movie_png = cgsnapshot(filename = outdir + 'Mercury Na.png', /nodialog)
        image = read_png(outdir + 'Mercury Na.png')
        void = video -> Put(stream, image)
;stop
        wset, 2
;        img_movie = bytscl(img_cube(img_extraction[0]:img_extraction[1],img_extraction[2]:img_extraction[3]), max_img*0.2, max_img)
        img_count_cut = smooth(img_count(img_extraction[0]:img_extraction[1],img_extraction[2]:img_extraction[3]),6,/edge_truncate)
;        img_movie = bytscl(img_cube_cut, 0.2*max(img_cube_cut), max(img_cube_cut))
;        img_movie = congrid(img_movie, xs_img*movie_scale, ys_img*movie_scale)
        img_movie = congrid(img_cube_cut * (count_factor*(max(img_count_cut)/img_count_cut)^2.), xs_img*movie_scale, ys_img*movie_scale)
        cgloadct, 0, /silent
        cgimage, img_movie
        cgcolorbar, /vertical, charsize=1.e-3
        cgtext, .65, .05, '# frames = ' + strtrim(ig+1), /data, /normal, charthick=0.375*movie_Scale, charsize=0.5*movie_scale
        cgtext, 0.02, 0.96, 'Solar continuum', /normal, charthick=0.75*movie_Scale, charsize=0.625*movie_Scale
        movie_png = cgsnapshot(filename = outdir + 'Mercury disk.png', /nodialog)
        image = read_png(outdir + 'Mercury disk.png')
        void = video2 -> Put(stream2, image)

        wset, 3
        cgimage, bytscl(big_frame_total/double(i),0,5000)
        skipbad: ;wait, 0.02   

;loadct, 20
;cgimage, reform(aligned_imaging_cube(*,*,i)), /axes, /keep_aspect_ratio
;x = spec_dwell*0.
;x(slit_columns,*) = 1.0
;cgcontour, x, /onimage, label=0
;print, slit_columns
;stop
   endfor ;ig
   
   ; LOOP over bad frames too for comparison figure ???
   dobad = 1
   if dobad then begin
     img_cube_bad = fltarr(s[0], s[1])
     for ib = 0, n_elements(bad_frames)-1 do begin                         ; loop over bad frames
        i = bad_frames(ib)
        if shift_array[i,0] eq -666 then goto, skipbad2

        big_frame = fltarr(s[0], 2*s[1])
        empty_frame = make_array(s[0], s[1], /float, value = !values.F_NaN)
        empty_frame[home[0]+shift_array[i,0]-1:home[0]+shift_array[i,0]+1,*] = [brightness[i, *], brightness[i, *], brightness[i, *]] ;Make it bigger?
        big_frame[*, s[1]:*] = rotate(empty_frame, 7) ;* 5. ;the factor of 5 is a HACK FOR NOW
        big_frame[*, 0:s[1]-1] = reform(aligned_imaging_cube[*,*,i < ngood_shifts-1])
        big_frame_total = big_frame_total + big_frame
;        big_cube[*,*,i] = big_frame
        img_cube_bad = img_cube_bad + reform(big_frame(*,0:s[1]))

        skipbad2: ;wait, 0.02
      endfor ;ib     
    img_cube_cut_bad = img_cube_bad(img_extraction[0]:img_extraction[1],img_extraction[2]:img_extraction[3])
   endif ;dobad==1
     
   
   wset, 3
;   test = median(big_cube, dimension = 3, /even)
   ;tv, bytscl(smooth(test,[3,3], /NaN), 0, 10000)
;   cgimage, bytscl(test, 0, 4500)
   spec_dwell_num = where(spec_dwell gt 0.)                               ; find pixels with spectral information included
   spec_cube(spec_dwell_num) = spec_cube(spec_dwell_num) / spec_dwell(spec_dwell_num)  ; normalize by # of frames summed
   
   ; Now add some blank frames to the Na image movie before plotting a smoothed frame
   ;spec_raw = spec_cube(spec_extraction[0]:spec_extraction[1],spec_extraction[2]:spec_extraction[3])
   wset, 1
   spec_raw = spec_cube(spec_extraction[0]:spec_extraction[1],spec_extraction[2]:spec_extraction[3])
   spec_raw = congrid(spec_raw, xs_spec*movie_scale, ys_spec*movie_scale) ; embiggen
;   spec_movie = bytscl(spec_raw, max_spec*0.05, max_spec*.8)
   loadct, ct, /silent
   if do_rotate then cgimage, rotate(spec_movie,5), /axes, charsize=1e-3, position=[0.0,0.0,1.,1.], axkeywords=axis_format $
                else cgimage, spec_movie, /axes, charsize=1e-3, position=[0.0,0.0,1.,1.], axkeywords=axis_format
   ;cgcontour, normalize(img_cube_cut), /onimage, levels=[0.2], label=0
   cgcontour, img_cube_cut / max(img_cube_cut, /NaN), /onimage, levels=[0.2], label=0
   cgcolorbar, /vertical, charsize=1.e-3, position=[0.9,0.15,0.95,0.85]
   cgtext, .65, .05, '# frames = ' + strtrim(ig), /normal, charthick=0.5*movie_Scale, font=-1, charsize=0.5*movie_Scale
   cgtext, 0.02, 0.96, 'Mercury Na', /normal, charthick=0.75*movie_Scale, charsize=0.625*movie_Scale, font=-1
   movie_png = cgsnapshot(filename = outdir + 'Mercury Na.png', /nodialog)
   image = read_png(outdir + 'Mercury Na.png')

   for i = 0, 30 do void = video -> Put(stream, image)
   smooth_width = 4
   if do_rotate then cgimage, rotate(smooth(spec_raw,smooth_width*movie_scale,/edge_truncate), 5) $              ; HACK
                else cgimage, smooth(spec_raw,smooth_width*movie_scale,/edge_truncate)
   cgcontour, img_cube_cut / max(img_cube_cut, /NaN), /onimage, levels=[0.2], label=0
   cgtext, .65, .05, '# frames = ' + strtrim(ig), /normal, charthick=0.5*movie_Scale, font=-1, charsize=0.5*movie_Scale
   cgtext, 0.02, 0.96, 'Mercury Na (smoothed=' + strtrim(smooth_width) + ')', /normal, charthick=0.75*movie_Scale, charsize=0.625*movie_Scale, font=-1
   movie_png = cgsnapshot(filename = outdir + 'Mercury Na - smoothed.png', /nodialog)
   image = read_png(outdir + 'Mercury Na - smoothed.png')
   for i = 0, 30 do void = video -> Put(stream, image)
      
   video -> Cleanup
   video2 -> Cleanup
   video3 -> Cleanup

   ; Overplot Na contours on top of Mercury disk
   loadct, 0, /silent
   cgimage, bytscl(img_movie), /axes
   spec_smoothed = bytscl(smooth(spec_raw,smooth_width*movie_scale,/edge_truncate))
;   cgcontour, spec_smoothed, /onimage, color='black', levels=[50,100,150,175,200,225,250], xthick=3 
;   cgcontour, spec_smoothed, /onimage, color='gold', levels=[50,100,150,175,200,225,250]
   cgcontour, spec_smoothed / max(spec_smoothed, /NaN), /onimage, color='gold', levels=[0.2,0.4,0.6,0.8]
   al_legend, ['Disk','Na'], box=0, /bottom, /right, charsize=0.5*movie_scale, line=[0,0], linsize=0.2, thick=[6,2], color=['gray','gold'], textcolor='white'
   combined = cgsnapshot(filename = outdir + 'Mercury Na + disk.png', /nodialog)

   ; Generate the "bad images" Mercury disk for comparison
   if dobad then begin
     bad_img = bytscl(img_cube_cut_bad, 0.2*max(img_cube_cut_bad), max(img_cube_cut_bad))
     cgimage, bad_img
;    cgimage, img_cube_cut_bad
     cgcolorbar, /vertical, charsize=1.e-3
     cgtext, .65, .05, '# frames = ' + strtrim(ib), /data, /normal, charthick=0.375*movie_Scale, charsize=0.5*movie_Scale
     cgtext, 0.02, 0.96, 'Solar continuum', /normal, charthick=0.75*movie_Scale, charsize=0.625*movie_Scale
     bad_png = cgsnapshot(filename = outdir + 'Mercury disk - bad frames.png', /nodialog)

     ; And might as well compare the good vs bad contours while we're here...
     cgimage, img_movie, /axes
;     cgcontour, img_movie, /onimage, color='gold'
     cgcontour, img_movie / max(img_movie, /NaN), /onimage, levels=[0.2,0.4,0.6,0.8], color='gold'
;     cgcontour, bad_img, /onimage, color='red'
     cgcontour, img_cube_cut_bad / max(img_cube_cut_bad, /NaN), /onimage, levels=[0.2,0.4,0.6,0.8], color='red'
     al_legend, ['good','bad'], box=0, /bottom, /right, charsize=2, line=[0,0], linsize=0.2, thick=[4,4], color=['orange','red'], textcolor=['orange','red']
     compare = cgsnapshot(filename = outdir + 'Mercury disk - good vs bad contours.png', /nodialog)
   endif ;dobad==1
   

   
   ; Create image showing the spatial coverage of the slit position
   wset, 2
;   dwell_img = spec_dwell(spec_extraction[0]:spec_extraction[1],spec_extraction[2]:spec_extraction[3])
   dwell_img = spec_dwell(img_extraction[0]:img_extraction[1],img_extraction[2]:img_extraction[3])
   loadct, 22, /silent
;   cgloadct, 22, ncolors=max(spec_dwell)/integration, bottom=0, /silent
   cgimage, bytscl(dwell_img, 0., max(dwell_img)), /axes, title='Slit Coverage'
   cgcontour, img_cube_cut / max(img_cube_cut, /NaN), /onimage, levels=[0.2]
   cgloadct, 22, ncolors=max(spec_dwell)/integration, bottom=0, /silent  
   cgcolorbar, /vertical, maxrange=max(spec_dwell), title='Total Na data coverage (sec)', /discrete, ncolors=max(spec_dwell)/integration, position=[0.9,0.2,0.95,0.8], tickinterval=integration*4.
   combined = cgsnapshot(filename = outdir + 'Mercury slit coverage.png', /nodialog)

   beep
endif

; cgPS2Raster, 'test_2.ps', /JPEG
;   IDL> cgPS_Open, 'test_2.ps'
;   IDL> Graphic_Display
;   IDL> cgPS_Close, /PNG
print, 'Run time = ',  (SYSTIME(/SECONDS) - Start_time) / 60, ' minutes'
end



