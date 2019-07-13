function Estimate_Seeing, FWHM_arcsec ; FWHM Arcseconds <----solve for this iteratively  
  Common Seeing_Estimate, A, B, Spectral_platescale, Imaging_platescale, x_fine, y_fine, Seeing_Est_Rebin_Factor

  if FWHM_arcsec gt 5.2/Seeing_Est_Rebin_Factor then return, 0.   ; prevents runaways and associated indexing issues
  seeing_sigma      = FWHM_arcsec / (2.0*sqrt(2.0*alog(2)))       ; convert FWHM to sigma
  seeing_sigma      = seeing_sigma/Imaging_platescale             ; convert to Gaussian sigma in pixels
  Blurred_B         = convol(B, GAUSSIAN_FUNCTION([seeing_sigma, seeing_sigma]), /EDGE_ZERO, /normalize, /NAN)
  
  ;Correlation based alignments don't like NaN Alignments 
  interp_slit_A     = A
  fill_missing, interp_slit_A, !values.F_Nan, 1
  
  zstack_align_images, interp_slit_A/max(interp_slit_A), Blurred_B/max(Blurred_B), x_fine, y_fine, corr_image, corr_dims
  ;Blurred_B = shift(Blurred_B, x_fine, y_fine)
  interp_slit_A = shift( interp_slit_A, [x_fine, y_fine] )
  A_shifted = shift( A, [x_fine, y_fine] )
  CORREL = CORREL_IMAGES( Blurred_B/max(Blurred_B), interp_slit_A/max(interp_slit_A), xoffset_b = 0, yoffset_B = 0, XSHIFT =0, ySHIFT = 0) ;2-D Correlation at this alignment
  
  ; Normalize both, excluding masked region of the slit
    junk  = where(finite(A), complement = mask)
    Blurred_B[mask]  = !values.F_NaN
    A_shifted        = A_shifted/total(A_shifted, /NAN)
    Blurred_B        = Blurred_B/total(Blurred_B, /NAN)

  ;Stat  = stddev(Blurred_B-A, /NaN)
  Stat  =  1./correl ;new test for 
  ;tv, bytscl(A_Shifted - Blurred_B)
  ;print, FWHM_arcsec[0]*Seeing_Est_Rebin_Factor, stat[0], x_fine, y_fine, correl
  ;stop
return, stat ;Scalar Value for AMOEBA TO MINIMIZE
end

pro RIPS_Mercury_Perkins, PART=part, NIGHT=night
  Common Seeing_Estimate, A, B, Spectral_platescale, Imaging_platescale, x_fine, y_fine, Seeing_Est_Rebin_Factor 
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
; Define all variables
; =====================================================================================================================
Body          = 'Mercury'                                                 ; For SPICE purposes
mercury_dir   = base_dir+'15\'                                            ; directory with Mercury FITS file(s)
dark_dir      = base_dir+'14\Carl_Keep\'                                  ; directory with dark file
flat_dir      = base_dir+'15\'                                            ; directory with flat file
sky_dir       = base_dir+'15\'                                            ; directory with sky file
outdir        =  Mercury_dir + 'Processed\'                               ; output directory
Arc_files     = 'Neon_NaSpectra_624slitwidth_NaND1imaging*.fits'          ; Arc frames 
Mercury_files = '*'+strcompress(789+indgen(17), /remove_all) + '.fits'    ; SUFFIX of the Mercury kinetic series to use (skipping 788 for now because of it's at a different rotiserizer angle)
sky_file      = 'RIPS1_Thu Mar 15 2018_01.12.12_785.fits'                 ; sky frame
dark_file     = 'Dark.fits'                                               ; dark frame
flat_files    = 'Flat_NaSpectra_624slitwidth_NaND1imaging*fits'           ; flat frames 
Na_D_rngs     = [235,255,475,495]                                         ; ranges of Na D Mercury emission (x pixels in spec domain for D1 and D2)
Na_D_rng      = 10                                                        ; number of pixels around Na D Mercury emission to avoid (in illum scaling for part 2)
do_realign    = 0                                                         ; 0=run part 1 as normal, 1=use previous cross-correlation as master template
frame_ref     = 203                                                       ; this number identifies a frame # in imaging_cube used as the reference key spatially aligning all other frames, works best if object is off-slit
offset_maxes  = [60,60]                                                   ; if the [x,y] spatial offset returned by image correlation is larger than these values, correlation fails and use the brightness centroid 
ims           = [283,759,39,374]                                          ; bounds of the imaging channel: x1,x2,y1,y2 (formerly "imaging_statsec")
sps           = [99,899,502,941]                                          ; bounds of the spectra channel: x1,x2,y1,y2 (formerly "spectra_statsec")
Spectral_platescale = 0.09669                                             ; SPECTRAL channel "/pix Jeff's Jun 5th 2018 email: Using the equatoial diameter of 44.48" on april 26, and the measured size on the spectral channel (460 pix) we get 0.09669 ["]/pix
imaging_platescale  = 0.114820                                            ; IMAGING channel "/pix from measureed Io - Europa seperation on 2018-04-27T07:29:33
effective_seeing  = 2.7                                                   ; Estimate of the effective Seeing FWHM in Arcseconds, This is later solved for iteratively  
Lucky_fraction    = 0.6                                                  ; Fraction of images within the cutoff we're calling "Lucky" 
Hapke_rotation_angle = 203.                                               ; Angle to align Hapke model with Spectral channel
Seeing_Est_Rebin_Factor = 2.                                              ; Rebinning during the seeing estimate & Hapke convolution makes things go alot faster, use 2 or 4. 
subframe      = 50                                                        ; half width of the subframe for the final images (pixels)

SetDefaultValue, night, '15'                                              ; the night of the run; default '15' (only kept for possible future runs with multiple nights)
SetDefaultValue, thresh, 12                                               ; hotpixel removal threshold (sigma above background)
SetDefaultValue, width, 5                                                 ; Median smoothing width to get local background reference (for pixels above thresh) 
SetDefaultValue, minfac, 0.9                                              ; maximum factor to apply to reference spectrum for Na emission extraction
SetDefaultValue, maxfac, 1.3                                              ; minimum factor to apply to reference spectrum for Na emission extraction
SetDefaultValue, dfac, 0.01                                               ; factor increment for minfac->maxfac
SetDefaultValue, ct, 22                                                   ; default color table
SetDefaultValue, smooth_width, 6                                          ; if do_smooth = 1, the pixel width of the smoothing kernel, HACK: CURRENTLY THE SAME FOR SPATIAL AND SPECTRAL. BE CAREFUL!
SetDefaultValue, stddev_cutoff, 1.                                        ; only use Mercury frames with stddev > [stddev_cutoff]_sigma, these are defined as "good" 
SetDefaultValue, Flux_Cal_Wavelength, 5893.                               ; Wavelength to extract and compare to Hapke code (Angstroms)
SetDefaultValue, Na_D1_WL, 5895.92424                                     ; Rest wavelength from NIST
SetDefaultValue, Na_D2_WL, 5889.95095                                     ; Rest wavelength from NIST

; COULD OMIT THIS AND CALCULATE IT IN PART 3 AND 4!
; Dispersion calculation: at y of 265, lines occur at 390 and 632 in x 
  dispersion          = (Na_D1_WL - Na_D2_WL) / (632. - 390.)  ; Angstroms per pixel
  dispersion_Velocity = cspice_clight()*dispersion / Na_D2_WL  ; Km/s per pixel at Na D2
  
; Imaging_platescale: NOTE at Perkins in the imaging Channel North is up (+x) with mrdfits.pro reader
  ;  Io - Europa Seperation on 2018-04-27T07:29:33 is...
  ;  pixel_sep = sqrt((794.-520.)^2+(422.-228.)^2)
  ;  angular_sep = 3600.*sqrt( ((227.76500 - 227.75408)*cos(-16.46847/!radeg))^2 + (16.46847-16.47068)^2 )
  ;  print, angular_sep / pixel_sep
  plate_scale_ratio  =  Spectral_platescale / imaging_platescale            ; Even if Jeff and Carl are both a few pixels off, this is accurate within ~ 3%
  
;; ======================================LOAD SPICE========================================================================
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
    
  ; Prepare the master dark 
    Dark = MRDFITS(dark_dir + dark_file, 0, header, /fscale, /silent )   ; read in dark
    dark = median(dark, dimension = 3)
    dark  = sigma_filter( dark, width, N_sigma=thresh)
    Imaging_dark = dark[ims[0]:ims[1],ims[2]:ims[3]]
    Spectra_Dark = dark[sps[0]:sps[1],sps[2]:sps[3]]
    MWRFITS, dark, outdir + 'MASTER_DARK_FULL_FRAME.fits', header, /create, /silent
    MWRFITS, Imaging_dark, outdir + 'Imaging_dark.fits', header, /create, /silent
    MWRFITS, Spectra_Dark, outdir + 'Spectra_Dark.fits', header, /create, /silent
    
  ; Prepare the Arc Spectrum  
    arc_frames = file_Search(Mercury_Dir + arc_files, count = n_arcs)
    big_array = fltarr(1024, 1024, n_arcs)
    for i = 0, n_arcs-1 do big_array[*,*,i] = MRDFITS(arc_frames[i], 0, header, /fscale, /silent )   ; read in arcs
    Neon = mean(big_array, dimension=3) - dark
    Neon = Neon[*,sps[2]:sps[3]]                                          ; Crop Arc Lamp to Just Spectral Channel
    MWRFITS, Neon, outdir + 'NEON_ARC_LAMP.fits', header, /create, /silent
    
  ; Prepare Imaging Channel Flats
    flat_frames = file_Search(Mercury_Dir + flat_files, count = n_flats)
    big_array = fltarr(1024, 1024, N_flats)
    for i = 0, N_flats-1 do big_array[*,*,i] = MRDFITS(flat_frames[i], 0, header, /fscale, /silent ) ; read in flats
    flat = median(big_array, dimension = 3) - dark                        ; combine and dark subtract the flat
    Imaging_Flat = flat[ims[0]:ims[1],ims[2]:ims[3]]                      ; Crop to just imaging portion
    Spectra_Flat = flat[sps[0]:sps[1],sps[2]:sps[3]]                      ; Crop to just Spectra portion
    Imaging_Flat = imaging_Flat / max(imaging_Flat)                       ; Unsure why I'm doing this 
    Imaging_Flat = Imaging_flat + (1. - median(Imaging_flat))             ; this normalizes such that the median of the flat is 1
    Imaging_flat[WHERE(Imaging_flat lt .01, /NULL)] = !values.f_Nan       ; reject unusual counts for centroid
    Imaging_flat[WHERE(Imaging_flat gt 2.0, /NULL)] = !values.f_Nan       ; reject unusual counts for centroid
    gooddata = where(Finite(Imaging_flat), ngooddata, comp=baddata, ncomp=nbaddata)
    if nbaddata gt 0 then Imaging_flat[baddata] = interpol(Imaging_flat[gooddata], gooddata, baddata) ; Interpolate over the slit
    MWRFITS, Imaging_flat, outdir + 'IMAGING_FLAT.fits', header, /create, /silent 
    
    ; Find indicies of all pixels under the slit in the imaging channel, and interpolate over them
      s = size(Imaging_flat, /dimensions)
      junk = min(total(Imaging_Flat, 2), slit_center)                             ; roughly find the slit center
      h    = histogram(Imaging_Flat, binsize = .05, REVERSE_INDICES=ri)           ; bin the Imaging_Flat into 0.05 bins
      slit_indices = ri[ri[0]:ri[11]-1]                                           ; find inidices of pixels in the lowest few bins of the Imaging_Flat's histogram
      slit_indices = array_indices(Imaging_Flat, slit_indices)                    ; convert into x & y indices
      keep = where(abs(slit_indices[0,*] - slit_center) lt 14, /Null)             ; keep only inidices within 14 pixels of slit center
      slit_indices = slit_indices[*,keep]
      dummy = Imaging_Flat
      dummy[ slit_indices[0,*], slit_indices[1,*] ] = max(Imaging_Flat)           ; show the location of pixels fully under the slit.
      window, 0, xs = s[0], ys = s[1], xpos=0, ypos=0, title = 'Imaging_Flat Field Inspection: Red = Behind Slit --> interpolate'
      cgimage, hist_equal(dummy)
      dummy[ slit_indices[0,*], slit_indices[1,*] ] = !Values.F_NaN
      fill_missing, dummy, !values.F_Nan, 1                                       ; interpolate over slit
      window, 1, xs = s[0], ys = s[1], title = 'Imaging_Flat Field Inspection: Slit be gone!'
      cgimage, hist_equal(dummy)
      Imaging_Flat_Slit_Interpolated = dummy 
      MWRFITS, Imaging_Flat_Slit_Interpolated, outdir + 'IMAGING_FLAT_SLIT_INTERPOLATED.fits', header, /create, /silent 
      
    ; Prepare a spectral flat 
      ispectra_flat     = spectra_flat / median(spectra_flat)
      acre, ispectra_flat, spectra_flat, thresh, width
      spectra_flat      = spectra_flat + (1. - median(spectra_flat))              ; this normalizes again such that median(flat) = 1
  
  ; Write the Mercury data cubes
    filenames = file_Search(Mercury_Dir + Mercury_files, count = n_flats)  
    header            = headfits(filenames[0])
    integration       = sxpar(header, 'EXPOSURE')                         ; integration time for individual frames
    n_frames_per_file = sxpar(header, 'NUMKIN')                           ; # of frames in the "kinetic series" for each "Mercury_file"
    nfiles            =  n_elements(Mercury_files)                         ; # of different Mercury kinetic frames
    s    = [1024, 1024, n_frames_per_file]
    cube = intarr(s[0], s[1], s[2]*nfiles)                                ; define a ***big** cube, where multiple kinetic series are stacked in the third dimension, INTEGER TYPE for now
    print, 'Found: ', N_elements(filenames), ' Mercury files'
    window, 0, xs = ims[1] - ims[0], ys = ims[3]-ims[2], title = 'Inspection: First frame, Imaging Channel' 
    for ifile = 0, nfiles-1 do begin                                      ; loop over kinetic series and add each to a stack
      icube = MRDFITS(filenames[ifile], 0, header, /unsigned, /silent )
      cube(*,*,s[2]*ifile:s[2]*ifile+n_frames_per_file-1) = icube
      integration       = sxpar(header, 'EXPOSURE')                       
      n_frames_per_file = sxpar(header, 'NUMKIN')   
      print, 'Loading: '+Mercury_files[ifile],' containing', strcompress(n_frames_per_file)+' frames of'+strcompress(integration)+' seconds each'    
      tv, bytscl(icube[ ims[0]:ims[1], ims[2]:ims[3], 0])
    endfor
    imaging_cube = float(reform(cube[ims[0]:ims[1],ims[2]:ims[3], *]))    ; extract the IMAGING portion of the kinetic series and combine
    spectra_cube = float(reform(cube[sps[0]:sps[1],sps[2]:sps[3], *]))    ; extract the SPECTRAL portion of the kinetic series and combine
    si        = size(imaging_cube, /dimensions)                           ; size of extracted IMAGING portion of frame
    ss        = size(spectra_cube, /dimensions)                           ; size of extracted SPECTRAL portion of frame
    imaging_cube = imaging_cube - rebin(Imaging_dark, si)                 ; Dark (Bias) subtract the data
    spectra_cube = spectra_cube - rebin(Spectra_dark, ss)                 ; Dark (Bias) subtract the data
    imaging_cube = imaging_cube / rebin(Imaging_flat, si)                 ; FLATTEN THE IMAGING DATA 
    spectra_cube = spectra_cube / rebin(Spectra_flat, ss)                 ; FLATTEN THE Spectral DATA
    MWRFITS, imaging_cube, outdir + 'imaging_cube.fits', header, /create, /silent 
    MWRFITS, spectra_cube, outdir + 'spectra_cube.fits', header, /create, /silent

  ; We'll want a "sky" cube as well  
    Sky_cube          = MRDFITS(sky_dir + sky_file, 0, header, /unsigned, /silent )
    imaging_sky_cube  = float(reform(sky_cube[ims[0]:ims[1],ims[2]:ims[3], *]))
    spectra_sky_cube  = float(reform(sky_cube[sps[0]:sps[1],sps[2]:sps[3], *]))
    imaging_sky_cube = imaging_sky_cube - Imaging_dark                    ; Dark (Bias) subtract the data
    spectra_sky_cube = spectra_sky_cube - Spectra_dark                    ; Dark (Bias) subtract the data
    imaging_sky_cube = imaging_sky_cube / rebin(Imaging_flat, si)         ; FLATTEN THE IMAGING DATA
    spectra_sky_cube = spectra_sky_cube / rebin(Spectra_flat, ss)         ; FLATTEN THE Spectral DATA
    MWRFITS, imaging_sky_cube, outdir + 'imaging_sky_cube.fits', header, /create, /silent
    MWRFITS, spectra_sky_cube, outdir + 'spectra_sky_cube.fits', header, /create, /silent
    beep
endif ; Note that above needs modification in case the integration times change between the fits files we're co-adding

; =====================================================================================================================
; Part 0.5 : Generate a Synthetic Mercury. Match it to the Platescale of each channel
; =====================================================================================================================
if Part eq 0.5 then begin
  if part ne 99 then imaging_cube = MRDFITS(outdir + 'imaging_cube.fits', 0, header, /fscale, /silent )
  if part ne 99 then spectra_cube = MRDFITS(outdir + 'spectra_cube.fits', 0, header, /fscale, /silent )
  ; Generate a Hapke Model of Mercury at the observation time/geometry. Match its plate scale to the RIPS spectral channel
    header = headfits(outdir + 'imaging_cube.fits')                       ; Grab a file header for the timestamp
    Surface_Flux_Calibration, Body = 'Mercury', Timestamp = sxpar(header, 'DATE'), Wavelength = Flux_Cal_Wavelength, Output_image = outdir + 'Hapke_For_Perkins.eps', /align_celestial_north, $
      Hapke_Platescale = Hapke_Platescale, MR_per_A = MR_per_A, Make_picture = Make_picture
    Flux_Cal_Img_Size = size(MR_per_A, /dimensions)
    Flux_Cal_Spectra  = congrid( MR_per_A, round(Flux_Cal_Img_Size[0]*Hapke_Platescale/Spectral_platescale), round(Flux_Cal_Img_Size[1]*Hapke_Platescale/Spectral_platescale), cubic = -0.5)
    Spectra_Hapke     = rot(Flux_Cal_Spectra, Hapke_rotation_angle, MISSING=0., cubic = -0.5)
    Flux_Cal_Imaging  = congrid( MR_per_A, round(Flux_Cal_Img_Size[0]*Hapke_Platescale/Imaging_platescale), round(Flux_Cal_Img_Size[1]*Hapke_Platescale/Imaging_platescale), cubic = -0.5)
    Imaging_Hapke     = rot(Flux_Cal_Imaging, Hapke_rotation_angle, MISSING=0., cubic = -0.5)
    Imaging_Hapke     = rotate(Imaging_Hapke, 7)                         ; mirror image of the spectral channel, X=X, Y=-Y
   
    Make_picture = rot(Make_picture, 28.7+9., cubic= -.5, /interp, MISSING=0.) ;HACK the angle to match the Imaging channel (Or Y flipped spectral channel), but we're certain north is truly north here

  ; Now pad them both with zeros, since we'll be blurring the hell out off them, match the size of each channel
    IH                = size(Imaging_Hapke, /dimensions)
    SH                = size(Spectra_Hapke, /dimensions)
    SI                = size(imaging_cube, /dimensions)
    SS                = size(Spectra_cube, /dimensions)
    padded_Imaging_Hapke = REPLICATE(0., SI[0], SI[1])
    padded_Spectra_Hapke = REPLICATE(0., SS[0], SS[1])
    padded_Imaging_Hapke[ (SI[0]- IH[0])/2., (SI[1]- IH[1])/2. ] = Imaging_Hapke
    padded_Spectra_Hapke[ (SS[0]- SH[0])/2., (SS[1]- SH[1])/2. ] = Spectra_Hapke
    Spectra_Hapke = padded_Spectra_Hapke
    Imaging_Hapke = padded_Imaging_Hapke
    save, Spectra_Hapke, Imaging_Hapke, Make_picture, filename = outdir + 'Hapke_Models.sav'
endif

; =====================================================================================================================
; Part 1 : find Mercury centroids and seeing by cross-correlation with a blurred Hapke image of matching platescale
; =====================================================================================================================
if part eq 1 or part eq 99 then begin ;Find the centroids by cross-correlation with the last image
  if part ne 99 then imaging_cube = MRDFITS(outdir + 'imaging_cube.fits', 0, header, /fscale, /silent )  ; don't need to read in if we're continuing
  if part ne 99 then spectra_cube = MRDFITS(outdir + 'spectra_cube.fits', 0, header, /fscale, /silent )
  if part ne 99 then restore, outdir + 'Hapke_Models.sav'
  
  s         = size(imaging_cube, /dimensions)                           ; size of extracted IMAGING portion of frame
  ss        = size(spectra_cube, /dimensions)                           ; size of extracted SPECTRAL portion of frame

  ; Define alignement/seeing reference and variables based on frame sizes
    reference           = Imaging_Hapke                                 ; Use the idealize Mercury disk as a reference
    Imaging_shift_array = intarr(s[2],2)                                ; array of x, y shift values for aligning images
    sharpness_metric    = fltarr(s[2])                                  ; standard deviations of imaging frames (TEMPORARILY used to define image quality)

  ; Set up windows based on image sizes
    window, 1, xs = 100, ys = 100, xpos=0,     ypos=0,                          title='correlation inspection'
    window, 2, xs = s[0], ys = s[1], xpos=winpos_x+s[0]+20,     ypos=0,         title='Raw Reference Frame for spatial alignment'
    window, 3, xs = s[0], ys = s[1], xpos=winpos_x+s[0]+20,     ypos=s[1]+40,   title='FRAME BY FRAME Comparison'
    window, 4, xs = s[0], ys = s[1], xpos=winpos_x+2*(s[0]+20), ypos=0,         title='Smoothed Raw Frame & Data Alignment'
    window, 5, xs = s[0], ys = s[1], xpos=winpos_x+2*(s[0]+20), ypos=s[1]+40,   title='Aligned Filtered/Smoothed Comparison'

  ; Find indicies of all pixels under the slit in the imaging channel
    flat = MRDFITS(outdir + 'IMAGING_FLAT.fits', 0, header, /fscale, /silent)
    junk = min(total(flat, 2), slit_center)                         ; roughly find the slit center
    h    = histogram(flat, binsize = .05, REVERSE_INDICES=ri)       ; bin the flat into 0.05 bins
    slit_indices = ri[ri[0]:ri[9]-1]                                ; find inidices of pixels in the lowest few bins of the flat's histogram
    slit_indices = array_indices(flat, slit_indices)                ; convert into x & y indices
    keep = where(abs(slit_indices[0,*] - slit_center) lt 10, /Null) ; keep only inidices within 10 pixels of slit center
    slit_indices = slit_indices[*,keep]

  ; Log the correlations to the reference image (correlation post-alignment)
    Correl     = fltarr(s[2])

  print, 'Calculating Frame-By-Frame Spatial Shifts for Co-alignment...'
  print, '-----Frame-#---X-------Y-----Seeing FWHM (")-----
  for i = 0, s[2]-1 do begin

    ; interpolate over slit 
      iframe = reform(imaging_cube[*,*,i])                            ; "raw" imaging frame bias and flat corrected
      iframe[ slit_indices[0,*], slit_indices[1,*] ] = !values.F_Nan  ; define slit pixels as NaN
      masked_frame = iframe
      fill_missing, iframe, !values.F_Nan, 1                          ; interpolate over slit    
      frame = iframe
  
    ;filter the frames and blur the Hapke reference by a fixed seeing estimate "effective_seeing" in arcsec FWHM
      img_bandpass      = float(bandpass_filter(frame, 0., 0.15, /butterworth))
      seeing_sigma      = effective_seeing / (2.0*sqrt(2.0*alog(2)))  ; convert FWHM to sigma
      seeing_sigma      = seeing_sigma/Imaging_platescale             ; convert to Gaussian sigma in pixels
      ref_bandpass      = convol(reference, GAUSSIAN_FUNCTION([seeing_sigma, seeing_sigma]), /EDGE_ZERO, /normalize, /NAN)
      ref_bandpass1     = ref_bandpass
  
      junk = max(total(ref_bandpass,2), ref_bandpass_centeroid_x) & junkx = max(total(img_bandpass,2), img_bandpass_centeroid_x)
      junk = max(total(ref_bandpass,1), ref_bandpass_centeroid_y) & junky = max(total(img_bandpass,1), img_bandpass_centeroid_y)

    ; Get the rough alignment
      zstack_align_images, img_bandpass, ref_bandpass, x_rough, y_rough
      ref_bandpass= ref_bandpass1
   
    ; Now estimate the seeing and find the optimal alignment
      subframe = 100
      B = float(reference[ s[0]/2. - subframe:s[0]/2. + subframe - 1, s[1]/2. - subframe:s[1]/2. + subframe - 1 ]) 
      A = shift(masked_frame, x_rough, y_rough)
      A = A[ s[0]/2. - subframe:s[0]/2. + subframe - 1, s[1]/2. - subframe:s[1]/2. + subframe - 1 ] 
      A = rebin(A, 2.*subframe / Seeing_Est_Rebin_Factor, 2.*subframe / Seeing_Est_Rebin_Factor)
      B = rebin(B, 2.*subframe / Seeing_Est_Rebin_Factor, 2.*subframe / Seeing_Est_Rebin_Factor)
  
      y_fine = 0. & x_fine = 0.                                       ; A default in case no shift fine tuning is returned in Amoeba
      wset, 1
      sharpness_metric[i] = AMOEBA(1.e-5, function_name='Estimate_Seeing', FUNCTION_VALUE = fval, P0 = [effective_seeing/Seeing_Est_Rebin_Factor], Scale = [1./Seeing_Est_Rebin_Factor], NCalls = 100 ) 
      if y_fine + x_fine gt 0. then print, 'Fine tuned Alignment' 
      Imaging_shift_array[i,*] = fix(round([x_rough + x_fine*Seeing_Est_Rebin_Factor, y_rough + y_fine*Seeing_Est_Rebin_Factor]))             ; Need to add fine tuned shift from amoeaba somehow
  
    ; If correlation fails then use the brightness centroid of the filtered/smoothed frames
      if abs(Imaging_shift_array[i,0]) gt offset_maxes[0] or abs(Imaging_shift_array[i,1]) gt offset_maxes[1] then begin ; If Failed Correlation
        Imaging_shift_array[i,*] = [ref_bandpass_centeroid_x - img_bandpass_centeroid_x, ref_bandpass_centeroid_y - img_bandpass_centeroid_y]
        print, i, Imaging_shift_array[i,0], Imaging_shift_array[i,1], sharpness_metric[i]*Seeing_Est_Rebin_Factor, '  <--- Correlation alignment failed, used brightness centroid'
      endif else print, i, Imaging_shift_array[i,0], Imaging_shift_array[i,1], sharpness_metric[i]*Seeing_Est_Rebin_Factor            ; update progress
  
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
  
  sharpness_metric = sharpness_metric*Seeing_Est_Rebin_Factor                       ; Seeing was determined using rebinned arrays, this is the true seeing if data unbinned
  Spectra_shift_array = fix(round(float(Imaging_shift_array) / plate_scale_ratio))  ; Spectral and imaging channels have different platescales, hence different shifts
  Spectra_shift_array[*,1] = -Spectra_shift_array[*,1]                              ; X_i = X_s, but +Y_i in the imaging channel is -Y_s in the spectral channel
  aligned_imaging_cube = fltarr(s)                                                  ; Imaging cube after everything is aligned
  aligned_spectra_cube = fltarr(ss)                                                 ; Spectral cube after alignment in Y only (for now) (note that

  ; Now loop over all frames to align both channels
    for i = 0, s[2]-1 do begin
      frame     =  reform(imaging_cube[*,*,i])                                      ; "raw" imaging frame
      specframe =  reform(spectra_cube[*,*,i])                                      ; "raw" spectral frame
      aligned_imaging_cube[*,*,i] = shift(frame, [Imaging_shift_array[i,*]])        ; Co-align imaging frame in X and Y
      aligned_spectra_cube[*,*,i] = shift(specframe, [0, Spectra_shift_array[i,1]]) ; Co-align Spectral frame **in Y only**. Y_Spec = -Y_img/PSR. Later on, we do the X_Spec = X_img/PSR.
    endfor

  save, Imaging_shift_array, Spectra_shift_array, correl, aligned_imaging_cube, aligned_spectra_cube, sharpness_metric, Slit_indices, filename = outdir + 'shift_array.sav'
  beep
endif

if part eq 1.5 then begin ;re-align now that we know the seeing in each
  if part ne 99 then restore, outdir + 'shift_array.sav'
  if part ne 99 then imaging_cube = MRDFITS(outdir + 'imaging_cube.fits', 0, header, /fscale, /silent )  ; don't need to read in if we're continuing
  if part ne 99 then spectra_cube = MRDFITS(outdir + 'spectra_cube.fits', 0, header, /fscale, /silent )
  if part ne 99 then restore, outdir + 'Hapke_Models.sav'
  
  s         = size(imaging_cube, /dimensions)                        ; Size of extracted IMAGING portion of frame
  reference = Imaging_Hapke                                          ; Use the idealized Mercury disk as a reference
  optimally_aligned_imaging_cube = fltarr(s)
  sharpness_metric  = Seeing_Est_Rebin_Factor*sharpness_metric       ;  HACK HACK HACK HACK HACK not needed once part 1 is run again
  
  for i = 0, s[2]-1 do begin

    ; interpolate over slit
      frame = reform(imaging_cube[*,*,i])                            ; "raw" imaging frame bias and flat corrected
      frame[ slit_indices[0,*], slit_indices[1,*] ] = !values.F_Nan  ; define slit pixels as NaN
      rough_aligned_frame = shift(frame, [Imaging_shift_array[i,*]])    
      rough_aligned_masked_frame = rough_aligned_frame               ; This one keeps slit as NaN
      fill_missing, rough_aligned_frame, !values.F_Nan, 1            ; rough_aligned_frame now has pixels interpolated over the slit

    ; filter the frames and blur the Hapke reference by a fixed seeing estimate "effective_seeing" in arcsec FWHM
      img_bandpass      = float(bandpass_filter(rough_aligned_frame, 0., 0.15, /butterworth))
      seeing_sigma      = sharpness_metric[i] / (2.0*sqrt(2.0*alog(2)))  ; convert FWHM to sigma
      seeing_sigma      = seeing_sigma/Imaging_platescale             ; convert to Gaussian sigma in pixels
      ref_bandpass      = convol(reference, GAUSSIAN_FUNCTION([seeing_sigma, seeing_sigma]), /EDGE_ZERO, /normalize, /NAN)
      ref_bandpass1     = ref_bandpass
      
      hapke_center      = fix(round(s/2.))
      sub_ref           = ref_bandpass[hapke_center[0]-subframe:hapke_center[0]+subframe-1,hapke_center[1]-subframe:hapke_center[1]+subframe-1]
      sub_img           = img_bandpass[hapke_center[0]-subframe:hapke_center[0]+subframe-1,hapke_center[1]-subframe:hapke_center[1]+subframe-1]

    ; Get the optimal alignment, this uses the instantantaeous seeing, where as the previous part used a guess for the seeing
      ; zstack_align_images, sub_img, sub_ref, x_optimum, y_optimum
      CORREL_OPTIMIZE, sub_img, sub_ref, x_optimum, y_optimum, magnification = 2
      print, 'Tweaking alignment of frame', i, ' by',  x_optimum, y_optimum
      
      optimally_aligned_imaging_cube[*,*,i] = smart_shift(reform(aligned_imaging_cube[*,*,i]), x_optimum, y_optimum, /interp)        ; Co-align imaging frame in X and Y
   endfor

  ranked = reverse(SORT(sharpness_metric))
  window, 0, xs = s[0], ys = s[1]
  cgimage, total(aligned_imaging_cube[*,*, ranked[1600:1699]], 3), /keep_aspect
  window, 1, xs = s[0], ys = s[1]
  cgimage, total(optimally_aligned_imaging_cube[*,*, ranked[1600:1699]], 3), /keep_aspect
  window, 2, xs = s[0], ys = s[1]
  cgimage, total(optimally_aligned_imaging_cube[*,*, ranked[500:1699]], 3) - total(aligned_imaging_cube[*,*, ranked[500:1699]], 3), /keep_aspect
  ;unclear if this is indeed any improvement, and there's some way off yonder---and its slow!
  save, optimally_aligned_imaging_cube, filename = outdir + 'optimal_shift_array.sav'
endif

;;******************************testing sharpness metric***********************
;if part eq 1.5 then begin ;testing for now! do after part 4
;  restore, outdir + 'shift_array.sav'  
;  restore, outdir + 'images_to_plot.sav'
;
;  s = size(aligned_imaging_cube, /dimensions)
;  sh = size(hapke, /dimensions)
;  hapke = congrid(hapke, sh[0]*plate_scale_ratio, sh[1]*plate_scale_ratio, cubic = -0.5)  ; Scale it
;  seeing_sigma      = 2.7 / (2.0*sqrt(2.0*alog(2)))             ; convert FWHM to sigma
;  seeing_sigma      = seeing_sigma/Spectral_platescale     ; convert to Gaussian sigma in pixels
;  Blurry_Hapke = convol(Hapke, GAUSSIAN_FUNCTION([seeing_sigma, seeing_sigma]), /EDGE_ZERO, /normalize)
;  Blurry_Hapke = rotate(Blurry_Hapke, 7)
;  hapke = rotate(hapke, 7)
;  ranked = reverse(SORT(correl))
;  find_shift_with_me = total(aligned_imaging_cube[*,*, ranked[0:99]], 3, /NaN)
;  CORREL_OPTIMIZE, Blurry_Hapke, find_shift_with_me, x, y, /NUMPIX
;
;  ; mask the slit:
;    for i = 0, s[2]-1 do begin
;      slit_mask = [ Slit_indices[0,*] + Imaging_shift_array[i,0], Slit_indices[1,*] + Imaging_shift_array[i,1] ]
;      frame = reform(aligned_imaging_cube[*,*,i])
;      frame[ slit_mask[0,*], slit_mask[1,*] ] = !values.F_NaN
;      aligned_imaging_cube[*,*,i] = frame
;      ;test = reform(aligned_imaging_cube[*,*,i])
;      ;tv, bytscl(test)
;    endfor
;  
;  ; now align the cube to hapke
;    corr_w_hapke = shift(aligned_imaging_cube, x, y, 0)
;
;  Hapke_corr = fltarr(s[2])
;  for i = 0, s[2]-1 do begin 
;    this_hapke = hapke[ 0:s[0]-1, 0:s[1]-1 ] ;crop it
;    junk = where(finite(corr_w_Hapke[*,*,i]), complement = mask)
;    ;this_hapke[mask] = !values.F_NaN ;don't maskk out hapke slit
;    if i eq 0 then begin
;      window, 0
;      tv, bytscl(corr_w_Hapke[*,*,i])
;      window, 1
;      tv, bytscl(this_Hapke)
;    endif  
;    B = this_Hapke
;    A = reform(corr_w_Hapke[*,*,i])
;    A = A[226:365, 130:289]
;    B = B[226:365, 130:289]
;;    A = rebin(a, 70, 80)
;;    B = rebin(b, 70, 80)
;    A = rebin(a, 35, 40)
;    B = rebin(b, 35, 40)
;
;    ;CORREL_OPTIMIZE, A, B, xoffset_optimum, yoffset_optimum, /NUMPIX
;
;   ;stop
;
;    ;Hapke_CORR[i] = AMOEBA(1.e-4, function_name='Estimate_Seeing', FUNCTION_VALUE = fval, P0 = [2.8], Scale = [1.], NCalls = 100 )
;    ;Hapke_CORR[i] = AMOEBA(1.e-4, function_name='Estimate_Seeing', FUNCTION_VALUE = fval, P0 = [1.4], Scale = [.5], NCalls = 100 )
;    Hapke_CORR[i] = AMOEBA(1.e-4, function_name='Estimate_Seeing', FUNCTION_VALUE = fval, P0 = [0.7], Scale = [.25], NCalls = 100 )
;    
;    ;stop
;    Print, i, ' seeing = ', 4.*Hapke_CORR[i]
;    ;indices = [ where( finite(A) eq 1 ), where( finite(B) eq 1 ) ]
;    ;common_indices = indices( UNIQ(indices, sort(indices)) )
;
;    ;Hapke_CORR[i] = CORRELATE( A[common_indices], B[common_indices] )
;    ;stop
;    ;Hapke_CORR[i] = CORREL_IMAGES( this_Hapke, corr_w_Hapke[*,*,i], xoffset_b = 0, yoffset_B = 0, XSHIFT =0, ySHIFT = 0) ; 2-D Correlation at this alignment
;  endfor
;
;  ranked = reverse(SORT(Hapke_corr))
;  window, 2
;  cgimage, total(corr_w_hapke[*,*, ranked[0:99]], 3, /NaN), /keep_aspect
;  window, 3
;  ;cgimage, total(aligned_imaging_cube[*,*, ranked[1500:1599]], 3), /keep_aspect
;  cgimage, total(corr_w_hapke[*,*,ranked[1500:1599]], 3, /NaN), /keep_aspect
;  ;
;  ;seems to work fairly well.
;  
;  save, Hapke_CORR, filename = outdir + 'Hapke_CORR.sav'
;endif
;;******************************testing sharpness metric***********************



; =====================================================================================================================
; Part 2 : Isolate the sodium emission in every spectral frame
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

  print, 'Applying a sky spectrum as the solar reflectance spectrum...'
  spectra_sky_cube = MRDFITS(outdir + 'spectra_sky_cube.fits', /silent ) 
  acre, spectra_sky_cube, reference, thresh, width                  ; Clean any remaining hot pixels

  ; Now find the xoffset in the reference (Sky, Moon, etc) reflectance spectrum (likely due to slightly different grating angles)
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
      frame       = reform(aligned_spectra_cube[*,*,i])               ; dark subtract each spectrum
  
      ; Now find best scaling for ref_spectrum and subtract (based on solar reflected emission between D lines        
        illum = frame / reference                                     ; scale the illumination against Mercury's disk
        for iD = 0, 2, 2 do illum[Na_D_rngs[iD+0]:Na_D_rngs[iD+1],*] = !values.F_NaN  ; carefully avoid sodium emissions when scaling to the disk
  
        illum_along_slit = MEDIAN(illum, dimension=1)
        illum_along_slit = MEDSMOOTH( illum_along_slit, 5 )           ; helps with hot pixels from bad flat-fielding where the slit has some dust
        scaled_reference = reference * rebin(transpose(illum_along_slit), s[0], s[1])
        Just_exosphere = Frame - scaled_reference * 1.

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
            ;if do_smooth then cgimage, bytscl(scaled_reference, 0, 200), /axes, title='smoothed and scaled by ' + strtrim(best_fac) $
            ;             else cgimage, bytscl(scaled_reference, 0, 200), /axes, title='scaled by ' + strtrim(best_fac)
            cgimage, bytscl(scaled_reference, 0, 200), /axes, title='scaled by ' + strtrim(best_fac)
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
    ;if do_smooth then img = smooth( exosphere_spectra_cube[*,*,n], [2, smooth_width] ) $
    ;             else img = exosphere_spectra_cube[*,*,n]
    img = exosphere_spectra_cube[*,*,n]
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
          D2_linewidth[n,i] = dispersion_Velocity*sqrt(((2.0*sqrt(2.0*alog(2)))*A[2])^2 - (2.45)^2) ;FWHM KM/S instrumental measured width is 2.45 pixels 
          D2_err_linewidth[n,i] = dispersion_Velocity*(2.0*sqrt(2.0*alog(2)))*Err_A[2]
       endif else begin
          D2_brightness[n,i] = !values.F_nan
          D2_linewidth[n,i] = !values.F_nan
       endelse 
    endfor
  endfor 
   
  Print, 'Fitting exosphere D1 in every spatial bin of every frame (slow, ~200 frames/min)...'
  A = [80.,D1_trace[0],1.5,0.] ; Initial guess for fitting
  for n = 0, s[2]-1 do begin
    ;if do_smooth then img = smooth( exosphere_spectra_cube[*,*,n], [2, smooth_width] ) $
    ;             else img = exosphere_spectra_cube[*,*,n]
    img = exosphere_spectra_cube[*,*,n]
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
        D1_linewidth[n,i] = dispersion_Velocity*sqrt(((2.0*sqrt(2.0*alog(2)))*A[2])^2 - (2.45)^2) ;FWHM KM/S instrumental measured width is 2.45 pixels
        D1_err_linewidth[n,i] = dispersion_Velocity*(2.0*sqrt(2.0*alog(2)))*Err_A[2]
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
  if part ne 99 then restore, outdir + 'Hapke_Models.sav'
  s  = size(surf_refl_spectra_cube, /dimensions)
  
  ; Reconstruct an image from the spectral channel at the calibration wavelength
    ; Find the pixel location of the calibration wavelength
      Cal_wavelength_trace    = fltarr(N_elements(D2_Trace))
      for i = 0, N_elements(D2_Trace)-1 do Cal_wavelength_trace[i] = interpol( [D2_Trace[i], D1_trace[i]], [Na_D2_Wl, Na_D1_Wl], Flux_Cal_Wavelength)
      Cal_wavelength_trace = round(Cal_wavelength_trace)
    
    ; As in part 3 for sodium, extact brightness slices. Use the surface reflectance cube, the exosphere cube already has the reflectance component removed
      DN_per_A  = fltarr(s[2], s[1])
      for n = 0, s[2]-1 do begin
        ;if do_smooth then img = smooth( surf_refl_spectra_cube[*,*,n], [2, smooth_width] ) $
                     ;else img = surf_refl_spectra_cube[*,*,n]
        img = surf_refl_spectra_cube[*,*,n]             
        for i = 0, s[1]-1 do begin              
          DN_per_A[n, i] = TOTAL(img[Cal_wavelength_trace[i]-(0.5/dispersion):Cal_wavelength_trace[i]+(0.5/dispersion), i]) 
        endfor
      endfor

  ; Reconstruct Calibration Image: Place this time series of 1D spectral slices along the slit into a (time averaged) 2D image
    flux_cube         = fltarr(s)   ; This datacube will hold the image of the slit for each frame, slit is extacted at the calibration wavelength over 1 Angstrom 
    flux_cube_weights = fltarr(s)   ; This datacube will hold the weights from the point spread function of the slit for each frame 
    Na_D1_cube        = fltarr(s)   ; This datacube will hold the image of the slit for each frame, slit is extacted at exosphere Na D1
    Na_D2_cube        = fltarr(s)   ; This datacube will hold the image of the slit for each frame, slit is extacted at exosphere Na D2
    
  ; Get the point spread function from Arc frames      
    Neon = MRDFITS(outdir + 'NEON_ARC_LAMP.fits', 0, header, /fscale, /silent )
    Neon = total(Neon[*,50:100], 2)          ; Always inspect that this this portion of the line is straight
    junk = max(Neon, line_center)
    Neon = Neon[line_center-5:line_center+5]
    PSF_1D          = Neon/total(Neon)       ; Point spread function accros the slit (at a 624 width setting the slit is surely resolved, and RIPS is acting like an imaging spectrograph) 
    PSF_Size        = N_elements(PSF_1D)
    spatial_weights = rebin(PSF_1D, psf_size, s[1])
    home = fix(round([s[0]/2.,s[1]/2.]))
        
    for i = 0, s[2]-1 do begin; loop over frames / time. Note that all spectral cubes have already been aligned in -Y/plate_scale_ratio. Only +X/plate_scale_ratio is needed  
      flux_cube[*,*,i]         = rebin(DN_per_A[i,*], s[0], s[1])
      Na_D2_cube[*,*,i]        = rebin(D2_Brightness[i,*], s[0], s[1])
      Na_D1_cube[*,*,i]        = rebin(D1_Brightness[i,*], s[0], s[1])
      weighting_frame          = fltarr(s[0], s[1]) 
      weighting_frame[home[0] + Spectra_shift_array[i,0] - PSF_Size/2.: home[0] + Spectra_shift_array[i,0] + PSF_Size/2.-1,*] = spatial_weights      
      flux_cube_weights[*,*,i] = weighting_frame
      ;cgimage, flux_cube_weights[*,*,i] * flux_cube[*,*,i]
      ;wait, 0.04
    endfor    

    ; now take the geometric mean using the weighting
      Calib_img    = total(flux_cube*flux_cube_weights, 3, /NAN)  / total(flux_cube_weights, 3)
      Na_D2_img    = total(Na_D2_cube*flux_cube_weights, 3, /NAN) / total(flux_cube_weights, 3)
      Na_D1_img    = total(Na_D1_cube*flux_cube_weights, 3, /NAN) / total(flux_cube_weights, 3)
      some_valid_data       = where(finite(Calib_img), complement = no_brightness_information)
      Calib_img[no_brightness_information] = 0.
      Na_D2_img[no_brightness_information] = 0. 
      Na_D1_img[no_brightness_information] = 0.
   
    ; Gather top 10%
      ranked = SORT(sharpness_metric)      ; Rank lowest to highest seeing
      cut_index = fix(s[2]*Lucky_fraction)

    ; now take the geometric mean of Just the lucky ones
      Lucky_Calib_img    = total(flux_cube[*,*,ranked[0:cut_index]]*flux_cube_weights[*,*,ranked[0:cut_index]], 3, /NAN)  / total(flux_cube_weights[*,*,ranked[0:cut_index]], 3)
      Lucky_Na_D2_img    = total(Na_D2_cube[*,*,ranked[0:cut_index]]*flux_cube_weights[*,*,ranked[0:cut_index]], 3, /NAN) / total(flux_cube_weights[*,*,ranked[0:cut_index]], 3)
      Lucky_Na_D1_img    = total(Na_D1_cube[*,*,ranked[0:cut_index]]*flux_cube_weights[*,*,ranked[0:cut_index]], 3, /NAN) / total(flux_cube_weights[*,*,ranked[0:cut_index]], 3)
      some_valid_data    = where(finite(Lucky_Calib_img), complement = no_brightness_information)
      Lucky_Calib_img[no_brightness_information] = 0.
      Lucky_Na_D2_img[no_brightness_information] = 0.
      Lucky_Na_D1_img[no_brightness_information] = 0.    

  ; Blur a Hapke Model of Mercury at the observation time/geometry. Match its plate scale to the RIPS spectral channel
    seeing_sigma      = effective_seeing / (2.0*sqrt(2.0*alog(2)))    ; convert FWHM to sigma
    seeing_sigma      = seeing_sigma/Spectral_platescale              ; convert to Gaussian sigma in pixels
    Blurred_Hapke = convol(Spectra_Hapke, GAUSSIAN_FUNCTION([seeing_sigma, seeing_sigma]), /EDGE_ZERO, /normalize)

    junk_s = max(Calib_img, loc_spec)
    spec_centroid = array_indices(Calib_img, loc_spec)
    junk_h = max(Blurred_Hapke, loc_hapke)
    hapke_centroid = array_indices(Blurred_Hapke, loc_hapke)
    
    sh = size(Blurred_Hapke, /dimensions)                            ; Hapke image dimensions
    si = size(aligned_imaging_cube, /dimensions)                     ; Imaging channel image dimensions
         
    ; Assemble the imaging channel results, expand them to match the spectral platescale 
      Imaging_Ch = total(aligned_imaging_cube, 3)                                                      & Lucky_Imaging_Ch = total(aligned_imaging_cube[*,*,ranked[0:cut_index]], 3)
      Imaging_Ch = rotate(Imaging_Ch, 7)                                                               & Lucky_Imaging_Ch = rotate(Lucky_Imaging_Ch, 7)
      Imaging_Ch = congrid(Imaging_Ch, si[0]/plate_scale_ratio, si[1]/plate_scale_ratio, cubic = -0.5) & Lucky_Imaging_Ch = congrid(Lucky_Imaging_Ch, si[0]/plate_scale_ratio, si[1]/plate_scale_ratio, cubic = -0.5)
      
    ; Crop and Co-align everything  
      Blurred_Hapke_align = Blurred_Hapke
     
      ; Crop everthing derived from the spectral channel and align it to the blurred Hapke Function
        ; crop calib image to only regions of valid data since all the zeros will throw off cross-correlation alignments
        valid = minmax( where(total(Calib_img, 2) ne 0.) )
        zstack_align_images, Calib_img[valid[0]:valid[1],*], Blurred_Hapke_align[valid[0]:valid[1],*], x, y
        Calib_img = shift(Calib_img, x, y)  
        Na_D2_img = shift(Na_D2_img, x, y)  
        Na_D1_img = shift(Na_D1_img, x, y)  

      ; Center of the planet's disk is at Hapke_Center, preserve this in the subframes
        hapke_center  = fix(round(sh/2.))    
        Blurred_Hapke = Blurred_Hapke[hapke_center[0]-subframe:hapke_center[0]+subframe-1,hapke_center[1]-subframe:hapke_center[1]+subframe-1]
        Calib_img     = Calib_img[hapke_center[0]-subframe:hapke_center[0]+subframe-1,hapke_center[1]-subframe:hapke_center[1]+subframe-1]
        Na_D1_img     = Na_D1_img[hapke_center[0]-subframe:hapke_center[0]+subframe-1,hapke_center[1]-subframe:hapke_center[1]+subframe-1]
        Na_D2_img     = Na_D2_img[hapke_center[0]-subframe:hapke_center[0]+subframe-1,hapke_center[1]-subframe:hapke_center[1]+subframe-1]
        
      ; And repeat for the Lucky ones in case that alignment differs
        valid = minmax( where(total(Lucky_Calib_img, 2) ne 0.) )
        zstack_align_images, Lucky_Calib_img[valid[0]:valid[1],*], Blurred_Hapke_align[valid[0]:valid[1],*], x, y
        Lucky_Calib_img = shift(Lucky_Calib_img, x, y)
        Lucky_Na_D2_img = shift(Lucky_Na_D2_img, x, y)
        Lucky_Na_D1_img = shift(Lucky_Na_D1_img, x, y)
        Lucky_Calib_img = Lucky_Calib_img[hapke_center[0]-subframe:hapke_center[0]+subframe-1,hapke_center[1]-subframe:hapke_center[1]+subframe-1]
        Lucky_Na_D1_img = Lucky_Na_D1_img[hapke_center[0]-subframe:hapke_center[0]+subframe-1,hapke_center[1]-subframe:hapke_center[1]+subframe-1]
        Lucky_Na_D2_img = Lucky_Na_D2_img[hapke_center[0]-subframe:hapke_center[0]+subframe-1,hapke_center[1]-subframe:hapke_center[1]+subframe-1]
        
      ; Now align the imaging channel derived quantitites
        sic = size(imaging_ch, /dimensions)
        zstack_align_images, Imaging_Ch, Blurred_Hapke_align[0:sic[0]-1, 0:sic[1]-1], x, y
        Imaging_Ch = shift(Imaging_Ch, x, y)
        Imaging_Ch = Imaging_Ch[hapke_center[0]-subframe:hapke_center[0]+subframe-1,hapke_center[1]-subframe:hapke_center[1]+subframe-1]

      ; And again repeat for the Lucky ones in case that alignment differs
        zstack_align_images, Lucky_Imaging_Ch, Blurred_Hapke_align[0:sic[0]-1, 0:sic[1]-1], x, y
        Lucky_Imaging_Ch = shift(Lucky_Imaging_Ch, x, y)
        Lucky_Imaging_Ch = Lucky_Imaging_Ch[hapke_center[0]-subframe:hapke_center[0]+subframe-1,hapke_center[1]-subframe:hapke_center[1]+subframe-1]

        ; Need to ratio the blurred Hapke image and the data-generated calibration image over the same to determine sensitivity
        MR_per_DN = TOTAL(Blurred_Hapke) / TOTAL(Calib_img) 

;        MR_per_DN = max(Blurred_Hapke[hapke_centroid[0]-30:hapke_centroid[0]+30, $
;          hapke_centroid[1]-40:hapke_centroid[1]+40]) $
;          / max(Calib_img[spec_centroid[0]-30:spec_centroid[0]+30, $
;          spec_centroid[1]-40:spec_centroid[1]+40]) ; ... and yet this looks better

  ; Calibrate the data in absolute flux, MegaRayleighs (line emission), MegaRayleighs/Angstrom (continuum) 
    Calib_img       = Calib_img * MR_per_DN
    Na_D2_img       = Na_D2_img * MR_per_DN
    Na_D1_img       = Na_D1_img * MR_per_DN
    Lucky_Calib_img = Lucky_Calib_img * MR_per_DN
    Lucky_Na_D2_img = Lucky_Na_D2_img * MR_per_DN
    Lucky_Na_D1_img = Lucky_Na_D1_img * MR_per_DN
    
  ; FLIP BACK All images to match the Imaging Channel / Plane of view on the sky.  
    Imaging_Ch       = rotate(Imaging_Ch, 7)  
    Lucky_Imaging_Ch = rotate(Lucky_Imaging_Ch, 7)  
    Calib_img        = rotate(Calib_img, 7)  
    Lucky_Calib_img  = rotate(Lucky_Calib_img, 7)                   
    Na_D1_img        = rotate(Na_D1_img, 7)
    Lucky_Na_D1_img  = rotate(Lucky_Na_D1_img, 7)
    Na_D2_img        = rotate(Na_D2_img, 7) 
    Lucky_Na_D2_img  = rotate(Lucky_Na_D2_img, 7) 
    Blurred_Hapke    = rotate(Blurred_Hapke, 7) 

  ; Setup inspection windows
    window, 0, xpos=winpos_x,              ypos=winpos_y,         xs=8*subframe, ys=8*subframe, title='EXTRACTED DN ALONG IN A 1 ANGSTROM BAND (DN/A)'
    window, 1, xpos=winpos_x+8*subframe+20,ypos=winpos_y,         xs=8*subframe, ys=8*subframe, title='EXTRACTED NA D2 DN'
    window, 2, xpos=winpos_x+8*subframe+20,ypos=winpos_y+8*subframe+40, xs=8*subframe, ys=8*subframe, title='EXTRACTED NA D1 DN'
    window, 3, xpos=winpos_x,              ypos=winpos_y+8*subframe+40, xs=8*subframe, ys=8*subframe, title='HAPKE MODEL (MR/A)'
    window, 4, xpos=winpos_x+8*subframe+20,ypos=winpos_y+8*subframe+40, xs=8*subframe, ys=8*subframe, title='ALIGNED & COADDED IMAGING CHANNEL, MATCHED TO SPECTRAL CHANNEL PLATESCALE'
    window, 5, xpos=winpos_x+2*(8*subframe+20),ypos=winpos_y, xs=8*subframe, ys=8*subframe, title='LUCKY EXTRACTED DN ALONG IN A 1 ANGSTROM BAND (DN/A)'
    
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
      cgimage, Blurred_Hapke, minvalue = 0, maxvalue = Continuum_colorbar_top ;Blurred_Hapke
    wset, 4
      cgimage, Imaging_Ch, /keep_aspect
      cgimage, Lucky_Imaging_Ch, /keep_aspect
    wset, 5
      cgimage, lucky_Calib_img, minvalue = 0, maxvalue = 40.
      cgColorbar, minrange = 0., maxrange = Continuum_colorbar_top, title = cgsymbol('Sigma')+'1'+cgsymbol('Angstrom') +' at '+ strcompress(Flux_Cal_Wavelength)+cgsymbol('Angstrom') + ' (MegaRayleighs/'+ cgsymbol('Angstrom')+')'   

    save, Lucky_Imaging_Ch, Imaging_Ch, Lucky_Calib_img, Calib_img, Lucky_Na_D2_img, Lucky_Na_D1_img, Na_D2_img, Na_D1_img, Blurred_Hapke, $
          filename = outdir + 'images_to_plot.sav' 
endif

; =====================================================================================================================
; Part 5 : Plot all images and and make movies of the exosphere.
; =====================================================================================================================
if part eq 5 or part eq 99 then begin
  if part ne 99 then restore, outdir + 'images_to_plot.sav'          ; contains all images to be plotted
  if part ne 99 then restore, outdir + 'Hapke_Models.sav'            ; contains all images to be plotted
  
  ; Find Mercury's True Anomaly Angle
    filenames         = FILE_SEARCH(Mercury_dir + '*_'+Mercury_files)
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
    print, 'Mercury''s True Anomaly Angle:', True_Anomaly

  ; Adjust the size of the Hapke image with albedo features
    cspice_bodvrd, Body, 'RADII', 3, radii                              ; Get Body shape constants
    cspice_spkpos, Body, ET, 'J2000', 'LT+S', 'Earth', ptarg, ltime
    R_M = 206264.806 * atan(radii[0] / norm(ptarg))                     ; Radius of the body in arcsec
    if body eq 'Mercury' then print, 'Angular diameter of ', body, ' =', 2.*R_M, ' arcsec'
    
    frame_width_arcsec = Spectral_platescale * subframe * 2. 
    sp = size(make_picture, /dimensions)
    picture_frame = fltarr(sp[0]*frame_width_arcsec / (2.*R_M), sp[0]*frame_width_arcsec / (2.*R_M))
    spf = size(picture_frame, /dimensions)
    picture_frame[(spf[0] - sp[0])/2., (spf[1] - sp[1])/2.] = make_picture
    
  cgPS_Open, filename = outdir+'Mercury_Perkins.eps', /ENCAPSULATED, xsize = 6, ysize = 6
    !P.font=1
    device, SET_FONT = 'Helvetica Bold', /TT_FONT
    !p.charsize = 1.2
    
    yspace = 4.0*(!D.Y_CH_SIZE)/!D.Y_SIZE
    xspace = 4.0*(!D.X_CH_SIZE)/!D.X_SIZE
    pos = cgLayout([2,2], xgap = 1., ygap = 1, oxmargin = 4, oymargin = 1)
    thick = 4.
    
    P = pos[*,0]
    loadct, 0
    cgimage, picture_frame, pos = P, /keep_aspect
    cgloadct, ct

    P = pos[*,1]
    p[0] = p[0] - .01
    p[2] = p[2] - .01
    levels = (findgen(5) + 1) /5.
    cgcontour, Blurred_Hapke / max(Blurred_Hapke, /NaN), levels=levels, label=0, /noerase, pos = P, color = 'red', aspect = 1, thick = 1., XTICKFORMAT = "(A1)", YTICKFORMAT= "(A1)"
    cgcontour, Imaging_Ch / max(Imaging_Ch, /NaN), levels=levels, label=0, pos = P, /noerase, color = 'blue', thick = 1. , XTICKFORMAT = "(A1)", YTICKFORMAT= "(A1)"
    cgcontour, Calib_img / max(Calib_img, /NaN), levels=levels, label=0, pos = P, /noerase, color = 'black', thick = 1. , XTICKFORMAT = "(A1)", YTICKFORMAT= "(A1)"
    ;cgcontour, Lucky_Calib_img / max(Lucky_Calib_img, /NaN), levels=levels, label=0, pos = P, /noerase, color = 'green', thick = 1. , XTICKFORMAT = "(A1)", YTICKFORMAT= "(A1)"

    P = pos[*,2]
    cgimage, Calib_img, /keep_aspect, pos = P, /noerase
    ;cgimage, Lucky_Calib_img, /keep_aspect, pos = P, /noerase
    p = [p[0] + xspace+0.01, p[1]+.02, p[0]+xspace+0.02, p[3]-.02]
    cgColorbar, pos = P, minrange = 0., maxrange = max(Calib_img),  /VERTICAL, $
       title = 'Continuum at '+ string(Flux_Cal_Wavelength, format = '(I5)') + cgsymbol('Angstrom') + ' (MR/'+ cgsymbol('Angstrom')+')'
    
    P = pos[*,3]
    p[0] = p[0] - .01
    p[2] = p[2] - .01
    loadct, 3
    axis_format = {XTICKFORMAT:"(A1)", YTICKFORMAT:"(A1)"}
    ;cgimage, smooth((Na_D1_img + Na_D2_img), [0, 2]), /keep_aspect, minvalue = 0, maxvalue = 2.5, pos = P, /noerase, /axes, axkeywords = axis_format
    cgimage, alog(smooth((Na_D1_img + Na_D2_img), [0, 2])), /keep_aspect, minvalue = 0, maxvalue = alog(max(smooth((Na_D1_img + Na_D2_img), [0, 2]))), pos = P, /noerase, /axes, axkeywords = axis_format
    ;cgimage, alog(smooth((Lucky_Na_D1_img + Lucky_Na_D2_img), [0, 2])), /keep_aspect, minvalue = 0, maxvalue = alog(2.5), pos = P, /noerase, /axes, axkeywords = axis_format
    cgcontour, smooth(Calib_img / max(Calib_img, /NaN), [0, 2]), /onimage, pos = P, levels=levels, C_CHARSIZE=2.4, color='green', $
    ;cgcontour, smooth(Lucky_Calib_img / max(Lucky_Calib_img, /NaN), [0, 2]), /onimage, color='green', pos = P, levels=levels, C_CHARSIZE=2.4, $
      C_annotation = string(levels)+' MR/'+cgsymbol('Angstrom'), C_LABELS=[1,1,1,1,1], thick = 0.5
    p = [p[0] + xspace+0.03, p[1]+.02, p[0]+xspace+0.04, p[3]-.02]
    cgColorbar, minrange = 0., maxrange = max(smooth((Na_D1_img + Na_D2_img), [0, 2])), pos = P, /VERTICAL, title = 'Sodium D1 + D2 (MR)', color = 'white'
    
    print, 'Hack! Color bar of sodiumhas wrong tick labels)
    
    cgText, 0.5, 0.9, 'RIPS --- ' + strmid(UTC,0,10) + ' ' + strmid(UTC,11,5) + ' --- True Anomaly '+ string(True_Anomaly, format = '(F4.1)')+ '!U'+cgsymbol('deg')+'!U' , ALIGNMENT=0.5, /NORMAL, charsize = 1.9
    cgText, 0.09, 0.845, 'Reflectance Model', /NORMAL, charsize = 0.9, color = 'white'
    cgText, 0.51, 0.845, 'Model w/ Effective Seeing = '+string(effective_seeing, format = '(F3.1)')+'"', /NORMAL, charsize = 0.9, color = 'red'
    cgText, 0.51, 0.8225, 'Co-Aligned Imaging', /NORMAL, charsize = 0.9, color = 'Blue'
    cgText, 0.51, 0.80, 'Co-Aligned Spectra', /NORMAL, charsize = 0.9, color = 'Black'
  cgPS_close
 
endif
stop
print, 'Run time = ',  (SYSTIME(/SECONDS) - Start_time) / 60, ' minutes'
end



