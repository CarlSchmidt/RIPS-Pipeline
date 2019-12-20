function Estimate_Seeing, FWHM_arcsec ; FWHM Arcseconds <----solve for this iteratively
  Common Seeing_Estimate, A, B, Spectral_platescale, Imaging_platescale, x_fine, y_fine, Seeing_Est_Rebin_Factor

  ; A and B images being compared have already be roughly shifted/aligned, but we want to know if that shift changes with the seeing value
  seeing_sigma      = FWHM_arcsec / (2.0*sqrt(2.0*alog(2.)))      ; convert FWHM to sigma
  seeing_sigma      = seeing_sigma/Imaging_platescale             ; convert to Gaussian sigma in pixels
  blur_kernel       = GAUSSIAN_FUNCTION([seeing_sigma, seeing_sigma])
  if N_elements(blur_kernel) gt N_elements(B) then begin          ; crop the kernel when seeing is really bad and its larger than the image
    SA = size(A, /dim)
    SK = size(blur_kernel, /dim)
    blur_kernel = blur_kernel[sk[0]/2.-SA[0]/2.:sk[0]/2.+SA[0]/2.-1, sk[1]/2.-SA[1]/2.:sk[1]/2.+SA[1]/2.-1]
  endif
  Blurred_B         = convol(B, blur_kernel, /EDGE_ZERO, /normalize, /NAN)
  zstack_align_images, A/max(A), Blurred_B/max(Blurred_B), x_fine, y_fine, corr_image, corr_dims
  CORREL            = CORREL_IMAGES( A/max(A), Blurred_B/max(Blurred_B), xoffset_b = x_fine, yoffset_B = y_fine, XSHIFT = 0, ySHIFT = 0) 

   ;debug = 1
    if keyword_set(debug) then begin
      s = size(A, /dimensions)
      window, 0, xs = s[0], ys = s[1], xpos = 0, ypos = 0
      tv, bytscl(A/max(A), 0, 1)
      window, 1, xs = s[0], ys = s[1], xpos = 0+s[0], ypos = 0
      tv, bytscl(Blurred_B/max(Blurred_B), 0, 1)
      print, [FWHM_arcsec*Seeing_Est_Rebin_Factor, 1./correl]
    endif

  Stat  =  1./correl ; Minimizing this statistic maximizes the correlation
  return, stat       ; Scalar Value for amoeba to minimize
end

FUNCTION Fit_Spec, X, P
  RETURN, P[0]*X + P[1]
END

;FUNCTION Match_Reference_Spectrum, X, P
;  ; Multiply P[0], Add P[1], Shift P[2] and Smooth P[3] a reference spectrum (X) until it best matches Y
;    shifted = interpolate(X, findgen(n_elements(x)) - P[2])      ; sub-pixel shifting
;    return, P[0]*gauss_smooth(shifted, p[3], /EDGE_TRUNCATE) + P[1] + P[4]*findgen(n_elements(x))
;end

pro RIPS_Mercury_Na, PART=part, NIGHT=night, Telescope = telescope
  Common Seeing_Estimate, A, B, Spectral_platescale, Imaging_platescale, x_fine, y_fine, Seeing_Est_Rebin_Factor
  ; *********************************************************************************************************************
  ; *********************************************************************************************************************
  ; Routine to extract Na emission from spectral scans over Mercury's disk with Perkins/RIPS in order to create a 2-D
  ; image of Na above and near Mercury's disk.
  ;
  ; THIS REDUCTION PROGRAM IS SPLIT INTO SEVERAL PARTS:
  ;  -1 = Calibrate the imaging plate scale and the spectral to imaging channel plateascale ratio 
  ;       input:   some Io and Europa frames from April 2018
  ;       output:  "/pixel for the Perkins telescope, AEOS PLATE SCALE IS YET UNKNOWN
  ;   0 = break kinetic series into manageable datacubes
  ;       input:   list of kinetic series FITS file(s), and calibration for a given night
  ;       output:  RIPS imaging (x,y,t) and spectral (wl,y,t) datacubes, also imaging with slit interpolated over (x,y,t)
  ;   1 = Make a theoretical image of Mercury using Hapke's fomulation 
  ;       input:   header of the imaging and spectra data cubes from part 0, used for timestamps and array sizes 
  ;       output:  a fake image of the planet in the plane-of-sky in MR/A units, matched to the plate scale of each channel 
  ;   2 = Find the rotation angle between the imaging channel and the plane-of-the-sky
  ;       input:   slitless imaging datcube from part 0 and a fake Mercury from Part 1
  ;       output:  An array of the rotational alignments in degrees vs frame number 
  ;   3 = find Mercury centroids by cross-corellation with previous image
  ;       input:   imaging datacube from Part 0
  ;       output:  co-aligned imaging datacube (x,y,t)
  ;   4 = isolate sodium emission in every frame, Do D1 & D2 separately
  ;       input:   spectral datacube from Part 0
  ;       output:  the calibrated spectral datacube (wl,y,t) for Na exosphere
  ;   5 = extract Na signal and place into 1D spectra
  ;       input:   spectral datacube from Part 2 (wl,y,t) with surface reflectance removed
  ;       output:  Na brightness (wl,y) and linewidth (wl,y)
  ;   6 = extract Surface reflectance and compare with a Hapke Model. Flux calibrate the data
  ;       input:   Na brightness and linewidth from Part 3 and imaging cube from Part 0
  ;       output:  Hapke Model, calibrated Na D images in Rayleighs and continuum in R/A
  ;   7 = plot everything
  ;       input:   save files from part 6
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
  ;           19 Aug 2019 - allow multiple rotation angles on the sky to be aligned and added (CS)

  Start_time = SYSTIME(/SECONDS)

  ; =====================================================================================================================
  ; Define all variables
  ; =====================================================================================================================
  Body          = 'Mercury'                                                 ; For SPICE purposes
 
  ;SetDefaultValue, Telescope,'AEOS'
  ;SetDefaultValue, Telescope,'Perkins'
  SetDefaultValue, thresh, 12                                               ; hotpixel removal threshold (sigma above background)
  SetDefaultValue, width, 5                                                 ; Median smoothing width to get local background reference (for pixels above thresh)
  SetDefaultValue, ct, 22                                                   ; default color table
  SetDefaultValue, Flux_Cal_Wavelength, 5893.                               ; Wavelength to extract and compare to Hapke code (Angstroms)
  SetDefaultValue, Na_D1_WL, 5895.92424                                     ; Rest wavelength from NIST
  SetDefaultValue, Na_D2_WL, 5889.95095                                     ; Rest wavelength from NIST
  SetDefaultValue, ims, [283,759,39,374]                                    ; bounds of the imaging channel: x1,x2,y1,y2 
  SetDefaultValue, sps, [99,899,502,941]                                    ; bounds of the spectra channel: x1,x2,y1,y2 
  
  case 1 of
    Telescope eq 'Perkins': begin
      Spectral_platescale     = 0.106723                                        ; SPECTRAL channel "/pix from measured imaging / spectra ratio in the lab in part eq '-1'
      imaging_platescale      = 0.119421                                        ; IMAGING channel "/pix from measured Io - Europa in part eq '-1'
      Seeing_Est_Rebin_Factor = 1.                                              ; Rebinning during the seeing estimate & Hapke convolution makes things go alot faster, use 2. or 4.
      ;Seeing_Est_Rebin_Factor = 4.                                              ; use for debug, speeds up part 3. 
      subframe                = 50                                              ; half width of the subframe for the final images (pixels)
      Local = FILE_SEARCH('D:\DATA\Perkins\Perkins RIPS - March 2018\*')        ; Carl's Computer
      Cloud = FILE_SEARCH('Y:\obs_18\Perkins_RIPS_March\*')                     ; Luke's Cloud
      if local[0] ne '' then base_dir = 'D:\DATA\Perkins\Perkins RIPS - March 2018\'
      if Cloud[0] ne '' then base_dir = 'Y:\obs_18\Perkins_RIPS_March\'
    end
    Telescope eq 'AEOS': begin
      Spectral_platescale     = 0.106723 / 2.3                                  ; 2.35 factor is only approximate, estimated visually 
      imaging_platescale      = 0.119421 / 2.3                                  ; 2.35 factor is only approximate, estimated visually             
      Seeing_Est_Rebin_Factor = 2.                                              ; Rebinning during the seeing estimate & Hapke convolution makes things go alot faster, use 2 or 4.
      subframe                = 120                                             ; half width of the subframe for the final images (pixels)
      Local = FILE_SEARCH('D:\DATA\AEOS\RIPS\*')                                ; Carl's Computer
      Cloud = FILE_SEARCH('?')                                                  ; Luke's Cloud
      if local[0] ne '' then base_dir = 'D:\DATA\AEOS\RIPS\'
      if Cloud[0] ne '' then base_dir = '?'
    end
  endcase
  
  if Telescope eq 'AEOS' then begin
    case 1 of 
      night eq '10': begin ;SPECTRAL CHANNEL OUT OF FOCUS, ghosting of an internal reflection visible in the spectrum, Never got this working, but maybe something useful
        Mercury_dir   = base_dir+'Dec 10\'                                        ; directory with Mercury FITS files
        dark_dir      = base_dir+'Dec 10\'                                        ; directory with dark file
        flat_dir      = base_dir+'Dec 11\'                                        ; directory with flat file
        sky_dir       = base_dir+'June 23 HST\'                                   ; directory with sky/solar file
        arc_dir       = base_dir+'Dec 11\'                                        ; directory with arc files
        Mercury_files = 'RIPS_setup_' + strcompress((indgen(8)+316), /remove_all) + '.fits'         ; Mercury kinetic series frames 
        sky_file      = 'RIPS_setup_174.fits'                                     ; sky/solar frame (164=sky Na, 174=lunar drift Na) ;HACK! this is with a different slitwidth, needs inspection at least!
        dark_file     = 'RIPS_setup_314.fits'                                     ; dark frame
        flat_files    = 'RIPS_setup_' + strcompress((indgen(3)+353), /remove_all) + '.fits'   
        Arc_files     = 'RIPS_setup_352.fits'                                     ; Arc frames
        Lucky_fraction= 0.25                                                      ; Fraction of images within the cutoff we're calling "Lucky"
        effective_seeing = 3.0                                                    ; Initial Estimate of the effective Seeing FWHM in Arcseconds, This is later solved for iteratively
        Smooth_Final  = [4,4]                                                     ; The final sodium image [x, y] smoothing, [0,0] = no smooth 
        plot_lucky_Na = 'true'                                                    ; Is the final sodium image showing just the lucky_fraction subset (= 'true') or everything (comment out)
        LSF = 'sum'                                                               ; How to integrated the line spread function, 'fit' = area under gaussian, 'sum' = total over search/2 bandpass from the trace 
        slit_interp_width = 14                                                    ; In part 0, the # of pixels used to interpolate over the slit, visually inspect this and make it as small as possible 
      end
      night eq '12': begin ;GREAT!
        Mercury_dir   = base_dir+'Dec 12\'                                        ; directory with Mercury FITS files
        dark_dir      = base_dir+'Dec 10\'                                        ; directory with dark file
        flat_dir      = base_dir+'Dec 11\'                                        ; directory with flat file
        sky_dir       = base_dir+'June 23 HST\'                                   ; directory with sky/solar file
        arc_dir       = base_dir+'Dec 11\'                                        ; directory with arc files
        Mercury_files = 'RIPS_setup_' + strcompress((indgen(7)+414), /remove_all) + '.fits'         ; Mercury kinetic series frames 
        sky_file      = 'RIPS_setup_174.fits'                                     ; Must use external solar spectrum for reflectance 
        dark_file     = 'RIPS_setup_314.fits'                                     ; dark frame
        flat_files    = 'RIPS_setup_' + strcompress((indgen(3)+353), /remove_all) + '.fits'    
        Arc_files     = 'RIPS_setup_352.fits'                                     ; Arc frames
        Lucky_fraction= 0.30                                                      ; Fraction of images within the cutoff we're calling "Lucky"
        effective_seeing = 2.2                                                    ; Initial Estimate of the effective Seeing FWHM in Arcseconds, This is later solved for iteratively
        Smooth_Final  = [4,4]                                                     ; The final sodium image [x, y] smoothing, [0,0] = no smooth 
        plot_lucky_Na = 'true'                                                    ; Is the final sodium image showing just the lucky_fraction subset (= 'true') or everything (comment out)
        LSF = 'fit'                                                               ; How to integrated the line spread function, 'fit' = area under gaussian, 'sum' = total over search/2 bandpass from the trace 
        slit_interp_width = 14                                                    ; In part 0, the # of pixels used to interpolate over the slit, visually inspect this and make it as small as possible 
      end
      night eq '13': begin ;Poor Seeing
        Mercury_dir   = base_dir+'Dec 13\'                                        ; directory with Mercury FITS files
        dark_dir      = base_dir+'Dec 10\'                                        ; directory with dark file
        flat_dir      = base_dir+'Dec 11\'                                        ; directory with flat file
        sky_dir       = base_dir+'June 23 HST\'                                   ; directory with sky/solar file
        arc_dir       = base_dir+'Dec 11\'                                        ; directory with arc files
        Mercury_files = 'RIPS_setup_' + strcompress((indgen(7)+467), /remove_all) + '.fits'         ; Mercury kinetic series frames
        sky_file      = 'RIPS_setup_174.fits'                                     ; Must use external solar spectrum for reflectance
        dark_file     = 'RIPS_setup_314.fits'                                     ; dark frame
        flat_files    = 'RIPS_setup_' + strcompress((indgen(3)+353), /remove_all) + '.fits'
        Arc_files     = 'RIPS_setup_352.fits'                                     ; Arc frames
        Lucky_fraction= 0.30                                                      ; Fraction of images within the cutoff we're calling "Lucky"
        effective_seeing = 2.2                                                    ; Initial Estimate of the effective Seeing FWHM in Arcseconds, This is later solved for iteratively
        Smooth_Final  = [4,4]                                                     ; The final sodium image [x, y] smoothing, [0,0] = no smooth
        plot_lucky_Na = 'true'                                                    ; Is the final sodium image showing just the lucky_fraction subset (= 'true') or everything (comment out)
        LSF = 'fit'                                                               ; How to integrated the line spread function, 'fit' = area under gaussian, 'sum' = total over search/2 bandpass from the trace
        slit_interp_width = 14                                                    ; In part 0, the # of pixels used to interpolate over the slit, visually inspect this and make it as small as possible
      end
      night eq '20': begin ;DECENT DATA, BUT NOT VERY MUCH COVERAGE
        mercury_dir = base_dir+'June 20 HST\'                                     ; directory with Mercury FITS files
        dark_dir    = base_dir+'June 21 HST\'                                     ; directory with dark file
        flat_dir    = base_dir+'June 21 HST\'                                     ; directory with flat file
        sky_dir     = base_dir+'June 23 HST\'                                     ; directory with sky file
        arc_dir     = base_dir+'June 21 HST\'                                     ; directory with arc files
        Mercury_files = ['Mercury - 500 kinetic 0.5 sec scan - Na 3A + ND imaging - Na order.fits',$
                         'Mercury - 500 kinetic 1 sec scan - Na 3A + ND imaging - Na order (2).fits'] ; Mercury kinetic series frames  ; 
        sky_file      = 'RIPS_setup_174.fits'                                     ; Lunar drift scan. Sky frames have too many absorption lines! 175 slit motor setting was unchanged this entire dataset
        dark_file     = 'Dark - 200 sec - Na 3A + ND imaging - Na order.fits'     ; dark frame
        flat_files    = 'Tungsten lamp flat - 200 sec - Na 3A + ND imaging - Na order_Autosave_recovered.fits'  ; flat frame
        Arc_files     = 'Ne lamp default motor settings - 200 sec - Na 3A + ND imaging - Na order.fits'         ; Arc frames
        Lucky_fraction= 0.30                                                      ; Fraction of images within the cutoff we're calling "Lucky"
        effective_seeing = 2.2                                                    ; Initial Estimate of the effective Seeing FWHM in Arcseconds, This is later solved for iteratively
        Smooth_Final  = [4,4]                                                     ; The final sodium image [x, y] smoothing, [0,0] = no smooth 
        plot_lucky_Na = 'true'                                                    ; Is the final sodium image showing just the lucky_fraction subset (= 'true') or everything (comment out)
        LSF = 'fit'                                                               ; How to integrated the line spread function, 'fit' = area under gaussian, 'sum' = total over search/2 bandpass from the trace 
        slit_interp_width = 14                                                    ; In part 0, the # of pixels used to interpolate over the slit, visually inspect this and make it as small as possible                                             ; Fraction of images within the cutoff we're calling "Lucky"
      end  
      night eq '22': begin ;DECENT DATA BUT FIXED POINTING, SO SPATIAL INFORMATION IS TELESCOPE DRIFT ONLY
        mercury_dir = base_dir+'June 22 HST\'                                     ; directory with Mercury FITS files
        dark_dir    = base_dir+'June 21 HST\'                                     ; directory with dark file
        flat_dir    = base_dir+'June 25 HST\'                                     ; directory with flat file
        sky_dir     = base_dir+'June 23 HST\'                                     ; directory with sky file
        arc_dir     = base_dir+'June 21 HST\'                                     ; directory with arc files
        Mercury_files = 'RIPS_setup_' + strcompress((indgen(7)+135), /remove_all) + '.fits'  ; Mercury kinetic series frames (135-141) no pushbroom scanning
        sky_file      = 'RIPS_setup_174.fits'                                     ; Lunar drift scan. Sky frames have too many absorption lines! 175 slit motor setting was unchanged this entire dataset
        dark_file     = 'Dark - 200 sec - Na 3A + ND imaging - Na order.fits'     ; dark frame
        flat_files    = 'RIPS_setup_205.fits'  ; flat frame
        Arc_files     = 'Ne lamp default motor settings - 200 sec - Na 3A + ND imaging - Na order.fits'         ; Arc frames
        Lucky_fraction= 0.30                                                      ; Fraction of images within the cutoff we're calling "Lucky"
        effective_seeing = 2.2                                                    ; Initial Estimate of the effective Seeing FWHM in Arcseconds, This is later solved for iteratively
        Smooth_Final  = [4,4]                                                     ; The final sodium image [x, y] smoothing, [0,0] = no smooth
        plot_lucky_Na = 'true'                                                    ; Is the final sodium image showing just the lucky_fraction subset (= 'true') or everything (comment out)
        LSF = 'fit'                                                               ; How to integrated the line spread function, 'fit' = area under gaussian, 'sum' = total over search/2 bandpass from the trace
        slit_interp_width = 14                                                    ; In part 0, the # of pixels used to interpolate over the slit, visually inspect this and make it as small as possible
      end
      night eq '25': begin ;LOOKS GREAT, BEST YET AEOS RESULTS!
        Mercury_dir   = base_dir+'June 25 HST\'                                   ; directory with Mercury FITS files
        dark_dir      = base_dir+'June 25 HST\'                                   ; directory with dark file
        flat_dir      = base_dir+'June 25 HST\'                                   ; directory with flat file
        sky_dir       = base_dir+'June 23 HST\'                                   ; directory with sky/solar file
        arc_dir       = base_dir+'June 21 HST\'                                   ; directory with arc files
        Mercury_files = 'RIPS_setup_' + strcompress((indgen(6)+210), /remove_all) + '.fits'         ; Mercury kinetic series frames (210-215)
        sky_file      = 'RIPS_setup_174.fits'                                     ; sky/solar frame (164=sky Na, 174=lunar drift Na)
        dark_file     = 'RIPS_setup_203.fits'                                     ; dark frame
        flat_files    = 'RIPS_setup_205.fits'                                     ; flat frame
        Arc_files     = 'Ne lamp default motor settings - 200 sec - Na 3A + ND imaging - Na order.fits'         ; Arc frames
        Lucky_fraction= 0.3                                                       ; Fraction of images within the cutoff we're calling "Lucky"
        effective_seeing = 2.2                                                    ; Initial Estimate of the effective Seeing FWHM in Arcseconds, This is later solved for iteratively
        Smooth_Final  = [4,4]                                                     ; The final sodium image [x, y] smoothing, [0,0] = no smooth 
        plot_lucky_Na = 'true'                                                    ; Is the final sodium image showing just the lucky_fraction subset (= 'true') or everything (comment out)
        LSF = 'fit'                                                               ; How to integrated the line spread function, 'fit' = area under gaussian, 'sum' = total over search/2 bandpass from the trace 
        slit_interp_width = 14                                                    ; In part 0, the # of pixels used to interpolate over the slit, visually inspect this and make it as small as possible
      end   
    endcase 
  endif
  if Telescope eq 'Perkins' then begin
    case 1 of
      night eq '15': begin
        mercury_dir   = base_dir+'15\'                                            ; directory with Mercury FITS file(s)
        dark_dir      = base_dir+'14\Carl_Keep\'                                  ; directory with dark file
        flat_dir      = base_dir+'15\'                                            ; directory with flat file
        sky_dir       = base_dir+'15\'                                            ; directory with Sky file
        Arc_dir       = base_dir+'15\'                                            ; directory with arc file
        Arc_files     = 'Neon_NaSpectra_624slitwidth_NaND1imaging*.fits'          ; Arc frames
        Mercury_files = '*'+strcompress(788+indgen(17), /remove_all) + '.fits'    ; SUFFIX of the Mercury kinetic series to use (note 788 is at a different rotiserizer angle)
        sky_file      = 'RIPS1_Thu Mar 15 2018_01.12.12_785.fits'                 ; sky frame
        dark_file     = 'Dark.fits'                                               ; dark frame
        flat_files    = 'Flat_NaSpectra_624slitwidth_NaND1imaging*fits'           ; flat frames
        effective_seeing= 2.7                                                     ; Inital Estimate of the effective Seeing FWHM in Arcseconds, This is later solved for iteratively
        Lucky_fraction= 0.75                                                      ; Fraction of images within the cutoff we're calling "Lucky" 
        plot_lucky_Na = 'true'                                                    ; Is the final sodium image showing just the lucky_fraction subset (= 'true') or everything (comment out)
        Smooth_Final  = [0,0]                                                     ; The final sodium image [x, y] smoothing, [0,0] = no smooth 
        LSF = 'fit'                                                               ; How to integrated the line spread function, 'fit' = area under gaussian, 'sum' = total over search/2 bandpass from the trace 
        slit_interp_width = 14                                                    ; In part 0, the # of pixels used to interpolate over the slit, visually inspect this and make it as small as possible
      end
    endcase  
  endif

  plate_scale_ratio  = Spectral_platescale / imaging_platescale            
  outdir             = Mercury_dir + 'Processed\'                                 ; output directory
  
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

; =====================================================================================================================
; Find display size, set window defaults, reset plotting variables just in case
; =====================================================================================================================
    cgloadct, ct
    CLEANPLOT                                                               ; Reset all plotting system variables (!P,!X,!Y,!Z) to their default values
    DEVICE, GET_SCREEN_SIZE = ss                                            ; find screen size
    winpos_x      = 0                                                       ; default initial x window location
    winpos_y      = 0                                                       ; default initial y window location

; =====================================================================================================================
; Part -1 : Calculate the imaging platescale, and get the spectral plate scale from their ratio in data of a fiducial
;           NOTE at Perkins in the imaging Channel North is up (+x), West is right (+y) with mrdfits.pro reader
; =====================================================================================================================
if part eq '-1' then begin ; This section calibrates angular plate scaling for each telescope. Only needs to be run once

  ; Determine an accurate plate scale at Perkins, based on the Io-Europa separation in the April 2018 data.
  Io_X     = float( [ 525, 520, 527, 519, 519, 519, 516, 525, 524, 519 ] )
  Europa_X = float( [ 609, 622, 659, 684, 681, 681, 676, 752, 770, 793 ] )
  Io_Y     = float( [ 282, 283, 269, 267, 263, 261, 259, 235, 232, 229 ] )
  Europa_Y = float( [ 46,   88, 148, 293, 343, 341, 389, 323, 364, 422 ] )
  Platescale_Filenames = ['Io_Jet_JovianDawnSide_K_13.fits'        , 'Io_Jet_JovianDawnSide_K_12.fits',       'Io_Jet_JovianDawnSide_Na_11.fits', $
    'Io_Equatorial_JovianDawnSide_Na_10.fits', 'Io_Equatorial_JovianDawnSide_K_9.fits', 'Io_Equatorial_JovianDawnSide_K_8.fits', $
    'Io_Jet_JovianDawnSide_K_7.fits'         , 'Io_Jet_JovianDawnSide_K_6.fits',        'Io_Jet_JovianDawnSide_Na_5.fits' ]
  Imaging_platescale_array = fltarr( N_elements( Platescale_Filenames ) )
  CSPICE_FURNSH, STRCOMPRESS(FILE_SEARCH(kernel_directory, 'Jup310.bsp'))   ; Jovian Moons Kernel
  for i = 0, N_elements( Platescale_Filenames )-1 do begin
    header = headfits( 'D:\DATA\Perkins\Perkins RIPS - April 2018\' + Platescale_Filenames[i] )
    timestamp = sxpar(header, 'DATE')
    Cspice_UTC2et, sxpar(header, 'DATE'), et
    et = et + float(sxpar(header, 'EXPOSURE')) * 0.5 * sxpar(header, 'NUMKIN')
    cspice_spkezr, 'Io', ET, 'J2000', 'LT+S', 'Earth', Io_Earth_State, ltime
    cspice_spkezr, 'Europa', ET, 'J2000', 'LT+S', 'Earth', Europa_Earth_State, ltime
    theta  = cspice_vsep(Europa_Earth_State[0:2], Io_Earth_State[0:2])
    Imaging_platescale_array[i] = theta*206265. / sqrt( (Io_X[i] - Europa_X[i])^2 + (Io_Y[i] - Europa_Y[i])^2 )
  endfor
  Imaging_platescale = mean(Imaging_platescale_array) ; "/pixel on the Perkins telescope
  Plate_Scale_Ratio  = (426.-73.) / (997.-602.)       ; 1/pixel ratio of Spectral_platescale / imaging_platescale using D:\DATA\___Calibration___\RIPS_Plate_Scale_Ratio.fits
  Spectral_platescale = Imaging_platescale * Plate_Scale_Ratio
  print, 'Perkins Telescope Spectral channel platescale: ', Spectral_platescale, '"/pixel'
  print, 'Perkins Telescope Imaging channel platescale: ', imaging_platescale, '"/pixel'
  print, 'Input these values above since they are universal'
endif

; =====================================================================================================================
; Part 0 : Perform basic corrections. Break kinetic series of mutliple frames into 2 manageable-sized datacubes 
;          for each channel. Interpolate over the slit and the write result to a seperate datacube. 
; =====================================================================================================================
if part eq 0 or part eq 99 then begin

  ; Prepare the master dark
    Dark = MRDFITS(dark_dir + dark_file, 0, header, /fscale, /silent )   ; read in dark
    if n_elements(size(dark, /dimensions)) gt 2 then dark = median(dark, dimension = 3)
    dark  = sigma_filter( dark, width, N_sigma=thresh)
    Imaging_dark = dark[ims[0]:ims[1],ims[2]:ims[3]]
    Spectra_Dark = dark[sps[0]:sps[1],sps[2]:sps[3]]
    MWRFITS, dark, outdir + 'MASTER_DARK_FULL_FRAME.fits', header, /create, /silent
    MWRFITS, Imaging_dark, outdir + 'Imaging_dark.fits', header, /create, /silent
    MWRFITS, Spectra_Dark, outdir + 'Spectra_Dark.fits', header, /create, /silent

  ; Prepare the Arc Spectrum
    arc_frames = file_Search(Arc_Dir + arc_files, count = n_arcs)
    big_array = fltarr(1024, 1024, n_arcs)
    for i = 0, n_arcs-1 do big_array[*,*,i] = MRDFITS(arc_frames[i], 0, header, /fscale, /silent )   ; read in arcs
    if n_arcs gt 1 then Neon = mean(big_array, dimension=3) else Neon = big_array
    Neon = Neon - dark
    Neon = Neon[*,sps[2]:sps[3]]                                          ; Crop Arc Lamp to just Spectral Channel
    MWRFITS, Neon, outdir + 'NEON_ARC_LAMP.fits', header, /create, /silent

  ; Prepare the Imaging Channel Flat
    flat_frames = file_Search(Flat_Dir + flat_files, count = n_flats)
    big_array = fltarr(1024, 1024, N_flats)
    for i = 0, N_flats-1 do big_array[*,*,i] = MRDFITS(flat_frames[i], 0, header, /fscale, /silent ) ; read in flats
    if n_flats gt 1 then begin                                            ; combine the flat (if there's multiple and that's needed)
      Imaging_Flat = median(big_array, dimension=3, /even) - dark         ; dark subtract the flat
      Spectra_Flat = median(big_array, dimension=3, /even) - dark         ; dark subtract the flat
    endif else begin
      Imaging_Flat = big_array - dark  
      Spectra_Flat = big_array - dark  
    endelse                             
    Imaging_Flat = Imaging_Flat[ims[0]:ims[1],ims[2]:ims[3]]              ; Crop to just imaging portion
    Spectra_Flat = Spectra_Flat[sps[0]:sps[1],sps[2]:sps[3]]              ; Crop to just Spectra portion
    Imaging_Flat = imaging_Flat / max(imaging_Flat)                       ; Unsure why I'm doing this
    Imaging_Flat = Imaging_flat + (1. - median(Imaging_flat))             ; this normalizes such that the median of the flat is 1
    Imaging_flat[WHERE(Imaging_flat lt .01, /NULL)] = !values.f_Nan       ; reject unusual counts for centroid
    Imaging_flat[WHERE(Imaging_flat gt 2.0, /NULL)] = !values.f_Nan       ; reject unusual counts for centroid
    gooddata = where(Finite(Imaging_flat), ngooddata, comp=baddata, ncomp=nbaddata)
    if nbaddata gt 0 then Imaging_flat[baddata] = interpol(Imaging_flat[gooddata], gooddata, baddata) ; Interpolate over the slit??????

      ; Find indicies of all pixels under the slit in the Imaging Flat. Interpolate over them.
        s = size(Imaging_flat, /dimensions)
        junk = min(total(Imaging_Flat, 2), slit_center)                             ; roughly find the slit center
        case telescope of 
          'Perkins' : begin
            h    = histogram(Imaging_Flat, binsize = .05, REVERSE_INDICES=ri)       ; bin the Imaging_Flat into 0.05 bins
            slit_indices_1D = ri[ri[0]:ri[8]-1]                                     ; find inidices of pixels in the lowest few bins of the Imaging_Flat's histogram ---> INSPECT slit interpolation in the data --> set to 8
            slit_indices = array_indices(Imaging_Flat, slit_indices_1D)             ; convert into x & y indices
            keep = where(abs(slit_indices[0,*] - slit_center) lt slit_interp_width, /Null)  ; keep only indices within 14 pixels of slit center
          end  
          'AEOS' : begin
            h    = histogram(Imaging_Flat, binsize = .05, REVERSE_INDICES=ri)      ; bin the Imaging_Flat into 0.025 bins
            slit_indices_1D = ri[ri[0]:ri[10]-1]                                    ; find inidices of pixels in the lowest few bins of the Imaging_Flat's histogram  ---> INSPECT slit interpolation in the data --> set to 28?
            slit_indices = array_indices(Imaging_Flat, slit_indices_1D)             ; convert into x & y indices
            keep = where(abs(slit_indices[0,*] - slit_center) lt slit_interp_width, /Null)         ; keep only indices within 8 pixels of slit center
          end   
        endcase
        slit_indices = slit_indices[*,keep] & slit_indices_1D = slit_indices_1D[keep]
        dummy = Imaging_Flat
        dummy[ slit_indices[0,*], slit_indices[1,*] ] = max(Imaging_Flat)           ; show the location of pixels fully under the slit.
        window, 0, xs = s[0], ys = s[1], xpos=0, ypos=0, title = 'Imaging_Flat Field Inspection: Red = Behind Slit --> interpolate'
        cgimage, hist_equal(dummy) 
        dummy[ slit_indices[0,*], slit_indices[1,*] ] = !Values.F_NaN
        sm = SMOOTH( dummy, slit_interp_width, /EDGE_TRUNCATE, /NAN)
        dummy[slit_indices_1D] = sm[slit_indices_1D]                      ; interpolate over slit                                  
        window, 1, xs = s[0], ys = s[1], title = 'Imaging_Flat Field Inspection: Slit be gone!'
        cgimage, hist_equal(dummy)
        Imaging_Flat_Slit_Interpolated = temporary(dummy)

  ; Prepare a spectral channel flat
    ispectra_flat     = spectra_flat / median(spectra_flat)
    acre, ispectra_flat, spectra_flat, thresh, width
    spectra_flat      = spectra_flat + (1. - median(spectra_flat))        ; this normalizes again such that median(flat) = 1

  ; Write the Mercury data cubes
    filenames = file_Search(Mercury_Dir + Mercury_files, count = n_flats)
    header            = headfits(filenames[0])
    integration       = sxpar(header, 'EXPOSURE')                         ; integration time for individual frames
    n_frames_per_file = sxpar(header, 'NUMKIN')                           ; # of frames in the "kinetic series" for each "Mercury_file"
    nfiles            =  n_elements(Mercury_files)                        ; # of different Mercury kinetic frames
    s    = [1024, 1024, n_frames_per_file]
    cube = intarr(s[0], s[1], s[2]*nfiles)                                ; define a **big** cube, where multiple kinetic series are stacked in the third dimension, INTEGER TYPE for now
    print, 'Found: ', N_elements(filenames), ' Mercury files'
    window, 0, xs = ims[1] - ims[0], ys = ims[3]-ims[2], title = 'Inspection: First frame, Imaging Channel'
    integration = fltarr(nfiles)
    for ifile = 0, nfiles-1 do begin                                      ; loop over kinetic series and add each to a stack
      icube = MRDFITS(filenames[ifile], 0, header, /unsigned, /silent )
      integration[ifile]       = sxpar(header, 'EXPOSURE')
      cube[*,*,s[2]*ifile:s[2]*ifile+n_frames_per_file-1] = icube   
      n_frames_per_file = sxpar(header, 'NUMKIN')
      print, 'Loading: '+Mercury_files[ifile],' containing', strcompress(n_frames_per_file)+' frames of'+strcompress(integration[ifile])+' seconds each'
      tv, bytscl(icube[ ims[0]:ims[1], ims[2]:ims[3], 0])
    endfor
    imaging_cube = float(reform(cube[ims[0]:ims[1],ims[2]:ims[3], *]))    ; extract the IMAGING portion of the kinetic series and combine
    spectra_cube = float(reform(cube[sps[0]:sps[1],sps[2]:sps[3], *]))    ; extract the SPECTRAL portion of the kinetic series and combine
    si        = size(imaging_cube, /dimensions)                           ; size of extracted IMAGING portion of frame
    ss        = size(spectra_cube, /dimensions)                           ; size of extracted SPECTRAL portion of frame
    imaging_cube = imaging_cube - rebin(Imaging_dark, si)                 ; Dark (Bias) subtract the data
    spectra_cube = spectra_cube - rebin(Spectra_dark, ss)                 ; Dark (Bias) subtract the data

    ; Now use the data to fine tune any pesky imperfections in the flat fielding... 

      ; Occasionally, if the flat is separated in time, there cann exist a few pixel shift between the flat and imaging channel science data. Try to correct for this:
        A         = total(imaging_cube, 3, /nan)
        B         = Imaging_flat
        x_crop    = minmax(slit_indices[0,*]) + [-100, 100]               ; region most sensitive to any potential shifts
        y_crop    = minmax(where(total(A,1) ge 0.2*max(total(A,1)),/Null)); Good S/N range 
        A         = bandpass_filter(A, 0.5, 1.)                           ; High pass filter things to isolate the small scale spatial structure
        B         = bandpass_filter(B, 0.5, 1.)                           ; High pass filter things to isolate the small scale spatial structure    
        A         = A[x_crop[0]:x_crop[1], y_crop[0]:y_crop[1]]
        B         = B[x_crop[0]:x_crop[1], y_crop[0]:y_crop[1]]
        zstack_align_images, B, A, shift_img_flat_x, shift_img_flat_y, corr_image, corr_dims
        if (shift_img_flat_x ne 0.) or (shift_img_flat_y ne 0.) then begin
          print, 'Shifting imaging channel''s flat field by:', fix([shift_img_flat_x, shift_img_flat_y])
          wset, 0
          tv, bytscl(total(imaging_cube, 3, /nan) / shift(Imaging_flat, 0, 0) )
          wset, 1
          tv, bytscl(total(imaging_cube, 3, /nan) / shift(Imaging_flat, shift_img_flat_x, shift_img_flat_y) )
          Imaging_flat    = shift(Imaging_flat, shift_img_flat_x, shift_img_flat_y)
          slit_indices    = [slit_indices[0,*]+shift_img_flat_x, slit_indices[1,*]+shift_img_flat_y]
          arr             = intarr(si[0:1])
          arr[slit_indices[0,*], slit_indices[1,*]] = 1
          slit_indices_1D = where(arr eq 1, /NULL)
        endif

        A         = total(spectra_cube, 3, /nan)
        B         = spectra_flat
        y_crop    = minmax(where(total(A,1) ge 0.1*max(total(A,1)),/Null)); Good S/N range 
        A         = total(A[*, y_crop[0]:y_crop[1]], 1) 
        B         = total(B[*, y_crop[0]:y_crop[1]], 1) 
        A         = A / FFT( FFT(A, -1) * BUTTERWORTH(N_elements(A), cutoff = 5), 1 ) ;use the high frequency structure to omit the planet
        B         = B / max(B)
        A         = A[5:-5]
        B         = B[5:-5]
        search_shifts = indgen(100) - 50
        junk          = max(C_CORRELATE(B, A, search_shifts), location) ;find any shift between the spatial dimension of the flat and the data
        if (search_shifts[location] ne 0) then begin
          print, 'Shifting spectral channel''s flat field by:'+string(search_shifts[location])+' pixels along the slit' 
          window, 2, xs = ss[0], ys = ss[1], title = 'Pre-spectral flat shifting' 
          tv, bytscl(total(spectra_cube, 3, /nan) / shift(spectra_flat, 0, 0) )
          window, 3, xs = ss[0], ys = ss[1], title = 'Post spectral flat shifting'
          tv, bytscl(total(spectra_cube, 3, /nan) / shift(spectra_flat, 0, search_shifts[location]) )
          Spectra_flat = shift(spectra_flat, 0, search_shifts[location])
        endif

  imaging_cube = imaging_cube / rebin(Imaging_flat, si)                 ; FLATTEN THE IMAGING DATA
  spectra_cube = spectra_cube / rebin(Spectra_flat, ss)                 ; FLATTEN THE Spectral DATA
  ; Normalize for integration time, critical if integration times changed between the fits files we're co-adding
    for ifile = 0, nfiles-1 do imaging_cube[*,*,s[2]*ifile:s[2]*ifile+n_frames_per_file-1] = imaging_cube[*,*,s[2]*ifile:s[2]*ifile+n_frames_per_file-1] / integration[ifile]   
    for ifile = 0, nfiles-1 do spectra_cube[*,*,s[2]*ifile:s[2]*ifile+n_frames_per_file-1] = spectra_cube[*,*,s[2]*ifile:s[2]*ifile+n_frames_per_file-1] / integration[ifile]         
  ; Make a slit-less imaging cube, interpolating over pixels near the slit location
    dummy                 = imaging_cube
    s                     = size(dummy, /dimensions)
    for i = 0, s[2]-1 do begin
      slice = reform(dummy[*,*,i])
      x = slice
      slice[slit_indices_1D] = !Values.F_NaN
      sm = smooth( slice, slit_interp_width, /EDGE_TRUNCATE, /NAN )
      slice[slit_indices_1D] = sm[slit_indices_1D]
      dummy[*,*,i] = slice
      ;cgimage, slice                                                    ; inspect the flat/dark/slit corrected imaging channel data
      ;wait, 0.05                                                        ; Set to zero if we're done inspecting 
    endfor
    slitless_imaging_cube = temporary(dummy)
  
  ; - Perkins and similar telescopes have the imaging channel in the plane of the sky
  ; - AEOS has it's imaging channel is flipped -X = X (East-West only) See the Saturn image #225 on June 25 2018: 
  ;   South is down and to the right, and Doppler shift indicates E is up and to the right, 
  ;   There's no way to rotate this into plane of sky coordinates (N-Up, West-Right) without a transpose
  case telescope of 
    'Perkins' : 
    'AEOS': begin
       imaging_cube          = reverse(imaging_cube, 1)          ; AEOS is flipped E-W (X = -X)
       slitless_imaging_cube = reverse(slitless_imaging_cube, 1) ; AEOS is flipped E-W (X = -X)
    end  
  endcase
  spectra_cube[where(spectra_cube gt mean(spectra_cube)+40.*stddev(spectra_cube), spec_rm, /Null)] = !values.F_NaN ; flag hot pixels
  imaging_cube[where(imaging_cube gt mean(imaging_cube)+30.*stddev(imaging_cube), imag_rm, /Null)] = !values.F_NaN ; flag hot pixels
  slitless_imaging_cube[where(spectra_cube gt mean(slitless_imaging_cube)+30.*stddev(slitless_imaging_cube), /Null)] = !values.F_NaN ; flag hot pixels
  print, 'Rejected ', spec_rm, ' spectral channel pixels and', imag_rm, ' imaging channel pixels'
  MWRFITS, imaging_cube, outdir + 'imaging_cube.fits', header, /create, /silent
  MWRFITS, spectra_cube, outdir + 'spectra_cube.fits', header, /create, /silent
  MWRFITS, slitless_imaging_cube, outdir + 'Slitless_imaging_cube.fits', header, /create, /silent

  ; We'll want a "sky" cube as well
    Sky_cube          = MRDFITS(sky_dir + sky_file, 0, header, /unsigned, /silent )
    imaging_sky_cube  = float(reform(sky_cube[ims[0]:ims[1],ims[2]:ims[3], *]))
    spectra_sky_cube  = float(reform(sky_cube[sps[0]:sps[1],sps[2]:sps[3], *]))
    si        = size(imaging_sky_cube, /dimensions)                           ; size of extracted IMAGING portion of frame
    ss        = size(spectra_sky_cube, /dimensions)                           ; size of extracted SPECTRAL portion of frame
    imaging_sky_cube = imaging_sky_cube - rebin(Imaging_dark, si)         ; Dark (Bias) subtract the data
    spectra_sky_cube = spectra_sky_cube - rebin(Spectra_dark, ss)         ; Dark (Bias) subtract the data
    imaging_sky_cube = imaging_sky_cube / rebin(Imaging_flat, si)         ; FLATTEN THE IMAGING DATA
    spectra_sky_cube = spectra_sky_cube / rebin(Spectra_flat, ss)         ; FLATTEN THE Spectral DATA
    if N_elements(si eq 3) then imaging_sky_cube = mean(imaging_sky_cube, dimension=3)
    if N_elements(ss eq 3) then Spectra_sky_cube = mean(Spectra_sky_cube, dimension=3)
    MWRFITS, imaging_sky_cube, outdir + 'imaging_sky_cube.fits', /create, /silent
    MWRFITS, spectra_sky_cube, outdir + 'spectra_sky_cube.fits', /create, /silent
    raw_data_cube = temporary(cube)
    MWRFITS, raw_data_cube, outdir + 'raw_data_cube.fits', /create, /silent
  beep
endif 

; =====================================================================================================================
; Part 1: Generate a Synthetic Image of Mercury using Hapke's fomulation. Match it to the platescale of each channel
; =====================================================================================================================
if Part eq 1 or part eq 99 then begin
  ; On Perkins, the imaging channel gives you a plane of sky orientation, which is then rotated at the rotiserizer's rotation angle.
  ; On AEOS, the imaging channel gives has a transpose of XXX relative to the plane of sky??? http://www.leibniz-kis.de/fileadmin/user_upload/observatorien/gre_docs/man/GRE-KIS-MAN-0008_v001_derotator.pdf
    imaging_header = headfits(outdir + 'imaging_cube.fits')               ; Get the timestamps
    spectra_header = headfits(outdir + 'spectra_cube.fits')

  ; Generate a Hapke Model of Mercury at the observation time/geometry. Match its plate scale to the RIPS spectral channel
    Surface_Flux_Calibration, Body = Body, Timestamp = sxpar(imaging_header, 'DATE'), Wavelength = Flux_Cal_Wavelength, Output_image = outdir + 'Hapke'+'_'+telescope+'_'+night+'.eps', /align_celestial_north, $
      Hapke_Platescale = Hapke_Platescale, MR_per_A = MR_per_A, Make_picture = Make_picture
    cgloadct, ct                                                          ; get the color table back, since the above subroutine uses a different one 

    Flux_Cal_Img_Size = size(MR_per_A, /dimensions)
    Flux_Cal_Spectra  = congrid(MR_per_A, round(Flux_Cal_Img_Size[0]*Hapke_Platescale/Spectral_platescale), round(Flux_Cal_Img_Size[1]*Hapke_Platescale/Spectral_platescale), cubic = -0.5)
    Spectra_Hapke     = rotate(Flux_Cal_Spectra, 7)                       ; Always a mirror image of the imaging channel, X=X, Y=-Y
    Flux_Cal_Imaging  = congrid(MR_per_A, round(Flux_Cal_Img_Size[0]*Hapke_Platescale/Imaging_platescale), round(Flux_Cal_Img_Size[1]*Hapke_Platescale/Imaging_platescale), cubic = -0.5)
    Imaging_hapke     = Flux_Cal_Imaging                                  ; No transpose needed, this in in the plane of the sky already
    Picture_Size      = size(Make_picture, /dimensions)
    Picture_Hapke     = congrid(Make_picture, round(Picture_Size[0]*Hapke_Platescale/Spectral_platescale), round(Picture_Size[1]*Hapke_Platescale/Spectral_platescale), cubic = -0.5)

  ; Now pad them both with zeros. This makes space since we'll be blurring the hell out off them, match the size of each channel
    IH                = size(Imaging_Hapke, /dimensions)
    SH                = size(Spectra_Hapke, /dimensions)
    SM                = size(Make_Picture,  /dimensions)
    SP                = size(Picture_Hapke, /dimensions)
    SI                = [sxpar(imaging_header, 'NAXIS1'), sxpar(imaging_header, 'NAXIS2')]
    SS                = [sxpar(Spectra_header, 'NAXIS1'), sxpar(Spectra_header, 'NAXIS2')]
    SF                = [Picture_Size[0]*2*subframe/sp[0], Picture_Size[1]*2*subframe/sp[1]]   
    padded_Imaging_Hapke = REPLICATE(0., SI[0], SI[1])
    padded_Spectra_Hapke = REPLICATE(0., SS[0], SS[1])
    padded_Picture       = REPLICATE(0., SS[0], SS[1])
    padded_Hi_Res        = REPLICATE(0., SF[0], SF[1])
    padded_Imaging_Hapke[ (SI[0]- IH[0])/2., (SI[1]- IH[1])/2. ] = Imaging_Hapke
    padded_Spectra_Hapke[ (SS[0]- SH[0])/2., (SS[1]- SH[1])/2. ] = Spectra_Hapke
    padded_Picture[ (SS[0]- SP[0])/2., (SS[1]- SP[1])/2. ]       = Picture_Hapke
    padded_Hi_Res[ (SF[0]- SM[0])/2., (SF[1]- SM[1])/2. ]        = Make_Picture      ; should be the same as Picture_Hapke but in high resolution (bigger size array)
    Spectra_Hapke = padded_Spectra_Hapke
    Imaging_Hapke = padded_Imaging_Hapke
    Picture_Hapke = padded_Picture
    Hi_Res_Hapke  = padded_Hi_Res
  save, Spectra_Hapke, Imaging_Hapke, Picture_Hapke, Hi_Res_Hapke, filename = outdir + 'Hapke_Models.sav'
endif

; =====================================================================================================================
; Part 2: Find the imaging channel's rotation angle from a "plane-of-sky" "north-is-up" orientation
; =====================================================================================================================
if Part eq 2 or part eq 99 then begin
  if part ne 99 then restore, outdir + 'Hapke_Models.sav'
  slitless_imaging_cube = MRDFITS(outdir + 'slitless_imaging_cube.fits', 0, header, /fscale, /silent )
  s                 = size(slitless_imaging_cube, /dimensions)    ; size of extracted IMAGING portion of frame
    window, 0, xs = 2*subframe, ys  = 2*subframe, xpos = 0, ypos = 0
    window, 1, xs = 2*subframe, ys  = 2*subframe, xpos = 2*subframe+30, ypos = 0
    window, 2, xs = s[0], ys = s[1] 
    window, 3, xs = s[0], ys = s[1] 
     
  ; blur the Hapke reference by a fixed seeing estimate "effective_seeing" in arcsec FWHM
    seeing_sigma      = effective_seeing / (2.0*sqrt(2.0*alog(2)))  ; convert FWHM to sigma
    seeing_sigma      = seeing_sigma/Imaging_platescale             ; convert to Gaussian sigma in pixels
    ref_bandpass      = convol(imaging_hapke, GAUSSIAN_FUNCTION([seeing_sigma, seeing_sigma]), /EDGE_ZERO, /normalize, /NAN)

  ; crop them to speed things up  
    B                 = float(ref_bandpass[ s[0]/2. - subframe:s[0]/2. + subframe - 1, s[1]/2. - subframe:s[1]/2. + subframe - 1 ])

  rough_rotational_alignments = intarr(s[2])
  xy_alignments               = intarr(s[2], 2)
  Print, '---File #----Rotation angle to align w/ Celestial North---X and Y Shift WRT Hapke Reference-----------------------'
  test_align = fltarr(s)
  for i = 0, s[2]-1 do begin
    
    ; filter the data, crop both images
      frame        = reform(slitless_imaging_cube[*,*,i])           ; slit-less imaging frame (bias and flat corrected)
      img_bandpass = float(bandpass_filter(frame, 0., 0.15, /butterworth))
      junk         = max(total(ref_bandpass,2), ref_bandpass_centeroid_x) & junkx = max(total(img_bandpass,2), img_bandpass_centeroid_x)
      junk         = max(total(ref_bandpass,1), ref_bandpass_centeroid_y) & junky = max(total(img_bandpass,1), img_bandpass_centeroid_y)
      subframe_x   = ref_bandpass_centeroid_x - img_bandpass_centeroid_x
      subframe_y   = ref_bandpass_centeroid_y - img_bandpass_centeroid_y
      A            = shift(img_bandpass, subframe_x, subframe_y)
      A            = A[ s[0]/2. - subframe:s[0]/2. + subframe - 1, s[1]/2. - subframe:s[1]/2. + subframe - 1 ]

    ; Get the rough angle that rotationally aligns the imaging channel with the blurred Hapke image 
      log = fltarr(360, 3)
      for ang = 0, 359 do begin
        rot_ref = rot(B, ang, missing = 0.)
        zstack_align_images, A, rot_ref, x, y, corr_image, corr_dims
        log[ang, *] = [ corr_image[ corr_dims[0], corr_dims[1] ], x, y ]
      endfor
      junk = max(log[*,0], angle)
      xy_alignments[i,*] = [ subframe_x + log[angle,1], subframe_y + log[angle,2] ] ; x, y alignments of a raw Imaging Ch frame to the Hapke into the hapke frame, apply BEFORE rotational alignement
      if angle gt 180 then angle = angle - 360
      rough_rotational_alignments[i] = angle

    ; For comparison with a full frame, centered & blurred hapke 
      test_align[*,*,i] = rot(shift( frame, xy_alignments[i,0], xy_alignments[i,1] ), -rough_rotational_alignments[i] )
      print, i, '          ', angle, ' degrees             ', xy_alignments[i,0], xy_alignments[i,1]
            
      wset, 2
      tv, bytscl(reform(test_align[*,*,i]))
      wset, 3
      tv, bytscl(ref_bandpass)
  endfor

  ; now the angles just calculated are a bit rough and shouldn't change over time, 
    nfiles = N_elements(Mercury_files)
    filenames = file_Search(Mercury_Dir + Mercury_files, count = n_flats)
    cumulative = 0L
    rotation_each_file = fltarr(nfiles)
    for i = 0, N_elements(Mercury_files)-1 do begin
      header            = headfits(filenames[i])
      n_frames_per_file = sxpar(header, 'NUMKIN')
      RESISTANT_Mean, rough_rotational_alignments[cumulative:cumulative + n_frames_per_file-1], 2, Mean_rotation_this_file, Sigma_Mean, Num_Rejected
      rotation_each_file[i] = Mean_rotation_this_file
      cumulative = cumulative + n_frames_per_file
    endfor
  
  Case telescope of
    'Perkins': rotational_alignments = congrid(rotation_each_file, nfiles*n_frames_per_file, /center)
    'AEOS': rotational_alignments = congrid(rotation_each_file, nfiles*n_frames_per_file, /center)     ;Hack, should include az-alt rotation on the sky, but doubt it makes any difference
  endcase
  cgplot, rough_rotational_alignments, psym = 4, /ynozero, title = 'Rotation angle on the Plane-of-Sky reference to match the imaging channel', xtitle = 'frame number #', ytitle = 'degrees
  cgplot, rotational_alignments, psym = 3, color = 'red', /overplot

  ; We'll use these to rotationally align data to a plane-of-sky refrence, but we just found rotations the other way (plane-of-sky reference to data) --> invert the angles
    rotational_alignments = -rotational_alignments
    rough_rotational_alignments = -rough_rotational_alignments ;unsure which is better!
  
  save, rotational_alignments, rough_rotational_alignments, xy_alignments, filename = outdir + 'rotational_alignments.sav' 
endif

; =====================================================================================================================
; Part 3 : find Mercury centroids and seeing by cross-correlation with a blurred Hapke image
; =====================================================================================================================
if part eq 3 or part eq 99 then begin ;Find the centroids by cross-correlation with the last image
  if part ne 99 then restore, outdir + 'Hapke_Models.sav'         ; don't need to read in if we're continuing
  if part ne 99 then restore, outdir + 'rotational_alignments.sav'
  imaging_cube          = MRDFITS(outdir + 'imaging_cube.fits', 0, header, /fscale, /silent ) 
  slitless_imaging_cube = MRDFITS(outdir + 'slitless_imaging_cube.fits', 0, header, /fscale, /silent )  
  spectra_cube          = MRDFITS(outdir + 'spectra_cube.fits', 0, header, /fscale, /silent )
  s         = size(imaging_cube, /dimensions)                     ; size of extracted IMAGING portion of frame
  ss        = size(spectra_cube, /dimensions)                     ; size of extracted SPECTRAL portion of frame

  ; Define alignement/seeing reference and variables based on frame sizes
    reference           = Imaging_Hapke                           ; Use the idealize Mercury disk as a reference
    Imaging_shift_array = intarr(s[2],2)                          ; array of x, y shift values for aligning images
    sharpness_metric    = fltarr(s[2])                            ; standard deviations of imaging frames (TEMPORARILY used to define image quality)
    aligned_imaging_cube= fltarr(s)                               ; Imaging cube after everything is aligned

  ; Set up windows based on image sizes
    window, 1, xs = 100, ys = 100, xpos=0,     ypos=0,                          title='correlation inspection'
    window, 2, xs = s[0], ys = s[1], xpos=winpos_x+s[0]+20,     ypos=0,         title='Raw Reference Frame for spatial alignment'
    window, 3, xs = s[0], ys = s[1], xpos=winpos_x+s[0]+20,     ypos=s[1]+40,   title='FRAME BY FRAME Comparison'
    window, 4, xs = s[0], ys = s[1], xpos=winpos_x+2*(s[0]+20), ypos=0,         title='Smoothed Raw Frame & Data Alignment'
    window, 5, xs = s[0], ys = s[1], xpos=winpos_x+2*(s[0]+20), ypos=s[1]+40,   title='Aligned Filtered/Smoothed Comparison'

  ; blur the Hapke reference by a fixed seeing estimate "effective_seeing" in arcsec FWHM
    seeing_sigma      = effective_seeing / (2.0*sqrt(2.0*alog(2)))  ; convert FWHM to sigma
    seeing_sigma      = seeing_sigma/Imaging_platescale             ; convert to Gaussian sigma in pixels
    ref_blurred       = convol(reference, GAUSSIAN_FUNCTION([seeing_sigma, seeing_sigma]), /EDGE_ZERO, /normalize, /NAN)

  ; Log the correlations to the reference image (correlation post-alignment)
    Correl     = fltarr(s[2])
    print, 'Calculating Frame-By-Frame Spatial Shifts for Co-alignment...'
    print, '-----Frame-#---X-------Y----Seeing FWHM (")-----Rotation to Plane-Of-Sky-----

  for i = 0, s[2]-1 do begin
    frame = reform(slitless_imaging_cube[*,*,i])

    ; filter the frames and blur the Hapke reference by a fixed seeing estimate "effective_seeing" in arcsec FWHM
      img_bandpass      = float(bandpass_filter(frame, 0., 0.15, /butterworth))
    
    ; Rotate & roughly align the plane-of-sky reference to match the imaging channel frame
      rot_ref = shift( rot( reference, -rotational_alignments[i] ), -xy_alignments[i,0], -xy_alignments[i,1] )

    ; Now estimate the seeing and find the optimal alignment
      case telescope of
        'Perkins': Crop2frame = 100
        'AEOS':    Crop2frame = 160
      endcase

      A = img_bandpass
      B = rot_ref
      A = A[ s[0]/2. - Crop2frame:s[0]/2. + Crop2frame - 1, s[1]/2. - Crop2frame:s[1]/2. + Crop2frame - 1 ]
      B = B[ s[0]/2. - Crop2frame:s[0]/2. + Crop2frame - 1, s[1]/2. - Crop2frame:s[1]/2. + Crop2frame - 1 ]
      A = rebin(A, 2.*Crop2frame / Seeing_Est_Rebin_Factor, 2.*Crop2frame / Seeing_Est_Rebin_Factor)
      B = rebin(B, 2.*Crop2frame / Seeing_Est_Rebin_Factor, 2.*Crop2frame / Seeing_Est_Rebin_Factor)

      y_fine = 0. & x_fine = 0.                                       ; A default in case no shift fine tuning is returned in Amoeba
      wset, 1
      sharpness_metric[i] = AMOEBAX(1.e-3, 1.e-3, function_name='Estimate_Seeing', FUNCTION_VALUE = fval, P0 = [effective_seeing/Seeing_Est_Rebin_Factor], Scale = [0.1], NCalls = 25 )
      if y_fine + x_fine gt 0. then print, 'Fine tuned Alignment'
      Imaging_shift_array[i,*]    = fix(round([xy_alignments[i,0] + x_fine*Seeing_Est_Rebin_Factor, xy_alignments[i,1] + y_fine*Seeing_Est_Rebin_Factor]))             ; add fine-tuned shift, if one is found
      Aligned_Frame               = rot( shift(frame, [Imaging_shift_array[i,*]]), rotational_alignments[i] )
      Aligned_Frame_Bandpass      = float(bandpass_filter(Aligned_Frame, 0., 0.15, /butterworth))
      aligned_imaging_cube[*,*,i] = aligned_Frame                                                                                         ; Co-align imaging frame in X and Y

    ; Plot & Inspect progress to track the quality of the alignments
      print, i, Imaging_shift_array[i,0], Imaging_shift_array[i,1], sharpness_metric[i]*Seeing_Est_Rebin_Factor, rotational_alignments[i] ; update progress
      wset, 2
      cgimage, reference, /keep_aspect
      wset, 4
      cgimage, Aligned_Frame, /keep_aspect
      cgcontour, ref_blurred/max(ref_blurred), levels=[0.2,0.4,0.6,0.8], position = [0,0,1,1], /noerase
      wset, 3
      cgcontour, ref_blurred/max(ref_blurred), levels=[0.2,0.4,0.6,0.8], position = [0,0,1,1], /fill
      cgcontour, ref_blurred/max(ref_blurred), levels=[0.2,0.4,0.6,0.8], position = [0,0,1,1], /noerase
      cgcontour, Aligned_Frame_Bandpass/max(Aligned_Frame_Bandpass), levels=[0.2,0.4,0.6,0.8], c_color='red', position = [0,0,1,1], /noerase
      wset, 5
      cgimage, Aligned_Frame, /keep_aspect
  endfor  ; end loop to determine Imaging_shift_array values to align images
  sharpness_metric[where(sharpness_metric lt 0., /NULL)] = max(sharpness_metric)        ; Hack! Hack! Hack! Hack! HACK! This shouldn't be needed, but it is!!!! 
  
  sharpness_metric         = sharpness_metric*Seeing_Est_Rebin_Factor                   ; Seeing was determined using rebinned arrays, this is the true seeing if data unbinned
  Spectra_shift_array      = fix(round(float(Imaging_shift_array) / plate_scale_ratio)) ; Spectral and imaging channels have different platescales, hence different shifts
  Spectra_shift_array[*,1] = -Spectra_shift_array[*,1]                                  ; X_i = X_s, but +Y_i in the imaging channel is -Y_s in the spectral channel
  aligned_spectra_cube     = fltarr(ss)                                                 ; Spectral cube after alignment in Y only (for now) (note that

  ; Now loop over all frames to align the x, y pixel shift. Ignore the rotation, and shift the s
  for i = 0, s[2]-1 do begin
    frame     =  reform(slitless_imaging_cube[*,*,i])                                   ; "raw" imaging frame
    specframe =  reform(spectra_cube[*,*,i])                                            ; "raw" spectral frame
    aligned_imaging_cube[*,*,i] = shift(frame, [Imaging_shift_array[i,*]])              ; Co-align imaging frame in X and Y
    aligned_spectra_cube[*,*,i] = shift(specframe, [0, Spectra_shift_array[i,1]])       ; Co-align Spectral frame **in Y only**. Y_Spec = -Y_img/PSR. Later on, we do the X_Spec = X_img/PSR.
  endfor

  ; inspect best / worse seeing and the lucky cut
  s = size(aligned_imaging_cube, /dimensions)
  cut_index = fix(s[2]*Lucky_fraction)
  ranked = SORT(sharpness_metric)
  window, 0, xs = s[0], ys = s[1]
  cgimage, total(aligned_imaging_cube[*,*,ranked[s[2]-100:s[2]-1]], 3), /keep_aspect  ; Worst seeing 100 frames
  window, 1, xs = s[0], ys = s[1]
  cgimage, total(aligned_imaging_cube[*,*,ranked[0:99]], 3), /keep_aspect             ; Best seeing 100 frames
  window, 2, xs = s[0], ys = s[1]
  cgimage, total(aligned_imaging_cube[*,*,ranked[0:cut_index]], 3), /keep_aspect      ; Frames flagged as the "lucky ones"
  window, 3

  save, Imaging_shift_array, Spectra_shift_array, correl, aligned_imaging_cube, aligned_spectra_cube, sharpness_metric, filename = outdir + 'shift_array.sav'
  beep
endif

;; test to see if we can smooth out some noise in the best shifts
;if part eq 3.5 then begin 
;  restore, outdir + 'shift_array.sav'
;  ;test_Imaging_shift_array = smooth(Imaging_shift_array, [5,1])
;  test_Imaging_shift_array = Imaging_shift_array
;  test_Imaging_shift_array[*,0] = medsmooth(test_Imaging_shift_array[*,0],7)
;  test_Imaging_shift_array[*,1] = medsmooth(test_Imaging_shift_array[*,1],7)
;  
;  slitless_imaging_cube = MRDFITS(outdir + 'slitless_imaging_cube.fits', 0, header, /fscale, /silent )
;  s         = size(slitless_imaging_cube, /dimensions)                     ; size of extracted IMAGING portion of frame  
;  aligned_imaging_cube     = fltarr(s)                                                 ; Spectral cube after alignment in Y only (for now) (note that
;  test_aligned_imaging_cube     = fltarr(s) 
;    ; Now loop over all frames to align the x, y pixel shift. Ignore the rotation, and shift the s
;  for i = 0, s[2]-1 do begin
;    frame     =  reform(slitless_imaging_cube[*,*,i])                                   ; "raw" imaging frame
;    aligned_imaging_cube[*,*,i] = shift(frame, [Imaging_shift_array[i,*]])              ; Co-align imaging frame in X and Y
;    test_aligned_imaging_cube[*,*,i] = shift(frame, [test_Imaging_shift_array[i,*]])              ; Co-align imaging frame in X and Y
;  endfor
;  cut_index = fix(s[2]*Lucky_fraction)
;  ranked = SORT(sharpness_metric)
;  window, 0, xs = 1500, ys=1000
;  cgimage, total(aligned_imaging_cube, 3), /keep_aspect, /axes, xr = [100,300], yr = [100, 300]
;  cgcontour, total(aligned_imaging_cube[*,*,ranked[0:99]], 3), /onimage, nlevels =4
;  cgcontour, total(test_aligned_imaging_cube[*,*,ranked[0:99]], 3), color ='red', /onimage, nlevels =4
;endif


; =====================================================================================================================
; Part 4 : Isolate the sodium emission in every spectral frame
; =====================================================================================================================
if part eq 4 or part eq 99 then begin
  if part ne 99 then restore, outdir + 'shift_array.sav'  ; contains Imaging_shift_array, Spectra_shift_array, Correl, sharpness_metric, aligned_imaging_cube and aligned_spectra_cube
  img_cube = MRDFITS(outdir + 'imaging_cube.fits', 0, header, /fscale, /silent )
  s         = size(aligned_spectra_cube, /dimensions)
  si        = size(img_cube, /dimensions)
  loadct, 0
  window, 0, xpos=winpos_x,         ypos=winpos_y,           xs=s[0],  ys=s[1],  title='IDL 0 - SPECTRAL FRAME'
  window, 1, xpos=winpos_x,         ypos=winpos_y+s[1]+20,   xs=s[0],  ys=s[1],  title='IDL 1 - SCALED SODIUM FREE REFERENCE FOR REMOVING REFLECTANCE (Sky, Moon etc...)'
  window, 2, xpos=winpos_x+0.9*s[0],ypos=winpos_y,           xs=s[0],  ys=s[1],  title='IDL 2 - EXOSPHERE'
  window, 3, xpos=winpos_x+1.8*s[0],ypos=winpos_y        ,   xs=si[0], ys=si[1], title='IDL 3 - IMAGING CHANNEL FRAME'
  window, 4, xpos=winpos_x+0.9*s[0],ypos=winpos_y+s[1]+20,   xs=1200,  ys=s[1],  title='IDL 4 - INSPECTION'

  ; Get the curvature of the line trace from an arc lamp frame
    search      = 10
    Neon        = MRDFITS(outdir + 'NEON_ARC_LAMP.fits', 0, arc_header, /fscale, /silent )
    spec        = total(Neon, 2)                    ; Always inspect that this portion of the line is straight
    junk        = max(spec, line_center)
    line_centers = fltarr(s[1])
    for i = 0, s[1]-1 do begin
      result = mpfitpeak(findgen(search*2. + 1), Neon[line_center-search:line_center+search, i], a, STATUS = STATUS, /positive)
      if STATUS ge 1 then line_centers[i] = A[1] else line_centers[i] = !values.F_nan
    endfor
    keep        = where(finite(line_centers), count, /Null)
    COEFF       = ROBUST_POLY_FIT(keep, line_centers[keep], 2)
    Neon_trace  = poly(findgen(s[1]), coeff) + line_center-search
    ;plot, findgen(s[1]), line_centers, psym = 3
    ;oplot, findgen(s[1]), Neon_trace - (line_center-search)

  ; Align the reference spectra that is being used fit Mercury's surface reflectance
    mean_spectral_ch = total(aligned_spectra_cube, 3, /NaN)
    y_Mercury        = where( total(mean_spectral_ch,1,/nan) ge 0.9*max(total(mean_spectral_ch,1,/nan)), /Null )   ; Good S/N range of Mercury over a small slit region w/ neglegible curvature
    Collapsed_spec   = total(mean_spectral_ch[*,y_mercury], 2, /nan)

  ; Read in a solar spectrum to use for reflectance substraction
    SOLAR_SPECTRUM_FILE = 'C:\IDL\Io\Kurucz_2005_irradthuwl.dat'  ; Solar Spectrum at 1AU in W/m2/nm
    READCOL, SOLAR_SPECTRUM_FILE, F='A,A', WL_nm, flux, STRINGSKIP = '#', /Silent ;flux is in W/m2/nm
    start = where(WL_nm eq '588.200')
    stop  = where(WL_nm eq '590.600')
    WL_nm = float(WL_nm[start:stop])
    flux  = float(flux[start:stop])
    ; change flux units from W/m^2/nm to photons / (cm^2 s A)
    ; multiply by ((lambda / hc) / 1 W)  * (1 m^2 / 1e4 cm^2) * (1 nm / 10 A)
    conversion = ((WL_nm*1.e-9)/(6.62606957e-34*299792458.D)) * (1./1.e4) * (1./10.)
    flux = flux * conversion                                                                    ; photons / (cm^2 s A)
    WL_A = temporary(WL_nm) * 10.                                                               ; Wavelength from nm into angstroms

  ; find the ***rough*** location of the Na D1 & D2 absorption wells in the collapsed datacube and stretch the solar reference to the same dispersion
    D2_abs_pixel     = (where(Collapsed_spec eq min(Collapsed_spec[0:400])))[0]
    D1_abs_pixel     = (where(Collapsed_spec eq min(Collapsed_spec[D2_abs_pixel+200:*])))[0]
    junk             = min(ABS(WL_A - Na_D2_WL), Ref_D2_abs_pixel)
    junk             = min(ABS(WL_A - Na_D1_WL), Ref_D1_abs_pixel)
    dispersion_ratio = float((D1_abs_pixel - D2_abs_pixel)) / float((Ref_D1_abs_pixel - Ref_D2_abs_pixel))
    ref = congrid(flux, dispersion_ratio*n_elements(WL_A))                                    ; stretch the solar reference, since it must have the same dispersion as the data

    ref = smooth(ref, 5)                                                                      ; Hack Hack Hack, this depends on slitwidth, but I don't know how to automate this!!!
    Collapsed_ref = ref[sps[0]:sps[1]]                                                        ; crop the solar reference, since it must be the same size the data for cross-correlation
    
;          p0 = [2.2e-7, 0., -40., 2., -6444.]   ; Initial guess (multiply, add)
;          p = mpfitfun('Match_Reference_Spectrum', float(Collapsed_ref), Collapsed_spec, replicate(stddev(Collapsed_spec[10:200])*.2, N_elements(Collapsed_spec)), p0, /NaN, status=status)
;          shifted = interpolate(Collapsed_ref, findgen(n_elements(Collapsed_ref)) - P[2])      ; sub-pixel shifting
;          test = P[0]*gauss_smooth(shifted, p[3], /EDGE_TRUNCATE) + P[1] + P[4]*findgen(n_elements(shifted))
;          cgplot, collapsed_spec
;          cgplot, test, /overplot, color = 'red'
;          stop

  search_shifts = indgen(512) - 256                                                             ; Array to search for shifting betweeen reference solar spectrum and the Mercury reflectance
  junk          = max(C_CORRELATE(Collapsed_ref, Collapsed_spec, search_shifts), location)      ; Find any shift between the reference solar spectrum and the Mercury reflectance
  ref_spectrum  = shift(Collapsed_ref, search_shifts[location-1])                               ; Align the reflectance spectrum with the data (assumes dispersion is unchanged)

  ; Re-find the ***rough*** location of the Na D1 & D2 absorption wells in the reflectance reference, now that it's been aligned
    D2_abs_pixel  = (where(ref_spectrum eq min(ref_spectrum[0:400])))[0]
    D1_abs_pixel  = (where(ref_spectrum eq min(ref_spectrum[D2_abs_pixel+200:*])))[0]

  ; Fine tune these locations---sodium emission can fill in the solar fraunhofer wells, so use the Nickel line about half way between them as a wavelength reference
    Ni_abs_pixel  = round(interpol([D2_abs_pixel,D1_abs_pixel], [Na_D2_WL, Na_D1_WL], 5892.88))
    search_shifts = indgen(11) - 5
    junk          = max(C_CORRELATE(ref_spectrum[Ni_abs_pixel-10:Ni_abs_pixel+10], Collapsed_spec[Ni_abs_pixel-10:Ni_abs_pixel+10], search_shifts), location)
    ref_spectrum  = shift(ref_spectrum, search_shifts[location])                                  ; Fine-tune alignment in the region of the Nickel line
    D2_abs_pixel  = (where(ref_spectrum eq min(ref_spectrum[0:400])))[0]
    D1_abs_pixel  = (where(ref_spectrum eq min(ref_spectrum[D2_abs_pixel+200:*])))[0]
    dispersion    = (Na_D1_WL - Na_D2_WL) / (D1_abs_pixel - D2_abs_pixel)                         ; Angstroms / pixel
    dispersion_v  = cspice_clight()*dispersion / mean([Na_D2_WL, Na_D1_WL])                       ; Km/s per pixel at Na
    reference     = fltarr(s[0:1])                                                                ; make the 2D reference
    x = lindgen(s[0])
    for i = 0, s[1]-1 do begin
      xshift = Neon_trace[i] - Neon_trace[mean(y_mercury)]
      reference[*,i] = interpolate(ref_spectrum, x-xshift)                                        ; match the spectral line's intrinsic curvature
    endfor
  ;loadct, 0
  ;window, 0
  ;tv, bytscl(spectra_sky_cube, 0, 5000)
  ;window, 1
  ;tv, bytscl(reference_test)

  ; Calculate the Doppler shift between absorption and emission line centers, i.e. Mercury's heliocentric velocity
    UTC_String = sxpar(header, 'DATE')
    cspice_UTC2ET, UTC_string, ET
    cspice_spkezr, Body, ET, 'J2000', 'NONE', 'Sun', Mercury_Sun_State, ltime
    theta  = cspice_vsep(Mercury_Sun_State[0:2], Mercury_Sun_State[3:5])
    Dopplershift = cos(theta) * norm(Mercury_Sun_State[3:5])                                      ; km/s heliocentric velocity: scalar projection of the relative velocity along the line of sight
    D1_exo_pixel = D1_abs_pixel+round(-Dopplershift / (dispersion*cspice_clight() / Na_D1_WL))
    D2_exo_pixel = D2_abs_pixel+round(-Dopplershift / (dispersion*cspice_clight() / Na_D2_WL))
    Print, Body+'''s heliocentric velocity =', Dopplershift, ' km/s'

  ; Identify the index ranges of pixels containing both D2 & D1 emissions
    Na_D_rngs    = fix(round([D2_exo_pixel+[-10.,10.], D1_exo_pixel+[-10.,10.]]))

  ; Roughly match the dynamic range of the reference to that of the collapsed spectrum
    ind = D2_abs_pixel-100 + indgen(D1_abs_pixel+100 - (D2_abs_pixel-100))                        ; indices to match the dynamic range over
    ind = cgSetDifference(ind, [na_D_rngs[0] + indgen(na_D_rngs[1]-na_D_rngs[0]), na_D_rngs[2] + indgen(na_D_rngs[3]-na_D_rngs[2])] ) ;exclude sodium exosphere
    dr_s         = minmax(Collapsed_spec[ind], /NaN)
    dr_r         = minmax(ref_spectrum[ind], /NaN)
    multiplier   = (dr_s[1]-dr_s[0]) / (dr_r[1]-dr_r[0])
    ref_spectrum = ref_spectrum * multiplier
    ref_spectrum = ref_spectrum + min(Collapsed_spec) - min(ref_spectrum)
    residual     = Collapsed_spec - ref_spectrum

  ; Fit D1 & D2 lines locally and individually rather than the whole spectrum together
    fit_range = 30  ; local fit range in pixels  
    omit_range= 10  ; the pixel width accross the exosphere's line to omit  
    lines = ['D2', 'D1'] 

  ;-------------------------------------------------Do the fits--------------------------------------------------------------
    print, 'Applying the solar reflectance spectrum and isolating exosphere emission...'
    Surf_refl_spectra_cube = fltarr( size(aligned_spectra_cube, /dim))                          ; variable to hold the raw Mercury spectra with exosphere + surface reflectance
    exosphere_spectra_cube = fltarr( size(aligned_spectra_cube, /dim))                          ; variable to hold the residual Mercury spectra with exosphere only
    plot_scaling           = 20.*mean(aligned_spectra_cube)
    
    For h = 0, N_elements(lines)-1 do begin
      wset, 4
      line = lines[h]
      Case line of
        'D2': begin
          incl_indicies = D2_exo_pixel+indgen(fit_range+1)-fit_range/2
          omit_indicies = D2_exo_pixel+indgen(omit_range+1)-omit_range/2
          fit_ind       = cgSetDifference(incl_indicies, omit_indicies) ;exclude sodium exosphere
          weights       = 1. / abs(fit_ind - D2_exo_pixel)^2
          wset, 4
          cgplot, Collapsed_spec, xr = [D2_abs_pixel-2*fit_range, D2_abs_pixel+2*fit_range],/nodata                  ; setup inspection
        end
        'D1': begin
          incl_indicies = D1_exo_pixel+indgen(fit_range+1)-fit_range/2
          omit_indicies = D1_exo_pixel+indgen(omit_range+1)-omit_range/2
          fit_ind       = cgSetDifference(incl_indicies, omit_indicies) ;exclude sodium exosphere
          weights       = 1. / abs(fit_ind - D1_exo_pixel)^2
          wset, 4
          cgplot, Collapsed_spec, xr = [D1_abs_pixel-2*fit_range, D1_abs_pixel+2*fit_range],/nodata                  ; setup inspection
        end
      endcase

      ; Inspection: Are the solar reflectance spectra and Mercury reflectance spectra aligned? Are the absorption and emission wavelengths correct?
        cgColorfill, [incl_indicies[0], incl_indicies[-1], incl_indicies[-1], incl_indicies[0]], [!y.crange[0], !y.crange[0], !y.crange[1], !y.crange[1]], Color='grey', /line_fill
        cgColorfill, [omit_indicies[0], omit_indicies[-1], omit_indicies[-1], omit_indicies[0]], [!y.crange[0], !y.crange[0], !y.crange[1], !y.crange[1]], Color='orange', /line_fill
        cgplot, Collapsed_spec, /overplot                                         ; Our data, collapsed in time and spatial dimensions
        cgplot, [Ni_abs_pixel,Ni_abs_pixel],[0, 1.e10], color = 'grey', /overplot ; Should match the Solar Nickel line
        cgplot, [D1_abs_pixel,D1_abs_pixel],[0, 1.e10], color = 'Green', /overplot; Should match reflectance absorption
        cgplot, [D2_abs_pixel,D2_abs_pixel],[0, 1.e10], color = 'Green', /overplot; Should match reflectance absorption
        cgplot, [D1_exo_pixel,D1_exo_pixel],[0, 1.e10], color = 'blue', /overplot ; Should match exosphere emission
        cgplot, [D2_exo_pixel,D2_exo_pixel],[0, 1.e10], color = 'blue', /overplot ; Should match exosphere emission
        cgplot, residual, color = 'orange', /overplot
        
        y = Collapsed_spec
        err_y = replicate(stddev(y[10:200]), N_elements(y))                                         ; Assume uniform weighting Hack must be a better way to weight things nearest the exosphere line centers
        x = ref_spectrum
        p0 = [1., 0.]   ; Initial guess (multiply, add)
        p = mpfitfun('Fit_Spec', x[fit_ind], y[fit_ind], err_y[fit_ind], p0, /quiet, /NaN, status=status, weights = weights) ; Fit a y = mx + b function to the spectrum at every spatial bin... SLOW
        cgplot, P[0]*ref_spectrum+P[1], color = 'orange', /overplot

      ; Look okay? If so, now scale and subtract... 
        for i = 0, s[2]-1 do begin                                                                  ; LOOP over each frame in the series
          ; find best scaling for ref_spectrum and subtract it (based on solar reflected emission between D lines)
          frame            = reform(aligned_spectra_cube[*,*,i])
          y_Mercury        = where( total(frame,1) ge 0.5*max(total(frame,1)), complement = not_mercury)                     ; FWHM range of Mercury
          illum = frame / reference                                                                 ; scale the illumination against Mercury's disk
          for iD = 0, 2, 2 do illum[Na_D_rngs[iD+0]:Na_D_rngs[iD+1],*] = !values.F_NaN              ; carefully avoid sodium emissions when scaling to the disk
          illum_along_slit = MEDIAN(illum, dimension=1)
          illum_along_slit = MEDSMOOTH( illum_along_slit, 5 )                                       ; helps with hot pixels from bad flat-fielding where the slit has some dust
          scaled_reference = reference * rebin(transpose(illum_along_slit), s[0], s[1])
          SMOOTHED_FRAME   = SMOOTH( frame, [0,5], /NAN, /EDGE_TRUNCATE )                           ; smoothing 5 pixels along slit helps get a decent spectral fit to noisy data
          for j = 0, s[1]-1 do begin                                                                ; Loop over each row and fit the spectra to the reference
            y = SMOOTHED_FRAME[*,J]
            sy = replicate(stddev(y[10:200]), N_elements(y))                                        ; Assume uniform weighting Hack must be a better way to weight things nearest the exosphere line centers
            x = scaled_reference[*,j]
            x = x[fit_ind]                                                                          ; hack! fit indices here are set as global parmaters but they should be local because of line curvature
            y = y[fit_ind]
            sy = sy[fit_ind]
            p0 = [1., 0.]   ; Initial guess (multiply, add)
            p = mpfitfun('Fit_Spec', x, y, sy, p0, /quiet, /NaN, status=status)                     ; Fit a y = mx + b function to the spectrum at every spatial bin... SLOW
            scaled_reference[*,j] = scaled_reference[*,j]*p[0] + p[1]
          endfor ;j row by row within frame # i

        ; Write the results of 2D scaled reference/reflection subrtraction for each line onto one half of the spectrum, using nickel 5892.88 as the split point
          Case line of                                                                              
            'D2': begin
                  surf_refl_spectra_cube[ 0:Ni_abs_pixel, *, i] = scaled_reference[0:Ni_abs_pixel,*]                                   ; write each each theortical reflectance spectra
                  exosphere_spectra_cube[ 0:Ni_abs_pixel, *, i] = Frame[0:Ni_abs_pixel,*] - scaled_reference[0:Ni_abs_pixel,*]           ; write each residual exosphere only spectra
            end
            'D1': begin
              surf_refl_spectra_cube[ Ni_abs_pixel+1:*, *, i] = scaled_reference[Ni_abs_pixel+1:*,*]                                   ; write each each theortical reflectance spectra
              exosphere_spectra_cube[ Ni_abs_pixel+1:*, *, i] = Frame[Ni_abs_pixel+1:*,*] - scaled_reference[Ni_abs_pixel+1:*,*]         ; write each residual exosphere only spectra
            end
          endcase
            
        ; Inspection
          wset, 0
          cgimage, bytscl(frame, 0, plot_scaling), /axes                                            ; plot the original spectral frame
          for iNa = 0, 3 do cgoplot, [1,1]*Na_D_rngs[iNa], !Y.Crange, color = 'white'               ; show the Na exclusion zones
          wset, 1
          cgimage, bytscl(scaled_reference, 0, plot_scaling), /axes, title=line+' Scaled Reference'
          for iNa = 0, 3 do cgoplot, [1,1]*Na_D_rngs[iNa], !Y.Crange, color = 'white'
          wset, 2
          cgimage, bytscl(Frame - scaled_reference, 0, plot_scaling), /axes                         ; plot the Mercury Na emission after subtracting disk
          for iNa = 0, 3 do cgoplot, [1,1]*Na_D_rngs[iNa], !Y.Crange, color = 'white'
          wset, 3
          img_frame = shift(reform(img_cube[*,*,i]), [Imaging_shift_array[i,*]])                    ; show frame just for visual reference
          cgimage, bytscl( img_frame )
          cgtext, 0, 0, 'Frame #'+strcompress(string(i, format = '(I)')), color = 'white'
          wset, 4
          cgplot, total(frame[*,y_Mercury],2, /nan)
          cgplot, total(scaled_reference[*,y_Mercury],2, /nan), /overplot, color = 'red'
          cgtext, .15, .8, line+' fit, ('+strcompress(n_elements(y_Mercury), /remove_all) +' spatial rows summed)', /normal        
      endfor ;i (loop over frames)
      print, 'Finished all '+line+' fitting
    endfor ;h which emission line was done  
  save, exosphere_spectra_cube, surf_refl_spectra_cube, Na_D_rngs, Neon_trace, dispersion, dispersion_v, filename = outdir + 'spectral_cubes.sav'
  cgloadct, ct
  beep
endif

; =====================================================================================================================
; Part 5 : extract the 1D exosphere spectra.
; =====================================================================================================================
if part eq 5 or part eq 99 then begin 
  if part ne 99 then restore, outdir + 'spectral_cubes.sav'
  if part ne 99 then restore, outdir + 'shift_array.sav'
  s = size(exosphere_spectra_cube, /dimensions) 
  window, 0, xpos=winpos_x,         ypos=winpos_y,           xs=s[0],  ys=s[1],  title='CO-ADDED RESIDUAL SPECTRUM (INSPECTION & LINE CENTERING)'
  window, 1, xpos=winpos_x,         ypos=winpos_y+s[1]+40,   xs=s[2],  ys=s[1],  title='EXTRACTED D2 BRIGHTNESS ALONG THE SLIT (Y) AS A TIME SERIES (X)'
  window, 2, xpos=winpos_x+s[2]+20, ypos=winpos_y+s[1]+40,   xs=s[2],  ys=s[1],  title='EXTRACTED D1 BRIGHTNESS ALONG THE SLIT (Y) AS A TIME SERIES (X)'

  ; Run MPFIT twice to trace the exosphere lines, first to find the D2 line center, and then the D1 line center 
    img       = median(exosphere_spectra_cube, dimension=3)    ; Co-add all the spectra that have surface reflectance removed 
    dummy     = img                                            ; A duplicate dummy image for inspection purposes
    search    = 8                                              ; Pixels to search over when finding the emission line center

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
      mean_spectral_ch = total(aligned_spectra_cube, 3, /NaN)
      real      = where( total(mean_spectral_ch,1) ge 0.4*max(total(mean_spectral_ch,1)), complement = no_data, /null )   ; range along slit where there's reliable exosphere signal
      D2_height[no_data] = 1.
      D2_trace  = neon_trace + mean(D2_center[real] - neon_trace[real]) + D2_expected_pixel-search 

    ; D1 Tracing
      D1_expected_pixel = mean(Na_D_rngs[2:3])             ; Rough D1 pixel location
      for i = 0, n_elements(img[0,*]) - 1 do begin
        result = mpfitpeak(findgen(search*2. + 1), img[D1_expected_pixel-search:D1_expected_pixel+search, i], a, STATUS = STATUS, /positive) 
        if STATUS ge 1 then D1_height[i] = A[0] else D1_height[i] = !values.F_nan
        if STATUS ge 1 then D1_center[i] = A[1] else D1_center[i] = !values.F_nan
        if STATUS ge 1 then D1_width[i] = A[2] else D1_width[i] = !values.F_nan
      endfor
      mean_spectral_ch = total(aligned_spectra_cube, 3, /NaN)
      real      = where( total(mean_spectral_ch,1) ge 0.4*max(total(mean_spectral_ch,1)), complement = no_data, /null )   ; range along slit where there's reliable exosphere signal
      D1_height[no_data] = 1.
      D1_trace  = neon_trace + mean(D1_center[real] - neon_trace[real]) + D1_expected_pixel-search 
      
    ; Get the point/line spread function from an Arc lamp frame
      Neon            = MRDFITS(outdir + 'NEON_ARC_LAMP.fits', 0, arc_header, /fscale, /silent )
      spec            = total(Neon, 2)                          ; Collapse the spatial dimension
      junk            = max(spec, rough_line_center)            ; Find the rough line center                      
      search_width    = 20.                                     ; Enough width to encompass both line curvature and the wings of the LSF
      PSF_2D          = fltarr(2.*search_width, s[1])           ; The PSF/LSF at every row
      for i = 0, s[1]-1 do begin                                ; Fit the arc lamp's LSF row by row 
        result        = mpfitpeak(findgen(search_width*2. + 1), Neon[rough_line_center-search_width:rough_line_center+search_width, i], A, $
                                   /positive, /NAN, STATUS = STATUS, NTERMS = 3)
        PSF_2D[*,i]   = interpolate(Neon[*, i], (rough_line_center-search_width) + (A[1]-search_width) + findgen(2.*search_width))               
      endfor
      PSF_1D          = median(PSF_2D, dimension = 2)           ; Collapse into 1 rejecting any outliers
      PSF_1D          = PSF_1D - min(PSF_1D)                    ; give it a floor of exactly zero
      PSF_1D          = PSF_1D / total(PSF_1D)                  ; normalize it to a total of 1
      PSF_1D          = psf_1D[where(PSF_1D ge max(PSF_1D)/s[2])];make it as small as possible, since a wider PSF produces wider *edge effects* in the final frame where only these weights exist
      B               = mpfitpeak(findgen(17), psf_1D, coeffs, /positive, /NAN, STATUS = STATUS, NTERMS = 3)

  ; Inspection
    dummy[round(D2_trace), indgen(s[1])] = max(dummy) 
    dummy[round(D1_trace), indgen(s[1])] = max(dummy) 
    wset, 0
    tv, bytscl(dummy)

  ; Setup MPFIT Gaussian parameters
    width              = median([D2_width[real], D1_width[real]])
    parinfo            = replicate( {value:0., fixed: 0b, limited: [0b,0b], limits: fltarr(2) }, 3)
    parinfo[0].limited = 1b                                 ; limit line amplitude 
    ;parinfo[0].limits  = [0., 8.*max(d2_height)]            ; positive amplitudes only
    parinfo[0].limits  = [-8.*max(d2_height), 8.*max(d2_height)]  
    parinfo[1].limited = 1b                                 ; limit center
    parinfo[1].limits  = [0.7*search+.5, 1.3*search+.5]     ; line center in pixels
    ;parinfo[2].limited = 1b                                 ; limit line sigma   
    ;parinfo[2].limits  = [.6*width, 1.4*width]              ; limit sigma width in pixels
    parinfo[2].fixed   = 1b                                 ; fix linewidth sigma
    parinfo[2].Value   = width                              ; Fix linewidth using the measured linewidth ----> hack: haven't yet made certain this is spatially invariant 
    ESTIMATES          = [parinfo[0].limits[0], search, width] ; A first guess at the three parameters above

  ; Define the arrays that we'll write to:
    D2_brightness = fltarr(s[2], s[1])     & D2_err_brightness = fltarr(s[2], s[1])
    D2_linewidth  = fltarr(s[2], s[1])     & D2_err_linewidth = fltarr(s[2], s[1])
    D1_brightness = fltarr(s[2], s[1])     & D1_err_brightness = fltarr(s[2], s[1])
    D1_linewidth  = fltarr(s[2], s[1])     & D1_err_linewidth = fltarr(s[2], s[1])
    D2_brightness_sum = fltarr(s[2], s[1]) & D1_brightness_sum = fltarr(s[2], s[1]) ; These are to compare the Gaussian fits of the line profile with a straight sum over all pixels within "search" of the trace 
  
  ; sum the LSF over 3 sigma 99.7% power  
    LSF_Sum_width = 3.*mean(D2_width[real])
    A_array = fltarr(s[2]*s[1],3)
    cound = 0
    
  Print, 'Fitting exosphere D2 in every spatial bin of every frame (slow, ~50 frames/min)...'
  for n = 0, s[2]-1 do begin
    img = exosphere_spectra_cube[*,*,n]
    err_img = sqrt( abs(img) )     
    for i = 0, s[1]-1 do begin
      if LSF eq 'sum' then D2_brightness_sum[n,i] = total(img[round(D2_trace[i]-LSF_Sum_width):round(D2_trace[i]+LSF_Sum_width), i], /NaN) ;sum of pixels to compare with the area under the Gaussian fit
      if LSF eq 'fit' then begin 
        ESTIMATES[0] = D2_height[i]
        result = mpfitpeak(findgen(search*2. + 1), img[D2_trace[i]-search:D2_trace[i]+search, i], A, PERROR = Err_A, ESTIMATES = ESTIMATES, $
                           error = err_img[D2_trace[i]-search:D2_trace[i]+search, i], /NAN, STATUS = STATUS, PARINFO = PARINFO, NTERMS = 3)
        ;if STATUS ge 1 then begin
         if (STATUS ge 1) and (Err_A[0] ne 0.) and (Err_A[1] ne 0.) then begin
          D2_brightness[n,i]     = A[0]*A[2]*SQRT(2*!DPI) 
          D2_err_brightness[n,i] = D2_brightness[n,i] * sqrt( (err_A[0]/A[0])^2. + (err_A[2]/A[2])^2. )
          D2_linewidth[n,i]      = Dispersion_V*sqrt(((2.0*sqrt(2.0*alog(2)))*A[2])^2 - (2.45)^2) ;FWHM KM/S instrumental measured width is 2.45 pixels 
          D2_err_linewidth[n,i]  = Dispersion_V*(2.0*sqrt(2.0*alog(2)))*Err_A[2]
          ;A_array[count,*]       = A

        endif else begin
          D2_brightness[n,i]     = !values.F_nan
          D2_linewidth[n,i]      = !values.F_nan
          ;A_array[count,*]       = !values.F_nan
       endelse
      endif ; fit lsf
      count = count+1
    endfor ; i row number
    print, '   Finished D2 extraction of frame', n 
  endfor ; n frame number 

  Print, 'Fitting exosphere D1 in every spatial bin of every frame (slow, ~50 frames/min)...'
  for n = 0, s[2]-1 do begin
    img = exosphere_spectra_cube[*,*,n]
    err_img = sqrt( abs(img) )
    for i = 0, s[1]-1 do begin
      if LSF eq 'sum' then D1_brightness_sum[n,i] = total(img[round(D1_trace[i]-LSF_Sum_width):round(D1_trace[i]+LSF_Sum_width), i], /NaN) ;sum of pixels to compare with the area under the Gaussian fit
        if LSF eq 'fit' then begin   
        ESTIMATES[0] = D1_height[i]
        result = mpfitpeak(findgen(search*2. + 1), img[D1_trace[i]-search:D1_trace[i]+search, i], A, PERROR = Err_A, ESTIMATES = ESTIMATES, $
          error = err_img[D1_trace[i]-search:D1_trace[i]+search, i], /NAN, STATUS = STATUS, PARINFO = PARINFO, NTERMS = 3)
        ;if STATUS ge 1 then begin
        if (STATUS ge 1) and (Err_A[0] ne 0.) and (Err_A[1] ne 0.) then begin
          D1_brightness[n,i]     = A[0]*A[2]*SQRT(2*!DPI)
          D1_err_brightness[n,i] = D1_brightness[n,i] * sqrt( (err_A[0]/A[0])^2. + (err_A[2]/A[2])^2. )
          D1_linewidth[n,i]      = Dispersion_V*sqrt(((2.0*sqrt(2.0*alog(2)))*A[2])^2 - (2.45)^2) ;FWHM KM/S instrumental measured width is 2.45 pixels
          D1_err_linewidth[n,i]  = Dispersion_V*(2.0*sqrt(2.0*alog(2)))*Err_A[2]
        endif else begin
          D1_brightness[n,i]     = !values.F_nan
          D1_linewidth[n,i]      = !values.F_nan
        endelse
      endif ; fit lsf
    endfor
    print, '   Finished D1 extraction of frame', n 
  endfor 
   
  if LSF eq 'sum' then D1_brightness = D1_brightness_sum
  if LSF eq 'sum' then D2_brightness = D2_brightness_sum
  
  loadct, 0 
  wset, 1
    cgimage, bytscl(D2_brightness, 0, 3000), /axes, /keep_aspect, xtitle = 'frame number', ytitle = 'D2 Brightness along slit' 
  wset, 2
    cgimage, bytscl(D1_brightness, 0, 3000), /axes, /keep_aspect, xtitle = 'frame number', ytitle = 'D1 Brightness along slit'
  cgloadct, ct

  ;D2_brightness = sigma_filter( D2_brightness, 15, N_sigma=6, /ALL,/MON ) ;hack hack hack dealing with a few hot pixels on Dec 12 frame 351
    
  save, D2_brightness, D2_brightness_sum, D2_linewidth, D1_brightness, D1_brightness_sum, D1_linewidth, D2_Trace, D1_Trace, PSF_1D, filename = outdir + 'brightness.sav'
  beep
endif

; =====================================================================================================================
; Part 6 : Extract the 1D surface reflectance, and flux calibrate the data
; =====================================================================================================================
if part eq 6 or part eq 99 then begin
  if part ne 99 then restore, outdir + 'spectral_cubes.sav' ; We'll need the cube of (co-aligned?) surface reflectance spectra
  if part ne 99 then restore, outdir + 'brightness.sav'     ; We'll need the D1 & D2 Traces to figure out which pixel the calibration wavelength is located at.
  if part ne 99 then restore, outdir + 'shift_array.sav'    ; Contains shift_arrays and aligned_cubes of both channels, sharpness_metric & correl
  if part ne 99 then restore, outdir + 'Hapke_Models.sav'
  if part ne 99 then restore, outdir + 'rotational_alignments.sav' 
  s  = size(aligned_spectra_cube, /dimensions)

  ; Reconstruct an image from the spectral channel at the calibration wavelength
    ; Find the pixel location of the calibration wavelength
      Cal_wavelength_trace    = fltarr(N_elements(D2_Trace))
      for i = 0, N_elements(D2_Trace)-1 do Cal_wavelength_trace[i] = interpol( [D2_Trace[i], D1_trace[i]], [Na_D2_Wl, Na_D1_Wl], Flux_Cal_Wavelength)
      Cal_wavelength_trace = round(Cal_wavelength_trace)

    ; As in part 5 for sodium, extract brightness slices from the "surface reflectance cube" (the "exosphere cube") already has the reflectance component removed
      DN_per_A  = fltarr(s[2], s[1])
      for n = 0, s[2]-1 do begin
        img = aligned_spectra_cube[*,*,n]
        for i = 0, s[1]-1 do begin
          DN_per_A[n, i] = TOTAL(img[Cal_wavelength_trace[i]-(0.5/dispersion):Cal_wavelength_trace[i]+(0.5/dispersion), i], /NaN) ;sum over 1 Angstrom
        endfor
      endfor

  ; Reconstruct Flux Calibration Image: Place this time series of 1D spectral slices along the slit into a (time averaged) 2D image
    flux_cube         = fltarr(s)   ; This datacube will hold the image of the slit for each frame, slit is extacted at the calibration wavelength over 1 Angstrom
    flux_cube_weights = fltarr(s)   ; This datacube will hold the weights from the point spread function of the slit for each frame
    Na_D1_cube        = fltarr(s)   ; This datacube will hold the image of the slit for each frame, slit is extacted at exosphere Na D1
    Na_D2_cube        = fltarr(s)   ; This datacube will hold the image of the slit for each frame, slit is extacted at exosphere Na D2

  ; Gather top % that are "lucky"
    ranked             = SORT(sharpness_metric)                      ; Rank lowest to highest seeing
    cut_index          = fix(s[2]*Lucky_fraction)
    mean_seeing        = mean(sharpness_metric[ranked[0:cut_index]])
    print, 'Lucky threshold effective seeing:', sharpness_metric[ranked[cut_index]]
    print, 'Mean effective seeing:', mean_seeing
    cghistoplot, sharpness_metric, title = 'Effective seeing (arcsec)', binsize = 0.1
    cgplot, [sharpness_metric[ranked[cut_index]],sharpness_metric[ranked[cut_index]]], [0, s[2]], /overplot

  ; Blur a Hapke Model of Mercury at the observation time/geometry. Match its plate scale to the RIPS spectral channel
    seeing_sigma      = mean_seeing / (2.0*sqrt(2.0*alog(2)))        ; convert FWHM to sigma
    seeing_sigma      = seeing_sigma/Spectral_platescale             ; convert to Gaussian sigma in pixels
    Blurred_Hapke = convol(Spectra_Hapke, GAUSSIAN_FUNCTION([seeing_sigma, seeing_sigma]), /EDGE_ZERO, /normalize)

    sh = size(Blurred_Hapke, /dimensions)                            ; Hapke image dimensions
    si = size(aligned_imaging_cube, /dimensions)                     ; Imaging channel image dimensions

  ; Put in the stripes long the slit with a width == line spread function from arc frames == spectral channel's "point spread function" 
    PSF_Size        = N_elements(PSF_1D)
    spatial_weights = rebin(PSF_1D, psf_size, s[1])           ; bin it into 2D to use in the weighting scheme of the next part 
    home            = fix(round([s[0]/2.,s[1]/2.]))
    
    ;window, 0

  ;    Spectra_shift_array[*,0] = MEDSMOOTH(Spectra_shift_array[*,0], 5 ) ; interesting, plot this up and check out the noise, smoothing the shift array seems to help *slightly*
  ;    Spectra_shift_array = smooth(Spectra_shift_array, [5,1])  
  for i = 0, s[2]-1 do begin                                  ; loop over frames / time. 
    
    ; Apply the imaging channel rotations from Celestial North measured earlier in Part 1
      aligned_imaging_cube[*,*,i] = rot(aligned_imaging_cube[*,*,i], rotational_alignments[i], /interp) ; The imaging channel has been aligned in 
    
    ; Spread the 1-D spatial information over the slit into 2-D at all X indices  
      flux_cube[*,*,i]            = rebin(DN_per_A[i,*], s[0], s[1])
      Na_D2_cube[*,*,i]           = rebin(D2_Brightness[i,*], s[0], s[1])
      Na_D1_cube[*,*,i]           = rebin(D1_Brightness[i,*], s[0], s[1])
 
    ; Weights will be at a certain X index where the slit was
    ; All spectral cubes have already been aligned in -Y/plate_scale_ratio. Only +X/plate_scale_ratio is needed now
      weighting_frame             = fltarr(s[0], s[1])        ; This will hold the weighting of the Gaussian LSF/PSF      
      weighting_frame[home[0] + Spectra_shift_array[i,0] - PSF_Size/2.: home[0] + Spectra_shift_array[i,0] + PSF_Size/2.-1,*] = spatial_weights
      flux_cube_weights[*,*,i]    = weighting_frame
  endfor

  ; Apply the rotations from Celestial north measured earlier
    for i = 0, s[2]-1 do begin      
      flux_cube[*,*,i]         = rot(flux_cube[*,*,i], -rotational_alignments[i], /interp)
      Na_D2_cube[*,*,i]        = rot(Na_D2_cube[*,*,i], -rotational_alignments[i], /interp)
      Na_D1_cube[*,*,i]        = rot(Na_D1_cube[*,*,i], -rotational_alignments[i], /interp)
      flux_cube_weights[*,*,i] = rot(flux_cube_weights[*,*,i], -rotational_alignments[i], /interp)
    endfor

  ; now take the geometric mean of the whole stack using the per-pixel weighting
    Calib_img    = total(flux_cube*flux_cube_weights, 3, /NAN)  / total(flux_cube_weights, 3) ; weighted average
    Na_D2_img    = total(Na_D2_cube*flux_cube_weights, 3, /NAN) / total(flux_cube_weights, 3) ; weighted average
    Na_D1_img    = total(Na_D1_cube*flux_cube_weights, 3, /NAN) / total(flux_cube_weights, 3) ; weighted average
    some_valid_data       = where(finite(Calib_img), complement = no_brightness_information)
    Calib_img[no_brightness_information] = 0.
    Na_D2_img[no_brightness_information] = 0.
    Na_D1_img[no_brightness_information] = 0.

  ; Similarly, take the geometric mean of just the lucky ones
    Lucky_Calib_img    = total(flux_cube[*,*,ranked[0:cut_index]]*flux_cube_weights[*,*,ranked[0:cut_index]], 3, /NAN)  / total(flux_cube_weights[*,*,ranked[0:cut_index]], 3, /NAN)
    Lucky_Na_D2_img    = total(Na_D2_cube[*,*,ranked[0:cut_index]]*flux_cube_weights[*,*,ranked[0:cut_index]], 3, /NAN) / total(flux_cube_weights[*,*,ranked[0:cut_index]], 3, /NAN)
    Lucky_Na_D1_img    = total(Na_D1_cube[*,*,ranked[0:cut_index]]*flux_cube_weights[*,*,ranked[0:cut_index]], 3, /NAN) / total(flux_cube_weights[*,*,ranked[0:cut_index]], 3, /NAN)
    some_valid_data    = where(finite(Lucky_Calib_img), complement = no_brightness_information)
    Lucky_Calib_img[no_brightness_information] = 0.
    Lucky_Na_D2_img[no_brightness_information] = 0.
    Lucky_Na_D1_img[no_brightness_information] = 0.

  ; Assemble the imaging channel results, expand them to match the spectral platescale
    Imaging_Ch = total(aligned_imaging_cube, 3, /Nan)                                                & Lucky_Imaging_Ch = total(aligned_imaging_cube[*,*,ranked[0:cut_index]], 3, /Nan)
    Imaging_Ch = rotate(Imaging_Ch, 7)                                                               & Lucky_Imaging_Ch = rotate(Lucky_Imaging_Ch, 7)
    Imaging_Ch = congrid(Imaging_Ch, si[0]/plate_scale_ratio, si[1]/plate_scale_ratio, cubic = -0.5) & Lucky_Imaging_Ch = congrid(Lucky_Imaging_Ch, si[0]/plate_scale_ratio, si[1]/plate_scale_ratio, cubic = -0.5)

  ; Crop and Co-align everything
    Blurred_Hapke_align = Blurred_Hapke

  ; Crop everthing derived from the spectral channel and align it to the blurred Hapke function
  ; crop calib image to only regions of valid data since all the zeros will throw off cross-correlation alignments
    valid = minmax( where(total(Calib_img, 2) ne 0.) )
    zstack_align_images, Calib_img[valid[0]:valid[1],*], Blurred_Hapke_align[valid[0]:valid[1],*], x, y
    Calib_img = shift(Calib_img, x, y)
    Na_D2_img = shift(Na_D2_img, x, y)
    Na_D1_img = shift(Na_D1_img, x, y)
    
;    ; inspection frame by frame
;      loadct, 0
;      flux_cube_weights = shift(flux_cube_weights, x, y, 0) 
;      flux_cube         = shift(flux_cube, x, y, 0) 
;      for i = 0, s[2]-1 do begin 
;        cgimage, flux_cube_weights[*,*,i] * flux_cube[*,*,i], /axes, minvalue = 0, maxvalue = .3*max(Calib_img)
;        cgcontour, Blurred_Hapke, /onimage, color = 'green', nlevels = 4, thick=.5
;        wait, 0.04
;      endfor  
;      cgloadct, ct
    
  ; Center of the planet's disk is at Hapke_Center, preserve this in the subframes
    hapke_center  = fix(round(sh/2.))
    Spectra_Hapke = Spectra_Hapke[hapke_center[0]-subframe:hapke_center[0]+subframe-1,hapke_center[1]-subframe:hapke_center[1]+subframe-1]
    Blurred_Hapke = Blurred_Hapke[hapke_center[0]-subframe:hapke_center[0]+subframe-1,hapke_center[1]-subframe:hapke_center[1]+subframe-1]
    Calib_img     = Calib_img[hapke_center[0]-subframe:hapke_center[0]+subframe-1,hapke_center[1]-subframe:hapke_center[1]+subframe-1]
    Na_D1_img     = Na_D1_img[hapke_center[0]-subframe:hapke_center[0]+subframe-1,hapke_center[1]-subframe:hapke_center[1]+subframe-1]
    Na_D2_img     = Na_D2_img[hapke_center[0]-subframe:hapke_center[0]+subframe-1,hapke_center[1]-subframe:hapke_center[1]+subframe-1]
    Abs_reference = Hi_Res_Hapke                                                                                                                 ;North-aligned imaging reference full resolution
    Abs_reference_small = Picture_Hapke[hapke_center[0]-subframe:hapke_center[0]+subframe-1,hapke_center[1]-subframe:hapke_center[1]+subframe-1] ; North-aligned imaging reference matched plate scale

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
    Spectra_Hapke    = rotate(spectra_Hapke, 7)

  ; Setup inspection windows and plot everything
    window, 0, xpos=winpos_x,                  ypos=winpos_y,               xs=8*subframe, ys=8*subframe, title='EXTRACTED DN ALONG IN A 1 ANGSTROM BAND (DN/A)'
    window, 1, xpos=winpos_x+8*subframe+20,    ypos=winpos_y,               xs=8*subframe, ys=8*subframe, title='EXTRACTED NA D2 DN'
    window, 2, xpos=winpos_x+8*subframe+20,    ypos=winpos_y+8*subframe+40, xs=8*subframe, ys=8*subframe, title='EXTRACTED NA D1 DN'
    window, 3, xpos=winpos_x,                  ypos=winpos_y+8*subframe+40, xs=8*subframe, ys=8*subframe, title='LUCKY EXTRACTED DN ALONG IN A 1 ANGSTROM BAND (DN/A)'
    window, 4, xpos=winpos_x+2*(8*subframe+20),ypos=winpos_y+8*subframe+40, xs=8*subframe, ys=8*subframe, title='ALIGNED & COADDED IMAGING CHANNEL, MATCHED TO SPECTRAL CHANNEL PLATESCALE'
    window, 5, xpos=winpos_x+2*(8*subframe+20),ypos=winpos_y,               xs=8*subframe, ys=8*subframe, title='HAPKE MODEL (MR/A)'   
    Continuum_colorbar_top = max(Blurred_Hapke) ; MR/A
    wset, 0
    cgimage, Calib_img, minvalue = 0, maxvalue = Continuum_colorbar_top
    cgColorbar, minrange = 0., maxrange = Continuum_colorbar_top, title = cgsymbol('Sigma')+'1'+cgsymbol('Angstrom') +' at '+ strcompress(Flux_Cal_Wavelength)+cgsymbol('Angstrom') + ' (MegaRayleighs /'+ cgsymbol('Angstrom')+')'
    wset, 1
    cgimage, Na_D2_img*1000.;, minvalue = 0, maxvalue = 2000.
    cgColorbar, minrange = 0, maxrange = max(Na_D2_img*1000.), title = 'Na D2 (KiloRayleighs)'
    wset, 2
    cgimage, Na_D1_img*1000.;, minvalue = 0, maxvalue = 1200.
    cgColorbar, minrange = 0, maxrange = max(Na_D1_img*1000.), title = 'Na D1 (KiloRayleighs)'
    wset, 3
    cgimage, lucky_Calib_img, minvalue = 0, maxvalue = Continuum_colorbar_top
    cgColorbar, minrange = 0., maxrange = Continuum_colorbar_top, title = cgsymbol('Sigma')+'1'+cgsymbol('Angstrom') +' at '+ strcompress(Flux_Cal_Wavelength)+cgsymbol('Angstrom') + ' (MegaRayleighs/'+ cgsymbol('Angstrom')+')'
    wset, 4
    cgimage, Imaging_Ch, /keep_aspect
    ;cgimage, Lucky_Imaging_Ch, /keep_aspect
    wset, 5
    cgimage, Blurred_Hapke, minvalue = 0, maxvalue = Continuum_colorbar_top ;Blurred_Hapke

  save, Lucky_Imaging_Ch, Imaging_Ch, Lucky_Calib_img, Calib_img, Lucky_Na_D2_img, Lucky_Na_D1_img, Na_D2_img, Na_D1_img, Blurred_Hapke, Spectra_Hapke, Abs_reference, Abs_reference_small, $
    mean_seeing, filename = outdir + 'images_to_plot.sav'
endif

; =====================================================================================================================
; Part 7 : Plot all images and and make movies of the exosphere.
; =====================================================================================================================
if part eq 7 or part eq 99 then begin
  if part ne 99 then restore, outdir + 'Hapke_Models.sav'                     ; contains all images to be plotted
  if part ne 99 then restore, outdir + 'images_to_plot.sav'                   ; contains all images to be plotted

  ; Find Mercury's True Anomaly Angle
    filenames = file_Search(Mercury_Dir + Mercury_files, count = n_flats)
    header = headfits(filenames[0])
    UTC = sxpar(header, 'DATE')
    cspice_str2et, UTC, et
    cspice_spkezr, Body, et, 'J2000', 'LT+S', 'Sun', state, ltime
    CSPICE_FURNSH, STRCOMPRESS(FILE_SEARCH(kernel_directory, 'Gravity.tpc'))  ; Solar Gravitational Constant  
    cspice_bodvrd, 'Sun', 'GM', 1, mu                                         ; The Keplerian GM of the Sun in units km^3/s^2
    cspice_oscelt, state, et, mu[0], elts                                     ; Comupute osculating orbital elements
    ecc = ELTS[1]                                                             ; Eccentricity
    MA  = ELTS[5]                                                             ; Mean Anomaly radians
    ; Now get true anomaly from Mean anomaly and Eccentricity, see https://en.wikipedia.org/wiki/True_anomaly
    ; Roy, A.E. (1988). Orbital Motion (1 ed.). Bristol, UK; Philadelphia, PA: A. Hilger. ISBN 0852743602.
      TAA =  MA + (2.*ecc - 0.25*ecc^3)*sin(MA) + 1.25*(ecc^2)*sin(2.*MA)+ (13./12.)*(ecc^3)*sin(3.*MA)
      True_Anomaly = TAA*!radeg ; True Anomaly in degrees
      print, 'Mercury''s True Anomaly Angle:', True_Anomaly

  ; Adjust the size of the Hapke image with its albedo features
    cspice_bodvrd, Body, 'RADII', 3, radii                                    ; Get Body shape constants
    cspice_spkpos, Body, ET, 'J2000', 'LT+S', 'Earth', ptarg, ltime
    R_M = 206264.806 * atan(radii[0] / norm(ptarg))                           ; Radius of the body in arcsec
    print, 'Angular diameter of ', body, ' =', 2.*R_M, ' arcsec'
    ;---------------------------------------------- Get a Mercury lat & lon grid ------------------------------------------------

    ; Find the position angle of the north pole vector, that is, the separation between Body's north pole and celestial north pole.
      CSPICE_SPKEZR, body, et, 'J2000', 'LT', 'Earth', BODY_state, ltime
      CSPICE_RECRAD, BODY_state[0:2], dist, ra, dec ;Convert rectangular coordinates to RA and Dec

    ; Find a vector from body's center to it's north pole in the J2000 frame at the ephemeris time, NEGLECT the planet's oblateness
      North_body_fixed = [0.,0.,1.]*radii[2]                    ; location of the North Pole in body-fixed coords
      cspice_pxform, 'IAU_'+Body, 'J2000', et - ltime, To_J2000 ; Find body-fixed coords to J2000 rotation matrix
      cspice_mxv, To_J2000, North_body_fixed, North_J2000       ; Rotate to J2000

    ; For comparison with JPL Horizons, get the RA and DEC of the pole direction
    ;CSPICE_RECRAD, North_J2000, dist, North_dir_ra, North_dir_dec ;Convert rectangular coordinates to RA and Dec

    ; Get the RA and Dec of the the body's north pole
      Pole_state = body_state + North_J2000
      CSPICE_RECRAD, Pole_state, dist, Pole_ra, Pole_dec ;Convert rectangular coordinates to RA and Dec

    ; Precess each RA and Dec to the Current Epoch. First we'll find the year and fraction of a year. we're interested in:
      Current_epoch = 2000. + et / 3.1556926e7                           ; J2000 + seconds past J2000 / seconds per year
      precess, Pole_ra, Pole_DEC, 2000, Current_epoch, /RADIAN           ; Precessed Ra and dec is now applied to Pole location
      precess, ra, DEC, 2000, Current_epoch, /RADIAN                     ; Precessed Ra and dec is now applied to body center location
      Delta_RA      = Pole_ra - ra
      Delta_Dec     = Pole_Dec - dec
      theta         = sqrt((Delta_RA*cos(Dec))^2.D + Delta_Dec^2.D)
      NPAng         = signum(Delta_RA)*Acos(Delta_Dec / theta) / cspice_rpd()
      print, 'North Pole Position Angle (CCW, E of N) = ', NPAng

    ; Find the sub-observer planetographic longitude and latitude
      cspice_subpnt, 'Near point: ellipsoid', Body, et, 'IAU_'+body, 'LT+S', 'Earth', $
        spoint, trgepc, srfvec
      f = ( radii[0]-radii[2] ) / radii[0]                               ; flatness parameter, mercury is round so this is zero
      cspice_recpgr, body, spoint, radii[0], f, spglon, spglat, spgalt
      lat_se = spglat * cspice_dpr()                                     ; planetographic latitude of the sub-observer point on the surface. (0 at equator, 90 at north pole, -90 at south pole)
      lon_se = spglon * cspice_dpr()                                     ; planetographic longitude of the sub-observer point on the surface. (0 to 360 west longitude).

    ; setup the x-y grid image
      s        = size(Calib_img, /dimensions)
      xdim     = s[0]*10
      ydim     = s[1]*10
      ctr_xpix = xdim / 2.
      ctr_ypix = ydim / 2.
      pix2km   = Spectral_platescale * (1./R_M) * radii[0]/10.               ; "/pix * radii/" * km/radii
      x        = (dindgen(xdim) - ctr_xpix) * pix2km
      y        = (dindgen(xdim) - ctr_ypix) * pix2km
      xsq      = dblarr(xdim, ydim)                                      ; xsq is an image the size of the calib-img where the values are the horizontal distance from the "center" in km
      ysq      = dblarr(xdim, ydim)                                      ; ysq is an image the size of the calib-img where the values are the vertical distance from the center.
      for i = 0, xdim-1 do begin
        xsq[*,i] = x
        ysq[i,*] = y
      endfor

    ; get the planetographic lon/lat at each pixel
      ob           = deprob(xsq, ysq, radii[0], radii[2], lat_se, lon_se, npang=NPAng)
      ob.lon[where(ob.lon eq -666, /Null)] = !values.F_nan
      ob.lat[where(ob.lat eq -666, /Null)] = !values.F_nan

    ; are plot annotations going on the left or rigth side? 
      CSPICE_SPKEZR, 'Sun', et, 'J2000', 'LT', 'Earth', sun_state, ltime
      CSPICE_RECRAD, sun_state[0:2], dist, ra_sun, dec_sun    ;Convert rectangular coordinates to RA and Dec
      if (ra - ra_sun)*!radeg lt -28. then ra     = ra + 2.*!pi       ; wrap zero crossings by 2pi
      if (ra - ra_sun)*!radeg gt  28. then ra_sun = ra_sun + 2.*!pi   ; wrap zero crossings by 2pi
      if (ra - ra_sun lt 0.) then side = 'right' else side = 'left'

  cgPS_Open, filename = outdir+body+'_'+Telescope+'_'+night+'.eps', /ENCAPSULATED, xsize = 6, ysize = 6
    !P.font=1
    device, SET_FONT = 'Helvetica Bold', /TT_FONT
    !p.charsize = 1.2
    yspace = 4.0*(!D.Y_CH_SIZE)/!D.Y_SIZE
    xspace = 4.0*(!D.X_CH_SIZE)/!D.X_SIZE
    pos = cgLayout([2,2], xgap = 1., ygap = 1, oxmargin = 4, oymargin = 1)
    thick = 4.
    
    ; Top Left
      P = pos[*,0]
      loadct, 0
      tvimage, bytscl(Abs_reference), position = P, /noframe
      ; Annotate with angular size of the disk
        s         = size(Abs_reference, /dimensions)
        collapsed = total(rot(Abs_reference, npang), 1)
        junk      = where(collapsed gt 0., count)
        halfwidth = (float(count) / float(N_elements(collapsed))) * (P[3] - P[1]) * 0.5
        case side of
          'left': begin
          cgarrow, p[0]+.1*(p[2]-P[0]), P[1]+(P[3]-P[1])/2., p[0]+.1*(p[2]-P[0]), P[1]+(P[3]-P[1])/2.+ halfwidth, /SOLID, color = 'white', /normal
          cgarrow, p[0]+.1*(p[2]-P[0]), P[1]+(P[3]-P[1])/2., p[0]+.1*(p[2]-P[0]), P[1]+(P[3]-P[1])/2.- halfwidth, /SOLID, color = 'white', /normal
          cgText,  p[0]+.09*(p[2]-P[0]), P[1]+(P[3]-P[1])/2. , string(2.*R_M, format = '(F4.2)')+'"', color = 'white', /NORMAL, ALIGNMENT=0.5, charsize = 1., orientation = 90  
          end
          'right': begin
            cgarrow, p[0]+.9*(p[2]-P[0]), P[1]+(P[3]-P[1])/2., p[0]+.9*(p[2]-P[0]), P[1]+(P[3]-P[1])/2.+ halfwidth, /SOLID, color = 'white', /normal
            cgarrow, p[0]+.9*(p[2]-P[0]), P[1]+(P[3]-P[1])/2., p[0]+.9*(p[2]-P[0]), P[1]+(P[3]-P[1])/2.- halfwidth, /SOLID, color = 'white', /normal
            cgText,  p[0]+.95*(p[2]-P[0]), P[1]+(P[3]-P[1])/2. , string(2.*R_M, format = '(F4.2)')+'"', color = 'white', /NORMAL, ALIGNMENT=0.5, charsize = 1., orientation = 90
          end
        endcase     
      cgloadct, ct
  
    ; Top Right
      P = pos[*,1]
      p[0] = p[0] - .01
      p[2] = p[2] - .01
      levels = (findgen(5) + 1) /5.
      cgcontour, Blurred_Hapke / max(Blurred_Hapke, /NaN), levels=levels, label=0, /noerase, pos = P, color = 'black', aspect = 1, thick = 2.0, XTICKFORMAT = "(A1)", YTICKFORMAT= "(A1)"
      cgcontour, Imaging_Ch / max(Imaging_Ch, /NaN), levels=levels, label=0, pos = P, /noerase, color = 'lime green', thick = 2.0, XTICKFORMAT = "(A1)", YTICKFORMAT= "(A1)", c_LINESTYLE = [1,1,1,1,1]     
      cgcontour, smooth(Calib_img / max(Calib_img, /NaN), smooth_Final, /edge_truncate), levels=levels, label=0, pos = P, /noerase, color = 'red', thick = 2.0, XTICKFORMAT = "(A1)", YTICKFORMAT= "(A1)", c_LINESTYLE = [1,1,1,1,1]
      ;cgcontour, spectra_hapke / max(spectra_hapke, /NaN), levels=levels, label=0, pos = P, /noerase, color = 'green', thick = 1., XTICKFORMAT = "(A1)", YTICKFORMAT= "(A1)"
      if keyword_set(plot_lucky_Na) then begin 
        cgcontour, Lucky_Imaging_Ch / max(Lucky_Imaging_Ch, /NaN), levels=levels, label=0, pos = P, /noerase, color = 'lime green', thick = 2.0, XTICKFORMAT = "(A1)", YTICKFORMAT= "(A1)" 
        cgcontour, smooth(Lucky_Calib_img / max(Lucky_Calib_img, /NaN), smooth_Final, /edge_truncate), levels=levels, label=0, pos = P, /noerase, color = 'red', thick = 2.0, XTICKFORMAT = "(A1)", YTICKFORMAT= "(A1)"
      endif

    if keyword_set(plot_lucky_Na) then Na_img = (Lucky_Na_D1_img + Lucky_Na_D2_img) else Na_img = (Na_D1_img + Na_D2_img)

    ; Bottom Left
      P = pos[*,2]
      axis_format = {XTICKFORMAT:"(A1)", YTICKFORMAT:"(A1)", Color:'white'} ;can't seem to get rid of the axis
      if keyword_set(plot_lucky_Na) then begin 
        cgimage, smooth(Lucky_Calib_img, smooth_Final, /edge_truncate, /Nan), /keep_aspect, pos = P, /noerase, minvalue = 0., maxvalue = max(smooth(Lucky_Calib_img, smooth_Final, /edge_truncate)), $
                 /axes, axkeywords = axis_format
        case side of 
          'left':  p_CB = [p[0] + xspace+0.02, p[1]+.02, p[0]+xspace+0.03, p[3]-.02]
          'right': p_CB = [p[2] - (xspace-0.03), p[1]+.02, p[2] - (xspace-0.04), p[3]-.02]
        endcase
        
        ; overlay the smoothed Na contours
        cgcontour, smooth(Na_img / max(Na_img, /NaN), smooth_Final, /edge_truncate), /onimage, color='black', pos = P, levels=levels, C_CHARSIZE=2.4, $
          C_annotation = string(levels)+' MR/'+cgsymbol('Angstrom'), C_LABELS=[1,1,1,1,1], thick = 0.5, axiscolor = 'white' ;can't seem to get rid of the axis, oh well!
        
        cgColorbar, pos = P_CB, minrange = 0., maxrange = max(Lucky_Calib_img),  /VERTICAL, $
          title = 'Continuum at '+ string(Flux_Cal_Wavelength, format = '(I5)') + cgsymbol('Angstrom') + ' (MR/'+ cgsymbol('Angstrom')+')'
      endif else begin   
        cgimage, smooth(Calib_img, smooth_Final, /edge_truncate), /keep_aspect, pos = P, /noerase, minvalue = 0., maxvalue = max(smooth(Calib_img, smooth_Final, /edge_truncate)) 
        case side of 
          'left':  p_CB = [p[0] + xspace+0.01, p[1]+.02, p[0]+xspace+0.02, p[3]-.02]
          'right': p_CB = [p[2] - (xspace-0.03), p[1]+.02, p[2] - (xspace-0.04), p[3]-.02]
        endcase
        cgColorbar, pos = P_CB, minrange = 0., maxrange = max(Calib_img),  /VERTICAL, $
          title = 'Continuum at '+ string(Flux_Cal_Wavelength, format = '(I5)') + cgsymbol('Angstrom') + ' (MR/'+ cgsymbol('Angstrom')+')'
      end
   
    ; Bottom Right
      P = pos[*,3]
      p[0] = p[0] - .01
      p[2] = p[2] - .01
      loadct, 3
      axis_format = {XTICKFORMAT:"(A1)", YTICKFORMAT:"(A1)"}
      
      Hapke_flat = Blurred_Hapke / Lucky_Calib_img
      Hapke_flat[where(finite(Hapke_flat, /INFINITY))] = 0.
      Hapke_flat[where(Hapke_flat gt 2.)] = 0.
      
      valid_pixels = where(Lucky_Calib_img ne 0.)
      Normalized_Imaging_Ch = Lucky_Imaging_Ch * total(Lucky_Calib_img[valid_pixels], /NAN) / total(Lucky_Imaging_Ch[valid_pixels], /NAN)
      Imaging_Ch_flat = Normalized_Imaging_Ch / Lucky_Calib_img
      Imaging_Ch_flat[where(finite(Imaging_Ch_flat, /INFINITY))] = 0.
      Imaging_Ch_flat[where(Imaging_Ch_flat gt 2.)] = 0.
      
        ;cgimage, smooth(Na_img, smooth_Final, /edge_truncate), /keep_aspect, minvalue = 0, maxvalue = max(smooth(Na_img, smooth_Final, /edge_truncate)), pos = P, /noerase, $
        ;cgimage, Na_img*Hapke_flat, /keep_aspect, minvalue = 0, maxvalue = max(Na_img*Hapke_flat), pos = P, /noerase, $
        ;cgimage, smooth(Na_img*Hapke_flat, smooth_Final, /edge_truncate), /keep_aspect, minvalue = 0, maxvalue = max(smooth(Na_img*Hapke_flat, smooth_Final, /edge_truncate)), pos = P, /noerase, $
        ;        /axes, axkeywords = axis_format
        cgimage, smooth(Na_img*Imaging_Ch_flat, smooth_Final, /edge_truncate), /keep_aspect, minvalue = 0, maxvalue = max(smooth(Na_img*Imaging_Ch_flat, smooth_Final, /edge_truncate)), pos = P, /noerase, $
                /axes, axkeywords = axis_format        
                
      if keyword_set(plot_lucky_Na) then begin
        ;cgcontour, smooth(Lucky_Calib_img / max(Lucky_Calib_img, /NaN), smooth_Final, /edge_truncate), /onimage, color='green', pos = P, levels=levels, C_CHARSIZE=2.4, $
        ;cgcontour, smooth(Blurred_Hapke / max(Blurred_Hapke, /NaN), smooth_Final, /edge_truncate), /onimage, color='lime green', pos = P, levels=levels, C_CHARSIZE=2.4, $
        cgcontour, smooth(Lucky_Imaging_Ch / max(Lucky_Imaging_Ch, /NaN), smooth_Final, /edge_truncate), /onimage, color='lime green', pos = P, levels=levels, C_CHARSIZE=2.4, $
          C_annotation = string(levels)+' MR/'+cgsymbol('Angstrom'), C_LABELS=[1,1,1,1,1], thick = 0.75
      endif else begin    
        cgcontour, smooth(Calib_img / max(Calib_img, /NaN), smooth_Final, /edge_truncate), /onimage, pos = P, levels=levels, C_CHARSIZE=2.4, color='lime green', $
        C_annotation = string(levels)+' MR/'+cgsymbol('Angstrom'), C_LABELS=[1,1,1,1,1], thick = 0.5
      endelse  
      case side of 
        'left':  p = [p[0] + xspace+0.02, p[1]+.02, p[0]+xspace+0.03, p[3]-.02]
        'right': p = [p[2] - (xspace-0.03), p[1]+.02, p[2] - (xspace-0.04), p[3]-.02]
      endcase
      cgColorbar, minrange = 0., maxrange = max(Na_img), pos = P, /VERTICAL, title = 'Sodium D1 + D2 (MR)', color = 'white'
  
    ; Annotations
      Case 1 of
        (True_Anomaly ge 100.)                           : format = '(F5.1)' 
        (True_Anomaly ge 10.) and (True_Anomaly lt 100.) : format = '(F4.1)'
        (True_Anomaly lt 10.)                            : format = '(F3.1)'
      Endcase
      cgText, 0.5, 0.9, 'RIPS --- ' + strmid(UTC,0,10) + ' ' + strmid(UTC,11,5) + ' --- True Anomaly '+ string(True_Anomaly, format=format)+ '!U'+cgsymbol('deg')+'!U', ALIGNMENT=0.5, /NORMAL, charsize = 1.9 
      cgText, 0.09, 0.845, 'Reflectance Model', /NORMAL, charsize = 0.9, color = 'white'
      cgText, 0.51, 0.845, 'Model w/ Effective Seeing = '+string(mean_seeing, format = '(F3.1)')+'"', /NORMAL, charsize = 0.9, color = 'black'
      cgText, 0.51, 0.8225, 'Co-Aligned Imaging       100%        '+strcompress(string(100*Lucky_fraction, format = '(I)'))+'%', /NORMAL, charsize = 0.9, color = 'lime green'
      cgText, 0.51, 0.80,   'Co-Aligned Spectra       100%        '+strcompress(string(100*Lucky_fraction, format = '(I)'))+'%', /NORMAL, charsize = 0.9, color = 'red'
      cgarrow, 0.68, 0.8075, 0.705, 0.8075, HSIZE = 0., color = 'red', /NORMAL, thick =2., linestyle = 1
      cgarrow, 0.68, 0.83, 0.705, 0.83, HSIZE = 0., color = 'lime green', /NORMAL, thick =2., linestyle = 1
      cgarrow, 0.76, 0.8075, 0.795, 0.8075, HSIZE = 0., color = 'red', /NORMAL, thick =2.
      cgarrow, 0.76, 0.83, 0.795, 0.83, HSIZE = 0., color = 'lime green', /NORMAL, thick =2.
  cgPS_close

  ;-------------------------------------------- Find the g-value (photon scattering rate) ----------------------------------------------
    READCOL, strcompress('C:\IDL\Generic Model V2\read_write\killen_5890_5900.txt'), Format='A,A', wavelength, intensity, /SILENT  ;sodium d lines
    start=where(wavelength eq '5890.000')
    finish=where(wavelength eq '5900.000')
    wavelength=double(wavelength(start:finish)) ;In angstroms
    intensity=double(intensity(start:finish))   ;Above the atmosphere at 1 AU, the differential flux in photons/cm2/s/Å
  
    CSPICE_SPKEZR, body, et, 'J2000', 'NONE', 'Sun', sun_state, ltime 
    sun_planet_vel = sun_state[3:5]
    sun_planet_pos = sun_state[0:2]
    r_sun = norm(sun_planet_pos)                         ; distance from the Sun (units of km for now)
    Vrad = (sun_planet_pos[0,*]*sun_planet_vel[0,*] + $  ; radial velocity in KM/S WRT the Sun 
            sun_planet_pos[1,*]*sun_planet_vel[1,*] + $
            sun_planet_pos[2,*]*sun_planet_vel[2,*]) / r_sun
    gvalue, 'Na-D', Vrad * 1000., r_sun / 149597871., wavelength, intensity, g
    print, 'Photon scattering rate ~'+string(g) +' photons/atom/s'    

  ;---------------------------------PSF deconvolution----------------------------------------------------------
    ;maxiter = 1 ;fast debug only
    ;if night eq '15' then maxIter           = 100 ; settled on this for simplicity after much fiddling
    ;if night eq '12' then maxIter           = 50 ; settled on this for simplicity after much fiddling
    if night eq '12' then maxIter           = 9 ; quick test
    if night eq '15' then maxIter           = 17 ; quick test
    if night eq '25' then maxIter           = 10 ; quick test
    
    seeing_sigma      = mean_seeing / (2.0*sqrt(2.0*alog(2.)))      ; convert FWHM to sigma
    seeing_sigma      = seeing_sigma/spectral_platescale            ; convert to Gaussian sigma in pixels
    PSF_2D            = GAUSSIAN_FUNCTION([seeing_sigma, seeing_sigma])
    PSF_2D            = PSF_2D/total(PSF_2D)
    PSF_size          = size(PSF_2D, /dim)
    s                 = size(Lucky_Imaging_Ch, /dimensions)
    
    ;orig = blurred_hapke
    orig = Lucky_Imaging_Ch
    error = fltarr(maxIter+1)
    multipliers = 0
    deconv_old = orig
    deconv_orig = orig
    error[0] = 1.d0
    j = 0
    window, 0, xs = 400, ys = 400 & loadct, 0
    while (error[j] gt 1.d-3 and j lt maxIter) do begin
      j = j + 1
      max_likelihood, orig, PSF_2D, deconv, /gaussian, ft_psf=psf_ft
      ;max_likelihood, orig, PSF_2D, deconv, ft_psf=psf_ft
      ;Max_Entropy, orig, PSF_2D, deconv, multipliers, FT_PSF=psf_ft

      res = convol(deconv, PSF_2D, /EDGE_ZERO, /normalize, /nan)
      error[j] = total((res-deconv_orig)^2) / (1.d0*100.*100.)
      print, 'Iter: ', j, ' - Error : ', error[j]
      if (j ne 0) then begin
        if (error[j] gt error[j-1]) then begin
          deconv = deconv_old
          print, 'Error is increasing. Going back to the previous image...'
        endif else begin
          deconv_old = deconv
          cgimage, deconv, /keep_aspect
          print, 'Correlation w/ Hapke = ', max(CORREL_IMAGES( Spectra_Hapke, deconv ))
        endelse
      endif
      if night eq '15' then wait, 0.2
    endwhile
    Lucky_Imaging_Ch_deconv = deconv
    ;test_PSF = GAUSSIAN_FUNCTION([(4.08/3.)*seeing_sigma, (4.08/3.)*seeing_sigma])
    ;test_PSF = test_PSF[0:99, 0:99]
    ;Calib_img_deconv = image_deconvolve(orig, test_psf, sqrt(abs(orig))+5., mask = where(orig le 0.5), /positive) ;, guess = Spectra_Hapke

    orig              = Na_img * Imaging_Ch_flat
    multipliers       = 0
    error      = fltarr(maxIter+1)
    deconv_old = orig & deconv_orig = orig
    error[0]   = 1.d0
    j = 0
    window, 1, xs = 400, ys = 400 & cgloadct, ct
    while (error[j] gt 1.d-3 and j lt maxIter) do begin
      j = j + 1
      max_likelihood, orig, PSF_2D, deconv, /gaussian, ft_psf=psf_ft ;best w/ 17 iter (perkins)
      ;max_likelihood, orig, PSF_2D, deconv, ft_psf=psf_ft            ;best w/ 4 iter
      ;Max_Entropy, orig, PSF_2D, deconv, multipliers, FT_PSF=psf_ft

      res = convol(deconv, PSF_2D, /EDGE_ZERO, /normalize, /nan)
      error[j] = total((res-deconv_orig)^2) / (1.d0*100.*100.)
      print, 'Iter: ', j, ' - Error : ', error[j]
      if (j ne 0) then begin
        if (error[j] gt error[j-1]) then begin
          deconv = deconv_old
          print, 'Error is increasing. Going back to the previous image...'
        endif else begin
          deconv_old = deconv
          cgimage, deconv, /keep_aspect
        endelse
      endif
      if night eq '15' then wait, 0.2
    endwhile
    Na_img_Deconv = deconv * 1.e12/g[0]                   ; convert to column density in atoms * cm^-2
    Na_img        = Na_img * 1.e12/g[0]                   ; convert to column density in atoms * cm^-2

    window, 0, xs = 400, ys = 400
    cgimage, Na_img*Imaging_Ch_flat, /keep_aspect
    window, 1, xs = 400, ys = 400
    cgimage, Na_img_deconv, /keep_aspect
    loadct, 0
    window, 2, xs = 400, ys = 400
    cgimage, Lucky_Imaging_Ch, /keep_aspect
    window, 3, xs = 400, ys = 400
    cgimage, Lucky_Imaging_Ch_deconv
    window, 4, xs = 400, ys = 400, xpos = 400, ypos = 00
    cgimage, Spectra_Hapke
    window, 5, xs = 400, ys = 400, xpos = 400, ypos = 600
    cgimage, Blurred_Hapke
    window, 6, xs = 400, ys = 400, xpos = 800, ypos = 00
    tv,bytscl(congrid(PSF_2D, 400*psf_size[0]/s[0], 400*psf_size[1]/s[1]))
    
    WRITEFITS, outdir + 'Deconv_PSF_2D'+'_'+night+'.fits', PSF_2D
    WRITEFITS, outdir + 'Deconv_Calib'+'_'+night+'.fits', Imaging_Ch
    WRITEFITS, outdir + 'Deconv_Na'+'_'+night+'.fits', Na_img

    cgPS_Open, filename = outdir+body+'_'+Telescope+'_'+night+'_Deconvolved.eps', /ENCAPSULATED, xsize = 6, ysize = 6
      !P.font=1
      device, SET_FONT = 'Helvetica Bold', /TT_FONT
      
      axis_format = {XTICKFORMAT:"(A1)", YTICKFORMAT:"(A1)"}   
      cgimage, Lucky_Imaging_Ch_deconv, /axes, /keep_aspect, axkeywords = axis_format, title = 'Na Column x 10!U10!N cm!U-2!N' + '     ' + strmid(UTC,0,10) + ' ' + strmid(UTC,11,5)
      ;cgimage, Spectra_Hapke, /axes, /keep_aspect, axkeywords = axis_format, title = strmid(UTC,0,10) + ' ' + strmid(UTC,11,5) + ' --- ' + string(mean_seeing, format = '(F3.1)')+'" PSF Deconvolved'
      ;cgimage, PSF_2D, /noerase, /keep_aspect, position =  [0.0, 0.0, float(N_elements(PSF_2D[0,*])) / float(N_elements(Calib_img_deconv[0,*])), float(N_elements(PSF_2D[1,*])) / float(N_elements(Calib_img_deconv[1,*]))]
      
      lat_contours = indgen(17)*10 - 80                                  ; every 10 deg lat
      lon_contours = indgen(24)*15                                       ; every 15 deg lon
      C_LABELS_LAT = replicate(1, N_elements(lat_contours))
      C_LABELS_LON = replicate(1, N_elements(lon_contours))
      C_LABELS_LAT[indgen(N_elements(lat_contours)/2)*2 -1] = 0
      C_LABELS_LON[indgen(N_elements(lon_contours)/2)*2] = 0
      C_LABELS_LON[where( (lon_contours gt lon_se + 90.) or (lon_contours lt lon_se - 90.), /NULL )] = 0
      if night eq '15' then begin
        C_LABELS_LON[[17, 19]] = 0 ;these contour labels don't look good on this date
      endif
      if night eq '12' then begin
        C_LABELS_LON[*] = 0 ;these contour labels don't look good on this date
        C_LABELS_LON[[23]] = 1 ;these contour labels look good on this date
      endif
           
      cgcontour, ob.lat, /onimage, levels = lat_contours, C_LABELS = C_LABELS_LAT, color = 'green', THICK = .5;, C_Spacing = 5
      cgcontour, ob.lon, /onimage, levels = lon_contours, C_LABELS = C_LABELS_LON, color = 'green', THICK = .5;, C_Spacing = 5 ;'snow'
      
      cold_pole_90 = ob.lon
      cold_pole_90[where((cold_pole_90 gt 100) or (cold_pole_90 lt 80))] = !values.F_NaN
      cold_pole_270 = ob.lon
      cold_pole_270[where((cold_pole_270 gt 280) or (cold_pole_270 lt 260))] = !values.F_NaN
      cgcontour, cold_pole_90, /onimage, levels = [90], color = 'Blue'                        ; marks "cold-pole" longitudes, Cassidy et al. 2016
      cgcontour, cold_pole_270, /onimage, levels = [270], color = 'Blue'                      ; marks "cold-pole" longitudes, Cassidy et al. 2016
      ;LEVELS = [4,8,12,16,20,24,28]
      LEVELS = (indgen(20)+1)*2.5
      
      if telescope eq 'AEOS' then begin
        cgcontour, smooth(Na_img_deconv,8)/1.e10, /onimage, color = 'orange', levels = LEVELS  ; plot column density x 10^10 cm^-2 SMOOTH SMOOTH SMOOTH!!
      endif else begin
        cgcontour, Na_img_deconv/1.e10, /onimage, color = 'orange', levels = LEVELS  ; plot column density x 10^10 cm^-2 
      endelse  
   cgPS_Close   
endif

; =====================================================================================================================
; Part 8 : Write an MPEG4 Video of the datacube
; =====================================================================================================================
if part eq 8 or part eq 99 then begin
  raw_data_cube = MRDFITS(outdir + 'raw_data_cube.fits', 0, header, /unsigned) 
  s = size(raw_data_cube, /dim)
  width = s[0] / 2    ;Bin by 2x to conserve file size for AAS Journal submission
  height = s[1] / 2   ;Bin by 2x
  fps = 20
  window, 0, xs = height, ys = width
  print, 'Writing video to '+outdir + 'Raw_data_movie_'+night+'.mp4'
  oVid = IDLffVideoWrite(outdir + 'Raw_data_movie_'+night+'.mp4', FORMAT='mp4')
  vidStream = oVid.AddVideoStream(width, height, fps, CODEC='mpeg4', PRESET='medium')

  FOR i = 0, N_elements(raw_data_cube[0,0,*])-1 do begin
    tv, bytscl( reform( rebin(raw_data_cube[*,*,i], width, height) ) )
    RGB = transpose( rebin(raw_data_cube[*,*,i], width, height, 3), [2,0,1] )
    !NULL = oVid.Put(vidStream, bytscl(RGB))
  ENDFOR

  oVid = 0
endif

print, 'Run time = ',  (SYSTIME(/SECONDS) - Start_time) / 60, ' minutes'
stop
end
