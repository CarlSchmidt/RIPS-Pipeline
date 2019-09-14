pro Surface_Flux_Calibration, Body = body, Timestamp = Timestamp, Wavelength = Wavelength, Hapke_Platescale = Hapke_Platescale, $
                              Output_image = Output_image, MR_per_A = MR_per_A, Make_picture = Make_picture, align_celestial_north = align_celestial_north

; NAME:
;    Surface_Flux_Calibration
; PURPOSE:
;    Given, a body, timestamp, and wavelengths, generate an image in units of 
;    Rayleighs / Angstrom for flux calibration
; EXPLANATION:
;    Use Kurucz (2005) Solar Spectrum
;    Use SPICE for viewing geometry and solar flux scaling. 
;    Use some flavor of Hapke Model for viewing geometry and solar flux scaling.
;    Assumes viewing from Earth center, but this could be modified.
;    No absorption---value is at the top of the Earth's atmosphere. 
; CALLING SEQUENCE:
;    Surface_Flux_Calibration, Body, Timestamp, Wavelength, Display_Outputs = Display_Outputs
; INPUTS:
;    Body                  - String, e.g. 'Mercury'
;    Timestamp             - String, UTC in any format SPICE can digest 
;    Wavelength            - Float, wavelength in ANGSTROMS
;
; OPTIONAL KEYWORDS     
;    Output_image          - Path and filename to output image, e.g., '.../desktop/pretty_picture.eps'
;    MR_per_A              - Variable name of a 3x "radius" pixel output array in MegaRayleigh / Angstrom Units 
;    Make_picture          - Variable name of a 2x "radius" pixel output image array with albedo features in ???? (I/F reflectance units i think CS) 
;    Align_Celestial_North - If set, this aligns the image so that celestial north is the vertical axis.
;                            Default is for image outputs with the body's north pole vertically aligned.
;    Hapke_Platescale      - Platescale of both the output arrays, with arcsec/pixel units                       
; OUTPUTS:
;    Image in R/A
; DEPENDENCIES:
;    ---A Kurucz Solar flux model 
;    ---A SPICE installation
;    ---A Surface map for returning pretty pictures centered at (0,0) lat, lon  
;    Hapke_cs.pro
;    LOAD_SPICE.pro
;    VACTOAIR.pro
; EXAMPLE:
;    Surface_Flux_Calibration, Body = 'Mercury', Timestamp = 'June 25 2018 06:03', Wavelength = 5889., Make_picture = pretty
; REVISION HISTORY:
;    Written: C. Schmidt & L. Moore, 2018 (Boston University)
;    Modified: 8/29/2019 - Fixed sign error in north pole position angle

Viewpoint = 'Earth'
 ;'1994-01-30T16:07:03'
Loadct, 0
radius     = 150.d0                                           ; PLATE SCALE: Radius of body in pixels
display_outputs = 1
Longitude_slice = -50.                                        ; Luminance longitude to display Hapke reflectance, 0 = central meridian.

;=====================================Define Paths==============================================================
kernel_directory    = 'C:\SPICE\'
SOLAR_SPECTRUM_FILE = 'C:\IDL\Io\Kurucz_2005_irradthuwl.dat'  ; Solar Spectrum at 1AU in W/m2/nm
inpdir              = 'C:\IDL\CS_lib\Surface_Maps\'

If Body eq 'Mercury' then begin
  mosaic_file       = 'Mercury_MESSENGER_MDIS_Basemap_LOI_Mosaic_Global_32ppd.jpg'
  ;'Monochrome_Morphology_20170512_PDS16_equirectangular_thumb-sm.png'
  ; Mosaics taken from: 
  ; https://astrogeology.usgs.gov/search/map/Mercury/Messenger/Global/Mercury_MESSENGER_MDIS_Basemap_LOI_Mosaic_Global_166m
endif 
If Body eq 'Moon' then begin
  mosaic_file       = 'Lunar_Clementine_UVVIS_750nm_Global_Mosaic_1.2km.jpg'
  ; Mosaics taken from:
  ; https://astrogeology.usgs.gov/search/map/Moon/Clementine/UVVIS/Lunar_Clementine_UVVIS_750nm_Global_Mosaic_118m_v2
endif 

;======================================LOAD SPICE===============================================================
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
    CSPICE_KTOTAL, 'all', count
    PRINT, STRCOMPRESS('Loaded ' + STRING(count) + ' new Spice kernel files . . .')

;===============================Determine the phase angle from the timestamp====================================
    cspice_utc2et, Timestamp, et
    phaseq        = cspice_phaseq(et, Body, 'Sun', viewpoint, 'LT')  ; phase angle [rad]
    npix          = fix(radius * 2. + 1.*radius)                     ; number of pixels for generated Mercury image
    rpd           = double(!dpi/180.d0)                              ; radian per degree
    Body_arr      = dblarr(npix,npix)                                ; create Body image array
    Lat_arr       = make_array(npix,npix,/Float,value=!values.F_NaN) ; create luminance lat array
    Lon_arr       = make_array(npix,npix,/Float,value=!values.F_NaN) ; create luminance lon array
    phase         = phaseq/rpd                                       ; phase angle [deg];0 is full; 180 is new

;================Find the sub-observer / subsolar lat & long in the body-fixed frame============================
    cspice_bodvrd, Body, 'RADII', 3, radii
    re = radii[0]
    rp = radii[2]
    f = (re-rp) / re
    cspice_subpnt, 'Near point: ellipsoid', body, et, strcompress('IAU_'+body), 'LT', $
                   viewpoint, sopoint, trgepc, srfvec
    ;------------'BEWARE! Positive longitude is EAST on the Moon, Earth & Sun-------------
    ;-------------------'Positive longitude is WEST on all bodies!'-----------------------
    if Body eq 'Moon' then cspice_recpgr, body, sopoint, re, f, solon, solat, soalt else $
                           cspice_recgeo, sopoint, re, f, solon, solat, soalt
    lon_targ = solon * cspice_dpr()          ; sub-observer body longitude (positive west)
    lat_targ = solat * cspice_dpr()          ; sub-observer body latitude
    cspice_subslr, 'Near point: ellipsoid', body, et, strcompress('IAU_'+body), 'LT', $
                   viewpoint, sspoint, trgepc, srfvec
    ;------------'BEWARE! Positive longitude is EAST on the Moon, Earth & Sun-------------
    ;-------------------'Positive longitude is WEST on all bodies!'-----------------------
    if Body eq 'Moon' then cspice_recpgr, body, sspoint, re, f, sslon, sslat, ssalt else $
                           cspice_recgeo, sspoint, re, f, sslon, sslat, ssalt
    lon_sun = sslon * cspice_dpr()           ; sub-solar body longitude (positive west)
    lat_sun = sslat * cspice_dpr()           ; sub-solar body latitude
    
;==========================Determine the Solar Irradiance at the Body======================================

    ; Find the Body-Earth light time & Velocity 
    ; TBD HACK!!!
;      cspice_subpnt, 'Near point: ellipsoid', 'Jupiter', ET, 'IAU_Jupiter', 'LT+S', viewpoint, Sub_Earth, trgepc, srfvec ;get sub observer point in catesian IAU Jupiter coords
;      cspice_bodvrd, 'Jupiter', 'RADII', 3, radii ;get Jupiter shape constants
;      cspice_spkpos, 'Jupiter', ET, 'J2000', 'LT+S', viewpoint, ptarg, ltime
;      R_j = 206264.806* atan(radii[0] / norm(ptarg)) 
    
    READCOL, SOLAR_SPECTRUM_FILE, F='A,A', WL_nm, flux, STRINGSKIP = '#', /Silent ;flux is in W/m2/nm
    start = where(WL_nm eq '299.100')
    WL_nm = float(WL_nm[start:*]) 
    flux  = float(flux[start:*])

    ; change flux units from W/m^2/nm to photons / (cm^2 s A)
    ; multiply by ((lambda / hc) / 1 W)  * (1 m^2 / 1e4 cm^2) * (1 nm / 10 A)    
    conversion = ((WL_nm*1.e-9)/(6.62606957e-34*299792458.D)) * (1./1.e4) * (1./10.)
    flux = flux * conversion                     ; Cross-checked this result against Huebner et al. (1992). 
    WL_A = temporary(WL_nm) * 10.                ; Wavelength from nm into angstroms 

    CSPICE_SPKEZR, body, et, 'J2000', 'LT+S', 'SUN', SUN_BODY, ltime
    helio_dist = norm(SUN_BODY[0:2]) / 149597870.7 
    flux = flux / helio_dist^2                   ; 'Flux' is in units of photons / (cm^2 s A)
    flux_at_wavelength = interpol(flux, WL_A, wavelength)
    
    
;    VACTOAIR, WL_A, WL_A_Air                    ; Vacuum to air wavelength conversion        
;    WL_A = temporary(WL_A_Air)    

;============================================Compute Reflectivity at each Pixel=======================================
for x = -radius, radius do begin                 ; loop through x positions
  for y = -radius, radius do begin               ; loop through y positions
    if (sqrt(double(x)^2+double(y)^2) le radius) then begin   ; if we're on the disk then continue
      ; Luminance lon and lat are the angular distance from central meridian and sub-observer point
      lat     = asin(double(y)/radius)           ; Luminance lat, radians
      lon     = asin(double(x)/(radius*cos(lat))); Luminance lon, radians
      cose    = cos(lon) * cos(lat)              ; cosine of the reflected angle    (passed to HAPKE function)
	   	cosi    = cos(lon + phase*rpd) * cos(lat)  ; cosine of the incident angle     (passed to HAPKE function)
	   	g       = phase*rpd                        ; the phase angle in radians       (passed to HAPKE function)
      RHAPKE  = HAPKE_CS(cose, cosi, g, $        ; Find the reflectance
      ;params  = 'Hapke_1984', wavelength=wavelength)          
      params  = 'Domingue_et_al_2016', wavelength=wavelength)
      body_arr[x+npix/2, y+npix/2]    = RHAPKE   ; Assign HAPKE returned value to Mercury I/F image array
      lat_arr[x+npix/2, y+npix/2]     = lat/rpd  ; log the luminance latitude 
      lon_arr[x+npix/2, y+npix/2]     = lon/rpd  ; log the luminance longitude
    endif
  endfor ;y
endfor ;x

; The above loop gives a NaN result at X = 0 & Y = 0, fix that
  junk = where(finite(body_arr), count, complement = badpix)
  if N_elements(badpix) eq 1 and (badpix ne -1) then begin
    here = array_indices(body_arr, badpix)
    body_arr[badpix] = mean(body_arr[here[0]-1:here[0]-1,here[1]-1:here[1]-1], /nan) 
  endif

; Write the output frames
flux_at_wavelength = flux_at_wavelength ; Convert to Rayleighs / Angstrom units
  R_per_A   = flux_at_wavelength*4.*body_arr /1.e6   ; Continuum Emission in R/A: HAPKE_CS.pro returns body_arr as pi*R_R. See Eqn 4 in Domingue, Sprague and Hunten (1996) Icarus
  MR_per_A  = R_per_A / 1.e6                          ; MR/A is often a more convenient unit 
  I_over_F  = body_arr                                ; Hapke Reflectivity pi*R_H

display_outputs = 1
if keyword_set(display_outputs) then begin        ; Plot things up
  DEVICE, GET_SCREEN_SIZE = ss                    ; Find screen size
  box           = npix > 200 < ss[1]/4            ; Still plot reasonbly-sized windows if npix is too small (but not too big!)
 
 ;=======================================PLOT AN IMAGE IN MR/A==============================================
  print, body+' at '+timestamp+' UTC has a phase angle of ', strcompress(phaseq/rpd)
  cspice_et2utc, et, 'C', 0, label_date  
  window, 0, xs=box*2, ys=box*2, xpos=ss[0]/2, ypos=box+40 
    cgimage, bytscl(MR_per_A, /nan, min=min(MR_per_A), max=max(MR_per_A)), ct=20
    cgcolorbar, ctindex=20, range=[min(MR_per_A),max(MR_per_A)], position=[0.15,0.84,0.9,0.87], font=0, $
      charthick = 1.6, charsize = 1.5, /top, title = strcompress(Body+ ' - '+label_date+' @ '+$
      string(fix(wavelength))+cgsymbol('Angstrom')+' (MR / '+cgsymbol('Angstrom')+')')

  ;================================PLOT THE I/F AS A FUNCTION OF LUMINANCE LONGITUDE=========================
  window, 1, xs=box*2, ys=box*2, xpos=ss[0]/2, ypos=0
    cgplot, lon_arr[npix/2.-radius:npix/2.+radius,npix/2+1], body_arr[npix/2.-radius:npix/2.+radius,npix/2+1], $
      xtitle = 'Luminance Longitude (deg)', ytitle = 'Reflectance '+cgsymbol('pi')+'r!DH!N', $ 
      title='phase = ' + strcompress(phase) + ' deg', xr = [-90.,90], psym = 4, xtickinterval = 30.
  
  ;=========PLOT THE I/F AS A FUNCTION OF LUMINANCE LATITUDE ALONG A GIVEN LONTITUDE=========================
  window, 2, xs=box*2, ys=box*2, xpos=0, ypos=0
    junk = min(abs(lon_arr - Longitude_slice), Fixed_longit_indicies, dimension = 1, /Nan) ; get the indices along a fixed longitude
    cgplot, lat_arr[Fixed_longit_indicies], body_arr[Fixed_longit_indicies], xtitle = 'Luminance Latitude (deg)',$
      title='phase = ' + strcompress(phase) + ' deg', xr = [0.,90], psym = 4, xtickinterval = 30.    
endif 
make_picture = 1
if keyword_set(make_picture) then begin
  ; Read in a body albedo map, and find parameter for longitude shift
    read_jpeg, inpdir + mosaic_file, m
    sm       = size(m, /dimensions)
    lons     = cgscalevector(indgen(sm[0]), -180., 180.); center of downloaded maps is 0 lon, by default
    lats     = cgscalevector(indgen(sm[1]), -90., 90.)  ; center of downloaded maps is 0 lat, by default
    if lon_targ gt 180. then lon_targ = lon_targ - 360. ; SPICE gives, 0 to 360, scale this -180 to 180 
    if lon_sun gt 180. then lon_sun = lon_sun - 360.    ; SPICE gives, 0 to 360, scale this -180 to 180 
    junk     = min(abs(lons), xlon_0)                   ; xlon_0 = approximate pixel location of zero longitude (center of downloaded maps by default)
    junk     = min(abs(lons - lon_targ), xlon_targ)     ; xlon_targ = approximate pixel location of target longitude
    junk     = min(abs(lons - lon_sun), xlon_sun)       ; xlon_sun = approximate pixel location of sub-solar longitude 
    junk     = min(abs(lats), xlat_0)                   ; xlat_0 = approximate pixel location of zero latitude (center of downloaded maps by default)
    junk     = min(abs(lats - lat_targ), xlat_targ)     ; xlat_targ = approximate pixel location of target latitude
    junk     = min(abs(lats - lat_sun), xlat_sun)       ; xlat_sun = approximate pixel location of sub-solar latitude
    
  ;===========================Plot the Albedo Map and mark the sub-observer lcation of the viewpoint==============================
  DEVICE, GET_SCREEN_SIZE = ss                               ; find screen size
  box      = sm[0] > 300 < ss[1]*.8 < ss[0]/2.               ; still plot reasonably-sized windows if npix is too small (but not too big!)
  window, 0, xs=box, ys=box/2, xpos=ss[0]/2, title = Body
    axis_format = {XTicks:12, Xminor:3, YTicks:6, Yminor:3}
    cgimage, bytscl(m,min=1,max=255), /axes, xr=[-180,180], yr=[-90,90], xtit='Longitude (deg)', ytit='Latitude (deg)', $
      title=mosaic_file, color='black', axkeywords=axis_format
    cgoplot, !X.Crange, [1,1]*lats[xlat_targ], col='red', line=2, thick=2
    cgoplot, [1,1]*lons[xlon_targ], !Y.Crange, col='red', line=2, thick=2
  
  ; Next, we want to keep north up but rotate the Hapke Model to the true illumination angle, 
  ; Use indices in a dummy array, post rotation & projection to find this rotation angle  
    Dummy = m 
    Dummy[xlon_targ-10:xlon_targ+10, xlat_targ-10:xlat_targ+10] = 255
    Dummy[xlon_sun-10:xlon_sun+10, xlat_sun-10:xlat_sun+10] = 255
  
  ; Generate image for pretty picture
    body_map = image(dummy, LIMIT=[-90,-180,90,180], GRID_UNITS=2, IMAGE_LOCATION=[-180, -90], COLOR='black', $
            IMAGE_DIMENSIONS=[360,180], MAP_PROJECTION='Orthographic', RGB_TABLE=0, LABEL_SHOW=0, $
            GRID_LATITUDE=180, GRID_LONGITUDE=360, DIMENSIONS=[1000,1000], BACKGROUND_COLOR='black', $
            CENTER_LATITUDE = lat_targ, CENTER_LONGITUDE = lon_targ, Location = [xlon_sun, xlat_sun], $
            SPHERE_RADIUS = re, /buffer)
        
  ; Find the rotation matrix to put sopoint at [0,0,r] and the body's north pole in the up direction
    cspice_twovec, sopoint, 3, [0.,0.,re], 2, mout
    cspice_mxv, mout, sopoint, projected_sopoint
    cspice_mxv, mout, sspoint, projected_sspoint      
    ss = projected_sspoint[0:1]
    if not finite(ss[0]) then print, 'Geometry error!' 
  
    ;body_map.rotate, -!radeg*atan(ss[1]/ss[0]) ; Rotate to align sub-solar with x
    body_map.save, inpdir + 'Mercury test.jpg', border =0, resolution=300
    body_map.close
  read_jpeg, inpdir + 'Mercury test.jpg', body_image, /grayscale 
  body_image       = congrid(float(body_image), 2.*radius, 2.*radius, /interp)
  angle            = cspice_vsep([-1.d,0.d,0.d], [ss,0.d])
  cropped_I_over_F = I_over_F[npix/2-radius:npix/2+radius-1, npix/2-radius:npix/2+radius-1] ;Hack, not sure about minus 1s!
  rotated_Hapke    = rot(cropped_I_over_F, angle*!radeg, MISSING = 0., /interp)
  MR_per_A         = rot(MR_per_A, angle*!radeg, MISSING = 0., /interp)
  R_per_A          = rot(R_per_A, angle*!radeg, MISSING = 0., /interp)
  
  window, 3, xs = 2*radius, ys = 2*radius
  tv, bytscl(rotated_Hapke)
  
  pretty           = rotated_Hapke * body_image
  Make_picture     = pretty/max(pretty)
  
  ;-----------Option: alignment to make celestial north up, default images have body's north pole up-------------------------------- 
  if keyword_set(align_celestial_north) then begin
    ; Find the position angle of the north pole vector, that is, the separation between Body's north pole and celestial north pole.
      CSPICE_SPKEZR, body, et, 'J2000', 'LT', viewpoint, BODY_state, ltime
      CSPICE_RECRAD, BODY_state[0:2], dist, ra, dec ;Convert rectangular coordinates to RA and Dec
    
    ; Find a vector from body's center to it's north pole in the J2000 frame at the ephemeris time, NEGLECT the planet's oblateness
      North_body_fixed = [0.,0.,1.]*rp ;location of the North Pole in body-fixed coords
      cspice_pxform, 'IAU_'+Body, 'J2000', et - ltime, To_J2000 ; Find body-fixed coords to J2000 rotation matrix
      cspice_mxv, To_J2000, North_body_fixed, North_J2000 ; Rotate to J2000
    
    ; For comparison with JPL Horizons, get the RA and DEC of the pole direction 
      ;CSPICE_RECRAD, North_J2000, dist, North_dir_ra, North_dir_dec ;Convert rectangular coordinates to RA and Dec
    
    ; Get the RA and Dec of the the body's north pole
      Pole_state = body_state + North_J2000
      CSPICE_RECRAD, Pole_state, dist, Pole_ra, Pole_dec ;Convert rectangular coordinates to RA and Dec
  
    ; Precess each RA and Dec to the Current Epoch. First we'll find the year and fraction of a year. we're interested in:
      Current_epoch = 2000. + et / 3.1556926e7 ;J2000 + seconds past J2000 / seconds per year 
      precess, Pole_ra, Pole_DEC, 2000, Current_epoch, /RADIAN ;Precessed Ra and dec is now applied to Pole location 
      precess, ra, DEC, 2000, Current_epoch, /RADIAN ;Precessed Ra and dec is now applied to body center location 
      Delta_RA      = Pole_ra - ra
      Delta_Dec     = Pole_Dec - dec
      theta         = sqrt((Delta_RA*cos(Dec))^2.D + Delta_Dec^2.D)
      PA            = signum(Delta_RA)*Acos(Delta_Dec / theta)
      print, 'Rotating by North Pole Position Angle (CCW, E of N) = ', PA / cspice_rpd()
      Make_picture  = rot(Make_picture, -PA / cspice_rpd(), /interp, missing = 0.) ; IDL rot is Clockwise ---> -PA
      MR_per_A      = rot(MR_per_A, -PA / cspice_rpd(), /interp, missing = 0.)     ; IDL rot is Clockwise ---> -PA
      R_per_A       = rot(R_per_A, -PA / cspice_rpd(), /interp, missing = 0.)      ; IDL rot is Clockwise ---> -PA
  endif
window, 4, xs = 2*radius, ys = 2*radius
tv, bytscl(Make_picture)
endif ;make_picture

; ==========================Determine the output image platescales in Arcsec / Pixel================================
  cspice_bodvrd, Body, 'RADII', 3, radii                            ; Get Body shape constants
  cspice_spkpos, Body, ET, 'J2000', 'LT+S', viewpoint, ptarg, ltime
  R_M = 206264.806 * atan(radii[0] / norm(ptarg))                     ; Radius of the body in arcsec
  Hapke_Platescale = R_M / radius
  if body eq 'Mercury' then print, 'Angular diameter of ', body, ' =', 2.*R_M, ' arcsec'
  if body eq 'Moon' then print, 'Angular diameter of ', body, ' =', 2.*R_M / 60., ' arcmin'

; Blur it to the desired PSF, assume Gausssian for now, sigma is in arcsec
;  PSF_sigma      = .25 ;arcseconds
;  angular_radius = 206265.*atan(re/norm(srfvec))
;  blur_function  = GAUSSIAN_FUNCTION([radius*(PSF_sigma/angular_radius),radius*(PSF_sigma/angular_radius)], WIDTH=radius, /normalize)
;  blur_function  = blur_function+.5*GAUSSIAN_FUNCTION([10*radius*(PSF_sigma/angular_radius),10*radius*(PSF_sigma/angular_radius)], WIDTH=radius, /normalize)
  cgPS_open, filename = strcompress(OUTPUT_IMAGE), /ENCAPSULATED
   Log_Display = 1
   !P.font=1
   device, SET_FONT = 'Helvetica Bold', /TT_FONT
  ;Get the convolution kernal PSF from a RIPS frame of a star.
;  star = mrdfits('C:\Users\schmidt\Desktop\RIPS\2018 - AEOS with RIPS\June 25 HST\RIPS_setup_217.fits', 0, header, /unsigned)
;  cropped_star  = star[619-49:619+50, 218-49:218+50] 
;  cropped_star  = congrid(cropped_star, 150, 150) 
  ;cropped_star  = cropped_star - min(cropped_star)
;  cropped_star  = cropped_star - 100.
;  blur_function = cropped_star / total(cropped_star)
  ;Blurred       = CONVOL(MR_per_A, blur_function, /edge_zero, /normalize)
  ;Blurred       = rot(Blurred, 7.06+180.)
  ;window, 1, xs = 3*radius, ys = 3*radius
  axis_format = {xticks:6, xtickname:['-1.5','-1','-.5','0','.5','1','1.5'], $
                 yticks:6, ytickname:['-1.5','-1','-.5','0','.5','1','1.5'], $
                 xtitle:body+' Radii',ytitle:body+' Radii'}
  ;cgimage, blurred, /scale, /axes, axkeywords = axis_format, title = body + ' - '+timestamp +' MR / A'
  ;Make_picture = rot(Make_picture, 204.4, /interp, missing = 0.0) ; Perfect for First RIPS Moon Data
  
  cgimage, Make_picture, /scale, /axes, axkeywords = axis_format, title = body + ' - '+timestamp +' MR / '+cgsymbol('Angstrom')
  
  make_contour = rotate(MR_per_A[npix/2-radius:npix/2+radius-1, npix/2-radius:npix/2+radius-1], 2)
  ;make_contour = rot(make_contour, 204.4, /interp, missing = 0.0) ; Perfect for First RIPS Moon Data

  ;levels = [indgen(max(blurred)/step)*step + step] 
  
  If Body eq 'Mercury' then step = 5. ;step for the contour in MR/A in the output image
  If Body eq 'Mercury' then format_Code = '(F4.1)'
  If Body eq 'Moon' then step = 0.5 ;step for the contour in MR/A in the output image
  If Body eq 'Moon' then format_Code = '(F3.1)'
  
  levels = [indgen(max(MR_per_A)/step)*step + step]
  ;cgcontour, blurred, /onimage, color='red', levels=levels, xthick=3, C_CHARSIZE=2.;, $
  cgcontour, make_contour, /onimage, color='green', levels=levels, xthick=3, C_CHARSIZE=2.4, $
  C_annotation = string(levels, Format = Format_code)+' MR/'+cgsymbol('Angstrom'), C_LABELS=[1,1,1,1,1]
cgPS_close
;--------------------------------------

end