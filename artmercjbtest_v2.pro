pro artmercjbtest_v2

;phase         = 10.D0                            ; phase angle [deg];0 is full; 180 is new
phase         = 105.9d0                            ; phase angle [deg];0 is full; 180 is new
;phase         = 1.D0                            ; phase angle [deg];0 is full; 180 is new
radius        = 150.d0                           ; radius of Mercury in pixels
;npix          = fix(radius * 2. + 0.5*radius)    ; number of pixels for generated Mercury image
npix          = fix(radius * 2. + 1.*radius)    ; number of pixels for generated Mercury image
rpd           = double(!dpi/180.d0)              ; radian per degree
Mercury_arr   = dblarr(npix,npix)                ; create Mercury image array
Lat_arr       = make_array(npix,npix,/Float,value=!values.F_NaN) ; create luminance lat array
Lon_arr       = make_array(npix,npix,/Float,value=!values.F_NaN) ; create luminance lon array

;--------for testing-----------
      Hapke_w         = 0.25d0           ; the Hapke (1984) w = 25, single scattering albedo parameter
      Hapke_h         = 0.00d0           ; the Hapke (1984) h, the angular width of the opposition surge 
      ;Hapke_theta_bar = 0.436d0          ; the Hapke (1984) theta-bar = 25 degrees IN RADIANS
      Hapke_theta_bar = 0.0d0          ; the Hapke (1984) theta-bar = 25 degrees IN RADIANS
      Hapke_b         = 0.0d0            ; the Hapke (1984) b coefficient, isotropically scattering
      Hapke_c         = 0.0d0            ; the Hapke (1984) c coefficient, isotropically scattering 
      Hapke_B_0       = 1.d0             ; the Hapke (1986) B_0, amplitude of the opposition effect
;------------------------------      

for x = -radius, radius do begin                 ; loop through x positions
  for y = -radius, radius do begin               ; loop through y positions
    if (sqrt(double(x)^2+double(y)^2) le radius) then begin   ; if we're on the disk then continue
      ; Luminance lon and lat are the angular distance from central meridian and sub-observer point
      lat     = asin(double(y)/radius)           ; Luminance lat, radians
      lon     = asin(double(x)/(radius*cos(lat))); Luminance lon, radians
      cose    = cos(lon) * cos(lat)              ; cosine of the reflected angle    (passed to HAPKE function)
	   	cosi    = cos(lon + phase*rpd) * cos(lat)  ; cosine of the incident angle     (passed to HAPKE function)
	   	g       = phase*rpd                        ; the phase angle in radians       (passed to HAPKE function)
      RHAPKE  = HAPKE_CS(cose, cosi, g)
      ;RHAPKE  = HAPKE_CS(cose < .9999999d, cosi < .9999999d, g)  ; the <.9998 prevents HAPKE blowing up for values of 1.0
      ;RHAPKE  = bidr2(w,emu,imu,g,holes,p,b0,theta)
      ;RHAPKE  = bidr2(Hapke_w,cose,cosi,g,Hapke_h,[Hapke_b,Hapke_c],Hapke_B_0,Hapke_theta_bar, Pparms = 2)
      Mercury_arr[x+npix/2, y+npix/2] = RHAPKE   ; Assign HAPKE returned value to Mercury I/F image array
      lat_arr[x+npix/2, y+npix/2]     = lat/rpd  ; log the luminance latitude 
      lon_arr[x+npix/2, y+npix/2]     = lon/rpd  ; log the luminance longitude
    endif
  endfor ;y
endfor ;x


; Plot things up
DEVICE, GET_SCREEN_SIZE = ss                     ; find screen size
box           = npix > 200 < ss[1]/4             ; still plot reasonbly-sized windows if npix is too small (but not too big!)
window, 0, xs=box, ys=box, xpos=ss[0]/2, title='0 - Mercury disk'
cgimage, bytscl(Mercury_arr,min=0,max=0.02), /axes, position=[0.,0.,1.,1.], color='black'   ; b/w image of Mercury
cgoplot, !X.Crange, [1,1]*npix/2, color='red'    ; show the line we are plotting
window, 1, xs=box*2, ys=box, xpos=ss[0]/2+box+15, ypos=0, title='1 - Reflectance cut'
cgplot, Mercury_arr(*,npix/2+1), xra=[npix/2-radius,npix/2+radius], psym=10, font=0, color='red', title='phase = ' + strcompress(fix(phase)) + ' deg'
window, 2, xs=box*2, ys=box*2, xpos=ss[0]/2, ypos=box+40, title='2 - Mercury disk (color)'
cgimage, bytscl(Mercury_arr,/nan,min=min(Mercury_arr),max=max(Mercury_arr)*1.), /axes, ct=20, font=0
cgcolorbar, ctindex=20, range=[min(Mercury_arr),max(Mercury_arr)*1.], /top, position=[0.15,0.9,0.9,0.95], font=0

window, 3
cgplot, lon_arr[npix/2.-radius:npix/2.+radius,npix/2+1], Mercury_arr[npix/2.-radius:npix/2.+radius,npix/2+1], xtitle = 'Luminance Longitude (deg)', $ 
  title='phase = ' + strcompress(fix(phase)) + ' deg', xr = [-90.,90], yr = [0,.06], psym = 4, xtickinterval = 30.

;================================Plot the I/F as a function of luminance latitude=========================
  window, 4
  Longitude_slice = -50. ; 0 = central meridian.
  junk = min(abs(lon_arr - Longitude_slice), Fixed_longit_indicies, dimension = 1, /Nan) ; get the indices along a fixed longitude
  cgplot, lat_arr[Fixed_longit_indicies], Mercury_arr[Fixed_longit_indicies], xtitle = 'Luminance Latitude (deg)',$
    title='phase = ' + strcompress(fix(phase)) + ' deg', xr = [0.,90], yr = [0,.047], psym = 4, xtickinterval = 30.


stop
save, Mercury_arr, filename = 'Mercury_arr_cos.sav'
end