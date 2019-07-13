
; ==========================================================================================================================
; Compute the scattering phase functions, P(g). For isotropic scattering P(g) = 1
; ==========================================================================================================================
FUNCTION P, x                     
  ; Here x is cos(g) = cos(phase angle)
  COMMON HAPKE_PARAMETERS, Phase_func, Hapke_w, Hapke_h, Hapke_theta_bar, Hapke_b, Hapke_c, Hapke_B_0
  b = Hapke_b & c = Hapke_c
    
  ; ========================================Using the first 3 Legendre Polynomials==========================================
  if Phase_func eq 'Legendre' then begin
    P         = LEGENDRE( x, [0,1,2], /DOUBLE )      ; Hapke, 1984 formulation
    P_of_x    = P[0] + b*P[1] + c*P[2]               ; Validated against Mercury figures in Hapke 1984 or 1986 (I forget)
  endif 
  
  ; =============================Using the Double-Lobbed Scattering Henyey-Greenstein Function==============================
  if Phase_func eq 'Henyey-Greenstein?' then begin
    P_of_x    = ((1.d0+c) / 2.d0) * (1.d0 - b^2) /  ( ((1.d0 - 2.d0*b*x + b^2)^1.5d) ) + $ ; Hapke (2012) formulation,
                ((1.d0-c) / 2.d0) * (1.d0 - b^2) /  ( ((1.d0 + 2.d0*b*x + b^2)^1.5d) )     ; *NOT YET VALIDATED*
  endif
  if Phase_func eq 'Henyey-Greenstein' then begin
    P_of_x    = (1.d0-c) * (1.d0 - b^2) /  ( ((1.d0 - 2.d0*b*x + b^2)^1.5d) ) + $    ; Domingue et al. (2016)
                     (c) * (1.d0 - b^2) /  ( ((1.d0 + 2.d0*b*x + b^2)^1.5d) )        ; Eqn. 3, validated against her Fig. 24
  endif
  return, P_of_x
END
; ==========================================================================================================================
; Compute the H functions, see Hapke (2012 Book) eqn. 8.70b 
; ==========================================================================================================================
FUNCTION H, x ;Approximation to Chandrasehkar's H function, See Hapke (2012) eq 8.70b
  COMMON HAPKE_PARAMETERS, Phase_func, Hapke_w, Hapke_h, Hapke_theta_bar, Hapke_b, Hapke_c, Hapke_B_0
  gamma     = sqrt(1.d0 - Hapke_w) 
  H_of_x    = (1.d0 + (2.d0*x)) / (1.d0 + 2.d0*gamma*x)   
  return, H_of_x
END
; ==========================================================================================================================
; Compute the backscatter shadow hiding opposition effect SHOE B_S(g) See Hapke (2012) 9.21 and Fig. 9.8 
; Assume h = h_S = h_C (not sure the validity of this, but most papers provide only the one "h" Hapke parameter)
; Assume B_0 = B_0S = B_0C (not sure the validity of this, but most papers provide only the one "B_0" Hapke parameter)
; ==========================================================================================================================
FUNCTION B_S, x 
  COMMON HAPKE_PARAMETERS, Phase_func, Hapke_w, Hapke_h, Hapke_theta_bar, Hapke_b, Hapke_c, Hapke_B_0 
  ; Amplitude of the opposition surge isn't needed/defined for crescent phases
  if (abs(x) lt !dpi/2.d0 and Hapke_B_0 ne 0.D0) then begin ; See Hapke (1984) Eq 54 
      y   = tan(x / 2.D0) / Hapke_h                         ; See Hapke (2012 Book) Eq 9.21b 
      ;B_s = sqrt(4.D0*!dpi/y) * exp(1.D0 / y) * (ERF(sqrt(4.D0 / y)) - ERF(sqrt(1.D0 / y))) + exp(-3.D0 / y) - 1.D0 ; Hapke (2012 Book) Eq 9.21b 
      B_S_of_x = (1.D0 + y)^(-1.d0) ; Approx solution, Hapke (2012 Book) Eq. 9.22 & Fig. 9.8, also Domingue et al. (2015) Eqn. 2      
  endif else B_S_of_x = 0.D0  ; Per Hapke (1984) Eq 54 
  return, B_S_of_x 
END

; ==========================================================================================================================
; Compute the coherent backscatter opposition effect CBOE B_C(g) See Hapke (2012) 9.43 and Fig. 9.12 
; Assume h = h_S = h_C (not sure the validity of this, but most papers provide only the one "h" Hapke parameter)
; Assume B_0 = B_0S = B_0C (not sure the validity of this, but most papers provide only the one "B_0" Hapke parameter)
; ==========================================================================================================================
FUNCTION B_C, x 
  COMMON HAPKE_PARAMETERS, Phase_func, Hapke_w, Hapke_h, Hapke_theta_bar, Hapke_b, Hapke_c, Hapke_B_0
  K     = 1.d0
  y     = tan(x / 2.d0) / Hapke_h 
  B_c_of_x    = (1.d0 / (1.42d0+K)) * (1.d0 / (1.d0+y^2)) * ( 1.d0 + ((1.d0 - exp(-1.42d*K*y))/y) )    ; Eqn. 9.41a
  ;B_c_of_x    = (1.d0 + (1.3d0+K)*(y+y^2) )^(-1.d0)   ;Approx solution Eqn. 9.43
  return, B_c_of_x 
END

; **************************************************************************************************************************
; **************************************************************************************************************************
FUNCTION HAPKE_CS, mu, mu_0, g, Params=Params, Wavelength=Wavelength
  COMMON HAPKE_PARAMETERS, Phase_func, Hapke_w, Hapke_h, Hapke_theta_bar, Hapke_b, Hapke_c, Hapke_B_0 

; INPUTS:
;    mu           - cosine of emission angle, radians
;    mu_0         - cosine of incidence angle, radians
;    g            - phase angle, radians 
;    params       - which parameters to use:
;                  'Hapke_test' = Placeholder for testing Hapke (1984)
;                  'Hapke_1984' = Mercury --- Verified against Hapke (1984) comparisons to Mariner 10.
;                  'Ververka_et_al_1988' = Mercury --- NOT WORKING, unable match to Veverka (1988)
;                  'Domingue_et_al_2016' = Mercury --- Wavelength dependent "Hapke Basic" parameters used in her Table 4.
; OPTIONAL INPUTS:
;    Wavelength   - Wavelength in ANGSTROMS; required whenever parameters are wavelength dependent
   
; OUTPUTS:
;    refl         - Hapke absolute relectivity: !pi * r_H 


; **************************************************************************************************************************
; **************************************************************************************************************************

mu      = double(mu)                          ; this is the cosine of the emergent (reflected) angle
mu_0    = double(mu_0)                        ; this is the cosine of the incident angle
g       = abs(g)                              ; the phase angle in radians

; ==========================================================================================================================
; Some Hapke Parameters from Various Publications
; ==========================================================================================================================
Case 1 of 
;  params eq 'Domingue_et_al_2015': begin ;(Quoting RMK Quoting Domingue et al. 2015---BEWARE!!!)
;      Hapke_w         = 0.25d0           ; the Hapke (1984) w, single scattering albedo parameter
;      Hapke_h         = 0.08d0           ; the Hapke (1984) h, the angular width of the opposition surge 
;      Hapke_theta_bar = 0.436d0          ; the Hapke (1984) theta-bar, an average topographic slope angle of IN RADIANS
;      Hapke_b         = 0.1365d0         ; the Hapke (1984) b coefficient for the H-G
;      Hapke_c         = 0.08d0           ; the Hapke (1984) c coefficient for the H-G
;      Hapke_B_0       = 2.9d0            ; the Hapke (1986) B_0, amplitude of the opposition effect
;  end  
;  params eq 'Ververka_et_al_1988': begin ;(Solution 1 in Table 1 of the Ververka et al 1988 chapter in the Arizona "Mercury" book)
;      Hapke_w         = 0.20d0           ; the Hapke (1984) w, single scattering albedo parameter
;      Hapke_h         = 0.11d0           ; the Hapke (1984) h, the angular width of the opposition surge 
;      Hapke_theta_bar = 0.367d0          ; the Hapke (1984) theta-bar, an average topographic slope angle of IN RADIANS
;      Hapke_b         = 0.20d0           ; the Hapke (1984) b coefficient for the 1st Lengendre term (0th term is just 1)
;      Hapke_c         = 0.18d0           ; the Hapke (1984) c coefficient for the 2nd Lengendre term (0th term is just 1)
;      Hapke_B_0       = 2.4d0            ; the Hapke (1986) B_0, amplitude of the opposition effect
;  end
  params eq 'Hapke_1984': begin    ;( From Hapke (1984) Figure 5 caption, Fit to Mariner 10 for g = 77 phase angle )
      Phase_func      = 'Legendre' 
      Hapke_w         = 0.20d0           ; the Hapke (1984) w, single scattering albedo parameter
      Hapke_h         = 0.40d0           ; the Hapke (1984) h, the angular width of the opposition surge 
      Hapke_theta_bar = 0.349d0          ; the Hapke (1984) theta-bar, an average topographic slope angle of IN RADIANS
      Hapke_b         = 0.579d0          ; the Hapke (1984) b coefficient for the 1st Lengendre term (0th term is just 1)
      Hapke_c         = 0.367d0          ; the Hapke (1984) c coefficient for the 2nd Lengendre term (0th term is just 1)
      Hapke_B_0       = 1.d0             ; Undefined in Hapke (1984), this term comes in in Hapke (1986)
  end 
  params eq 'Hapke_test': begin          ;Test for Hapke (2012) Fig. 12.5, or Hapke (1984) where !pi*r is plotted in Fig. 3
      Phase_func      = 'Legendre'
      Hapke_w         = 0.25d0           ; the Hapke (1984) w = 25, single scattering albedo parameter
      Hapke_h         = 0.00d0           ; the Hapke (1984) h, the angular width of the opposition surge 
      ;Hapke_theta_bar = 0.436d0         ; the Hapke (1984) theta-bar = 25 degrees IN RADIANS
      Hapke_theta_bar = 0.0d0            ; the Hapke (1984) theta-bar = 0 degrees IN RADIANS
      Hapke_b         = 0.0d0            ; the Hapke (1984) b coefficient, isotropically scattering
      Hapke_c         = 0.0d0            ; the Hapke (1984) c coefficient, isotropically scattering 
      Hapke_B_0       = 1.d0             ; the Hapke (1986) B_0, amplitude of the opposition effect
  end 
  params eq 'Domingue_et_al_2016': begin ;See Domingue et al. 2016 Table 4. 
    ; This uses her Eqn. 38, which she terms the "Hapke Basic" formulation, K=1, B_C0 = 0. 
      Phase_func      = 'Henyey-Greenstein'
      Table_4         = double([[433.2, 0.2111, 0.1514, 0.3341, 0.1551, 0.6248, 0.1261, 26.4272, 14.6013],$
                                [479.9, 0.2308, 0.1697, 0.3248, 0.1474, 0.6135, 0.1062, 26.0306, 14.7452],$ 
                                [558.9, 0.2589, 0.1974, 0.3147, 0.1365, 0.6025, 0.0818, 25.7300, 14.7801],$
                                [628.8, 0.2809, 0.2187, 0.3094, 0.1292, 0.5983, 0.0704, 25.7043, 14.6528],$
                                [748.7, 0.3129, 0.2491, 0.3057, 0.1223, 0.5990, 0.0729, 25.9386, 14.2707],$
                                [828.4, 0.3299, 0.2655, 0.3051, 0.1220, 0.6021, 0.0902, 26.1120, 14.0295],$
                                [898.8, 0.3417, 0.2778, 0.3045, 0.1248, 0.6046, 0.1160, 26.1607, 13.9099],$
                                [996.2, 0.3542, 0.2921, 0.3018, 0.1339, 0.6052, 0.1679, 25.8967, 14.0090]])   
      Hapke_w         = interpol(Table_4[2,*], Table_4[0,*], wavelength / 10.)      
      Hapke_h         = 0.090d0           
      Hapke_b         = interpol(Table_4[4,*], Table_4[0,*], wavelength / 10.)  
      Hapke_c         = interpol(Table_4[6,*], Table_4[0,*], wavelength / 10.)  
      Hapke_theta_bar = interpol(Table_4[8,*], Table_4[0,*], wavelength / 10.) / !radeg
      Hapke_B_0      = 3.086d0        
  end    
  params eq 'Sato_et_al_2014': begin ;See Sato et al. 2014 Figure 17 & Domingue et al. 2016 Table 11 
    ; This uses her Eqn. 38, which she terms the "Hapke Basic" formulation, K=1, B_C0 = 0. 
      Phase_func      = 'Henyey-Greenstein'
      ;------------------------lambda,    w,   h_s,     b,    c, B_S0]
      Mare_Table      = double([[321., 0.09, 0.065, 0.232, 0.39, 0.27],$
                                [360., 0.12, 0.061, 0.250, 0.20, 0.26],$ 
                                [415., 0.14, 0.056, 0.254, 0.17, 0.25],$
                                [566., 0.20, 0.053, 0.261, 0.10, 0.21],$
                                [604., 0.22, 0.051, 0.262, 0.09, 0.21],$
                                [643., 0.24, 0.049, 0.264, 0.08, 0.20],$
                                [689., 0.26, 0.049, 0.265, 0.07, 0.19]])  
                                
      Highlands_Table = double([[321., 0.18, 0.080, 0.237, 0.33, 0.24],$
                                [360., 0.21, 0.075, 0.237, 0.33, 0.23],$ 
                                [415., 0.25, 0.081, 0.239, 0.31, 0.20],$
                                [566., 0.36, 0.079, 0.235, 0.35, 0.17],$
                                [604., 0.39, 0.076, 0.233, 0.37, 0.17],$
                                [643., 0.41, 0.074, 0.232, 0.38, 0.17],$
                                [689., 0.44, 0.077, 0.235, 0.35, 0.16]]) 
      
      All_Table       = double([[321., 0.16, 0.076, 0.237, 0.33, 0.25],$ ;integrated +/- 30 deg lat
                                [360., 0.19, 0.071, 0.241, 0.29, 0.23],$ 
                                [415., 0.21, 0.074, 0.243, 0.27, 0.21],$
                                [566., 0.33, 0.072, 0.239, 0.31, 0.18],$
                                [604., 0.36, 0.068, 0.237, 0.33, 0.17],$
                                [643., 0.39, 0.067, 0.238, 0.33, 0.17],$
                                [689., 0.41, 0.069, 0.240, 0.30, 0.16]]) 
                                 
      Hapke_w         = interpol(Table_4[1,*], Table_4[0,*], wavelength / 10.)      
      Hapke_h         = interpol(Table_4[2,*], Table_4[0,*], wavelength / 10.)          
      Hapke_b         = interpol(Table_4[3,*], Table_4[0,*], wavelength / 10.)  
      Hapke_c         = interpol(Table_4[4,*], Table_4[0,*], wavelength / 10.)  
      Hapke_theta_bar = 23.4 / !radeg ;Sato et al. 2014 Table 1
      Hapke_B_S0      = interpol(Table_4[5,*], Table_4[0,*], wavelength / 10.) 
      Hapke_B_C0      = 0.d ;No CBOE      
  end
Endcase

; First some workarounds and explicit definitions
  if (mu*mu_0 le 1.e-8) then return, 0     ; For very small values of mu*mu_0, exit & return 0 reflectance
  if not keyword_set(Hapke_B_S0) then Hapke_B_S0 = Hapke_B_0
  if not keyword_set(Hapke_h_s) then Hapke_h_s = Hapke_h 
  if not keyword_set(Hapke_B_C0) then Hapke_B_C0 = Hapke_B_0
  if not keyword_set(Hapke_h_c) then Hapke_h_c = Hapke_h 
  if not keyword_set(Hapke_K) then Hapke_K = 1. 
  
;--------------------------All the following computes Equation 12.55 in the Hapke (2012) Book-------------------------------

K             = 1.                             ; Absorption coefficient of the medium. Formally appears in 12.55, not needed here.
e             = acos(mu)                       ; This is the emergent (reflected) angle in radians
i             = acos(mu_0)                     ; This is the incident angle in radians
coti          = 1.d0/tan(i)                    ; cotangent of i = 1/tangent(i)
cote          = 1.d0/tan(e)                    ; cotangent of e = 1/tangent(e)
cottheta_bar  = 1.D0/tan(Hapke_theta_bar)
coti2         = coti^2
cote2         = cote^2
cottheta_bar2 = cottheta_bar^2
      
; ==========================================================================================================================
; Some basic definitions are now needed in order to get the mu_0e and mu_e variables used in 12.55
; These new 'mu sub e' variables satisfy Hapke (2012) eqn. 12.10. 
; They're found via 12.46 & 12.47 when i < e & via 12.52 & 12.53 when e < i or e = i
; ==========================================================================================================================

  if (i eq 0.d) or (e eq 0.d) then psi = 0.d0                      ; Azimuthal angle "psi" between incident and emergent planes
  cospsi   = (cos(g) - (cos(i)*cos(e))) / (sin(i)*sin(e))          ; Hapke (1984) eqn 3, cospsi is cosine of "psi"
  if (abs(cospsi) le 1.d0) then psi = acos(cospsi) else psi = 0.d0 ; Psi defined
  f_of_psi = exp(-2.d0*tan(psi/2.d0))                              ; Hapke (2012) Book eqn 12.29 
  chi = 1.d0 / sqrt(1.d0 + !dpi * (tan(Hapke_theta_bar))^2)        ; Hapke (2012) Book eqn 12.45a 
  
  ;Now we can calculate E1 and E2 per eqn. 12.45b & 12.45c First, get all the cotangent terms for convenience

      e1i = exp( (-2.D0/!dpi) * cottheta_bar * coti   )     ; eqn. 12.45b incident
      e2i = exp( (-1.D0/!dpi) * cottheta_bar2 * coti2 )     ; eqn. 12.45c incident
      e1e = exp( (-2.D0/!dpi) * cottheta_bar * cote   )     ; eqn. 12.45b emergent
      e2e = exp( (-1.D0/!dpi) * cottheta_bar2 * cote2 )     ; eqn. 12.45c emergent

  ;Now we have everything needed to define 'mu sub e' variables
      eta_e = chi * ( mu + sin(e) * tan(Hapke_theta_bar) * e2e / ( 2.d0 - e1e ) )    ;eqn. 12.48 
      eta_0e = chi * ( mu_0 + sin(i) * tan(Hapke_theta_bar) * e2i / ( 2.d0 - e1i ) ) ;eqn. 12.49
      
      ; e >= i
      if e ge i then begin
         denom = 2.d0 - e1e - (psi/!dpi)*e1i                        
         mu_0e = chi * ( mu_0 + sin(i) * tan(Hapke_theta_bar) * $ 
                       (( cospsi*e2e + ((sin(psi/2.d0))^2) * e2i ) / denom ))   ; eqn. 12.46
         mu_e = chi * ( mu + sin(e) * tan(Hapke_theta_bar) * $
                      (( e2e - ((sin(psi/2.d0))^2)*e2i ) / denom ))             ; eqn. 12.47
         sfun = (mu_e/eta_e) * (mu_0/eta_0e) * $
                chi / (1.d0 - f_of_psi + f_of_psi*chi*mu_0/eta_0e )             ; eqn. 12.50
      endif

      ; e < i
      if e lt i then begin
         denom = 2.d0 - e1i - (psi/!dpi)*e1e
         mu_0e = chi * ( mu_0 + sin(i) * tan(Hapke_theta_bar) * $
                       (( e2i - ((sin(psi/2.d0))^2)*e2e ) / denom ))            ; eqn. 12.52
         mu_e = chi * ( mu + sin(e) * tan(Hapke_theta_bar) * $
                       (( cospsi*e2i + ((sin(psi/2.d0))^2)*e2e ) / denom ))     ; eqn. 12.53
         sfun = (mu_e/eta_e) * (mu_0/eta_0e) * $
                chi / ( 1.d0 - f_of_psi + f_of_psi*chi*mu/eta_e )               ; eqn. 12.54
      endif

; Compute the Hapke r_R in Equation 12.55
  Lambertian      = (mu_0 / !dpi)        ; mind the w and 4pi
  Lommel_Seeliger = (mu_0 / (mu + mu_0)) ; mind the w and 4pi
  test            = (Hapke_w / (4.d0 *!dpi)) * (mu_0 / (mu_0 + mu)) * H(mu_0) * H(mu)  
  Hapke_1981      = (Hapke_w / (4.d0 *!dpi)) * (mu_0 / (mu_0 + mu)) * ( P(cos(g)) * (1.d0+Hapke_B_0*B_S(g)) + (H(mu_0)*H(mu) - 1.d0) ) 
  Hapke_1986      = (Hapke_w / (4.d0 *!dpi)) * (mu_0e / (mu_0e + mu_e)) * ( P(cos(g)) * (1.d0+Hapke_B_0*B_S(g)) + (H(mu_0e)*H(mu_e) - 1.d0) ) * Sfun
  ;Hapke_2012      = (Hapke_w / (4.d0 *!dpi)) * (mu_0e / (mu_0e + mu_e)) * ( P(cos(g)) * (1.d0+Hapke_B_S0*B_S(g)) + $
  ;                  (H(mu_0e)*H(mu_e) - 1.d0) ) * (1.d0 + Hapke_B_C0*B_C(g)) * Sfun   ;Hapke (2012) Eqn. 12.55
  Basic_Hapke     = (Hapke_w / (4.d0 *!dpi)) * (mu_0e / (mu_0e + mu_e)) * ( P(cos(g)) * (1.d0+Hapke_B_0*B_S(g)) + (H(mu_0e)*H(mu_e) - 1.d0) ) * Sfun  
                      ; Domingue et al. (2016) Eqn. 38
                    
; Note that Hapke (2012) Figure 12.5 gives the reflectance as r_R, but
; Hapke (1984) Figure 3 gives the reflectance as pi*r_R, and that's correct for most practical definitions
; See Domingue Sprague and Hunten (1996) Icarus ---> (4*pi*r_R) * F_solar/10^6 is the brightness in Rayleighs/Angstrom
if Params eq 'Hapke_1984' then return, Hapke_1986 * !dpi 
if Params eq 'Domingue_et_al_2016' then return, Basic_Hapke * !dpi 
end
