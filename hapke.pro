FUNCTION HAPKE, AMU, AMUNOT, GIR

; **************************************************************************************************************************
; **************************************************************************************************************************
; HAPKE function based on code from Rosemary Killen, edited for clarity on 5 Nov 2018 by LM
; 
; Input:
;    amu           - cosine of emission angle, radian
;    amunot        - cosine of incidence angle, radian
;    gir           - phase angle, radian (???? NOTE Jeff, original said "in deg")
;     
; Output:
;    refl          - output Hapke reflectance
;
; Original comments remaining (even if I don't understand them...):
;
;	NEW VERSION 3-1-89 TO MATCH HAPKE '86 AND HELFENSTEIN PARAMS
;	modified 11/13/2008 to match Mallama et al., 2002 by RMK
;
;      implicit double precision(a-h,o-z)
;      double precision munpn,mupnot,munotp,mupr,ff
;	common /stuff/ corr
;      data Hapke_theta_bar,Hapke_h/0.349d0,0.11d0/
;	Mallama et al Icarus 2002
;      data Hapke_theta_bar = 0.279d0, Hapke_h = 0.065d0/
;	Warell, 2004
;      data Hapke_w = 0.15d0, thetab = 0.175d0, Hapke_h = 0.08d0
;	Domingue et al., 2015, 558 nm values
; 	test for Hapke paper w=0.25, Hapke_theta_bar=0.436

;CS: Domingue et al. Icarus 2015 gives:
    ;lambda = 559.2 nm 
    ;w = 0.18509333 
    ;A = 0.16461793 
    ;B = 0.07097318 
    ;theta = 8.6778805
    ;theta_b = 

params = 'Ververka_et_al_1988'

; **************************************************************************************************************************
; **************************************************************************************************************************

amu      = double(amu)                           ; this is the cosine of the emergent (reflected) angle
amunot   = double(amunot)                        ; this is the cosine of the incident angle
g        = abs(GIR)                              ; the phase angle in radians

; ==========================================================================================================================
; Some Hapke variable initializations,  
; ==========================================================================================================================
Case 1 of 
  params eq 'Domingue_et_al_2015': begin ;(Quoting RMK Quoting Domingue et al. 2015---BEWARE!!!)
      Hapke_w         = 0.25d0           ; the Hapke (1984) w, single scattering albedo parameter
      Hapke_h         = 0.08d0           ; the Hapke (1984) h, the angular width of the opposition surge 
      Hapke_theta_bar = 0.436d0          ; the Hapke (1984) theta-bar, an average topographic slope angle of IN RADIANS
      Hapke_b         = 0.1365d0         ; the Hapke (1984) b coefficient for the 1st Lengendre term (0th term is just 1)
      Hapke_c         = 0.08d0           ; the Hapke (1984) c coefficient for the 2nd Lengendre term (0th term is just 1)
      Hapke_B_0       = 2.9d0            ; the Hapke (1986) B_0, amplitude of the opposition effect
  end  
  params eq 'Ververka_et_al_1988': begin ;(Solution 1 in Table 1 of the Ververka et al 1988 chapter in the Arizona "Mercury" book)
      Hapke_w         = 0.20d0           ; the Hapke (1984) w, single scattering albedo parameter
      Hapke_h         = 0.11d0           ; the Hapke (1984) h, the angular width of the opposition surge 
      Hapke_theta_bar = 0.367d0          ; the Hapke (1984) theta-bar, an average topographic slope angle of IN RADIANS
      Hapke_b         = 0.20d0           ; the Hapke (1984) b coefficient for the 1st Lengendre term (0th term is just 1)
      Hapke_c         = 0.18d0           ; the Hapke (1984) c coefficient for the 2nd Lengendre term (0th term is just 1)
      Hapke_B_0       = 2.4d0            ; the Hapke (1986) B_0, amplitude of the opposition effect
  end    
Endcase

reflh           = 0.
if(amu*amunot le 3.d-4) then return, reflh              ; leave without doing anything for such small valus of amu*amunot (reflh = 0 in this case)

ae       = acos(amu)                             ; this is the emergent (reflected) angle in radians
ai       = acos(amunot)                          ; this is the incident angle in radians
cotai    = 1.d0/tan(ai)                          ; cotangent of ai = 1/tangent(ai)
cotae    = 1.d0/tan(ae)                          ; cotangent of ae = 1/tangent(ae)

; ==========================================================================================================================
; Compute the B function
; ==========================================================================================================================
if (g gt !dpi/2.d0) then Hapke_B_0 = 0.d0             ; Amplitude of the opposition surge isn't needed/defined for crescent phases (RMK true?)
;bg       = Hapke_B_0 / (1.D0+TAN(G/2.D0)/Hapke_h)     ; See domingue et al. (2015) eq 2
bg       = 1.D0 + (Hapke_B_0 / (1.D0+TAN(G/2.D0)/Hapke_h))    ; See Shkurtov et al. (2012)

; ==========================================================================================================================
; COMPUTE THE ANGLES
; ==========================================================================================================================
tanthb   = tan(Hapke_theta_bar)
cotthb   = 1.d0/tanthb
t1       = 1.d0/sqrt(1.d0+!dpi*tanthb^2)
t2       = sin(ai)*tanthb
if(amunot eq 1.d0) then begin
  t3     = 0.d0
  t4     = 2.d0
endif else begin
  t3     = exp(-(cotthb^2)/((tan(ai)^2)*!dpi))
  t4     = 2.d0-exp(-2.d0*cotthb/(tan(ai)*!dpi))
endelse

MUNPN    = t1*(amunot+(t2*t3)/t4)
munpn    = abs(munpn)
t22      = sin(ae)*tanthb
if(amu eq 1.d0) then begin
  t32    = 0.d0
  t42    = 2.d0
endif else begin
  t32    = exp(-(cotthb^2)/((tan(ae)^2)*!dpi))
  t42    = 2.d0-exp(-2.d0*cotthb/(tan(ae)*!dpi))
endelse
mupnot   = abs(t1*(amu+(t22*t32)/t42))
if(ai eq 0) or (ae eq 0) then psi = 0.d0         ; azimuthal angle "psi" between incident and emergent planes
cpsi     = (cos(g)-amu*amunot)/(sin(ai)*sin(ae)) ; cpsi is cosine of "psi"
if(abs(cpsi) gt 1.d0) then psi = 0.d0 $
                      else psi = acos(cpsi)

if (AI GE AE) then begin
  af1    = exp(-(cotthb^2)*(cotai^2)/!dpi)
  af2    = (sin(psi/2.d0)^2)*exp(-(cotthb^2)*cotae^2/!dpi)
  af3    = 2.d0-exp(-2.d0*cotthb*cotai/!dpi)
  af4    = psi*exp(-2.d0*cotthb*cotae/!dpi)/!dpi
  af     = (af1-af2)/(af3-af4)
  munotp = t1*(amunot+sin(ai)*tanthb*af)
  af1    = cos(psi)*af1
  af     = (af1-af2)/(af3-af4)
  mupr   = t1*(amu+sin(ae)*tanthb*af)
ENDIF ELSE begin
  af1    = exp(-(cotthb^2)*(cotae^2)/!dpi)
  af2    = (sin(psi/2.d0)^2)*exp(-(cotthb^2)*cotai^2/!dpi)
  af3    = 2.d0-exp(-2.d0*cotthb*cotae/!dpi)
  af4    = psi*exp(-2.d0*cotthb*cotai/!dpi)/!dpi
  af     = (af1-af2)/(af3-af4)
  MUPR   = t1*(amu+sin(ae)*tanthb*af)
  af1    = cos(psi)*af1
  af     = (af1-af2)/(af3-af4)
  MUNOTP = t1*(amunot+sin(ai)*tanthb*af)
ENDelse
mupr     = abs(mupr)
munotp   = abs(munotp)

; ==========================================================================================================================
; Compute ff
; ==========================================================================================================================
ff       = exp(-2.d0*tan(psi/2.d0))

; ==========================================================================================================================
; Compute pg: Using the first 3 Legendre Polynomials, values from Domingue et al., 2015
; ==========================================================================================================================
p0       = 1.d 
p1       = cos(g)
p2       = (-1.d0 +3.d0*cos(g)^2) / 2.d0
pg       = P0 + Hapke_b*p1 + Hapke_c*p2

; ==========================================================================================================================
; Compute the H functions
; ==========================================================================================================================
den      = 1.d0+2.d0*munotp*sqrt(1.d0-Hapke_w)
hmunot   = (1.d0+2.d0*munotp)/den
den      = 1.d0+2.d0*mupr*sqrt(1.d0-Hapke_w)
hmu      = (1.d0+2.d0*mupr)/den
if(hmu lt 0.d0)then hmu = 1.d0
if(hmunot lt 0.d0)then hmunot = 1.d0
if(ai lt ae) then begin
  ang    = amunot/munpn
endif else begin
  ang    = amu/mupnot
  ang2   = amunot/(amu+amunot)
endelse
root     = sqrt(1.d0+!dpi*(tan(Hapke_theta_bar)^2))
geom1    = munotp/(munotp+mupr)
geom2    = mupr*amunot/(mupnot*munpn)
geom     = geom1*geom2
geom3    = amunot /(amu + amunot)

den      = root*(1.d0-ff+(ff*ang*root))
reflh    = Hapke_w*((1.d0+bg)*pg-1.d0+hmu*hmunot)/!dpi/4.d0
reflh    = reflh * geom3

return, REFLH

end
