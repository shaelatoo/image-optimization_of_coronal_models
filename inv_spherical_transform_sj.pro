
function inv_spherical_transform_sj,B,cp,lmax=lmax, $
  mrange=mrange,period=period,phirange=phirange,cprange= $
  cprange,calc_coeffs=calc_coeffs
  
;+
;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
;
;  Purpose: This routine performs an inverse spherical harmonic
;             transform on a 2-D array.  Heavily based on the 
;             routine with nearly the same name by Marc DeRosa.
;
;  Inputs: B - complex array to be transformed, ordered (l,m)
;          cp - cosine of theta collocation points for theta grid
;
;  Returns: inv_trans - transformed array, ordered (phi,theta)
;
;  Keywords:
;             period = periodicity factor in phi, assumes input array
;                      contains m values which are integral multiples
;                      of period; ignored unless calc_coeffs keyword is set
;             lmax = set to max l value we want to use
;             mrange = set to be range of m values we want to use
;             phirange = (optional) a one- or two-element array containing the
;                       range of phi to return, in radians.  This option is
;                       useful for high-resolution transforms where the region
;                       of interest is bounded in longitude.  If phirange is
;                       one element, range is [0,phirange].  If two elements,
;                       range is [phirange(0),phirange(1)].  If not specified,
;                       the range is set to is [0,2*pi/period]; ignored unless
;                       calc_coeffs keyword is set
;             cprange = (optional) a one- or two-element array containing the
;                       range of cp to return.  This option is useful for
;                       high-resolution transforms where the region of
;                       interest is bounded in latitude.  If cprange is one
;                       element, range is [-cprange,cprange].  If two
;                       elements, range is [min(cprange),max(cprange)].  If
;                       not specified, default range is [-1,1] is used;
;                       ignored unless calc_coeffs keyword is set
;             thindex = (optional) custom array of theta (colatitude)
;                         coordinates, in radians, for output grid.  If not
;                         specified, the argument cp and optional keyword
;                         cprange are used to determine the theta grid;
;                         ignored unless calc_coeffs keyword is set
;             phindex = (optional) custom array of phi (longitude)
;                       coordinates, in radians, for output grid.  If not
;                       specified, an equally spaced grid with
;                       2*n_elements(cp) elements is used;  ignored unless
;                       calc_coeffs keyword is set
;             calc_coeffs - if set, normalization and legendre coefficients
;                       are calculated from scratch based on theta grid; 
;                       otherwise, they are retrieved from the common block
;                       variable coeffs, stored in transform_coefficients.pro,
;                       as are theta and phi grids
;
;  Marc's notes: - All calculations are done in double precision.
;         - Default is to return an array of size
;           (2*n_elements(cp),n_elements(cp)) unless limits are put on cp
;           and/or phi ranges via the cprange and phirange keywords, or a
;           custom grid is specified using thindex and phindex (in which case
;           phirange and cprange are ignored).
;         - Routine is increasingly less accurate for higher l.  To see why,
;           look at table of Legendre functions (m=0 example) in Arfken - 
;           they are alternating series with increasingly larger numbers being
;           added and subtracted from each other.
;
;  Shaela's notes:
;         - Marc's on-the-fly algorithm should be preferred for very low or very
;           high resolution applications; the use of stored values of the mult-
;           iplicative coefficients is more efficient in the mid-resolution 
;           range but can become too memory-intensive as the size of the
;           input array increases
;
;  M.DeRosa - 12 Sep 2000 - created
;             24 Oct 2001 - fixed nasty bug related to sign of m=0 components
;              2 Apr 2003 - added phirange and cprange keywords
;              2 Apr 2007 - now calculates phases of B in a cleaner fashion
;              2 Apr 2007 - added thindex and phindex keywords
;              
;  S.Jones -   1 July 2017 - created
;             10 July 2017 - changed so that phiix and thetaindex
;                             are stored rather than re-calculated
;
;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
;-

; common blocks
@transform_coefficients


;  preliminaries
if KEYWORD_SET(period) then period=ROUND(period(0)) else period=1
if KEYWORD_SET(lmax) then lmax=ROUND(LONG(lmax(0))) $
  else lmax=N_ELEMENTS(B(*,0))-1
case N_ELEMENTS(mrange) of
  0:  mrange=[0,lmax]
  1:  mrange=[0,(ROUND(mrange(0))<lmax)]
  else:  mrange=[ROUND(mrange(0)),(ROUND(mrange(1))<lmax)]
endcase



; check reasonableness of coefficients array
  if lmax gt 513 then begin
    PRINT,'Coefficients array will likely reduce efficiency of this routine.  Consider using original version of this routine, inv_spherical_transform.pro'
    RETURN,-1
  endif


; if coefficient array needs to be calculated, do so
if KEYWORD_SET(calc_coeffs) then begin
 
  ;  determine output (co-)latitude grid
  if N_ELEMENTS(thindex) gt 0 then begin
    if MAX(thindex,min=minthindex) gt !dpi then begin
      PRINT,'  ERROR in inv_spherical_transform: thindex out of range'
      RETURN,-1
    endif else if minthindex lt 0d0 then begin
      print,'  ERROR in inv_spherical_transform: thindex out of range'
      return,-1
    endif 
    ntheta=N_ELEMENTS(thindex)
    costheta=COS(DOUBLE(thindex))
    sintheta=SQRT(1-costheta*costheta)
  endif else begin
    case N_ELEMENTS(cprange) of
      0: begin
        cp1i=-1.0
        cp2i=1.0
        end
      1: begin
        cp2i=ABS(cprange(0))<1
        cp1i=-cp2i
        end
      else: begin
        cp1i=MIN(cprange)>(-1)
        cp2i=MAX(cprange)<1
        end
    endcase
    wh=WHERE((cp ge cp1i) and (cp le cp2i),ntheta)
    if ntheta eq 0 then begin
      PRINT,'  ERROR in inv_spherical_transform: invalid cprange'
      RETURN,-1
    endif 
    costheta=DOUBLE(cp(wh))
    sintheta=SQRT(1-costheta*costheta)
    thindex=ACOS(costheta)
  endelse

;  determine output longitude grid
  if N_ELEMENTS(phindex) gt 0 then begin
    nphi=N_ELEMENTS(phindex)
    phiix=double(phindex)
  endif else begin
    nphi=2*N_ELEMENTS(cp)
    phiix=2*!dpi*DINDGEN(nphi/period)/nphi
    case N_ELEMENTS(phirange) of
      0: begin
        ph1i=0.0
        ph2i=2*!dpi/period
        end
      1: begin
        ph1i=0.0 
        ph2i=(2*!dpi/period)<phirange(0)
        end
      else: begin
        ph1i=0.0>MIN(phirange)
        ph2i=(2*!dpi/period)<MAX(phirange)
        end 
    endcase
    wh=WHERE((phiix ge ph1i) and (phiix le ph2i),nphi)
    if nphi eq 0 then begin
      print,'  ERROR in inv_spherical_transform: invalid phirange'
      return,-1
    endif
    phiix=phiix(wh)
  endelse


  ;  calculate array of amplitudes and phases of B
  Bamp=ABS(B)
  phase=ATAN(B,/phase)


  ;  set up output array inv_trans, coefficient array
  nr=N_ELEMENTS(rix)
  inv_trans=DBLARR(nphi,ntheta)
  coeffs=DBLARR(lmax+1,lmax+1,ntheta)
  

 ;  take care of modes where m=0
  CP_0_0=1/SQRT(4*!dpi)
  coeffs[0,0,*]=CP_0_0
  if mrange[0] eq 0 then begin   ; m's start at 0

    ;  start with m=l=0 mode
    inv_trans=inv_trans+REPLICATE(Bamp(0,0)*COS(phase(0,0))*CP_0_0,[nphi,ntheta])
    
    ;  now do l=1 m=0 mode
    CP_1_0=sqrt(3d0)*costheta*CP_0_0
    Y=replicate(cos(phase(1,0)),nphi) # CP_1_0
    inv_trans=inv_trans+(Bamp(1,0)*Y)
    coeffs[1,0,*]=CP_1_0
    
    ;  do other l modes for which m=0
    if lmax gt 1 then begin
      CP_lm1_0=CP_0_0
      CP_l_0=CP_1_0
      for l=2,lmax do begin
        ld=double(l)
        CP_lm2_0=CP_lm1_0
        CP_lm1_0=CP_l_0
        c1=sqrt(4*ld^2-1)/ld
        c2=sqrt((2*ld+1)/(2*ld-3))*((ld-1)/ld)
        CP_l_0=c1*costheta*CP_lm1_0-c2*CP_lm2_0
        coeffs[l,0,*]=CP_l_0
        Y=replicate(cos(phase(l,0)),nphi) # CP_l_0
        inv_trans=inv_trans+(Bamp(l,0)*Y)
      endfor
    endif
  endif 
  
  ;  loop through m's for m>0, and then loop through l's for each m
  CP_m_m=CP_0_0
  for m=1,mrange(1) do begin
    
    md=double(m)
      
    ;  do l=m mode first
    CP_mm1_mm1=CP_m_m
    CP_m_m=-sqrt(1+1/(2*md))*sintheta*CP_mm1_mm1
    coeffs[m,m,*]=CP_m_m
    if (mrange(0) le m) and ((m mod period) eq 0) then begin
      
      angpart=cos(md*phiix + phase(m,m/period))
      inv_trans=inv_trans+Bamp(m,m/period)*(angpart#CP_m_m)

      ;  now do l=m+1 mode
      if lmax ge m+1 then begin
        CP_mp1_m=sqrt(2*md+3)*costheta*CP_m_m
        coeffs[m+1,m,*]=CP_mp1_m
        angpart=cos(md*phiix+phase(m+1,m/period))
        inv_trans=inv_trans+Bamp(m+1,m/period)*(angpart#CP_mp1_m)

      endif
        
      ;  now do other l's
      if lmax ge m+2 then begin
        CP_lm1_m=CP_m_m
        CP_l_m=CP_mp1_m
        for l=m+2,lmax do begin
          ld=double(l)
          CP_lm2_m=CP_lm1_m
          CP_lm1_m=CP_l_m
          c1=sqrt((4*ld^2-1)/(ld^2-md^2))
          c2=sqrt(((2*ld+1)*((ld-1)^2-md^2))/((2*ld-3)*(ld^2-md^2)))
          CP_l_m=c1*costheta*CP_lm1_m-c2*CP_lm2_m
          coeffs[l,m,*]=CP_l_m
          angpart=cos(md*phiix+phase(l,m/period))
          inv_trans=inv_trans+Bamp(l,m/period)*(angpart#CP_l_m)

        endfor
       endif
     endif
  endfor
  
  ; store theta and phi arrays
  thindex_stored=thindex
  phiix_stored=phiix
  
    
endif else begin       ; if coeffs array populated

  ; check size of coefficient array
  szcoeffs=SIZE(coeffs)
  if szcoeffs[0] ne 3 then begin
    print,'Error.  Coefficient array must be provided via transport_coeffs common variable, or keyword "calc_coeffs" must be set.'
    return,-1
  endif
  
  ;  calculate array of amplitudes and phases of B
  Bamp=ABS(B)
  phase=atan(IMAGINARY(B),DOUBLE(B))

  ;  set up output array inv_trans
  phiix=phiix_stored
  thindex=thindex_stored
  ntheta=N_ELEMENTS(thindex)
  nphi=N_ELEMENTS(phiix)
  inv_trans=DBLARR(nphi,ntheta)

  ;  take care of modes where m=0
  if (mrange[0] eq 0) then begin

    ;  start with m=l=0 mode
    inv_trans=inv_trans+REPLICATE(Bamp[0,0]*COS(phase[0,0])*coeffs[0], $
         [nphi,ntheta])

    ;  now do l=1 m=0 mode
    Y=(REPLICATE(COS(phase[1,0]),nphi)) # REFORM(coeffs[1,0,*])
    inv_trans=inv_trans+(Bamp[1,0]*Y)

    ;  do other l modes for which m=0
    if lmax gt 1 then begin
      for l=2,lmax do begin
        Y=REPLICATE(COS(phase[l,0]),nphi) # REFORM(coeffs[l,0,*])
        inv_trans=inv_trans+(Bamp[l,0]*Y)
      endfor
    endif
  endif

  ;  loop through m's for m>0, and then loop through l's for each m
  for m=1,mrange(1) do begin

    if (mrange(0) le m) and ((m mod period) eq 0) then begin
      md=DOUBLE(m)
      
      ;  do l=m mode first
      angpart=COS(md*phiix + phase(m,m/period))
      inv_trans=inv_trans+Bamp(m,m/period)*(angpart#REFORM(coeffs[m,m,*]))

      ;  now do l=m+1 mode
      if lmax ge m+1 then begin
        angpart=COS(md*phiix+phase(m+1,m/period))
        inv_trans=inv_trans+Bamp(m+1,m/period)*(angpart#REFORM(coeffs[m+1,m,*]))
      endif

      ;  now do other l's
      if lmax ge m+2 then begin
        for l=m+2,lmax do begin
          angpart=COS(md*phiix+phase(l,m/period))
          inv_trans=inv_trans+Bamp(l,m/period)*(angpart#REFORM(coeffs[l,m,*]))
        endfor
      endif
      
    endif 
  endfor          ; loop over m values
  
endelse          ; multiplications only

RETURN,inv_trans

end
