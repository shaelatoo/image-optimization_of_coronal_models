function optimization_inv_transform,transform,altitude, $
         component,lmaxi=lmaxi,recursive=recursive

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                                                             ;
; Purpose: Re-calculates the inverse spherical harmonic trans-;
;            form.  For use with pfss_potl_field_sj, where    ;
;            some components (those with l > maxlvar) are     ;
;            unchanged and therefor their contributions to    ; 
;            the field can be stored and re-used for multiple ;
;            inverse transform calculations.                  ;
;                                                             ;
; Inputs: transform - spherical harmonic transform to be inv- ;
;           erted                                             ;
;         altitude - index giving the height in the PFSS grid ;
;           at which the inversion is being calculated; in    ;
;           the SolarSoft PFSS library the field is calculated;
;           by repeated calculation of a inverse transform    ;
;           at each radius and for each component             ;
;           (r,phi,theta)                                     ;
;         component - field component being calculated (to    ;
;           reference stablefieldcomps), 1=b_r,2=b_ph, 3=b_th ;
;                                                             ;
; Keywords: lmaxi - if being called from pfss_potl_field_sj   ;
;             and 'trunc' keyword was set, indicates the max- ;
;             imum l value to include in the summation; for   ;
;             models with high lmax it can be impractical     ;
;             and/or wasteful to use all the l values at      ;
;             higher altitudes                                ;
;           recursive - if set, the routine is being called by;
;             itself in order to initialize the               ;
;             stablefieldcomps variale; it will proceed even  ;
;             though the values in this variable are zero     ;
;                                                             ;
; Notes: Note that this routine cannot be run unless the      ;
;     coeffs array in the transform_coefficients block has    ;
;     been set.                                               ;
;                                                             ;
; Created: July 2017                                          ;
;                                                             ;
; Created by: Shaela Jones                                    ;
;                                                             ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


; parameters
maxmaxlvar=8


; common block variables
@pfss_opt_parameters
@transform_coefficients
@pfss_data_block
; contains:
;   normalization and Legendre multipliers for each term in the
;      inverse transform
;   theta and phi grid values
;   saved values for the portion of the field that comes from 
;      components with l greather than maxlvar (the components
;      not being optimized)



; check size of coefficient array, print warning if necessary
szcoeffs=SIZE(coeffs)
if szcoeffs[0] ne 3 then begin
  print,'Error.  Coefficient array must be provided via transport_coeffs common variable, or keyword "calc_coeffs" must be set.'
  return,-1
endif


; check size of saved field components array
szstablefieldcomps=SIZE(stablefieldcomps)
if szstablefieldcomps[0] ne 4 then begin
  print,'Error.  Variable stablefieldcomps has not been created.  To create it, run inv_spherical_transform_sj.pro'
  return,-1
endif


; check that saved field components has been calculated for this altitude
if NOT(fieldcomps_set[altitude,component]) then begin
  if NOT(KEYWORD_SET(recursive)) then begin
    full_trans=INV_SPHERICAL_TRANSFORM_SJ(transform,cth,lmax= $
         lmaxi)
    small_l_trans=OPTIMIZATION_INV_TRANSFORM(transform, $
         altitude,component,lmaxi=lmaxi,/recursive)
    stablefieldcomps[*,*,altitude,component]=full_trans- $
         small_l_trans
    fieldcomps_set[altitude,component]=1
    RETURN,full_trans
  endif else begin
    if NOT(MAX(ABS(stablefieldcomps[*,*,altitude, $
         component]))) eq 0 then begin
      print,'optimization_inverse_transform:  Variable stablefieldcomps non-zero while fieldcomps_set=0.  Re-calculate or halt?'
      print,'1=recalculate'
      print,'other=halt execution'
      read,user_input
      case user_input of
        1: stablefieldcomps[*,*,altitude,component]=0d
        else:stop
      endcase
    endif
  endelse
endif


; initializations
phiix=phiix_stored
thindex=thindex_stored
nphi=N_ELEMENTS(phiix)
ntheta=N_ELEMENTS(thindex)
field=REFORM(stablefieldcomps[*,*,altitude,component])
if keyword_set(lmaxi) then lmaxi=round(long(lmaxi(0))) $
     else lmaxi=n_elements(transform(*,0))-1




; calculate field at specified altitude
localmarr=INDGEN(maxlvar+1)
lmaxbymmax=(maxlvar+1)*(maxlvar+1)
if maxlvar le maxmaxlvar then begin    ; if arrays can be kept small enough, do more efficient calculation
;   lousy condition - should be based on total size of array - fix later

  ; calculate amplitude and phase arrays
  subtrans=transform[0:maxlvar,0:maxlvar]
  amp=ABS(subtrans)
  phase=ATAN(imaginary(subtrans),double(subtrans))

  ; calculate argument of trig function
  mbyphi=phiix#DOUBLE(localmarr)
  mbyphi=REBIN(REFORM(mbyphi,nphi,maxlvar+1,1), $
         [nphi,maxlvar+1,maxlvar+1])
  mbyphi=TRANSPOSE(mbyphi,[0,2,1])
;  angpart=mbyphi+REBIN(REFORM(phase,1,lmaxi+1,maxlvar+1), $
;       [nphi,lmaxi+1,maxlvar+1])
  angpart=mbyphi+REBIN(REFORM(phase,1,maxlvar+1,maxlvar+1), $
       [nphi,maxlvar+1,maxlvar+1])
  angpart=REFORM(angpart,[nphi,lmaxbymmax])

  ; calculate summand, use matrix multiplication to sum over m,l
  lpoly=REBIN(REFORM(amp,1,lmaxbymmax),[nphi,lmaxbymmax])* $
       COS(REFORM(angpart,nphi,lmaxbymmax))
;  fieldcomps=TRANSPOSE(lpoly)## $
;       TRANSPOSE(REFORM(coeffs[0:lmaxi,0:maxlvar,*], $
;       lmaxbymmax,ntheta))
  fieldcomps=TRANSPOSE(lpoly)## $
       TRANSPOSE(REFORM(coeffs[0:maxlvar,0:maxlvar,*], $
       lmaxbymmax,ntheta))

  ; transpose and add to stable field components
  fieldcomps=TRANSPOSE(TEMPORARY(fieldcomps))
  field=field+fieldcomps
  
endif else begin      ; otherwise revert to nested loop solution
  
  ; calculate amplitude, phase
  amp=ABS(transform)
  phase=ATAN(imaginary(transform),double(transform))

  ;  start with m=l=0 mode
  field=field+REPLICATE(amp[0,0]*cos(phase[0,0])*coeffs[0], $
       [nphi,ntheta])
    
  ;  do other l modes for which m=0
  for l=1,maxlvar do begin
    trigfunc=REPLICATE(COS(phase[l,0]),nphi)# $
      REFORM(coeffs[l,0,*])
    field=field+(amp[l,0]*trigfunc)
  endfor
  
;  loop through m's for m>0, looping through l's > 0 for each m
  mbyphi=phiix#DOUBLE(localmarr)
  for m=1,maxlvar do begin
    for l=m,maxlvar do begin
      angpart=cos(mbyphi[*,m]+phase(l,m))
      field=field+amp(l,m)*(angpart#REFORM(coeffs[l,m,*]))
    endfor
  endfor  
  
endelse


return,field
end



; obsolete

;saveall=DBLARR(nphi,ntheta,lmax+1,lmax+1)
;saveall[*,*,0,0]=REPLICATE(amp[0,0]*cos(phase[0,0])*coeffs[0], $
;       [nphi,ntheta])
;saveall[*,*,l,0]=amp[l,0]*trigfunc
;saveall[*,*,l,m]=amp(l,m)*(angpart#REFORM(coeffs[l,m,*]))




