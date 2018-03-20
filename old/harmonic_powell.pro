function harmonic_powell,magfile,ftol,function_value=y, $
  maxits=maxits,initmag=initmag,inittrans=inittrans, $
  _extra=extraKeys

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; NAME: harmonic_powell
;;
;; PURPOSE: Optimizes agreement between a PFSS extrapolation of the
;;            magnetogram in magfile and the constraints specified in
;;            angles,coords by varying the elements of the transform
;;            of the original magnetgram.
;;
;; INPUTS: magfile - filename containing magnetogram to be optimized
;;         ftol - the fractional decrease in the penalty function value
;;           in the terminating step.  This should never be less than the
;;           machine precision.
;;         from the (powell_variables) common block:
;;         angles - n-element list of angles describing the approximate
;;           orientation of the B field at the locations indicated by
;;           coords
;;         coords - 3xn array of coordinates for the locations at which
;;           we have B field constraints; three elements in first dimension
;;           are (longitude,co-latitude,radius) in units of
;;           (radians,radians,solar_radii)
;;         spcCoords - 2xn array giving the latitude,        
;;           longitude from which the constraint angles are to 
;;           be measured, in degrees                                    
;;
;; KEYWORD PARAMETERS:
;;    maxits - maximum number of improvement steps routine should attempt
;;    maxlvar: largest angular index l to optimize over; restricts 
;;      optimization to lower-order spherical harmonics
;;    weights - array of n values that quantify the reliability of the
;;      angle measurements in variable angles
;;    initmag - if set, magfile will not be processed and initmag will be
;;      used as the initial magnetogram; can save time for high-res. re-runs
;;    inittrans - if set, should be a transform of a magnetogram for 
;;      use in developing the initial state; this transform is used in place
;;      of the calculated value of the inital transform;primarily for use 
;;      in artificial test cases
;;
;; Returns:
;;   Result: if the optimization converged, the magnetogram corresponding 
;;     to the phase values that minimize the objective function, else -1.
;;
;; Dependencies: powell_objective (called by idl POWELL routine), $
;;     find_angle_values2, pfss_mag_create_sj,calc_phis, extract_transform, $
;;     PFSS solarsoft package
;;
;; Created: 01/12/17 by Shaela Jones, based on harmonic_amoeba
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 

 
; initialize pfss common block
@pfss_data_block
@pfss_opt_parameters
@powell_variables
  
  
; parameters
defaultmaxits=5000L
berr=0.1
pen_mult=10.
noreset=0
  
  
  
; check that pfss_opt_parameters have been defined
if N_ELEMENTS(rss) eq 0 or N_ELEMENTS(rgrid) eq 0 or $
     N_ELEMENTS(magtype) eq 0 or N_ELEMENTS(nlat0) eq 0 or $
     N_ELEMENTS(noreset) eq 0 then begin
  print,'Important variables from pfss_opt_parameters not defined'
  return,-1
endif
  
  
; initializations
if n_elements(maxits) eq 0 then maxits = defaultmaxits
if N_ELEMENTS(maxlvar) eq 0 then maxlvar=nlat0
lmax=nlat0    ;; make lmax a keyword?
nangles=N_ELEMENTS(angles)
if N_ELEMENTS(weights) eq 0 then weights=FLTARR(nangles)+1
  


; read magfile and initialize some elements of data block
if NOT(KEYWORD_SET(initmag)) then begin
  pfss_mag_create_sj,magnetogram,magtype,nlat0,file=magfile,/quiet
  magnetogram=magnetogram-mean(magnetogram)
  omag=magnetogram
  nlat=n_elements(theta)
  nlon=2*nlat
endif else begin
  szinitmag=SIZE(initmag)
  nlat=szinitmag[2]
  nlon=2*nlat
  lat=linrange(nlat+1,-90,90)
  lat=(lat(0:nlat-1)+lat(1:nlat))/2
  theta=(90-lat)*!dpi/180
  lon=(linrange(nlon+1,0,360))(0:nlon-1)
  phi=lon*!dtor
  magnetogram=initmag
endelse
cth=cos(theta)


  
; initialize vector to be optimized & magt
if KEYWORD_SET(inittrans) then begin
  magt=inittrans
endif else magt=SPHERICAL_TRANSFORM(magnetogram,cth,lmax=lmax)
realmagt=REAL_PART(magt)
imagmagt=IMAGINARY(magt)
sim0=REFORM(realmagt[1,0:1])
for i=2,maxlvar do sim0=[sim0,REFORM(realmagt[i,0:i])]
for i=1,maxlvar do sim0=[sim0,REFORM(imagmagt[i,1:i])]
nvert=N_ELEMENTS(sim0)  
  
  
; initialize vector matrix
vectors=FLTARR(nvert,nvert)
for i=0,nvert-1 do vectors[i,i]=1.      ;;; CHECK THIS ;;;;;

 
  
; call IDL's POWELL function to perform optimization
TIC
POWELL,sim0,vectors,ftol,y,'powell_objective',itmax=maxits,iter=iter
convtime=TOC()
newmagt=EXTRACT_TRANSFORM(sim0,magt,maxlvar)
result=INV_SPHERICAL_TRANSFORM(newmagt,cth) 
noreset=0     
print,'Convergence Time: ',convtime/60.,' minutes
print,'Objective Function at Final Iteration: ',y
return, result    


end
