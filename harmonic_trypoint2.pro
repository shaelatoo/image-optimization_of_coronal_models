function harmonic_trypoint2,simplex,y,psum,windex,fac, $
  angles,coords,spcCoords,penalty_only=penalty_only, $
  bounds=bounds,weights=weights,penalize_magchange= $
  penalize_magchange,magweights=magweights,diff_power= $
  diff_power,rindex=rindex,omag=omag
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;                                                            ;;
  ;; Purpose: This routine is called by harmonic_amoeba to try  ;;
  ;;            altering one of the vertices in the simplex var-;;
  ;;            iable.  If the altered vertex gives a lower pen-;;
  ;;            alty function value it replaces the old vertex  ;;
  ;;            in the simplex.  Alternatively, this routine can;;
  ;;            be used to simply evaluate the penatly function ;;
  ;;            at the specified vertex. Harmonic_trypoint2 and ;;
  ;;            harmonic_amoeba2 were updates to the originals, ;;
  ;;            designed to use the more efficient transform    ;;
  ;;            computation algorithms in pfss_potl_field_sj.pro;;
  ;;                                                            ;;
  ;; Inputs: simplex - nvert x nvert+1 array of vertices in the ;;
  ;;           solution space of the harmonic_amoeba optimizat- ;;
  ;;           tion problem                                     ;;
  ;;         y - value of penalty function corresponding to each;;
  ;;           set of phase values
  ;;         psum - "center of mass" of the phase states in the ;;
  ;;           simplex array = total(simplex,2) before perturb. ;;
  ;;         windex - index of the worst vertex, or index of    ;;
  ;;           interest if penalty_only keyword is set          ;;
  ;;         fac - factor by which to alter position of worst   ;;
  ;;           vertex; or, set to 1 if penalty_only keyword is  ;;
  ;;           to be used                                       ;;
  ;;         angles - set of constraints indicating the angular ;;
  ;;           separation between the horizontal and the orient-;;
  ;;           ation of the B field                             ;;
  ;;         coords - heliocentric coordinates for the points   ;;
  ;;           where the elements of the angles variable were   ;;
  ;;           measured, in the form (longitude,co-latitude,radius);;
  ;;           with units (radians,radians,solar radii)         ;;
  ;;         spcCoords - 2xn array giving the latitude,        ;;
  ;;           longitude from which the constraint angles are to ;;
  ;;           be measured, in degrees                          ;;
  ;;         magweights - an optional array of weights indicat- ;;
  ;;           ing the relative reliability of each pixel in the;;
  ;;           original synoptic magnetogram                    ;;
  ;;                                                            ;;
  ;; Returns: the value of the penalty function at the point    ;;
  ;;;           being tested, as indicated by windex            ;;
  ;;                                                            ;;
  ;; Keywords: penalty_only - this keyword allows the function  ;;
  ;;           to be called to calculate the penalty function   ;;
  ;;           without altering any of the input variables (see ;;
  ;;           explanation above)                               ;;
  ;;         bounds: array giving upper,lower bounds for the    ;;
  ;;           theoretically possible values along each dim-    ;;
  ;;           ension of the state space; if set, penalty func- ;;
  ;;           tion is increased by a factor of ten whenever any;;
  ;;           magnetogram value is outside the specificied     ;;
  ;;           bounds                                           ;;
  ;;         weights - array of values that quantify the reliab-;;
  ;;           ility of the constraints in the variable angles; ;;
  ;;           if specified, each element of the fidelity pen-  ;;
  ;;           alty is multiplied by the corresponding element  ;;
  ;;           of weights; angles values believed to be more    ;;
  ;;           accurate should have higher weighting and there- ;;
  ;;           for a greater importance in the optimization     ;;
  ;;         penalize_magchange - if set, a term is added to   ;;
  ;;           the penalty funciton that penalizes changes to   ;;
  ;;           the magnetogram; the penalty due to each pixel   ;;
  ;;           can be individually weighted according to how re-;;
  ;;           liable the magnetogram pixel is believed to be   ;; 
  ;;         diff_power - by default the penalty function is pro-;;
  ;;           portional to the square of the residuals; set    ;;
  ;;           diff_power if any value other than two is desired;;
  ;;                                                            ;;
  ;; Dependencies: PFSS software in Solarsoft,                  ;;
  ;;                 pfss_opt_parameters                        ;;
  ;;                 extract_transform                          ;;
  ;;                 calc_phis                                  ;;
  ;;                 find_angle_values2                         ;;
  ;;                                                            ;;
  ;; Created: 07/17                                             ;;
  ;;                                                            ;;
  ;; Created by: Shaela Jones                                   ;;
  ;;                                                            ;;
  ;; Modified: 8/1/17 - changed penalty term implemented when   ;;
  ;;             penalize_magchange is set, so that inverse     ;;
  ;;             transform is not computed; B_r is used instead ;;
  ;;           8/7/17 - added omag keyword                      ;;
  ;;                                                            ;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  
  ; initializations
  @pfss_data_block
  @pfss_opt_parameters
  nlat=N_ELEMENTS(lat)
  nrix=N_ELEMENTS(rix)
  nconstraints=N_ELEMENTS(angles)
  szsimplex=SIZE(simplex)
  if KEYWORD_SET(rindex) then rgrid=3
  
  
  ; set keyword defaults - only necessary when routin is run stand-alone, otherwise defaults set by harmonic_amoeba.pro
  if N_ELEMENTS(weights) eq 0. then weights=FLTARR(nconstraints)+1.
  if NOT(KEYWORD_SET(diff_power)) then diff_power=2.
  if N_ELEMENTS(penalize_magchange) eq 0 then begin
    penalize_magchange=0
  endif else begin
    if N_ELEMENTS(magweights) eq 0 then magweights=1.
  endelse
  
  
  
  ; find new vertex
  fac1=(1-fac)/N_ELEMENTS(psum)
  fac2=fac1-fac
  new_vertex=psum*fac1-simplex[*,windex]*fac2
  
  
  ; re-configure transform, calculate phiat,phibt in pfss common block
  newmagt=EXTRACT_TRANSFORM(new_vertex,magt,maxlvar)
  CALC_PHIS,newmagt

  
  ; extrapolate magnetic field of vertex
  PFSS_POTL_FIELD_SJ,/trunc,/quiet


  ; find angle values to compare to constraints
  newangles=FIND_ANGLE_VALUES2(coords,spcCoords)
  
  
  ; calculate fidelity term of penalty function
  diffs=ABS(angles-newangles)
  list=WHERE(diffs ge !dpi/2.,cnt)
  if cnt ne 0 then diffs[list]=ABS(diffs[list]-!dpi)
  diffs=diffs^diff_power
  normal=nconstraints*MEDIAN(weights)*(berr*MEDIAN(angles))^diff_power
  penalty=TOTAL(weights*diffs)/normal


  ; add additional penalty term if vertex is outside of permitted solution space
  if KEYWORD_SET(bounds) then begin
    updiff=bounds-newmagt
    lowdiff=newmagt-bounds
    if MIN(updiff) lt 0 or MIN(lowdiff) lt 0 then penalty= $
         pen_mult*penalty
  endif
  
  
  
  ; add optional penalty terms for net flux or magnetogram changes
  if penalize_netflux ne 0. then begin
    if KEYWORD_SET(omag) then oldmag=omag else oldmag= $
      INV_SPHERICAL_TRANSFORM_SJ(magt,cth,lmax=lmax)
    newmag=OPTIMIZATION_INV_TRANSFORM(newmagt,nrix,0,lmax=lmax)
    penalty=penalty+penalize_netflux*TOTAL(newmag)^2.
    if penalize_magchange ne 0. then penalty=penalty+ $
        penalize_magchange*TOTAL(magweights*ABS(newmag-mag0))

  endif else if penalize_magchange ne 0. then begin
    if KEYWORD_SET(omag) then oldmag=omag else oldmag= $
         INV_SPHERICAL_TRANSFORM_SJ(magt,cth,lmax=lmax)
    penalty=penalty+penalize_magchange*TOTAL(magweights* $
         ABS(oldmag-OPTIMIZATION_INV_TRANSFORM(newmagt,nrix, $
         0,lmax=lmax)))
  endif
  
  
  ; if new vertex is an improvement, switch it into simplex
  if NOT(KEYWORD_SET(penalty_only)) then begin
    if penalty lt y[windex] then begin
      y[windex]=penalty
      simplex[*,windex]=new_vertex
      psum=TOTAL(simplex,2)
    endif
  endif
  
  
  RETURN,penalty
end
