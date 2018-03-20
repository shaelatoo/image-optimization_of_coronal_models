function powell_objective,new_vertex
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                            ;;
;; Purpose: Calculates the value of the objective function    ;;
;;            for use by the IDL POWELL function to optimize  ;;
;;            the vector x.                                   ;;
;;                                                            ;;
;; Inputs: x - location in the solution space at which to     ;;
;;             calculate the objective funtion value          ;;
;;                                                            ;;
;;         via common block:                                  ;;
;;         angles - set of constraints indicating the angular ;;
;;           separation between the horizontal and the orient-;;
;;           ation of the B field                             ;;
;;         coords - heliocentric coordinates for the points   ;;
;;           where the elements of the angles variable were   ;;
;;           measured, in the form (longitude,latitude,radius);;
;;           with units (radians,radians,solar radii)         ;;
;;         spcCoords - 2xn array giving the longitude,        ;;
;;           latitude from which the constraint angles are to ;;
;;           be measured, in degrees                          ;;
;;                                                            ;;
;; Returns: the value of the penalty function at the point    ;;
;;            being tested                                    ;;
;;                                                            ;;
;; Keywords: via common block:                                ;; 
;;         weights - array of values that quantify the reliab-;;
;;           ility of the constraints in the variable angles; ;;
;;           if specified, each element of the fidelity pen-  ;;
;;           alty is multiplied by the corresponding element  ;;
;;           of weights; angles values believed to be more    ;;
;;           accurate should have higher weighting and there- ;;
;;           for a greater importance in the optimization     ;;
;;                                                            ;;
;; Dependencies: PFSS software in Solarsoft,                  ;;
;;                 pfss_opt_parameters                        ;;
;;                 extract_transform                          ;;
;;                 calc_phis                                  ;;
;;                 find_angle_values2                         ;;
;;                                                            ;;
;; Created: 01/12/17                                          ;;
;;                                                            ;;
;; Created by: Shaela Jones                                   ;;
;;                                                            ;;
;;                                                            ;;
 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  
; access common blocks
@pfss_data_block
@pfss_opt_parameters
@powell_variables


; parameters
diff_power=2.


; initializations 
nlat=N_ELEMENTS(lat)
nrix=N_ELEMENTS(rix)
nconstraints=N_ELEMENTS(angles)
if N_ELEMENTS(weights) eq 0. then weights=FLTARR(nconstraints)+1.
  
  
; calculate new transform, phiat,phibt in pfss common block
newmagt=EXTRACT_TRANSFORM(new_vertex,magt,maxlvar)
CALC_PHIS,newmagt

  
; extrapolate magnetic field
PFSS_POTL_FIELD,rss,rgrid,/trunc,/quiet


; find angle values to compare to constraints
newangles=FIND_ANGLE_VALUES2(coords,spcCoords)
  
  
; calculate fidelity term of penalty function
diffs=ABS(angles-newangles)
list=WHERE(diffs ge !dpi/2.,cnt)
if cnt ne 0 then diffs[list]=ABS(diffs[list]-!dpi)
diffs=diffs^diff_power
normal=nconstraints*MEDIAN(weights)*(berr*MEDIAN(angles))^diff_power
objective=TOTAL(weights*diffs)/normal

  
RETURN,objective
end
