function objective_function,newmagt,angles,coords,spccoords, $
      weights=weights,diff_power=diff_power,spherical_field= $
      spherical_field,calc_spherical=calc_spherical,diffs=diffs, $
      newangles=newangles

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                            ;;
;; Purpose: To calculate the value of the objective function  ;;
;;            for a given spherical harmonic transform, given ;;
;;            the constraints provided by the angles, coords, ;;
;;            and spccoords variables.                        ;;
;;                                                            ;;
;; Inputs: newmagt - the spherical harmonic transform of a    ;;
;;           synoptic magnetogram; this transform is used to  ;;
;;           calculate the field for which the value of the   ;;
;;           objective function is desired; can be produced   ;;
;;           using spherical_transform(magnetogram,cos(theta),;;
;;           lmax=lmax                                        ;;
;;                                                            ;;
;; Returns: the value of the objective function               ;;
;;                                                            ;;
;; Keywords: weights - numberical weights to apply to resid-  ;;
;;             als, one for each value of angles; good choices;;
;;             might be to heavily weight constraints that are;;
;;             part of a long feature, or those that are not  ;;
;;             close to another, discrepant feature; defaults ;;
;;             to one for all values                          ;;
;;           diff_power - power to which residuals should be  ;;
;;             raised in the objective function; default value;;
;;             is two
;;                                                            ;;
;; Dependencies: ;;
;;                                                            ;;
;; Created: 04/03/17                                          ;;
;;                                                            ;;
;; Created by: Shaela Jones                                   ;;
;;                                                            ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; get common blocks
@pfss_data_block
@pfss_opt_parameters


; initializations
nconstraints=N_ELEMENTS(angles)
if NOT(KEYWORD_SET(diff_power)) then diff_power=2.
if NOT(KEYWORD_SET(weights)) then weights=REPLICATE(1.,nconstraints)


; calculate field
CALC_PHIS,newmagt
PFSS_POTL_FIELD,rss,rgrid,/trunc,/quiet
if KEYWORD_SET(calc_spherical) then PFSS_TO_SPHERICAL,spherical_field

; find angle values to compare to constraints
newangles=FIND_ANGLE_VALUES2(coords,spcCoords)


; calculate fidelity term of penalty function
diffs=ABS(angles-newangles)
list=WHERE(diffs ge !dpi/2.,cnt)
if cnt ne 0 then diffs[list]=ABS(diffs[list]-!dpi)
diffs_sq=diffs^diff_power
normal=nconstraints*MEDIAN(weights)*(berr*MEDIAN(angles))^diff_power
penalty=TOTAL(weights*diffs_sq)/normal



return,penalty
end
