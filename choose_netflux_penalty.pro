function choose_netflux_penalty,magfile,angles,coords,spccoords

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                            ;;
;; Purpose: Calculates a value for the parameter              ;;
;;             'netfluxpenalty' input to harmonic_amoeba2     ;;
;;             based on the initial value of the objective    ;;
;;             function given by harmonic_trypoint2, such that;;
;;             the net flux penalty term of the objective     ;;
;;             function would be 20% of the initial discre-   ;;
;;             pancy term if the mean flux density in the in- ;;
;;             put magnetogram were 0.01G.                    ;;               
;;                                                            ;;
;; Inputs: magfile - synoptic magnetogram file used by        ;;
;;           harmonic_amoeba2 to initialize coronal magnetic  ;;
;;           field model                                      ;;
;;         angles,coords,spccoords - input variables to       ;;
;;           harmonic_trypoint, used to calculate value of    ;;
;;           objective function
;;                                                            ;;
;; Returns: Returns the suggested value of the parameter      ;;
;;            that balances the importance of the fidelity    ;;
;;            and net flux penalty terms in the objective     ;;
;;            function used to optimize the coronal magnetic  ;;
;;            field model.                                    ;;
;;                                                            ;;
;; Keywords: none yet                                         ;;
;;                                                            ;;
;; Dependencies: pfss solarsoft library, form_initial_vertex, ;;
;;                 build_rgrid,inv_spherical_transform_sj,    ;;
;;                 optimization_inv_transform                 ;;
;;                                                            ;;
;; Created: 08/21/18                                          ;;
;;                                                            ;;
;; Side effects: will alter some variables in pfss-related    ;;
;;                 common blocks                              ;;
;;                                                            ;;
;; Created by: Shaela Jones                                   ;;
;;                                                            ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; common block variables
@pfss_opt_parameters
@pfss_data_block


; check that pfss_opt_parameters have been defined
if N_ELEMENTS(rss) eq 0 or N_ELEMENTS(rgrid) eq 0 or $
  N_ELEMENTS(magtype) eq 0 or N_ELEMENTS(nlat0) eq 0 or $
  N_ELEMENTS(noreset) eq 0 or N_ELEMENTS(maxlvar) eq 0 $
  then begin
  print,'Important variables from pfss_opt_parameters not defined'
  print,'Please provide values for variables:'
  print,'rss (suggestion: 2.5)'
  print, 'rgrid (suggestion: 1-equal radial spacing, 2-grows as r^2'
  print, 'magtype (suggestion: use find_magtype(<magfile>)'
  print, 'nlat0 (suggestion: 180)
  print, 'noreset (suggestion: initialize to 0)'
  print, 'maxlvar (suggestion: 6)'
  return,-1
endif


; initialize penalize_netflux to get a penalty-free objective
;   function value from harmonic_amoeba2
penalize_netflux=0.


; process magnetogram, create transform
PFSS_MAG_CREATE_SJ,magnetogram,magtype,nlat0,file=magfile, $
  /quiet
nlat=N_ELEMENTS(theta)
nlon=2*nlat
cth=COS(theta)
magt=SPHERICAL_TRANSFORM(magnetogram,cth,lmax=lmax)


; create initial simplex, radial grid
sim0=FORM_INITIAL_VERTEX(magt,maxlvar)
nvert=N_ELEMENTS(sim0)
simplex=sim0#REPLICATE(1.0,nvert+1)
BUILD_RGRID
nr=N_ELEMENTS(rix)


; initialize field components common block
phiix_stored=phi
thindex_stored=theta
stablefieldcomps=DBLARR(2*nlat0,nlat0,nr+1,3)
fieldcomps_set=BYTARR(nr+1,3)
PFSS_GET_POTL_COEFFS,magnetogram,rtop=rss
invmagt=INV_SPHERICAL_TRANSFORM_SJ(magt,cth,/calc_coeffs, $
  lmax=lmax)   ; populates coeffs array
invmagt=OPTIMIZATION_INV_TRANSFORM(magt,nr,0,lmaxi=lmax)  ; populates transform_coefficients array for input magnetogram
fieldcomps_set[nr,0]=1
PFSS_POTL_FIELD_SJ,/trunc,/quiet   ; populates fieldcomps array
; calculate value of objective function with no net flux penalty


; calculate initial objective function value
y=FLTARR(nvert+1)
psum=TOTAL(simplex,2)
obj_fcn=HARMONIC_TRYPOINT2(simplex,y,psum,0,1., $
  angles,coords,spcCoords,/penalty_only)


; create net flux penalty that balances initial penalty function
netfluxpen=20.*obj_fcn



return,netfluxpen
end
