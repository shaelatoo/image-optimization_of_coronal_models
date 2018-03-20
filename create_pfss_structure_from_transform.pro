function create_PFSS_structure_from_transform,transform, $
      magfile=magfile,quiet=quiet

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                            ;;
;; Purpose: Creates the PFSS data structure resulting from    ;;
;;           the input transform                              ;;
;;                                                            ;;
;; Inputs: transform - an array of complex double-precision   ;;
;;           values describing the spherical harmonic trans-  ;;
;;           form of the magnetogram on which to base the     ;;
;;           PFSS structure                                   ;;
;;                                                            ;;
;; Returns: a "spherical" PFSS field data structure           ;;
;;                                                            ;;
;; Keywords: magfile - optional initialize several variables  ;;
;;             set by the pfss_mag_create routine             ;;
;;           quiet - if set, messages are suppressed          ;;
;;                                                            ;;
;; Dependencies: PFSS package from SolarSoft                  ;;
;;               CALC_PHIS                                    ;;
;;                                                            ;;
;; Warning: alters pfss common block values                   ;;
;;                                                            ;;
;; Created: 05 Aug 2016                                       ;;
;;                                                            ;;
;; Created by: Shaela Jones                                   ;;
;;                                                            ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; pfss common blocks
@pfss_data_block
@pfss_opt_parameters


; calculate magnetogram if keyword set
if KEYWORD_SET(magfile) then begin
  PFSS_MAG_CREATE_SJ,foo,magtype,nlat0,file=magfile,/quiet
  if NOT(KEYWORD_SET(quiet)) then print,'nlat= ',nlat
endif else begin
  if NOT(KEYWORD_SET(quiet)) then begin
    print,'using pre-existing values for variables set by PFSS_MAG_CREATE'
    print,'nlat= ',N_ELEMENTS(lat)
  endif
endelse


; fill in phiat and phibt values
CALC_PHIS,transform


; calculate field in the volume
PFSS_POTL_FIELD,rss,rgrid,/quiet,lmax=lmax
if NOT(KEYWORD_SET(quiet)) then begin
  print,'rss= ',rss
  print,'rgrid= ',rgrid
endif
if KEYWORD_SET(trunc) then print,'Truncating field calculation.'


; save data in a field data structure
PFSS_TO_SPHERICAL,output_pfss


return,output_pfss
end
