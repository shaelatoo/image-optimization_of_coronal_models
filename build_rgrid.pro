pro build_rgrid

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                            ;;
;; Purpose: Calculates the radial grid points for the PFSS    ;;
;;            model.  Heavily based on code from Marc DeRosa's;;
;;            pfss_potl_field.pro, but I wanted to be able to ;;
;;            just calculate the radial grid without doing the;;
;;            field calculation.                              ;;
;;                                                            ;;
;; Inputs: none - all input variables are parameters from     ;;
;;           pfss_opt_parameters common block                 ;;
;;                                                            ;;
;; Outputs: none - establishes values for rix variable in     ;;
;;            pfss_data_block common block                    ;;
;;                                                            ;;
;; Keywords: none                                             ;;
;;                                                            ;;
;; Dependencies: pfss_data_block.pro, pfss_opt_parameters.pro ;;
;;                                                            ;;
;; Created: 07/12/17                                          ;;
;;                                                            ;;
;; Created by: Shaela Jones                                   ;;
;;                                                            ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; declare common blocks
@pfss_opt_parameters
@pfss_data_block



; calc radial resoution
delta_rad=(!dpi/nlat0)  ;  r grid spacing at r=1, make it half avg lat grid spacing
r_range=[1,DOUBLE(rss)]  ;  range of r


; calculate grid
case rgrid of
  2: begin  ;  radial gridpoint separation is proportional to r^2
    rix=[r_range[0]]
    lastr=r_range[0]
    repeat begin
      nextr=lastr+delta_rad*(lastr/r_range[0])^2
      rix=[rix,nextr]
      lastr=nextr
    endrep until nextr ge r_range[1]
    rix2=rix/((max(rix)-r_range[0])/(r_range[1]-r_range[0]))
    rix=rix2+(r_range[0]-rix2[0])
  end
  3: begin  ;  custom radial grid
    if N_ELEMENTS(rindex) eq 0 then begin
      print,'  ERROR in pfss_potl_field: rindex must be set if rgrid=3'
      return
    endif
    rix=rindex
  end
  else: begin  ;  radial gridpoints uniformly spaced
    nr=ROUND((r_range[1]-r_range[0])/delta_rad)
    rix=LINRANGE(nr,r_range[0],r_range[1])
  end
endcase


end

