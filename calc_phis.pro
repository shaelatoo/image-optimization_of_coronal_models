pro calc_phis,newmagt

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                            ;;
;; Purpose: Calculates the values of phiat and phibt, which   ;;
;;            the PFSS method uses to calculate the magnetic  ;;
;;            field in the extrapolation volume.  Normally    ;;
;;            they are calculated using the routine           ;;
;;            pfss_get_potl_coeffs                            ;;
;;                                                            ;;
;; Inputs: newmagt - an alter spherical harmonic transform of ;;
;;           the magnetogram being optimized                  ;;
;;                                                            ;;
;; Outputs: (in common block) phiat,phibt - see PFSS document-;;
;;            ation                                           ;;
;;                                                            ;;
;; Keywords: none                                             ;;
;;                                                            ;;
;; Dependencies: pfss solarsoft package,pfss_opt_parameters   ;;
;;                                                            ;;
;; Created: 03/19/15                                          ;;
;;                                                            ;;
;; Created by: Shaela Jones                                   ;;
;;                                                            ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; common blocks
@pfss_data_block
@pfss_opt_parameters


; get l and m index arrays, if necessary
szlarr=SIZE(larr)
if szlarr[1] ne lmax+1 then begin
  lix=lindgen(lmax+1)
  mix=lix
  larr=lix#replicate(1,lmax+1)
  marr=replicate(1,lmax+1)#mix
  wh=where(marr gt larr)
  larr(wh)=0  &  marr(wh)=0
endif


;  determine coefficients
phibt=-newmagt/(1+larr*(1+rss^(-2.*larr-1)))
phiat=-phibt/(rss^(2.*larr+1))
wh=where(finite(phiat) eq 0b,nwh)
if nwh gt 0 then phiat(wh)=complex(0,0)


end
