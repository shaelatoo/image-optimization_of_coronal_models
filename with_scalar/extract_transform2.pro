function extract_transform2,vertex,magt,maxlvar

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                            ;;
;; Purpose: Copies elements of a spherical harmonic transform ;;
;;            from a compressed vector into the full trans-   ;;
;;            form.  The size of the transform and            ;;
;;            any elements of the transform not included in   ;;
;;            the vector are copied from the variable magt.   ;;
;;            Modified from original extract_transform.pro    ;;
;;            by including the l=0,m=0 element, which was     ;;
;;            previously assumed excluded.                    ;;
;;                                                            ;;
;; Inputs: vertex - compressed vector, a vertex from the simp-;;
;;           lex created by harmonic amoeba                   ;;
;;         magt - spherical transform into which vertex is to ;;
;;           be transferred                                   ;;
;;         maxlvar - largest l value for which vertex contains;;
;;           values                                           ;;
;;                                                            ;;
;; Returns: an array of the same size as magt, which the first;;
;;            maxlvar columns comprised of values from vertex ;;
;;            and the remainder (if any) of values from magt  ;;
;;                                                            ;;
;; Keywords: none                                             ;;
;;                                                            ;;
;; Dependencies: none                                         ;;
;;                                                            ;;
;; Created: 10/19/17                                          ;;
;;                                                            ;;
;; Created by: Shaela Jones                                   ;;
;;                                                            ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; initializations
szmagt=SIZE(magt)
lmax=szmagt[1]


; set new transform equal to original transform
newmagt=magt

; replace altered values
numalteredreals=0.5*((maxlvar+1)^2.+maxlvar+1)  
realvert=vertex[0:numalteredreals-1]    ; "real" transform coefficients
newmagt[0,0]=vertex[0]
imagvert=0D
realreadpt=1
imagreadpt=numalteredreals    ; starting position for imaginary transform coefficients in vertex
for i=1,maxlvar do begin
  imagvert=[0.,vertex[imagreadpt:imagreadpt+i-1]]   ; imaginary coefficients with m=0 always 0, not included in vertex
  newmagt[i,0:i]=COMPLEX(vertex[realreadpt:realreadpt+i],imagvert)
  realreadpt=realreadpt+i+1
  imagreadpt=imagreadpt+i
endfor


return,newmagt
end
