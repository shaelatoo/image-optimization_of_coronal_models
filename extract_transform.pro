function extract_transform,vertex,magt,maxlvar

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                            ;;
;; Purpose: Called by harmonic_trypoint to copy the non-zero  ;;
;;            elements of a harmonic transform from a comp-   ;;
;;            ressed vector.  The size of the transform and   ;;
;;            any elements of the transform not included in   ;;
;;            the vector are copied from the variable magt    ;;
;;                                                            ;;
;; Inputs: vertex - compressed vector, a vertex from the simp-;;
;;           lex created by harmonic amoeba                   ;;
;;         magt - spherical transform on which vertex is based;;
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
;; Created: 03/20/15                                          ;;
;;                                                            ;;
;; Created by: Shaela Jones                                   ;;
;;                                                            ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; initializations
;halfway=N_ELEMENTS(vertex)/2
szmagt=SIZE(magt)
lmax=szmagt[1]


; set equal to original transform
newmagt=magt

; replace altered values
numalteredreals=0.5*((maxlvar+1)^2.+maxlvar+1)-1     ; -1 is for [0,0] value, which is left out

realvert=vertex[0:numalteredreals-1]
imagvert=0D
realreadpt=0
imagreadpt=numalteredreals
for i=1,maxlvar do begin
  imagvert=[0.,vertex[imagreadpt:imagreadpt+i-1]]
  newmagt[i,0:i]=COMPLEX(vertex[realreadpt:realreadpt+i],imagvert)
  realreadpt=realreadpt+i+1
  imagreadpt=imagreadpt+i
endfor


;realreadpt=0
;imagreadpt=numalteredreals
;for i=1,maxlvar do begin
;  realpart=vertex[realreadpt:realreadpt+1]
;  imagpart=[0.,vertex[imagreadpt:imagreadpt+i-1]]
;  newmagt[i,0:i]=COMPLEX(realpart,imagpart)
;endfor

; re-combine real and imaginary parts
;newmagvec=COMPLEX(vertex[0:halfway-1],vertex[halfway:*])


;; sort columns 
;newmagt=magt
;readpt=0
;for i=1,maxlvar do begin
;  newmagt[i,0:i]=newmagvec[readpt:readpt+i]
;  readpt=readpt+i+1
;endfor
;if maxlvar ne lmax then newmagt[maxlvar:*,*]=magt[maxlvar:*,*]


return,newmagt
end
