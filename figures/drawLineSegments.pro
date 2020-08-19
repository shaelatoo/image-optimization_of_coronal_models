function drawLineSegments,angles,coordstruc,img,strength, $
     angles2,coords2,strength2,_extra=ex,xsz=xsz, $
     residual=residual

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                            ;;
;; Purpose: Plots a coronal image and then plots a vector     ;;
;;            field over it (sort of).  For use in comparing  ;;
;;            Vadim's image-based constraints to model field  ;;
;;            projections.                                    ;;
;;                                                            ;;
;; Inputs: angles - a set of angles defining the slopes of the ;;
;;           vectors to be drawn                              ;;
;;         coordstruc - structure containing x and y locations    ;;
;;           (in pixels) where the vectors should be centered ;;
;;         img - a 2D image array (or fits file name) to be   ;;
;;           plotted                                          ;;
;;                                                            ;;
;; Returns: an image object                                   ;;
;;                                                            ;;
;; Keywords: xsz - dimension of output image plot window      ;;
;;                                                            ;;
;; Dependencies: none                                         ;;
;;                                                            ;;
;; Created: 02/03/16                                          ;;
;;                                                            ;;
;; Created by: Shaela Jones                                   ;;
;;                                                            ;;
;; Modifications: (3/28/17) modified to work with Vadim's     ;;
;;                  new coordinate structure by checking the  ;;
;;                  tag names associated with x,y coords      ;;
;;                (4/21/17) added xsz keyword to allow plott- ;;
;;                  ing smaller images                        ;;
;;                (4/19/18) added residual keyword to color-  ;;
;;                  code the output based on the magnitude of ;;
;;                  a second vector argument; intended for use;;
;;                  in identifying areas where the agreement  ;;
;;                  with a model is poor, either before or    ;;
;;                  after optimization                        ;;
;;                                                            ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; initializations
szimg=SIZE(img)
if szimg[0] ne 2 then begin
  print,'Reading FITS file'
  newimg=READFITS(img)
endif else newimg=img


; plot image
if NOT(KEYWORD_SET(xsz)) then xsz=800
imgObj=IMAGE(newimg,dimensions=[xsz,xsz],_extra=ex,/buffer)


; check which coordstruc tags were provided
if TAG_EXIST(coordstruc,'x') then begin
  xs=coordstruc.x
  ys=coordstruc.y
  if N_ELEMENTS(coords2) ne 0 then begin
    xs2=coords2.x
    ys2=coords2.y
  endif
endif else begin
  xs=coordstruc.angles_xx_r
  ys=coordstruc.angles_yy_r
  if N_ELEMENTS(coords2) ne 0 then begin
    xs2=coords2.angles_xx_r
    ys2=coords2.angles_yy_r
  endif
endelse



; overplot vector field
if N_ELEMENTS(strength) eq 0 then strength=185
u=strength*COS(angles)
v=strength*SIN(angles)
if KEYWORD_SET(residual) then begin
  vecs=VECTOR(u,v,xs,ys,/overplot,rgb_table=56,vector_colors=residual) 
endif else begin
  vecs=VECTOR(u,v,xs,ys,/overplot)
  ;if max(strength) eq 255 then vecs.COLOR='white'
  vecs.COLOR='white'
endelse
vecs.data_location=1
vecs.length_scale=0.1
vecs.head_size=0



; if desired, change vector colors to reflect some additional property of the data
;if KEYWORD_SET(residual) then begin
;  vecs.rgb_table=56
;  vecs.vector_colors=BYTSCL(residual)
;endif else if max(strength) eq 255 then vecs.COLOR='white'


; optionally, overplot second vector field
if N_ELEMENTS(angles2) ne 0 then begin
  if N_ELEMENTS(strength2) eq 0 then strength2=255
  u2=strength2*COS(angles2)
  v2=strength2*SIN(angles2)
  ;vecs2=VECTOR(u2,v2,xs2,ys2,/overplot,auto_color=1)
  vecs2=VECTOR(u2,v2,xs2,ys2,/overplot)
  vecs2.data_location=2
  vecs2.length_scale=0.1
  vecs2.head_size=0
  vecs2.color='red'
endif



return,imgObj
end
