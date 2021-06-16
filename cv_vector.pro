function cv_vector,vector,coords

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                            ;;
;; Purpose: Currently, converts a vector or array of vectors  ;;
;;            from spherical coordinates to cartesian coor- ;;
;;            dinates.  Later, may be enhanced to be something;;
;;            more general like IDL's CV_COORD.               ;;
;;                                                            ;;
;; Inputs: vector - 3xn array of spherical vector field comp-;;
;;            onents (radial, polar, and azimuthal)           ;;
;;         coords - 3xn array of spherical coordinate locations;;
;;           at which vectors are given (longitude, co-latit- ;;
;;           ude, radial), with angles given in radians; co-  ;;
;;           latitude is between 0 and pi                     ;;
;;                                                            ;;
;; Returns: 3xn array of vectors in terms of cartesian coord- ;;
;;;           inates  (x,y,z)                                 ;;
;;                                                            ;;
;; Keywords: none                                             ;;
;;                                                            ;;
;; Dependencies: none                                         ;;
;;                                                            ;;
;; Created: 11/18/15                                          ;;
;;                                                            ;;
;; Modified: 2/12/16                                          ;;
;;                                                            ;;
;; Created by: Shaela Jones                                   ;;
;;                                                            ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; initializations
szcoords=SIZE(coords)
cmatrix=FLTARR(3,3,szcoords[2])
sinthetas=REFORM(SIN(coords[1,*]))
sinphis=REFORM(SIN(coords[0,*]))
costhetas=REFORM(COS(coords[1,*]))
cosphis=REFORM(COS(coords[0,*]))
cartvec=FLTARR(3,szcoords[2])


; define matrix elements
cmatrix[*,0,*]=TRANSPOSE([[sinthetas*cosphis],[costhetas*cosphis], $
      [-sinphis]])
cmatrix[*,1,*]=TRANSPOSE([[sinthetas*sinphis],[costhetas*sinphis], $
      [cosphis]])
cmatrix[*,2,*]=TRANSPOSE([[costhetas],[-sinthetas],[REPLICATE(0, $
      szcoords[2])]])


; multiply vectors by conversion matrix
for i=0,szcoords[2]-1 do cartvec[*,i]=cmatrix[*,*,i]##vector[*,i]




return,cartvec
end



