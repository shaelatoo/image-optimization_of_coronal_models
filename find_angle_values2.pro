function find_angle_values2,coords,spcCoords,projected_strength= $
     projected_strength,direction=direction
     

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                            ;;
;; Purpose: Calculates the orientation of the magnetic field  ;;
;;            stored in the PFSS common block in a plane perp-;;
;;            endicular to a specified line of sight, at the  ;;
;;            specified coordinates.  Angles are returned in  ;;
;;            radians between 0 and pi, where 0 indicates the ;;
;;            magnetic field in the plane is parallel to the  ;;
;;            horizontal midline of the image, and pi/2 indic-;;
;;            ates it is perpendicular to it.                 ;;
;;            Note: this differs from the previous version of ;;
;;            this function in that it takes the solar B angle;;
;;            into account when calculating the orientation of;;
;;            the B field.  The previous version assumed the  ;;
;;            observing satellite was in the solar equatorial ;;
;;            plane.                                          ;;
;;                                                            ;;
;; Inputs: coords - 3xn array of heliocentric coordinates for ;;
;;           the points where the angles should be measured   ;;
;;         spcCoords - 2xn array giving the latitude,        ;;
;;           longitude from which the constraint angles are to ;;
;;           be measured, in degrees                          ;;
;;                                                            ;;
;; Returns: a 1D array of angles (in radians) between 0 and pi;;
;;                                                            ;;
;; Keywords: none                                             ;;
;;                                                            ;;
;; Dependencies: pfss solarsoft package                       ;;
;;                                                            ;;
;; Created: 11/04/15                                          ;;
;;                                                            ;;
;; Created by: Shaela Jones                                   ;;
;;                                                            ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; common block
@pfss_data_block
@pfss_opt_parameters


; initializations
latvec=90.-REFORM(coords[1,*])*180./!dpi
szcoords=SIZE(coords)
case szcoords[0] of
   2: nangs=szcoords[2] 
   1: nangs=1
   else: begin
       print,'something is wrong with coords array'
     end
endcase



; interpolate b components at given coords
phinds=GET_INTERPOLATION_INDEX(phi,REFORM(coords[0,*]))
thetainds=GET_INTERPOLATION_INDEX(lat,latvec)
rinds=GET_INTERPOLATION_INDEX(rix,REFORM(coords[2,*]))
bthetas=INTERPOLATE(bth,phinds,thetainds,rinds)
brads=INTERPOLATE(br,phinds,thetainds,rinds)
bphis=INTERPOLATE(bph,phinds,thetainds,rinds)


;; Find vectors that define horizontal and vertical image direction
;;       in the spherical coordinate grid on which B is defined,
;;       for each location in coords.  Calculate dot product of
;;       B with the image vectors at each location, and use the 
;;       ratio of these dot products to find the orientation angle
;;       of the model field in that image plane

if (noreset eq 0) then begin  ; if this is the first time the routine has been called
  coscrlt=REFORM(COS(spcCoords[0,*]*!dtor))
  sincrlt=REFORM(SIN(spcCoords[0,*]*!dtor))
  coscrln=REFORM(COS(spcCoords[1,*]*!dtor))
  sincrln=REFORM(SIN(spcCoords[1,*]*!dtor))
  imageNormals=TRANSPOSE([[coscrlt*coscrln],[coscrlt*sincrln], $
      [sincrlt]])
  imageY=TRANSPOSE([[-sincrln],[coscrln],[REPLICATE(0,nangs)]])
  imageZ=CROSSPN(imageNormals,imageY)
endif

bFieldSpherical=TRANSPOSE([[brads],[bthetas],[bphis]])
bFieldCart=CV_VECTOR(bFieldSpherical,coords)       ;CV_COORD(from_sphere=bFieldSpherical,/to_rect)
bDotIy=DOT_PRODUCT(bFieldCart,imageY)
bDotIz=DOT_PRODUCT(bFieldCart,imageZ)
extrap_angles=ATAN(bDotIz,bDotIy)
projected_strength=SQRT(bDotIz^2.+bDotIy^2.)



; get angle between 0 and pi
list=WHERE(extrap_angles lt 0,nlt)
direction=INTARR(nangs)+1   ; useful if using for line-tracing
direction[list]=-1
if nlt ne 0 then extrap_angles[list]=extrap_angles[list]+!dpi


return,extrap_angles
end



; saved garbage just in case it's useful

;I_y=TRANSPOSE([[SIN(spcCoords[0,*])],[COS(spcCoords[0,*])], $
;        [REPLICATE(0,nangs)]])
;I_z=TRANSPOSE([[-COS(spcCoords[0,*]*SIN(spcCoords[1,*])], $
;        [SIN(spcCoords[0,*]*SIN(spcCoords[1,*])], $
;        [COS(spcCoords[1,*])]])
;; transform to spherical coords  I_y,I_z(r,theta,phi) - different for each coordinate location
;I_ySph=CV_COORD(from_rect=I_y,/to_sphere)
;I_zSph=CV_COORD(from_rect=I_z,/to_sphere)
;; perform dot products with B at each coordinate location
;bvec=TRANSPOSE([[brads],[bphis],[bthetas]])
;bDotIy=DOT_PRODUCT(bvec,I_ySph)
;bDotIz=DOT_PRODUCT(bvec,I_zSph)
;; calc orientation angle from atan(B dot I_z, B dot I_x)
;extrap_angles=ATAN(bDotIz,bDotIy)
