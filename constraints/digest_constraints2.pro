pro digest_constraints2,constptrs,hdrs,img_enhs,angles,coords, $
      lengths,weights,spcCoords,filter=filter

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                            ;;
;; Purpose: Reads in sample data from Vadim Uritsky's image   ;;
;;            processing software and uses it to create       ;;
;;            variables useful to harmonic_amoeba; this update;;
;;            to diget_constraints accepts data structures    ;;
;;            with updated tag names                          ;;
;;                                                            ;;
;; Inputs: constptrs - an array of pointers to constraint     ;;
;;           structures; each pointer points to a structure   ;;
;;           with 10 tags; the important ones for this routine;;
;;           are:                                             ;;
;;               angles_p - vector of angles indicating       ;;
;;                 feature orientation                        ;;
;;               angles_xx_r,angles_yy_r - vectors giving loc-;;
;;                 ations where angles were measured, in pixel;;
;;                 coordinates                                ;;
;;               l - lengths of detected features             ;;
;;         hdrs - array of header files corresponding to each ;;
;;           structure in constptrs                           ;;
;;         img_enhs - set of enhanced coronagraph images from ;;
;;           which features contained in constptrs were       ;;
;;           obtained                                         ;;
;;                                                            ;;
;; Outputs: angles - the measured POS B field angles          ;;
;;          coords - 3xn array of coordinates where angles    ;;
;;            were measured (longitude,latitude,radius) in    ;;
;;            units of (radians,radians,solar radii)          ;;
;;          spcCoords - 2xn array giving the latitude,       ;;
;;            longitude from which the constraint angles are to;;
;;            be measured, in degrees                         ;;
;;                                                            ;;
;; Keywords: filter - if set, calls filter_constraints.pro to ;;
;;             allow user to manually go through the constr-  ;;
;;             aints and remove those that are "pathological" ;;
;;             or combine duplicate featuresinto one, either  ;;
;;             by removing one or by picking a representative ;;
;;             set from the two                               ;;
;;                                                            ;;
;; Dependencies: none                                         ;;
;;                                                            ;;
;; Created: 04/05/16                                          ;;
;;                                                            ;;
;; Modifications: 02/20/18 - added filter keyword             ;;
;;                                                            ;;
;; Created by: Shaela Jones                                   ;;
;;                                                            ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; initializations
nconstptrs=N_ELEMENTS(constptrs)


; read in data structures, concatenate
for k=0,nconstptrs-1 do begin   ; read in data structures
  ; for each image, find feature and image parameters
  constraints=*constptrs[k]
;  angs=REFORM(constraints.angles_p, $
;      N_ELEMENTS(constraints.angles_p))
  hdr=REFORM(hdrs[k,*])
  wcs=FITSHEAD2WCS(hdr)
  crlnobs=wcs.position.crln_obs*!dtor
  ctr=wcs.crpix
;  xs=REFORM(constraints.angles_xx_r, $
;       N_ELEMENTS(constraints.angles_xx_r))
;  ys=REFORM(constraints.angles_yy_r, $
;       N_ELEMENTS(constraints.angles_yy_r))
  xs=constraints.angles_xx_r
  ys=constraints.angles_yy_r
  angs=constraints.angles_p
  ls=REFORM(constraints.l)

  ; if desired, manually filter out bad and/or duplicate features
  if KEYWORD_SET(filter) then begin
    filter_constraints,REFORM(img_enhs[k,*,*]),ctr,xs,ys, $
         angs,keepers
    keeperlist=WHERE(keepers eq 1,keepercnt)
    if keepercnt eq 0 then begin
      print,'Something is wrong with the filter results.'
      return
    endif
    xs=REFORM(xs[*,keeperlist],5*keepercnt)
    ys=REFORM(ys[*,keeperlist],5*keepercnt)
    angs=REFORM(angs[*,keeperlist],5*keepercnt)
    ls=ls[keeperlist]
  endif else begin
    nfeatures=N_ELEMENTS(constraints.angles_xx_r)
    xs=REFORM(xs,nfeatures)
    ys=REFORM(ys,nfeatures)
    angs=REFORM(angs,nfeatures)
  endelse
  
  ; some initializations
  nls=N_ELEMENTS(ls)
  lens=FLTARR(5L*nls)
  convert=wcs.position.dsun_obs/206265.0/1000.0/695500
  nanglesi=N_ELEMENTS(angs)
  coordsi=FLTARR(3,nanglesi)

;  transform from image coordinates to angular sep from sun
  imcoords=transpose([[xs],[ys]])
  cartcoords=WCS_GET_COORD(wcs,imcoords)

;  find plane-of-sky longitudes
  listgt=WHERE(cartcoords[0,*] ge 0,ngth)
  listlt=WHERE(cartcoords[0,*] lt 0,nlt)
  if ngth gt 0 then coordsi[0,listgt]=crlnobs+!dpi/2.
  if nlt gt 0 then coordsi[0,listlt]=crlnobs-!dpi/2.


; calc radial coordinates
  coordsi[2,*]=REFORM(SQRT(cartcoords[0,*]^2.+cartcoords[1,*]^2.))* $
       convert


;  calc latitudes
  thets=ATAN(cartcoords[1,*],cartcoords[0,*])
  list=WHERE(thets gt !dpi/2.,cnt)        ; second quadrant
  if cnt ne 0 then coordsi[1,list]=thets[list]-!dpi/2.
  list=WHERE(thets le !dpi/2. and thets ge -!dpi/2,cnt)  ; first & third quadrants
  if cnt ne 0 then coordsi[1,list]=!dpi/2.-thets[list]
  list=WHERE(thets lt -!dpi/2.,cnt)         ; fourth quadrant
  if cnt ne 0 then coordsi[1,list]=3*!dpi/2.+thets[list]


;  check angles for out of bounds values
  listgt=WHERE(coordsi[0,*] ge 2*!dpi,ngth)
  if ngth ne 0 then coordsi[0,listgt]=coordsi[0,listgt]-2*!dpi
  listlt=WHERE(coordsi[0,*] lt 0,nlt)
  if nlt ne 0 then coordsi[0,listlt]=coordsi[0,listlt]+2*!dpi


;  correct for solar B angle
  solarb=wcs.position.solar_b0*!dtor
  csolb=COS(solarb)
  mcsolb=1-csolb
  ssolb=SIN(solarb)
  ux=COS(crlnobs+!dpi/2.)
  uy=SIN(crlnobs+!dpi/2.)
; adjust lats to be between -pi/2 and pi/2
  tempcoords=coordsi
  tempcoords[1,*]=!pi/2.-tempcoords[1,*]
  rectcoords=CV_COORD(from_sphere=tempcoords,/to_rect)
  rot_matrix=[[csolb+ux^2.*mcsolb, ux*uy*mcsolb, -uy*ssolb], $
            [uy*ux*mcsolb, csolb+uy^2.*mcsolb, ux*ssolb], $
            [uy*ssolb, -ux*ssolb, csolb]]   ; matrix is transposed from usual sense due to IDL's strange way of multiplying matrices
;  newrectcoords=TRANSPOSE(rot_matrix##TRANSPOSE(rectcoords))
  newrectcoords=TRANSPOSE(MATRIX_MULTIPLY(rectcoords,rot_matrix,/atranspose))
  coordsi=CV_COORD(from_rect=newrectcoords,/to_sphere)
  coordsi[1,*]=!dpi/2.-coordsi[1,*]
  list=WHERE(coordsi[0,*] lt 0.,ncnt)
  if ncnt ne 0 then coordsi[0,list]=coordsi[0,list]+2*!dpi
  list=WHERE(coordsi[0,*] ge 2*!dpi,cnt)
  if cnt ne 0 then coordsi[0,list]=coordsi[0,list]-2*!dpi
  

  ;   create spcCoords array
  crlts=REPLICATE(wcs.position.crlt_obs,nanglesi)
  crlns=REPLICATE(wcs.position.crln_obs,nanglesi)
  spccs=TRANSPOSE([[crlts],[crlns]])
    
  
;  concatenate coordinate arrays
  for i=0L,nls-1 do lens[i*5:(i+1)*5-1]=REPLICATE(ls[i],5)
  if k eq 0 then begin
    angles=angs
    lengths=lens
    coords=coordsi
    spcCoords=spccs
  endif else begin
    angles=[angles,angs]
     lengths=[lengths,lens]
    coords=[[coords],[coordsi]]
    spcCoords=[[spcCoords],[spccs]]
  endelse


endfor


; calculate weights
weights=lengths/TOTAL(lengths)


; check for angle values greater than pi or less than 0
list=WHERE(angles ge !dpi,cnt)
if cnt ne 0 then begin
  angles[list]=angles[list]-!dpi
  print,'Note: some angle values were greater than pi, correcting...'
endif
list=WHERE(angles lt 0,cnt)
if cnt ne 0 then begin
  angles[list]=angles[list]+!dpi
  print,'Note: some angle values were less than zero, correcting...'
endif


end
