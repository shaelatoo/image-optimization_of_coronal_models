function calc_realworld_coords,xs,ys,wcs
  ; calculates the real world coordinates from the image Cartesian
  ;    coordinates and the wcs header

; initializations
  crlnobs=wcs.position.crln_obs*!dtor
  ctr=wcs.crpix
  convert=wcs.position.dsun_obs/206265.0/1000.0/695500
  nangs = n_elements(xs)
  coords=FLTARR(3,nangs)

;  transform from image coordinates to angular sep from sun
  imcoords=transpose([[xs],[ys]])
  cartcoords=WCS_GET_COORD(wcs,imcoords)

;  find plane-of-sky longitudes
  listgt=WHERE(cartcoords[0,*] ge 0,ngth)
  listlt=WHERE(cartcoords[0,*] lt 0,nlt)
  if ngth gt 0 then coords[0,listgt]=crlnobs+!dpi/2.
  if nlt gt 0 then coords[0,listlt]=crlnobs-!dpi/2.

; calc radial coordinates
  coords[2,*]=REFORM(SQRT(cartcoords[0,*]^2.+cartcoords[1,*]^2.))* $
       convert

;  calc latitudes
  thets=ATAN(cartcoords[1,*],cartcoords[0,*])
  list=WHERE(thets gt !dpi/2.,cnt)        ; second quadrant
  if cnt ne 0 then coords[1,list]=thets[list]-!dpi/2.
  list=WHERE(thets le !dpi/2. and thets ge -!dpi/2,cnt)  ; first & third quadrants
  if cnt ne 0 then coords[1,list]=!dpi/2.-thets[list]
  list=WHERE(thets lt -!dpi/2.,cnt)         ; fourth quadrant
  if cnt ne 0 then coords[1,list]=3*!dpi/2.+thets[list]

;  check angles for out of bounds values
  listgt=WHERE(coords[0,*] ge 2*!dpi,ngth)
  if ngth ne 0 then coords[0,listgt]=coords[0,listgt]-2*!dpi
  listlt=WHERE(coords[0,*] lt 0,nlt)
  if nlt ne 0 then coords[0,listlt]=coords[0,listlt]+2*!dpi

;  correct for solar B angle
  solarb=wcs.position.solar_b0*!dtor
  csolb=COS(solarb)
  mcsolb=1-csolb
  ssolb=SIN(solarb)
  ux=COS(crlnobs+!dpi/2.)
  uy=SIN(crlnobs+!dpi/2.)
; adjust lats to be between -pi/2 and pi/2 - expected by CV_COORD
  tempcoords=coords
  tempcoords[1,*]=!pi/2.-tempcoords[1,*]
  rectcoords=CV_COORD(from_sphere=tempcoords,/to_rect)
  rot_matrix=[[csolb+ux^2.*mcsolb, ux*uy*mcsolb, -uy*ssolb], $
            [uy*ux*mcsolb, csolb+uy^2.*mcsolb, ux*ssolb], $
            [uy*ssolb, -ux*ssolb, csolb]]   ; matrix is transposed from usual sense due to IDL's strange way of multiplying matrices
  newrectcoords=TRANSPOSE(MATRIX_MULTIPLY(rectcoords,rot_matrix,/atranspose))
  coords=CV_COORD(from_rect=newrectcoords,/to_sphere)
  coords[1,*]=!dpi/2.-coords[1,*]
  list=WHERE(coords[0,*] lt 0.,ncnt)
  if ncnt ne 0 then coords[0,list]=coords[0,list]+2*!dpi
  list=WHERE(coords[0,*] ge 2*!dpi,cnt)
  if cnt ne 0 then coords[0,list]=coords[0,list]-2*!dpi

return, coords
end
