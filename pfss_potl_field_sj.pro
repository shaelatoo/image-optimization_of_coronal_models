
pro pfss_potl_field_sj,trunc=trunc,quiet=quiet
  
; modified version of pfss_potl_field.pro by Marc DeRosa
;  see original docs below


;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
;
;  pfss_potl_field.pro - This procedure computes the potential magnetic field
;                        B(r,theta,phi) given the potential Phi(l,m,r)
;
;  usage: pfss_potl_field,rtop,rgrid,rindex=rindex,thindex=thindex,
;           phindex=phindex,lmax=lmax,/trunc,potl=potl,/quiet
;         where rtop=radius of uppermost gridpoint
;               rgrid=sets radial gridpoint spacing:
;                      1 = equally spaced (default)
;                      2 = grid spacing varies with r^2
;                      3 = custom radial grid given by the rindex keyword
;               rindex = custom array of radial coordinates for output grid
;               thindex = (optional) custom array of theta (colatitude)
;                         coordinates, in radians, for output grid.  If not
;                         specified existing latitudinal grid is used.
;               phindex = (optional) custom array of phi (longitude)
;                         coordinates, in radians, for output grid.  If not
;                         specified, existing longitudinal grid is used.
;               lmax=if set, only use this number of spherical harmonics in
;                    constructing the potential (and thus the field)
;               trunc=set to use fewer spherical harmonics when
;                     reconstructing B as you get farther out in radius
;               quiet = set for minimal screen output
;
;         and in the common block we have:
;               phiat=on input, (l,m) array of dcomplex coeffs,
;                     corresponding to r^l eigenfunction
;               phibt=on input, (l,m) array of dcomplex coeffs,
;                     corresponding to 1/r^(l+1) eigenfunction
;               (br,bth,bph)=in output, (r,theta,phi) components of B-field
;
;  Notes: -If thindex and phindex are set, the variables
;          nlon,nlat,lon,lat,theta,phi in the common block are reset to the
;          values commensurate with the arrays given in thindex and phindex.
;          The old values (which correspond to the magnetic field used in the
;          potential extrapolation) are lost forever!
;         -Be careful when using thindex to put points at theta=0 and/or at
;          theta=!dpi.  Usually this routine returns Inf's at these points in
;          Bth and Bph, which might cause problems down the road.
;
;  M.DeRosa - 30 Jan 2002 - converted from earlier script
;              8 Feb 2002 - added lmax keyword
;              2 Jul 2002 - added quiet keyword
;             30 Apr 2003 - now utilizes memory more efficiently, based on a
;                           suggestion from Bart De Pontieu
;             12 May 2003 - converted common block to PFSS package format
;             27 Mar 2007 - added rgrid option 3 and rindex,thindex,phindex
;                           keywords
;             23 May 2007 - fixed bug with nlat when thindex is specified
;             30 Nov 2007 - fixed bug with truncated l arrays
;             15 Apr 2009 - added check to see if rgrid is provided
;
;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
;-
  
  ;  include common block
  @pfss_data_block
  @pfss_opt_parameters
  
  
  ;  check input parameters
  if N_ELEMENTS(rix) eq 0 then begin
    PRINT,'Radial grid must be set in common variable pfss_data_block.'
    RETURN 
  endif
  if N_ELEMENTS(phiat) eq 0 or N_ELEMENTS(phibt) eq 0 then begin
    PRINT,'Variables phiat and phibt must be set in pfss_data_block.  Run pfss_get_potl_coeffs.pro first.'
  endif
  
  
  ;  initializations
  ntheta=N_ELEMENTS(theta)
  nphi=N_ELEMENTS(phi)

  
  
  ;  get l and m index arrays of transform
  phisiz=size(phiat,/dim)
  lix=lindgen(phisiz(0))
  mix=lindgen(phisiz(1))
  larr=lix#replicate(1,phisiz(1))
  marr=replicate(1,phisiz(0))#mix
  wh=where(marr gt larr)
  larr(wh)=0  &  marr(wh)=0
  
  
  ;  set up planar sin(theta) array
  thindex=theta
  stharr=replicate(1,nphi)#sqrt(1-(cos(thindex))^2)
  
  
  ;  compute lmax for each radius
  lmaxarr=lonarr(nr,/noz)
  if keyword_set(trunc) then begin  ;  include fewer l modes as you get higher up
    lmaxarr=lonarr(nr,/noz)
    lmaxarr[0]=nlat0<lmax
    for i=1,nr-1 do begin
      wh=where(rix(i)^lindgen(nlat0+1) gt 1e6,nwh)
      if nwh eq 0 then lmaxarr(i)=lmax else lmaxarr(i)=(wh(0)<lmax)
    endfor
  endif else lmaxarr=REPLICATE(nlat0<lmax,nr)  
  

  ;  compute Br in (r,l,m)-space
  bt=make_array(dim=[phisiz,nr],/noz,/dcomplex)
  for i=0,nr-1 do $
    bt(*,*,i)= phiat*larr*rix(i)^(larr-1.) - phibt*(larr+1)*rix(i)^(-larr-2.)
    
  ;  ...and then transform to (r,theta,phi)-space
  br=make_array(dim=[nphi,ntheta,nr],/float,/noz)
  for i=0,nr-1 do begin
    if not keyword_set(quiet) then $
         pfss_print_time,'  pfss_potl_field: computing Br:  ',i+1,nr,tst,slen,/elap
    br(*,*,i)=OPTIMIZATION_INV_TRANSFORM(bt(*,*,i),i,0, $
         lmaxi=lmaxarr(i))    
  endfor
  
  ;  compute sin(theta) * Bth in (r,l,m)-space...
  factor=sqrt(double(larr^2-marr^2)/double(4*larr^2-1))
  for i=0,nr-1 do begin
    bt(*,*,i)=(larr-1)*factor* $
      (shift(phiat,1,0)*rix(i)^(larr-2) + shift(phibt,1,0)*rix(i)^(-larr-1)) $
      - (larr+2)*shift(factor,-1,0)* $
      (shift(phiat,-1,0)*rix(i)^larr + shift(phibt,-1,0)*rix(i)^(-larr-3))
    bt(0,0,i)=-2*factor(1,0)*(phiat(1,0) + phibt(1,0)*rix(i)^(-3))
    bt(lmax,*,i)=(lmax-1)*factor(lmax,*)* $
      (phiat(lmax-1,*)*rix(i)^(lmax-2) + phibt(lmax-1,*)*rix(i)^(-lmax-1))
  endfor
  
  ;  ...and then compute Bth in (r,theta,phi)-space
  bth=make_array(dim=[nphi,ntheta,nr],/float,/noz)
  for i=0,nr-1 do begin
    if not keyword_set(quiet) then pfss_print_time, $
         '  pfss_potl_field: computing Bth:  ',i+1,nr,tst, $
         slen,/elap
    bth(*,*,i)=OPTIMIZATION_INV_TRANSFORM(bt(*,*,i),i,2,lmaxi= $
         lmaxarr[i])/stharr
  endfor
  
  ;  compute sin(theta) * Bph in (r,l,m)-space...
  for i=0,nr-1 do bt(*,*,i)=complex(0,1)*marr* $
    (phiat*rix(i)^(larr-1) + phibt*rix(i)^(-larr-2))
    
  ;  ...and then compute Bph in (r,theta,phi)-space
  bph=make_array(dim=[nphi,ntheta,nr],/float,/noz)
  for i=0,nr-1 do begin
    if not keyword_set(quiet) then $
      pfss_print_time,'  pfss_potl_field: computing Bph:  ',i+1,nr,tst,slen,/elap
    bph(*,*,i)=OPTIMIZATION_INV_TRANSFORM(bt(*,*,i),i,1,lmaxi= $
         lmaxarr[i])/stharr
  endfor
  
end


; obsolete

;  ;  now transform the field potential to (r,theta,phi)-space
;  if n_elements(potl) gt 0 then begin
;    potl=make_array(dim=[nlon,nlat,nr],/float,/noz)
;    for i=0,nr-1 do begin
;      if not keyword_set(quiet) then $
;        pfss_print_time,'  pfss_potl_field: computing the field potential:  ',$
;        i+1,nr,tst,slen,/elap
;      potl(*,*,i)=OPTIMIZATION_INV_TRANSFORM(phibt*rix(i)^ $
;           (-larr-1)+phiat*rix(i)^larr,i,lmaxi=lmaxarr[i])
;    endfor
;  endif

