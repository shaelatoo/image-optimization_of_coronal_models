function harmonic_amoeba2,magfile,angles,coords,ftol,scale, $
  spcCoords,magmov,function_value=y,maxits=maxits,simplex= $
  simplex,tims=tims,min_values=min_values,bounds=bounds, $
  time_cutoff=time_cutoff,weights=weights, $
  magmovset=magmovset,initmag=initmag,inittrans=inittrans, $
  netfluxpenalty=netfluxpenalty,magchangepenalty= $
  magchangepenalty,magweights=magweights,rindex=rindex, $
  _extra=extraKeys,nums=nums,adapt_ind=adapt_ind
  
  
  
  ; Reference: Numerical Recipes, 2nd Edition, Page 411.
  ;
  ; NAME: harmonic_amoeba
  ;
  ; PURPOSE: Optimizes agreement between a PFSS extrapolation of the
  ;            magnetogram in magfile and the constraints specified in
  ;            angles,coords by varying the elements of the transform
  ;            of the original magnetgram.
  ;
  ; INPUTS: magfile - filename containing magnetogram to be optimized
  ;         angles - n-element list of angles describing the approximate 
  ;           orientation of the B field at the locations indicated by 
  ;           coords
  ;         coords - 3xn array of coordinates for the locations at which
  ;           we have B field constraints; three elements in first dimension
  ;           are (longitude,co-latitude,radius) in units of 
  ;           (radians,radians,solar_radii) 
  ;         ftol - the fractional decrease in the penalty function value
  ;           in the terminating step.  This should never be less than the
  ;           machine precision.
  ;         scale - scalar or nvert-element array specifying the size of 
  ;           the initial perturbations to the spherical transform components
  ;           (see simplex keyword for description of nvert)
  ;         spcCoords - 2xn array giving the latitude,        ;;
  ;           longitude from which the constraint angles are to ;;
  ;           be measured, in degrees                          ;;
  ;           
  ; OPTIONAL OUTPUT: magmov - image series showing progress of the mag-
  ;   netograms as the optimization progresses
  ;
  ; KEYWORD PARAMETERS:
  ;    function_value: on exit, an Nvert+1 element vector containing the
  ;      function values at the simplex vertices; if the program converged,
  ;      the minimum penalty function value should be in first element
  ;    maxits - maximum number of improvement steps routine should attempt
  ;    simplex: (nvert x nvert+1) array of vertices.  each vertex 
  ;      represents a point in the solution space at which the penalty 
  ;      function was evaluated.  if the optimization converged, the first 
  ;      image, simplex[*,*,0], will contain the solution corresponding to 
  ;      the minimum penalty function.  otherwise, it gives some idea where 
  ;      in the solution space the algorithm was looking.  nvert is deter-
  ;      mined by the maxlvar keyword, nvert=(maxlvert+1)*(maxlvert+2);
  ;      as simplex grows the computation becames more time-consuming and 
  ;      less accurate
  ;    tims - array giving the time (in seconds) since beginning optimization,
  ;      with a time corresponding to each element of min_values
  ;    min_values - array giving the lowest value in the function_value 
  ;      array at the end of each run through the optimization loop
  ;    bounds: nvert array giving upper bounds for the transform magnitudes;
  ;      if set, the penalty funciton is increased by a factor of ten 
  ;      whenever any transform component is outside the specificied bounds; 
  ;      nvert may also be a scalar, indicating the bound is the same
  ;      for every element of the transform 
  ;    time_cutoff: the time, in seconds, after which the optimization should
  ;      be terminated and the best known solution should be returned; default
  ;      is two weeks
  ;    maxlvar: largest angular index l to optimize over; restricts 
  ;      optimization to lower-order spherical harmonics
  ;    weights - array of n values that quantify the reliability of the
  ;      angle measurements in variable angles (see harmonic_trypoint)
  ;    magmovset - if set, magmov variable will be calculated; this may take
  ;      a significant amount of memory/time, depending on optimization
  ;      parameters
  ;    initmag - if set, should be a magnetogram for use in developing the
  ;      initial simplex; in this case the magfile will not be processed
  ;    inittrans - if set, should be a transform of a magnetograms for 
  ;      use in developing the initial simplex; in this case the magfile
  ;      will not be processed and initmag, if set, will be ignored; 
  ;      primarily for use in artificial test cases
  ;    netfluxpenalty - if set, a term is added to the penalty function
  ;      proportional to the net flux in the magnetogram; may be set to a 
  ;      constant value (not equal to one), or if set to one will be auto-
  ;      matically calculated at 20 times the initial value of the 
  ;      fidelity term
  ;    radial_grid - if set, should be a custom arrangement of radial
  ;      grid points at which to calculate the model; added for work
  ;      with Leon, who wanted a linear grid at high resolution to
  ;      facilitate interpolation to a higher-res cartesian grid
  ;
  ; Returns:
  ;   Result: if the optimization converged, the magnetogram corresponding 
  ;     to the phase values that minimize the penalty function, else -1.
  ;
  ; Dependencies: harmonic_trypoint, find_angle_values2, pfss_mag_create_sj,
  ;     calc_phis, extract_transform, PFSS solarsoft package
  ;
  ; Created: 03/18/15 by Shaela Jones
  ;
  ; Modified: 11/30/15 - modified to accept spcCoords variable and
  ;     pass it to harmonic_trypoint  (sij)
  ;            2/28/17 - removed maxlvar input and added to pfss_opt_parameters
  ;            5/30/17 - modified to accept a custom radial model grid
  ;            7/27/17 - removed mean subraction from magnetogram, now performed 
  ;                        in pfss_mag_create_sj before re-gridding
  ;            8/21/18 - incorporated caculation of net flux penalty
  ;                        parameter based on initial fidelity term, changed 
  ;                        usage of netfluxpenalty keyword
  
  
  
  ; initialize pfss common block
  @pfss_data_block
  @pfss_opt_parameters
  @transform_coefficients
  
  
  ; parameters
  defaultscale=3.
  defaultmaxits=5000L
  expansion=2.0
  reflection=-1.0
  contraction=0.5
  reduction=0.5
  minrange=0.001     ; minimum allowed value of the range variable, y(max)+y(min)
  tc_default=2.*7.*24.*60.*60.   ; default time cutoff: two weeks!
  berr=0.1
  pen_mult=10.
  noreset=0
  if maxlvar ge 4 then printcnt=100 else printcnt=30
  
  
  ; begin timing
  TIC
  
  
  ; if user supplied a radial grid, set rgrid to custom
  if N_ELEMENTS(rindex) ne 0 then rgrid=3
  

  ; check that pfss_opt_parameters have been defined
  if N_ELEMENTS(rss) eq 0 or N_ELEMENTS(rgrid) eq 0 or $
    N_ELEMENTS(magtype) eq 0 or N_ELEMENTS(nlat0) eq 0 or $
    N_ELEMENTS(noreset) eq 0 or N_ELEMENTS(maxlvar) eq 0 $
    then begin
    print,'Important variables from pfss_opt_parameters not defined'
    print,'Please provide values for variables:'
    print,'rss (suggestion: 2.5)'
    print, 'rgrid (suggestion: 1-equal radial spacing, 2-grows as r^2'
    print, 'magtype (suggestion: use find_magtype(<magfile>)'
    print, 'nlat0 (suggestion: 180)'
    print, 'noreset (suggestion: initialize to 0)'
    print, 'maxlvar (suggestion: 6)'
    return,-1
  endif
  
  
  ; initializations
  min_values=DBLARR(maxits)
  tims=DBLARR(maxits)
  nums=DBLARR(maxits)
  nangles=N_ELEMENTS(angles)
  if N_ELEMENTS(lmax) eq 0 then lmax=2/3.*nlat0
  if n_elements(maxits) eq 0 then maxits = defaultmaxits
  if N_ELEMENTS(scale) eq 0 then scale=defaultscale
  if NOT(KEYWORD_SET(time_cutoff)) then time_cutoff=tc_default
  if N_ELEMENTS(weights) eq 0 then weights=FLTARR(nangles)+1
  if KEYWORD_SET(magchangepenalty) then begin
    penalize_magchange=magchangepenalty
    if N_ELEMENTS(magweights) eq 0 then begin
      magweights=1.
      print,'Magweights not set; assuming uniform reliability.'
    endif
  endif else penalize_magchange=0.
  if n_elements(rgrid) eq 0 then rgrid=1

  
  ; read magfile and initialize some elements of data block
  if NOT(KEYWORD_SET(initmag)) then begin
    PFSS_MAG_CREATE_SJ,magnetogram,magtype,nlat0,file=magfile, $
         /quiet,adapt_ind=adapt_ind     ; keyword irrelevant if not using magfile=5
    nlat=N_ELEMENTS(theta)
    nlon=2*nlat
  endif else begin
    szinitmag=SIZE(initmag)
    nlat=szinitmag[2]
    nlon=2*nlat
    lat=LINRANGE(nlat+1,-90,90)
    lat=(lat(0:nlat-1)+lat(1:nlat))/2
    theta=(90-lat)*!dpi/180
    lon=(LINRANGE(nlon+1,0,360))(0:nlon-1)
    phi=lon*!dtor
    magnetogram=initmag
  endelse
  cth=COS(theta)
  magmov=magnetogram

  
  ; double-check magweights size
  szmagweights=SIZE(magweights)
  szmagnetogram=SIZE(magnetogram)
  if szmagweights[1] ne szmagnetogram[1] and szmagweights[0] ne 0 then begin
    print,'Variable magweights does not match magnetogram size.'
    return,-1
  endif
  
  
  ; initialize simplex & magt
  if KEYWORD_SET(inittrans) then begin
    magt=inittrans
  endif else begin
    magt=SPHERICAL_TRANSFORM(magnetogram,cth,lmax=lmax)
  endelse
  if maxlvar gt 1 then begin
    sim0=FORM_INITIAL_VERTEX(magt,maxlvar)
    nvert=N_ELEMENTS(sim0)
    simplex=sim0#REPLICATE(1.0,nvert+1)
    for i=0,nvert-1 do simplex[i,i+1]=sim0[i]+ $
        scale[i<(N_ELEMENTS(scale)-1)]
  endif else begin
    simplex=FLTARR(1,2)
    simplex[0,0]=realmagt[1,0]
    simplex[0,1]=realmagt[1,0]+scale
  endelse
  
  
  ; establish stablefieldcomps, save phi and theta grids
  if N_ELEMENTS(rindex) eq 0 then BUILD_RGRID else rix=rindex
  nr=N_ELEMENTS(rix)
  phiix_stored=phi
  thindex_stored=theta
  stablefieldcomps=DBLARR(2*nlat0,nlat0,nr+1,3)
  fieldcomps_set=BYTARR(nr+1,3)
  PFSS_GET_POTL_COEFFS,magnetogram,rtop=rss  
  invmagt=INV_SPHERICAL_TRANSFORM_SJ(magt,cth,/calc_coeffs, $
       lmax=lmax)   ; populates coeffs array
  invmagt=OPTIMIZATION_INV_TRANSFORM(magt,nr,0,lmaxi=lmax)  ; populates transform_coefficients array for input magnetogram
  fieldcomps_set[nr,0]=1
  PFSS_POTL_FIELD_SJ,/trunc,/quiet   ; populates fieldcomps array

  
  ; establish net flux penalty proportionality constant
  if NOT(KEYWORD_SET(netfluxpenalty)) then begin
    penalize_netflux=0.
  endif else begin
    if netfluxpenalty eq 1 then begin
      y=FLTARR(nvert+1)
      psum=TOTAL(simplex,2)
      penalize_netflux=0.
      obj_fcn=HARMONIC_TRYPOINT2(simplex,y,psum,0,1., $
           angles,coords,spcCoords,/penalty_only)
      penalize_netflux=20.*obj_fcn
    endif else penalize_netflux=netfluxpenalty
  endelse

  
  
  
  
  ; initialize vector of penalty function values
  psum=0
  y=FLTARR(nvert+1)
  for i=0,nvert do y[i]=HARMONIC_TRYPOINT2(simplex,y,psum,i,1., $
    angles,coords,spcCoords,/penalty_only,bounds=bounds, $
    weights=weights,magweights=magweights,penalize_magchange= $
    penalize_magchange,rindex=rindex,_extra=extraKeys,omag= $
    magnetogram)
  ncalls=LONG(nvert)    ; number of penalty function calculations
  noreset=1

   
  
  ; minimization loop
  cnt=0L
  psum = TOTAL(simplex,2)
    while cnt lt maxits do begin   ;Each iteration
    ord = SORT(y)
    lowest = ord[0]    ;Lowest point
    min_values[cnt]=y[lowest]   ; keeps track of min value each iteration
    highest = ord[nvert]   ;Highest point
    next_highest = ord[nvert-1]  ;Next highest point
    range = abs(y[highest]) + abs(y[lowest]) ;Denominator = interval
    if range ge minrange then rtol = 2.0 * abs(y[highest]-y[lowest])/range $
    else rtol = ftol / 2.   ;Terminate if interval is 0
    
    if rtol lt ftol then begin ;Done?
      t = y[0] & y[0] = y[lowest] & y[lowest] = t ;Sort so fcn min is 0th elem
      t = simplex[*,lowest] & simplex[*,lowest] = simplex[*,0] & simplex[*,0] = t
      print,'Simplex minimum size has been reached'
      newmagt=EXTRACT_TRANSFORM(t,magt,maxlvar)
      result=INV_SPHERICAL_TRANSFORM(newmagt,cth)
;    result=OPTIMIZATION_INV_TRANSFORM(newmagt,0,0)
      noreset=0
      return,result                
    endif
    
    if TOC() gt time_cutoff then begin   ; time cutoff reached?
      t = y[0] & y[0] = y[lowest] & y[lowest] = t ;Sort so fcn min is 0th elem
      t = simplex[*,lowest] & simplex[*,lowest] = simplex[*,0] & simplex[*,0] = t
      print,'Optimization time cutoff has been reached'
      newmagt=EXTRACT_TRANSFORM(REFORM(t),magt,maxlvar)
      result=INV_SPHERICAL_TRANSFORM(newmagt,cth)
;      result=OPTIMIZATION_INV_TRANSFORM(newmagt,0) 

      noreset=0     
      return, result
    endif
    
    ; try a reflection
;    print,'reflect'
    ytry=HARMONIC_TRYPOINT2(simplex,y,psum,highest,reflection, $
      angles,coords,spcCoords,bounds=bounds,weights= $
      weights,penalize_magchange=penalize_magchange,magweights= $
      magweights,rindex=rindex,_extra=extraKeys,omag=magnetogram)
    ncalls++
    
    ; if ytry is better than the best point, expand
;    print,'expand'
    if ytry le y[lowest] then begin
      ytry=HARMONIC_TRYPOINT2(simplex,y,psum,highest,expansion, $
        angles,coords,spcCoords,bounds=bounds,weights= $
        weights,penalize_magchange=penalize_magchange,magweights= $
        magweights,rindex=rindex,_extra=extraKeys,omag=magnetogram)
      ncalls++
    endif else if ytry ge y[next_highest] then begin
      ; if ytry is the new worst point, contract
;      print,'contract'
      ysave = y[highest]
      ytry=HARMONIC_TRYPOINT2(simplex,y,psum,highest,contraction, $
        angles,coords,spcCoords,bounds=bounds,weights= $
        weights,penalize_magchange=penalize_magchange, $
        magweights=magweights,rindex=rindex,_extra=extraKeys, $
        omag=magnetogram)
      ncalls++
      if ytry ge ysave then begin
        ; if the contracted point is still the new worst, reduce
        for i=0, nvert do begin
          if i ne lowest then begin
            psum = reduction * (simplex[*,i] + simplex[*,lowest])
            simplex[*,i] = psum
            y[i]=HARMONIC_TRYPOINT2(simplex,y,psum,i,1.0,angles, $
               coords,spcCoords,/penalty_only,bounds= $
               bounds,weights=weights,penalize_magchange= $
               penalize_magchange,magweights=magweights, $
               rindex=rindex,_extra=extraKeys,omag=magnetogram)
          endif
        endfor
        ncalls = ncalls + nvert
        psum = TOTAL(simplex,2)
      endif   ;ytry ge ysave
    endif
    tims[cnt]=TOC()
    nums[cnt]=cnt
    if (cnt mod printcnt) eq 0 then begin
      print,"Elapsed time = ",tims[cnt]/60.,' mins'
      print,'Penalty function at best vertex = ',min_values[cnt]
      if KEYWORD_SET(magmovset) then begin
        bestmagt=EXTRACT_TRANSFORM(REFORM(simplex[*,lowest]),magt,maxlvar)
;        newmag=OPTIMIZATION_INV_TRANSFORM(bestmagt,0)
        magmov=[[[magmov]],[[newmag]]]
      endif
    endif
    cnt++
  endwhile
  print,'phase_varying_amoeba failed to converge to a solution in ', $
    'the specified number of iterations.'
  t = y[0] & y[0] = y[lowest] & y[lowest] = t ;Sort so fcn min is 0th elem
  t = simplex[*,lowest] & simplex[*,lowest] = simplex[*,0] & simplex[*,0] = t
  newmagt=EXTRACT_TRANSFORM(REFORM(t),magt,maxlvar)
  result=INV_SPHERICAL_TRANSFORM(newmagt,cth)
;  result=OPTIMIZATION_INV_TRANSFORM(newmagt,0)  
  noreset=0
  return, result    ; the function failed to converge.
end
